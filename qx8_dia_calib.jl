#addprocs(15);
proj_info = (id = "cov2",
             data_ver = "20200411",
             fit_ver = "20200411",
             modelobj = :protgroup,
             msinstrument = "QX8",
             quantobj = :protgroup,
             quanttype = :intensity,
             quantcol = :intensity,
             ms_folder = "spectronaut_oeproteome_20200411")
using Pkg
Pkg.activate(@__DIR__)
using Revise
using StatsBase, DataFrames, BlackBoxOptim, CodecZlib, JSON, CSV, FileIO

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl")
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

@info "Calibrating instrument $(proj_info.msinstrument) for $(proj_info.id) project (data v$(proj_info.data_ver), fit v$(proj_info.fit_ver))"

include(joinpath(misc_scripts_path, "bbo_utils.jl"))
include(joinpath(misc_scripts_path, "msinstrument.jl"))
Revise.includet(joinpath(misc_scripts_path, "msinstrument_calib.jl"))
include(joinpath(misc_scripts_path, "frame_utils.jl"))
include(joinpath(misc_scripts_path, "delimdata_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "spectronaut_utils.jl"));

data_path = joinpath(base_analysis_path, proj_info.id, "data")
scratch_path = joinpath(base_analysis_path, proj_info.id, "scratch")

protgroup_data, report_colinfo = SpectronautUtils.read_proteins_report(joinpath(data_path, proj_info.ms_folder, "20200410_COV2_B1_DIA_DDA library_proteinreport.csv"),
                import_data=[:quantity])
quantity_colinfo_df = SpectronautUtils.metrics_colinfo(report_colinfo[:quantity])
msruns_df = unique!(select(quantity_colinfo_df, [:msrun_ix, :rawfile]))
dsname = "20200326_QX8_OzKa_SA_COV2_FPMS"
rawfilename_matches = match.(Ref(Regex(string("^", dsname, "_(?<bait_code>[^_]+)_(?<replicate>\\d+)(?:_\\d+)?\\.raw\$"))),
       levels(msruns_df.rawfile))
msruns_df[!, :bait_code] .= string.(getindex.(rawfilename_matches, :bait_code))
msruns_df[!, :replicate] .= parse.(Int, getindex.(rawfilename_matches, :replicate))
msruns_df[!, :mstag] .= "Sum"
protgroups_df = select(protgroup_data, report_colinfo[:protgroup])
protgroups_df.is_contaminant = occursin.(Ref(r"(?:^|;)CON__"), protgroups_df.protgroup_sn_id)
protgroups_df.use_for_calib = .!protgroups_df.is_contaminant

pg_intensities_df = FrameUtils.pivot_longer(protgroup_data, :protgroup_id,
                        measure_vars_regex=r"^(?<value>[^.]+)\.\[(?<var>\d+)\]\s(?<rawfile>.+)?$",
                        var_col=:msrun_ix, value_col=:metric)
pg_intensities_df[!, :msrun_ix] = parse.(Int, string.(pg_intensities_df[!, :msrun_ix]))
pg_used_intensities_df = pg_intensities_df[pg_intensities_df.protgroup_id .∈ Ref(Set(sample(protgroups_df[protgroups_df.use_for_calib, :protgroup_id], 1500))), :]
mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(pg_used_intensities_df, protgroups_df, msruns_df,
        mschannel_col=:msrun_ix, exp_col=:bait_code, tech_repl_col=:replicate)

instr_calib_filename = "instr_$(proj_info.msinstrument)_$(proj_info.quanttype)_$(proj_info.quantobj)_calib_$(proj_info.id)_$(proj_info.data_ver)"

using JLD2, JSON

function save_results(calib_problem, calib_result, jld_file, json_file)
    calib_model = MSInstrumentCalibration.params2model(calib_problem.factory, best_candidate(calib_result))
    @info "Best instrument params: $calib_model"

    @info "Saving optimization results into $jld_file..."
    @save(jld_file, mscalib_data, calib_result, calib_model)

    @info "Saving JSON $json_file..."
    open(json_file, "w") do io JSON.print(io, Dict{Symbol,Any}(
        :raw_params => best_candidate(calib_result),
        :instr_calib => calib_model,
        :fitness => best_fitness(calib_result)
        ))
    end
end

res_jld_file = joinpath(scratch_path, "$(instr_calib_filename)_borg.jld2")
res_json_file = joinpath(scratch_path, "$(instr_calib_filename)_borg.json")

ini_population = Matrix{Float64}(undef, 0, 0)

MaxTime = 600.0
min_ini_sigma = 1E-2
ini_sigma = 1.0
ini_lnB = nothing

using BlackBoxOptim

@info "Running Borg..."
flush(STDERR)
flush(STDOUT)
borg_opt = bbsetup(mscalib_data; PopulationSize=1000,
                   MaxTime=MaxTime, MaxSteps=10^8, NThreads=Threads.nthreads()-1,
                   TraceMode=:verbose, TraceInterval=5.0,
                   ϵ=[1E-7, 1E-7], Population=ini_population);
borg_res = bboptimize(borg_opt)
ini_population = BBOUtils.popmatrix(borg_res)
@info "Done Borg..."
save_results(BlackBoxOptim.problem(borg_opt), borg_res, res_jld_file, res_json_file)
flush(STDERR)
flush(STDOUT)

@info "Done"

using PlotlyJS
instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(borg_opt).factory, best_candidate(borg_res))

log2_intensity_range = 2:0.01:35.0
intensity_range = 2 .^ log2_intensity_range
plot([scatter(x = intensity_range, name="rel_sd",
              y = MSInstrument.signal_std.(Ref(instr_calib_model), intensity_range) ./ intensity_range),
      scatter(x = intensity_range, y = fill(1.0, length(intensity_range)), name="1", mode=:lines, line_dash=:dash)],
     Layout(xaxis=attr(type="log"), yaxis=attr(type="log")))
plot([scatter(x = intensity_range, name="sd", legendgroup="sd", showlegend=true, line_color=:blue,
              y = intensity_range .+ MSInstrument.signal_std.(Ref(instr_calib_model), intensity_range)),
      scatter(x = intensity_range, name="sd", legendgroup="sd", showlegend=false, line_color=:blue,
              y = intensity_range .- MSInstrument.signal_std.(Ref(instr_calib_model), intensity_range)),
      scatter(x = intensity_range, y = intensity_range, name="I", mode=:lines, line_color=:orangered, line_dash=:dash),
      scatter(x = intensity_range, y = 2intensity_range, name="I+I", mode=:lines, line_color=:orangered, line_dash=:dot),
      scatter(x = intensity_range, y = 0.01intensity_range, name="I+I", mode=:lines, line_color=:orangered, line_dash=:dot)],
      Layout(xaxis=attr(type="log"), yaxis=attr(type="log")))

plot([scatter(x = intensity_range, name="detected", line_color=:red,
              y = exp.(MSInstrument.detection_likelihood_log.(Ref(instr_calib_model), true, log.(intensity_range)))),
      scatter(x = intensity_range, name="missed", line_color=:blue,
              y = exp.(MSInstrument.detection_likelihood_log.(Ref(instr_calib_model), false, log.(intensity_range))))],
    Layout(xaxis=attr(type="log")))

plot_intensities_df = @where(instr_calib_data.intensities_df, !isnull(:expected_log) & !isnull(:error));
plot_intensities_df[:log_error] = log10.(abs.(dropnull(plot_intensities_df[:error])));
scatter(dropnull(plot_intensities_df[:expected_log]),
        dropnull(plot_intensities_df[:log_error]))

#=
using Cairo, Gadfly

draw(PDF( joinpath(scratch_path,"$(instr_calib_filename)_intensities.pdf"), 16inch, 12inch),
plot( instr_calib_data.intensities_df, x="expected_log", y = "error",
      Geom.point ) )

PCP.signal_precision(instrument, 40.0)

draw(PDF( joinpath(scratch_path,"$(instr_calib_filename)_std.pdf"), 16inch, 12inch),
plot( [x->x, x -> PCP.signal_std(instrument, x)],
#plot( [x->x, x -> x + 1.0/PCP.signal_precision(instrument, x), x -> x - 1.0/PCP.signal_precision(instrument, x) ],
      1,1E+7 ) )

names(instr_calib_jlser["features"])

keys(instr_calib_jlser)

using ProfileView

ProfileView.view()
=#
