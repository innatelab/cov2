#addprocs(15);
proj_info = (id = "cov2",
             data_ver = "20200423",
             fit_ver = "20200423",
             modelobj = :ptmgroup,
             msinstrument = "QX7",
             quantobj = :ptmgroup,
             quanttype = :intensity,
             quantcol = :intensity,
             ms_folder = "phospho_qx7_calib_20200423")
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

calib_data_df = CSV.read(joinpath(data_path, proj_info.ms_folder, "COV2_DIA_phospho_control samples.txt"),
                         delim='\t', truestrings=["True"], falsestrings=["False"])
dsname = "20200417_QX7_MaTa_SA_A549"
msruns_df = DataFrame(report_colname = names(calib_data_df)[occursin.(dsname, string.(names(calib_data_df)))])
msruns_df.rawfile = string.(msruns_df.report_colname)
rawfilename_matches = match.(Ref(Regex(string("^", dsname, "_(?<msrun>phospho_test_(?<replicate>\\d+))\$"))),
       msruns_df.rawfile)
msruns_df.msrun = string.(getindex.(rawfilename_matches, :msrun))
msruns_df.replicate = parse.(Int, getindex.(rawfilename_matches, :replicate))
msruns_df[!, :condition] .= dsname
msruns_df[!, :mstag] .= "F"
ptms_df = select(calib_data_df, Not(msruns_df.report_colname))
rename!(n -> replace(replace(string(n), r"^[NT]:\s" => ""), "." => "_"), ptms_df)
ptms_df.use_for_calib = .!ptms_df.EG_IsDecoy
ptms_df.ptmgroup_id = categorical(ptms_df.PTM_collapse_key)
calib_data_df.ptmgroup_id = copy(ptms_df.ptmgroup_id)
ptmgroup_intensities_df = FrameUtils.pivot_longer(calib_data_df, :ptmgroup_id,
                        measure_vars_regex=Regex("^(?<var>$(dsname)\\w+)\$"),
                        var_col=:rawfile, value_col=:log2_intensity)
filter!(r -> isfinite(r.log2_intensity), ptmgroup_intensities_df)
ptmgroup_intensities_df.intensity = 2 .^ ptmgroup_intensities_df.log2_intensity
ptmgroup_intensities_df = join(ptmgroup_intensities_df, select(msruns_df, [:report_colname, :msrun, :replicate, :condition, :mstag]),
                          on=:rawfile=>:report_colname, kind=:left)
select!(ptmgroup_intensities_df, Not(:rawfile))
ptm_used_intensities_df = ptmgroup_intensities_df[ptmgroup_intensities_df.ptmgroup_id .∈ Ref(Set(sample(ptms_df[ptms_df.use_for_calib, :ptmgroup_id], 10000))), :]
mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(ptm_used_intensities_df, ptms_df, msruns_df,
        object_col=:ptmgroup_id, mschannel_col=:msrun, exp_col=:condition, tech_repl_col=:replicate)

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
borg_opt = bbsetup(mscalib_data; PopulationSize=1000,
                   MaxTime=MaxTime, MaxSteps=10^8, NThreads=Threads.nthreads()-1,
                   MaxStepsWithoutProgress=50000,
                   TraceMode=:verbose, TraceInterval=5.0,
                   ϵ=[1E-7, 1E-7], Population=ini_population);
borg_res = bboptimize(borg_opt)
ini_population = BBOUtils.popmatrix(borg_res)
@info "Done Borg..."
save_results(BlackBoxOptim.problem(borg_opt), borg_res, res_jld_file, res_json_file)

@info "Done"

using PlotlyJS
instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(borg_opt).factory, best_candidate(borg_res))

log2_intensity_range = 2:0.01:25.0
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
      scatter(x = intensity_range, y = 1.5intensity_range, name="I+I", mode=:lines, line_color=:orangered, line_dash=:dot),
      scatter(x = intensity_range, y = 0.66intensity_range, name="I+I", mode=:lines, line_color=:orangered, line_dash=:dot)],
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
