proj_info = (id = "cov2",
             data_ver = "20200425",
             fit_ver = "20200425",
             modelobj = :protgroup,
             msinstrument = "QEP5",
             quantobj = :pepmodstate,
             quanttype = :intensity,
             quantcol = :intensity,
             msfolder = "qep5_apms_calib_20200425")
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
include(joinpath(misc_scripts_path, "mq_utils.jl"))

data_path = joinpath(base_analysis_path, proj_info.id, "data")
scratch_path = joinpath(base_analysis_path, proj_info.id, "scratch")

using CSV

msruns_df = CSV.read(joinpath(data_path, proj_info.msfolder, "combined", "experimentalDesign.txt"), delim='\t')
rename!(msruns_df, :Experiment => :msrun)
msruns_df.experiment = replace.(msruns_df.msrun, r"_\d+ul_\d+$" => "")
msruns_df[!, :tech_replicate] .= 1
msruns_df.mschannel = copy(msruns_df.msrun)
msruns_df[!, :mstag] .= "F"
for techrepls_df in groupby(msruns_df, :experiment)
    techrepls_df[:, :tech_replicate] .= 1:nrow(techrepls_df)
end
mqevidence = MaxQuantUtils.read_evidence(joinpath(data_path, proj_info.msfolder, "combined", "txt"))
mqevidence.observations.mschannel = copy(mqevidence.observations.msrun)
#protgroups_df, protgroups_colgroups = MaxQuantUtils.read_protgroups(joinpath(data_path, "combined", "txt"))

mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(mqevidence, msruns_df,
                                                              object=proj_info.quantobj, intensity_col=proj_info.quantcol,
                                                              max_objects = 10000)

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

MaxTime = 3600.0
min_ini_sigma = 1E-2
ini_sigma = 1.0
ini_lnB = nothing

using BlackBoxOptim

@info "Running Borg..."
borg_opt = bbsetup(mscalib_data; PopulationSize=200,
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

log2_intensity_range = 12:0.01:35.0
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
