proj_info = (id = "cov2",
             data_ver = "20200428",
             fit_ver = "20200428",
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
include(joinpath(misc_scripts_path, "msinstrument_calib.jl"))
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
                                                              allowed_ident_types = ["ISO-MSMS", "MULTI-MSMS"],
                                                              max_objects = 3000)

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

instr_calib_filename = "instr_$(proj_info.msinstrument)_$(proj_info.quanttype)_$(proj_info.quantobj)_calib_$(proj_info.id)_$(proj_info.data_ver)"
res_jld_file = joinpath(scratch_path, "$(instr_calib_filename)_borg.jld2")
res_json_file = joinpath(scratch_path, "$(instr_calib_filename)_borg.json")

ini_population = Matrix{Float64}(undef, 0, 0)

MaxTime = 600.0
min_ini_sigma = 1E-2
ini_sigma = 1.0
ini_lnB = nothing

using BlackBoxOptim

@info "Noise Model Calibration..."
bbox_opt = bbsetup(mscalib_data; PopulationSize=1000,
                   MaxTime=MaxTime, MaxSteps=10^8, NThreads=Threads.nthreads()-1,
                   #MaxStepsWithoutProgress=50000,
                   TraceMode=:verbose, TraceInterval=5.0,
                   #fitness_scheme=ScalarFitnessScheme{false}(),
                   detectionMaxMin=0.98,
                   ϵ=[1E-2, 1E-1], Population=ini_population);
bbox_res = bboptimize(bbox_opt)
ini_population = BBOUtils.popmatrix(bbox_res)
@info "Done Calibration..."
save_results(BlackBoxOptim.problem(bbox_opt), bbox_res, res_jld_file, res_json_file)

@info "Done"

include(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))

@load("/pool/analysis/astukalov/cov2/scratch/instr_QEP5_intensity_pepmodstate_calib_cov2_20200425_borg.jld2", calib_model)
instr_calib_model = calib_model
instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(bbox_opt).factory, best_candidate(bbox_res))

PlotlyUtils.MSInstrument.plot_rel_sd(instr_calib_model, intensity_range = 2 .^ (15:0.01:35))
PlotlyUtils.MSInstrument.plot_sd(instr_calib_model, intensity_range = 2 .^ (15:0.01:35))
PlotlyUtils.MSInstrument.plot_detection(instr_calib_model, intensity_range = 2 .^ (15:0.01:35))
