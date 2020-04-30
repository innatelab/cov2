proj_info = (id = "cov2",
             data_ver = "20200423",
             calib_ver = "20200430",
             modelobj = :ptmgroup,
             msinstrument = "QX7",
             quantobj = :ptmgroup,
             quanttype = :intensity,
             quantcol = :intensity,
             ms_folder = "phospho_qx7_calib_20200423",
             #ms_folder = "cov2timecourse_phospho_dia_20200423"
             )
using Pkg
Pkg.activate(@__DIR__)
using Revise
using StatsBase, DataFrames, BlackBoxOptim, CodecZlib, JSON, CSV, FileIO

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl")
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

@info "Calibrating instrument $(proj_info.msinstrument) for $(proj_info.id) project (data v$(proj_info.data_ver)), v$(proj_info.calib_ver)"

Revise.includet(joinpath(misc_scripts_path, "bbo_utils.jl"))
include(joinpath(misc_scripts_path, "msinstrument.jl"))
Revise.includet(joinpath(misc_scripts_path, "msinstrument_calib.jl"))
include(joinpath(misc_scripts_path, "frame_utils.jl"))
include(joinpath(misc_scripts_path, "delimdata_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "spectronaut_utils.jl"));

data_path = joinpath(base_analysis_path, proj_info.id, "data")
scratch_path = joinpath(base_analysis_path, proj_info.id, "scratch")

calib_dsname = "20200417_QX7_MaTa_SA_A549"
calib_filename = "COV2_DIA_phospho_control samples.txt"
#calib_filename = "COV2_DIA_phospho_0.75probablity_no normalization.txt"
calib_data_df = CSV.read(joinpath(data_path, proj_info.ms_folder, calib_filename),
                         delim='\t', truestrings=["True"], falsestrings=["False"])
msruns_df = DataFrame(report_colname = names(calib_data_df)[occursin.(calib_dsname, string.(names(calib_data_df)))])
msruns_df.rawfile = string.(msruns_df.report_colname)
rawfilename_matches = match.(Ref(Regex(string("^", calib_dsname, "_(?<msrun>(?<condition>phospho_test)_(?<replicate>\\d+))\$"))),
       msruns_df.rawfile)
#rawfilename_matches = match.(Ref(Regex(string("^", calib_dsname, "_(?<msrun>(?<condition>(?:SARS_COV2|mock)_\\d+h(?:pi)?)_(?<replicate>\\d+))\$"))),
#       msruns_df.rawfile)
msruns_df.msrun = string.(getindex.(rawfilename_matches, :msrun))
msruns_df.replicate = parse.(Int, getindex.(rawfilename_matches, :replicate))
msruns_df.condition = string.(getindex.(rawfilename_matches, :condition))
msruns_df[!, :mstag] .= "F"
ptms_df = select(calib_data_df, Not(msruns_df.report_colname))
rename!(n -> replace(replace(string(n), r"^[NT]:\s" => ""), "." => "_"), ptms_df)
ptms_df.use_for_calib = .!ptms_df.EG_IsDecoy
ptms_df.ptmgroup_id = categorical(ptms_df.PTM_collapse_key)
calib_data_df.ptmgroup_id = copy(ptms_df.ptmgroup_id)
ptmgroup_intensities_df = FrameUtils.pivot_longer(calib_data_df, :ptmgroup_id,
                        measure_vars_regex=Regex("^(?<var>$(calib_dsname)\\w+)\$"),
                        var_col=:rawfile, value_col=:log2_intensity)
filter!(r -> isfinite(r.log2_intensity), ptmgroup_intensities_df)
ptmgroup_intensities_df.intensity = 2 .^ ptmgroup_intensities_df.log2_intensity
ptmgroup_intensities_df = join(ptmgroup_intensities_df, select(msruns_df, [:report_colname, :msrun, :replicate, :condition, :mstag]),
                          on=:rawfile=>:report_colname, kind=:left)
select!(ptmgroup_intensities_df, Not(:rawfile))
mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(ptmgroup_intensities_df, ptms_df, msruns_df,
        object_col=:ptmgroup_id, msrun_col=:msrun, exp_col=:condition,
        missed_intensity_factor = 0.9, # assume intensity was missed in replicate because it was really less abundant
        bin_oversize=5, nbins=40, error_scale=0.8) # assume 80% of variation comes from measurement

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

instr_calib_filename = "instr_$(proj_info.msinstrument)_$(proj_info.quanttype)_$(proj_info.quantobj)_calib_$(proj_info.id)_$(proj_info.calib_ver)"
res_jld_file = joinpath(scratch_path, "$(instr_calib_filename).jld2")
res_json_file = joinpath(scratch_path, "$(instr_calib_filename).json")

ini_population = Matrix{Float64}(undef, 0, 0)

MaxTime = 300.0
min_ini_sigma = 1E-2
ini_sigma = 1.0
ini_lnB = nothing

using BlackBoxOptim

@info "Noise Model Calibration..."
bbox_opt = bbsetup(mscalib_data; PopulationSize=500,
                   MaxTime=MaxTime, MaxSteps=10^8, NThreads=Threads.nthreads()-1,
                   MaxStepsWithoutProgress=50000,
                   TraceMode=:verbose, TraceInterval=5.0,
                   #fitness_scheme=ScalarFitnessScheme{false}(),
                   detectionMaxMin=0.98,
                   ϵ=[1.0, 1.0], Population=ini_population);
bbox_res = bboptimize(bbox_opt)
ini_population = BBOUtils.popmatrix(bbox_res)
@info "Done Calibration..."
save_results(BlackBoxOptim.problem(bbox_opt), bbox_res, res_jld_file, res_json_file)

@info "Done"

Revise.includet(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))

@load("/pool/analysis/astukalov/cov2/scratch/instr_QX7_intensity_ptmgroup_calib_cov2_20200428_borg.jld2", calib_model)
instr_calib_model = calib_model
instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(bbox_opt).factory, best_candidate(bbox_res))

ref_intensities = MSInstrumentCalibration.reference_intensities(mscalib_data)
plot_intensities = exp.((log(first(ref_intensities))-3):0.01:(log(last(ref_intensities))+3))
PlotlyUtils.MSInstrument.plot_rel_sd(instr_calib_model, intensity_range = plot_intensities)
PlotlyUtils.MSInstrument.plot_sd(instr_calib_model, intensity_range = plot_intensities)
PlotlyUtils.MSInstrument.plot_detection(instr_calib_model, intensity_range = plot_intensities)
