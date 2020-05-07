proj_info = (id = "cov2",
             data_ver = "20200506",
             calib_ver = "20200506",
             modelobj = :protgroup,
             msinstrument = "QX7",
             quantobj = :protgroup,
             quanttype = :intensity,
             quantcol = :intensity,
             msfolder = "cov2timecourse_dia_20200423")
using Pkg
Pkg.activate(@__DIR__)
using Revise
using RData, StatsBase, DataFrames, BlackBoxOptim, CodecZlib, JSON, CSV, FileIO

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl")
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

@info "Calibrating instrument $(proj_info.msinstrument) for $(proj_info.id) project (data v$(proj_info.data_ver)), v$(proj_info.calib_ver)"

include(joinpath(misc_scripts_path, "bbo_utils.jl"))
include(joinpath(misc_scripts_path, "msinstrument.jl"))
Revise.includet(joinpath(misc_scripts_path, "msinstrument_calib.jl"))
include(joinpath(misc_scripts_path, "frame_utils.jl"))
include(joinpath(misc_scripts_path, "delimdata_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "spectronaut_utils.jl"));

data_path = joinpath(base_analysis_path, proj_info.id, "data")
scratch_path = joinpath(base_analysis_path, proj_info.id, "scratch")

protgroup_data, report_colinfo = SpectronautUtils.read_proteins_report(joinpath(data_path, proj_info.msfolder, "COVID proteome NO normalization_Report.txt"),
                import_data=[:quantity], delim='\t')
quantity_colinfo_df = SpectronautUtils.metrics_colinfo(report_colinfo[:quantity])
msruns_df = unique!(select(quantity_colinfo_df, [:msrun_ix, :rawfile]))
calib_dsname = "20200418_QX7_MaTa_SA_proteome_A549"
rawfilename_matches = match.(Ref(Regex(string("^", calib_dsname, "_(?<msrun>(?<condition>.+)_(?<replicate>\\d+))(?:_\\d+)?\\.raw\$"))),
       levels(msruns_df.rawfile))
msruns_df.msrun = string.(getindex.(rawfilename_matches, :msrun))
msruns_df.replicate = parse.(Int, getindex.(rawfilename_matches, :replicate))
msruns_df.condition = string.(getindex.(rawfilename_matches, :condition))
protgroups_df = select(protgroup_data, report_colinfo[:protgroup])
protgroups_df.is_contaminant = occursin.(Ref(r"(?:^|;)CON__"), protgroups_df.protgroup_sn_id)
protgroups_df.use_for_calib = .!protgroups_df.is_contaminant

pg_intensities_df = FrameUtils.pivot_longer(protgroup_data, :protgroup_id,
                        measure_vars_regex=r"^(?<value>[^.]+)\.\[(?<var>\d+)\]\s(?<rawfile>.+)?$",
                        var_col=:msrun_ix, value_col=:metric)
pg_intensities_df[!, :msrun_ix] = parse.(Int, string.(pg_intensities_df[!, :msrun_ix]))
mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(pg_intensities_df, protgroups_df, msruns_df,
        msrun_col=:msrun_ix, exp_col=:condition,
        missed_intensity_factor = 1.0, # assume intensity was missed in replicate because it was really less abundant
        bin_oversize=2, nbins=40, max_objXexp=5000, error_scale=0.8) # assume 80% of variation comes from measurement

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
MaxTime = 200.0
min_ini_sigma = 1E-2
ini_sigma = 1.0
ini_lnB = nothing

using BlackBoxOptim

@info "Noise Model Calibration..."
bbox_opt = bbsetup(mscalib_data; PopulationSize=1000,
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

msdata_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.msfolder)_$(proj_info.data_ver).RData"))
msdata_dict = msdata_rdata["msdata_full"]
pepmods_df = copy(msdata_rdata["msdata_full"]["pepmods"])
pepmods_df[!, :use_for_calib] .= true # FIXME exclude contaminants

pepmod_mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(
    msdata_dict["pepmod_intensities"], pepmods_df, msdata_dict["msruns"],
    object_col=:pepmod_id, msrun_col=:msrun, exp_col=:condition,
    missed_intensity_factor = 1.0, # assume intensity was missed in replicate because it was really less abundant
    bin_oversize=2, nbins=40, max_objXexp=5000, error_scale=0.33) # assume 1/3 of variation comes from measurement

ini_population = Matrix{Float64}(undef, 0, 0)
MaxTime = 200.0
min_ini_sigma = 1E-2
ini_sigma = 1.0
ini_lnB = nothing

using BlackBoxOptim

@info "Noise Model Calibration..."
bbox_opt = bbsetup(pepmod_mscalib_data; PopulationSize=1000,
                   MaxTime=MaxTime, MaxSteps=10^8, NThreads=Threads.nthreads()-1,
                   MaxStepsWithoutProgress=50000,
                   TraceMode=:verbose, TraceInterval=5.0,
                   #fitness_scheme=ScalarFitnessScheme{false}(),
                   detectionMaxMin=0.98,
                   ϵ=[1.0, 1.0], Population=ini_population);
bbox_res = bboptimize(bbox_opt)
ini_population = BBOUtils.popmatrix(bbox_res)
@info "Done Calibration..."

instr_calib_filename = "instr_$(proj_info.msinstrument)_$(proj_info.quanttype)_pepmod_calib_$(proj_info.id)_$(proj_info.calib_ver)"
save_results(BlackBoxOptim.problem(bbox_opt), bbox_res,
    joinpath(scratch_path, "$(instr_calib_filename).jld2"),
    joinpath(scratch_path, "$(instr_calib_filename).json"))

instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(bbox_opt).factory, best_candidate(bbox_res))

ref_intensities = MSInstrumentCalibration.reference_intensities(mscalib_data)
plot_intensities = exp.((log(first(ref_intensities))-3):0.01:(log(last(ref_intensities))+3))
PlotlyUtils.MSInstrument.plot_rel_sd(instr_calib_model, intensity_range = plot_intensities)
PlotlyUtils.MSInstrument.plot_sd(instr_calib_model, intensity_range = plot_intensities)
PlotlyUtils.MSInstrument.plot_detection(instr_calib_model, intensity_range = plot_intensities)
