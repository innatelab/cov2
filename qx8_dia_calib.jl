#addprocs(15);
proj_info = (id = "cov2",
             data_ver = "20200923",
             calib_ver = "20200923",
             modelobj = :pepmodstate,
             msinstrument = "QX8",
             quantobj = :pepmodstate,
             quanttype = :intensity,
             quantcol = :intensity,
             msfolder = "snaut_oefp_20200923",
             qvalue_max = 1E-3)
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

include(joinpath(misc_scripts_path, "bbo_utils.jl"))
include(joinpath(misc_scripts_path, "msinstrument.jl"))
includet(joinpath(misc_scripts_path, "msinstrument_calib.jl"))
include(joinpath(misc_scripts_path, "frame_utils.jl"))
include(joinpath(misc_scripts_path, "delimdata_utils.jl"))
include(joinpath(misc_scripts_path, "spectronaut_utils.jl"));

data_path = joinpath(base_analysis_path, proj_info.id, "data")
scratch_path = joinpath(base_analysis_path, proj_info.id, "scratch")

pms_data, report_colinfo = SpectronautUtils.read_report(joinpath(data_path, proj_info.msfolder, "EG_Report.xls"),
                import_data=[:peptide, :pepmodstate, :quantity])
pms_data.pepmodstate_id = 1:nrow(pms_data)
quantity_colinfo_df = SpectronautUtils.metrics_colinfo(last.(report_colinfo[:quantity]))
msruns_df = unique!(select(quantity_colinfo_df, [:rawfile_ix, :rawfile]))
calib_dsname = "20200426_QX8_OzKa_SA_COV2_FPMS"
rawfilename_matches = match.(Ref(Regex(string("^", calib_dsname, "_(?<bait_code>LMix_|B5_[ABCD])(?<replicate>\\d+)\\.raw\$"))),
       levels(msruns_df.rawfile))
msruns_df.is_used = .!isnothing.(rawfilename_matches)
sum(msruns_df.is_used)
msruns_df.bait_code = [isnothing(m) ? missing : string(m[:bait_code]) for m in rawfilename_matches]
msruns_df.replicate = [isnothing(m) ? missing : parse(Int, m[:replicate]) for m in rawfilename_matches]

pepmodstates_df = select(pms_data, vcat(last.(report_colinfo[:peptide]), last.(report_colinfo[:pepmodstate])));
pepmodstates_df.is_contaminant = occursin.(Ref(r"(?:^|;)CON__"), pepmodstates_df.peptide_protein_acs);
pepmodstates_df.use_for_calib = .!pepmodstates_df.is_contaminant;
pepmodstates_df

pms_intensities_df = SpectronautUtils.pivot_longer(pms_data, :pepmodstate_id)
rename!(pms_intensities_df, :EG_intensity => :intensity_norm)
pms_intensities_df.intensity = pms_intensities_df.intensity_norm ./ pms_intensities_df.EG_normfactor
rawfile_normfactors_df = filter!(r -> !ismissing(r.EG_normfactor), select(pms_intensities_df, [:rawfile_ix, :EG_normfactor])) |> unique!
msruns_df = leftjoin(msruns_df, rawfile_normfactors_df, on=:rawfile_ix)
msruns_df.msrun_mult = inv.(msruns_df.EG_normfactor)
used_pepmodstate_ids = Set(pepmodstates_df[pepmodstates_df.use_for_calib, :pepmodstate_id])
used_rawfile_ixs = Set(msruns_df[msruns_df.is_used, :rawfile_ix])
pms_used_intensities_df = filter(r -> (r.pepmodstate_id ∈ used_pepmodstate_ids) &&
                                      (r.rawfile_ix ∈ used_rawfile_ixs) &&
                                      (coalesce(r.EG_qvalue, 1) <= proj_info.qvalue_max), pms_intensities_df)
mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(pms_used_intensities_df, pepmodstates_df, msruns_df,
        msrun_col=:rawfile_ix, exp_col=:bait_code, intensity_col=:intensity, object_col=:pepmodstate_id,
        missed_intensity_factor = 1.0, # assume intensity was missed in replicate because it was really less abundant
        bin_oversize=3, nbins=30, max_objXexp=5000, error_scale=0.8) # assume 80% of variation comes from measurement

ini_population = Matrix{Float64}(undef, 0, 0)

MaxTime = 600.0
min_ini_sigma = 1E-2
ini_sigma = 1.0
ini_lnB = nothing

using BlackBoxOptim

mscalibprefix = MSInstrumentCalibration.calibfileprefix(proj_info)

@info "Noise Model Calibration..."
bbox_opt = bbsetup(mscalib_data; PopulationSize=1000,
                   MaxTime=MaxTime, MaxSteps=10^8, NThreads=Threads.nthreads()-1,
                   MaxStepsWithoutProgress=50000,
                   TraceMode=:verbose, TraceInterval=5.0,
                   #fitness_scheme=ScalarFitnessScheme{false}(),
                   detectionMaxMin=0.98,
                   ϵ=[0.3, 1.0], Population=ini_population);
bbox_res = bboptimize(bbox_opt)
ini_population = BBOUtils.popmatrix(bbox_res)
@info "Done Calibration..."
MSInstrumentCalibration.save_result(joinpath(scratch_path, "$(mscalibprefix).json"),
                                    BlackBoxOptim.problem(bbox_opt), bbox_res, verbose=true,
                                    info=proj_info, data=mscalib_data)
MSInstrumentCalibration.save_result(joinpath(scratch_path, "$(mscalibprefix).jld2"),
                                    BlackBoxOptim.problem(bbox_opt), bbox_res, verbose=true,
                                    info=proj_info, data=mscalib_data)

@info "Done"

using PlotlyBase
Revise.includet(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))
@load("/pool/analysis/astukalov/cov2/scratch/instr_QX8_intensity_protgroup_calib_cov2_20200411_borg.jld2", calib_model)
@load("/pool/analysis/astukalov/cov2/scratch/instr_QX7_intensity_protgroup_calib_cov2_20200430.jld2", calib_model)
instr_calib_model = calib_model

instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(bbox_opt).factory, best_candidate(bbox_res))

ref_intensities = MSInstrumentCalibration.reference_intensities(mscalib_data)
plot_intensities = exp.((log(first(ref_intensities))-3):0.01:(log(last(ref_intensities))+3))

calib_plot_prefix = "instr_QX7_intensity_protgroup_calib_cov2_20200430"
calib_plot_prefix = mscalibprefix
relsd_plot = PlotlyUtils.MSInstrument.plot_rel_sd(instr_calib_model, intensity_range = plot_intensities)
savefig(relsd_plot.plot, joinpath(analysis_path, "plots", "$(calib_plot_prefix)_relsd.svg"))

detect_plot = PlotlyUtils.MSInstrument.plot_detection(instr_calib_model, intensity_range = plot_intensities)
savefig(detect_plot.plot, joinpath(analysis_path, "plots", "$(calib_plot_prefix)_detect.svg"))

PlotlyUtils.MSInstrument.plot_sd(instr_calib_model, intensity_range = plot_intensities)
PlotlyUtils.MSInstrument.plot_detection(instr_calib_model, intensity_range = plot_intensities)
