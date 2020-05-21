#addprocs(15);
proj_info = (id = "cov2",
             data_ver = "20200519",
             calib_ver = "20200519",
             modelobj = :protgroup,
             msinstrument = "QX8",
             quantobj = :protgroup,
             quanttype = :intensity,
             quantcol = :intensity,
             msfolder = "spectronaut_oeproteome_20200519")
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
Revise.includet(joinpath(misc_scripts_path, "msinstrument_calib.jl"))
include(joinpath(misc_scripts_path, "frame_utils.jl"))
include(joinpath(misc_scripts_path, "delimdata_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "spectronaut_utils.jl"));

data_path = joinpath(base_analysis_path, proj_info.id, "data")
scratch_path = joinpath(base_analysis_path, proj_info.id, "scratch")

protgroup_data, report_colinfo = SpectronautUtils.read_proteins_report(joinpath(data_path, proj_info.msfolder, "oe_proteome_unnormalized_20200519_dedup.txt"),
                import_data=[:quantity])
quantity_colinfo_df = SpectronautUtils.metrics_colinfo(report_colinfo[:quantity])
msruns_df = unique!(select(quantity_colinfo_df, [:msrun_ix, :rawfile]))
calib_dsname = "20200426_QX8_OzKa_SA_COV2_FPMS"
rawfilename_matches = match.(Ref(Regex(string("^", calib_dsname, "_(?<bait_code>LMix_|B5_[ABCD])(?<replicate>\\d+)\\.raw\$"))),
       levels(msruns_df.rawfile))
msruns_df.is_used = .!isnothing.(rawfilename_matches)
sum(msruns_df.is_used)
msruns_df[!, :bait_code] .= [isnothing(m) ? missing : string(getindex(m, :bait_code)) for m in rawfilename_matches]
msruns_df[!, :replicate] .= [isnothing(m) ? missing : parse(Int, getindex(m, :replicate)) for m in rawfilename_matches]
protgroups_df = select(protgroup_data, report_colinfo[:protgroup])
protgroups_df.is_contaminant = occursin.(Ref(r"(?:^|;)CON__"), protgroups_df.protgroup_sn_id)
protgroups_df.use_for_calib = .!protgroups_df.is_contaminant

pg_intensities_df = FrameUtils.pivot_longer(protgroup_data, :protgroup_id,
                        measure_vars_regex=r"^(?<value>[^.]+)\.\[(?<var>\d+)\]\s(?<rawfile>.+)?$",
                        var_col=:msrun_ix, value_col=:metric)
pg_intensities_df[!, :msrun_ix] = parse.(Int, string.(pg_intensities_df[!, :msrun_ix]))
pg_used_intensities_df = pg_intensities_df[(pg_intensities_df.protgroup_id .∈ Ref(Set(protgroups_df[protgroups_df.use_for_calib, :protgroup_id]))) .&
                                           (pg_intensities_df.msrun_ix .∈ Ref(Set(msruns_df[msruns_df.is_used, :msrun_ix]))), :]
mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(pg_used_intensities_df, protgroups_df, msruns_df,
        msrun_col=:msrun_ix, exp_col=:bait_code,
        missed_intensity_factor = 1.0, # assume intensity was missed in replicate because it was really less abundant
        bin_oversize=3, nbins=40, max_objXexp=5000, error_scale=0.8) # assume 80% of variation comes from measurement

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
@load("/pool/analysis/astukalov/cov2/scratch/instr_QX8_intensity_protgroup_calib_cov2_20200411_borg.jld2", calib_model)
@load("/pool/analysis/astukalov/cov2/scratch/instr_QX7_intensity_protgroup_calib_cov2_20200430.jld2", calib_model)
instr_calib_model = calib_model

instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(bbox_opt).factory, best_candidate(bbox_res))

ref_intensities = MSInstrumentCalibration.reference_intensities(mscalib_data)
plot_intensities = exp.((log(first(ref_intensities))-3):0.01:(log(last(ref_intensities))+3))
PlotlyUtils.MSInstrument.plot_rel_sd(instr_calib_model, intensity_range = plot_intensities)
PlotlyUtils.MSInstrument.plot_sd(instr_calib_model, intensity_range = plot_intensities)
PlotlyUtils.MSInstrument.plot_detection(instr_calib_model, intensity_range = plot_intensities)
