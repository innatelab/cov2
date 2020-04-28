proj_info = (id = "cov2",
             data_ver = "20200428",
             fit_ver = "20200428",
             modelobj = :protgroup,
             msinstrument = "QX7",
             quantobj = :protgroup,
             quanttype = :intensity,
             quantcol = :intensity,
             ms_folder = "cov2timecourse_dia_20200423")
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

protgroup_data, report_colinfo = SpectronautUtils.read_proteins_report(joinpath(data_path, proj_info.ms_folder, "COVID proteome NO normalization_Report.txt"),
                import_data=[:quantity], delim='\t')
quantity_colinfo_df = SpectronautUtils.metrics_colinfo(report_colinfo[:quantity])
msruns_df = unique!(select(quantity_colinfo_df, [:msrun_ix, :rawfile]))
calib_dsname = "20200418_QX7_MaTa_SA_proteome_A549"
rawfilename_matches = match.(Ref(Regex(string("^", calib_dsname, "_(?<msrun>(?<condition>.+)_(?<replicate>\\d+))(?:_\\d+)?\\.raw\$"))),
       levels(msruns_df.rawfile))
msruns_df.msrun = string.(getindex.(rawfilename_matches, :msrun))
msruns_df.replicate = parse.(Int, getindex.(rawfilename_matches, :replicate))
msruns_df.condition = string.(getindex.(rawfilename_matches, :condition))
msruns_df[!, :mstag] .= "F"
protgroups_df = select(protgroup_data, report_colinfo[:protgroup])
protgroups_df.is_contaminant = occursin.(Ref(r"(?:^|;)CON__"), protgroups_df.protgroup_sn_id)
protgroups_df.use_for_calib = .!protgroups_df.is_contaminant

pg_intensities_df = FrameUtils.pivot_longer(protgroup_data, :protgroup_id,
                        measure_vars_regex=r"^(?<value>[^.]+)\.\[(?<var>\d+)\]\s(?<rawfile>.+)?$",
                        var_col=:msrun_ix, value_col=:metric)
pg_intensities_df[!, :msrun_ix] = parse.(Int, string.(pg_intensities_df[!, :msrun_ix]))
pg_used_intensities_df = pg_intensities_df[pg_intensities_df.protgroup_id .∈ Ref(Set(sample(protgroups_df[protgroups_df.use_for_calib, :protgroup_id], 3000))), :]
mscalib_data = MSInstrumentCalibration.MSErrorCalibrationData(pg_used_intensities_df, protgroups_df, msruns_df,
        mschannel_col=:msrun_ix, exp_col=:condition, tech_repl_col=:replicate,
        error_scale=0.8) # assume 80% of variation comes from measurement

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
bbox_opt = bbsetup(mscalib_data; PopulationSize=300,
                   MaxTime=MaxTime, MaxSteps=10^8, #NThreads=Threads.nthreads()-1,
                   MaxStepsWithoutProgress=50000,
                   TraceMode=:verbose, TraceInterval=5.0,
                   fitness_scheme=ScalarFitnessScheme{false}(),
                   detectionMaxMin=0.98,
                   ϵ=[1E-7, 1E-7], Population=ini_population);
bbox_res = bboptimize(bbox_opt)
ini_population = BBOUtils.popmatrix(bbox_res)
@info "Done Calibration..."
save_results(BlackBoxOptim.problem(bbox_opt), bbox_res, res_jld_file, res_json_file)

@info "Done"

Revise.includet(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))

instr_calib_model = MSInstrumentCalibration.params2model(BlackBoxOptim.problem(bbox_opt).factory, best_candidate(bbox_res))

PlotlyUtils.MSInstrument.plot_rel_sd(instr_calib_model)
PlotlyUtils.MSInstrument.plot_sd(instr_calib_model)
PlotlyUtils.MSInstrument.plot_detection(instr_calib_model, intensity_range = 2 .^ (2:0.01:20))
