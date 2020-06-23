proj_info = (id = "cov2",
             data_ver = "20200525",
             fit_ver = "20200525",
             ms_folder = "mq_apms_20200525")
using Pkg
Pkg.activate(@__DIR__)

using Revise
using CSV, RData, DataFrames

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "protgroup_assembly.jl"));
Revise.includet(joinpath(misc_scripts_path, "protgroup_crossmatch.jl"));

pepmods_rdata = load(joinpath(data_path, proj_info.ms_folder,
                              "$(proj_info.id)_$(proj_info.ms_folder)_$(proj_info.data_ver)_pepmods.RData"))
pepmods_df = pepmods_rdata["pepmods.df"]
proteins_df = pepmods_rdata["proteins.df"]
pepmod2protgroups = Dict(r.pepmod_id => (Set(parse.(Int, split(r.protgroup_ids, ';'))), r.pepmod_rank)
                         for r in eachrow(pepmods_df))
pepmod2acs = Dict(r.pepmod_id => (Set{String}(split(coalesce(r.protein_acs, r.lead_protein_acs), ';')), r.pepmod_rank)
                  for r in eachrow(pepmods_df))
protregroups = ProtgroupAssembly.assemble_protgroups(pepmod2protgroups, verbose=true,
                                                     nspec_peptides=2, rspec_peptides=0.25)
protregroups_df = ProtgroupAssembly.dataframe(protregroups,
                            protgroup_col=:protregroup_id, protein_col=:protgroup_id)

protregroups_acs = ProtgroupAssembly.assemble_protgroups(pepmod2acs, verbose=true,
                                                         nspec_peptides=2, rspec_peptides=0.25)
protregroups_acs_df = ProtgroupAssembly.dataframe(protregroups_acs,
                            protein_ranks=Dict(r.protein_ac => ProtgroupXMatch.rank_uniprot(r) for r in eachrow(proteins_df)),
                            protgroup_col=:protregroup_id, protein_col=:protein_ac)

CSV.write(joinpath(data_path, proj_info.ms_folder,
                   "$(proj_info.id)_$(proj_info.ms_folder)_$(proj_info.data_ver)_protregroups.txt"),
          protregroups_df, delim='\t', missingstring="")

CSV.write(joinpath(data_path, proj_info.ms_folder,
                   "$(proj_info.id)_$(proj_info.ms_folder)_$(proj_info.data_ver)_protregroups_acs.txt"),
          protregroups_acs_df, delim='\t', missingstring="")
