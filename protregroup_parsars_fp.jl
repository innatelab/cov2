proj_info = (id = "cov2",
             data_ver = "20200829",
             fit_ver = "20200829",
             msfolder = "snaut_parsars_fp_20200829")
using Pkg
Pkg.activate(@__DIR__)

using Revise
using CSV, RData, DataFrames

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const scratch_path = joinpath(analysis_path, "scratch")

Revise.includet(joinpath(misc_scripts_path, "protgroup_assembly.jl"));
Revise.includet(joinpath(misc_scripts_path, "protgroup_crossmatch.jl"));

pepmods_rdata = load(joinpath(data_path, proj_info.msfolder,
                              "$(proj_info.id)_$(proj_info.msfolder)_$(proj_info.data_ver)_pepmods.RData"))
pepmods_df = pepmods_rdata["pepmods.df"]
proteins_df = pepmods_rdata["all_proteins.df"]
pepmod2acs = Dict(r.pepmod_id => (Set{String}(split(r.protein_acs, ';')), r.pepmod_rank)
                  for r in eachrow(pepmods_df))

protregroups_acs = ProtgroupAssembly.assemble_protgroups(pepmod2acs, verbose=true,
                                                         nspec_peptides=2, rspec_peptides=0.25)
protregroups_acs_df = ProtgroupAssembly.dataframe(protregroups_acs,
                            protein_ranks=Dict(r.protein_ac => ProtgroupXMatch.rank_uniprot(r) for r in eachrow(proteins_df)),
                            protgroup_col=:protregroup_id, protein_col=:protein_ac)

CSV.write(joinpath(data_path, proj_info.msfolder,
                   "$(proj_info.id)_$(proj_info.msfolder)_$(proj_info.data_ver)_protregroups_acs.txt"),
          protregroups_acs_df, delim='\t', missingstring="")
