proj_info = (id = "cov2",
             data_ver = "20200923",
             fit_ver = "20200923",
             msfolder = "snaut_oefp_20200923")
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

peptides_rdata = load(joinpath(data_path, proj_info.msfolder,
                              "$(proj_info.id)_$(proj_info.msfolder)_$(proj_info.data_ver)_peptides.RData"))
peptides_df = peptides_rdata["peptides.df"]
proteins_df = peptides_rdata["all_proteins.df"]
peptide2acs = Dict(r.peptide_id => (Set{String}(split(r.protein_acs, ';')), r.peptide_rank)
                  for r in eachrow(peptides_df))

protregroups_acs = ProtgroupAssembly.assemble_protgroups(peptide2acs, verbose=true,
                                                         nspec_peptides=2, rspec_peptides=0.25)
protregroups_acs_df = ProtgroupAssembly.dataframe(protregroups_acs,
                            protein_ranks=Dict(r.protein_ac => ProtgroupXMatch.rank_uniprot(r) for r in eachrow(proteins_df)),
                            protgroup_col=:protregroup_id, protein_col=:protein_ac)

CSV.write(joinpath(data_path, proj_info.msfolder,
                   "$(proj_info.id)_$(proj_info.msfolder)_$(proj_info.data_ver)_protregroups_acs.txt"),
          protregroups_acs_df, delim='\t', missingstring="")
