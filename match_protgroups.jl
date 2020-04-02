proj_info = (id = "cov2",
             data_ver = "20200331")
using Pkg
Pkg.activate(@__DIR__)
using Revise
using DataFrames, CSV, Statistics, RData

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const networks_path = joinpath(analysis_path, "networks")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(base_scripts_path, "misc_jl", "protgroup_crossmatch.jl"))

rdataset_names = [
    "APMS" => "$(proj_info.id)_msglm_data_mq_apms_20200329_20200329",
    "OeProteome" => "$(proj_info.id)_msdata_full_spectronaut_qc_20200331_20200331"
]
protgroups_dfs = [k => begin
        rdata = load(joinpath(scratch_path, "$v.RData"), convert=true)
        if haskey(rdata, "msdata") && haskey(rdata["msdata"], "protgroups")
            rdata["msdata"]["protgroups"]
        elseif haskey(rdata, "msdata_full") && haskey(rdata["msdata_full"], "protgroups")
            rdata["msdata_full"]["protgroups"]
        else
            @warn "No protgroups dataframe found for $k"
            nothing
        end
    end for (k, v) in rdataset_names]

proteins_df = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_spectronaut_qc_20200331_20200331.RData"), convert=true)["msdata_full"]["proteins"]
for noiso_group_df in groupby(proteins_df, :protein_ac_noiso)
    pe = noiso_group_df.protein_existence
    if any(ismissing, pe) && any(!ismissing, pe)
        noiso_group_df.protein_existence .= coalesce.(pe, minimum(skipmissing(pe)))
    elseif any(ismissing, pe) && any(noiso_group_df.is_viral)
        noiso_group_df.protein_existence .= 1
    end
end
# calculate AC ranks - the more canonical proptein is, the better
proteins_df[!, :ac_rank] .= [
    ifelse(ismissing(r.gene_name), 10, 0) +
    ifelse(r.protein_ac == r.protein_ac_noiso, 0, 1) +
    ifelse(r.src_db == "sp", 0, 3) +
    coalesce(r.protein_existence, 7) - 1
    for r in eachrow(proteins_df)]

# dataset names ordered by priority
ds_names = ["SKN_phospho", "ZIKV_apms", "hNPC_proteome", "candidates_apms"]
pg_dfs = [ds_name => protgroup_acs_rdata["all_protgroups"][ds_name] for ds_name in ds_names];
protein_ac_ranks = Dict(r.protein_ac => r.ac_rank for r in eachrow(proteins_df))

# save the matches
pg_matches_df = ProtgroupXMatch.match_protgroups(protgroups_dfs, protein_ac_ranks);
pg_matches_df[coalesce.(pg_matches_df.gene_names, "") .== "DPYSL2",
              [#:gene_names, :best_match_rank, :worst_match_rank,
               Symbol("row_ix.SKN_phospho"), #Symbol("row_ix.ZIKV_apms"),
               Symbol("row_ix.hNPC_proteome"), Symbol("row_ix.candidates_apms")]]
CSV.write(joinpath(analysis_path, "$(project_id)_all_protgroups_matches_$(analysis_version).txt"),
          pg_matches_df, delim='\t', null="NULL");
#output_filepath = joinpath(analysis_path, "$(project_id)_all_protgroups_matches_$(analysis_version).RData");
#@rput(output_filepath, output_filepath)
#@rput(combined_pgs_df, combined_pgs_df)
#R"save(combined_pgs_df, file=output_filepath)";
