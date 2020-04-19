proj_info = (id = "cov2",
             match_ver = "20200417")
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
protein_ac_ranks = Dict(r.protein_ac => ProtgroupXMatch.rank_uniprot(r) for r in eachrow(proteins_df))

# save the matches
pg_matches_df = ProtgroupXMatch.match_protgroups(protgroups_dfs, protein_ac_ranks);
for (dsname, pg_df) in protgroups_dfs
    pg_col = Symbol("protgroup_id_", dsname)
    rowix_col = Symbol("rowix_", dsname)
    pg_matches_df[!, pg_col] = missings(nonmissingtype(eltype(pg_df.protgroup_id)), nrow(pg_matches_df))
    for (i, rowix) in enumerate(pg_matches_df[!, rowix_col])
        if !ismissing(rowix)
            pg_matches_df[i, pg_col] = pg_df[rowix, :protgroup_id]
        end
    end
    select!(pg_matches_df, Not(rowix_col))
end

nonunique_matches_dfs = [dsname => begin
    pg_col = Symbol("protgroup_id_", dsname)
    by(filter(r -> !ismissing(r[pg_col]), pg_matches_df), pg_col) do df
        return nrow(df) > 1 ? df : df[1:0, :]
    end
end for (dsname, _) in protgroups_dfs]

nonunique_matches_dfs[1][2].protgroup_id_OeProteome
CSV.write(joinpath(analysis_path, "$(proj_info.id)_protgroup_matches_$(proj_info.match_ver).txt"),
          pg_matches_df, delim='\t', missingstring="");
#output_filepath = joinpath(analysis_path, "$(project_id)_all_protgroups_matches_$(analysis_version).RData");
#@rput(output_filepath, output_filepath)
#@rput(combined_pgs_df, combined_pgs_df)
#R"save(combined_pgs_df, file=output_filepath)";
