proj_info = (id = "cov2",
             oesc_ver = "20200606")

datasets = Dict(
    :cov2ts_proteome => (msfolder="cov2timecourse_dia_20200423",
                         fit_ver="20200429",
                         analysis=:msglm, ptm=nothing,
                         datatype=:proteome,
                         label="Proteome (DIA)",
                         color="#226430"),
    :cov2ts_phospho => (msfolder="cov2timecourse_phospho_dia_20200423",
                        fit_ver="20200428",
                        analysis=:msglm, ptm=:phospho,
                        datatype=:phospho,
                        label="Phosphoproteome (DIA)",
                        color="#5e268f"),
    :cov2ts_rnaseq => (msfolder="cov2ts_tx",
                       fit_ver="20200518",
                       analysis=:deseq, ptm=:nothing,
                       datatype=:rnaseq,
                       label="RNA-Seq",
                       color="#383c9b"),
    :cov2el_proteome => (msfolder="cov2earlylate_fp_phos_ubi_dda_20200429",
                         fit_ver="20200514",
                         analysis=:msglm, ptm=nothing,
                         datatype=:proteome,
                         label="Proteome (DDA)",
                         color="#226430"),
    :cov2el_phospho => (msfolder="cov2earlylate_fp_phos_ubi_dda_20200601",
                        fit_ver="20200510",
                        analysis=:perseus, ptm=:phospho, prefix=:pho,
                        datatype=:phospho,
                        label="Phosphoproteome (DDA)",
                        filename="output_Phospho (STY)Sites",
                        color="#5e268f"),
    :cov2el_ubi => (msfolder="cov2earlylate_fp_phos_ubi_dda_20200601",
                    fit_ver="20200510",
                    analysis=:perseus, ptm=:ubiquitin, prefix=:ubi,
                    datatype=:ubi,
                    label="Ubiquitinome (DDA)",
                    filename="output_S0-01_GlyGly (K)Sites",
                    color="#bf1c2c"),
)

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, CSV, DataFrames, FastaIO
using JLD2
using StatsBase

@info "Project '$(proj_info.id)' analysis version=$(proj_info.oesc_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "mq_utils.jl"));

objid_col = :protgroup_id_united;

msglm_rdata = Dict(begin
    @info "Loading $ds analysis (fit_ver=$(dsinfo.fit_ver))..."
    ds => (input_data = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(dsinfo.msfolder)_$(dsinfo.fit_ver).RData"), convert=true),
           full_data = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(dsinfo.msfolder)_$(dsinfo.fit_ver).RData"), convert=true),
           fit = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(dsinfo.msfolder)_$(dsinfo.fit_ver).RData"), convert=true))
end for (ds, dsinfo) in pairs(datasets) if dsinfo.analysis == :msglm)

rnaseq_rdata = load(joinpath(scratch_path, "$(datasets[:cov2ts_rnaseq].msfolder)_pairwise_comp_uniprot_mapped_$(datasets[:cov2ts_rnaseq].fit_ver).RData"), convert=true)
rnaseq_contrasts_df = rnaseq_rdata["all_pairwise_comp.df"]
# recode ComparisonName into the same format as contrast
comparison_matches = match.(Ref(r"(?:(?<tp_lhs>\d+)h|(?<trt_lhs>SARSCoV2))_vs_(?:(?<tp_rhs>\d+)h|(?<trt_rhs>MOCK))within(?:(?<tp>\d+)h|(?<trt>SARSCoV2|MOCK))"),
                            rnaseq_contrasts_df.ComparisonName)
getmatch(regmatch, m, mfallback) = !isnothing(regmatch[m]) ? regmatch[m] : regmatch[mfallback]
rnaseq_contrasts_df.ComparisonName[isnothing.(comparison_matches)]
rnaseq_contrasts_df.dataset = :cov2ts_rnaseq
rnaseq_contrasts_df.timepoint_lhs = [parse(Int, getmatch(m, :tp_lhs, :tp)) for m in comparison_matches]
rnaseq_contrasts_df.timepoint_rhs = [parse(Int, getmatch(m, :tp_rhs, :tp)) for m in comparison_matches]
rnaseq_contrasts_df.treatment_lhs = [replace(replace(getmatch(m, :trt_lhs, :trt), "SARSCoV2" => "SARS_CoV2"), "MOCK" => "mock")
                                     for m in comparison_matches]
rnaseq_contrasts_df.treatment_rhs = [replace(replace(getmatch(m, :trt_rhs, :trt), "SARSCoV2" => "SARS_CoV2"), "MOCK" => "mock")
                                     for m in comparison_matches]
rnaseq_contrasts_df.contrast = [r.treatment_lhs * "@" * string(r.timepoint_lhs) * "h_vs_" *
                                r.treatment_rhs * "@" * string(r.timepoint_rhs) * "h"
                                for r in eachrow(rnaseq_contrasts_df)]
rnaseq_protgroups_df = combine(groupby(
    filter(r -> !ismissing(r.protein_ac), unique!(select(rnaseq_contrasts_df, [:GeneName, :protein_ac]))),
    :GeneName)) do df
    res = DataFrame(protein_acs = join(df.protein_ac, ";"),
                    gene_names = df.GeneName[1:1])
    res.majority_protein_acs = copy(res.protein_acs)
    res[!, :is_contaminant] .= false
    res[!, :is_reverse] .= false
    res[!, :organism] .= false
    return res
end
rnaseq_protgroups_df.protgroup_id = 1:nrow(rnaseq_protgroups_df)


using BioAlignments, BioSequences
@load(joinpath(scratch_path, "$(datasets[:cov2ts_phospho].msfolder)_data_$(datasets[:cov2ts_phospho].fit_ver).jld2"),
      ptmcollapsed2psitep_best_df)
ptmcollapsed2psitep_best_df

perseus_analysis_orig_df = vcat([begin
    @info "Loading $ds analysis..."
    wide_df = CSV.read(joinpath(data_path, dsinfo.msfolder,
            "curban_analysis", "$(dsinfo.filename).txt"),
            header=1, datarow=3, comment="#", delim='\t')
    long_df = FrameUtils.pivot_longer(wide_df,
        Symbol.(["Positions within proteins", "Leading proteins",
         "Protein", "Protein names", "Gene names"]),
        measure_vars_regex=Regex("^(?<value>[^.]+) (?<var>$(dsinfo.prefix)_.+)\$"),
        var_col=:comparison)
    comparison_matches = match.(Ref(Regex("$(dsinfo.prefix)_SARS_COV2_(\\d+)h_$(dsinfo.prefix)_mock_(?:\\d+)h")),
                                string.(long_df.comparison))
    long_df[!, :contrast_lhs] .= "SARS_COV2"
    long_df[!, :contrast_rhs] .= "mock"
    long_df[!, :ptm_type] .= dsinfo.ptm
    long_df[!, :dataset] .= ds
    long_df.timepoint = parse.(Int, getindex.(comparison_matches, 1))
    long_df.contrast = ["$(r.contrast_lhs)@$(r.timepoint)h_vs_$(r.contrast_rhs)@$(r.timepoint)h" for r in eachrow(long_df)]
    rename!(long_df,
        "Positions within proteins" => :ptm_poses,
        "Leading proteins" => :lead_protein_acs,
        "Protein" => :protein_acs,
        "Protein names" => :protein_descriptions,
        "Gene names" => :gene_names,
        "Student's T-test Significant" => :is_signif,
        "-Log Student's T-test p-value" => :mlog10_pvalue,
        "Student's T-test Test statistic" => :ttest_stat,
        "Student's T-test q-value" => :qvalue,
        "Student's T-test Difference" => :mean)
    long_df.is_signif = coalesce.(long_df.is_signif, ".") .== "+"
    long_df.change = ifelse.(long_df.is_signif,
        ifelse.(coalesce(long_df.mean .> 0), "+", "-"),
        ".")
    long_df
end for (ds, dsinfo) in datasets if dsinfo.analysis == :perseus]...)
countmap(perseus_analysis_orig_df.contrast)
cov2el_new_protgroups_df = rename!(select(rename!(x -> replace(x, "fp_" => ""), CSV.read(joinpath(data_path, datasets[:cov2el_ubi].msfolder,
        "curban_analysis", "output_proteinGroups.txt"),
        header=1, datarow=3, comment="#", delim='\t')),
        ["Protein IDs", "Majority protein IDs", "Protein names", "Gene names", "id"]),
        ["Protein IDs" => :protein_acs, "Majority protein IDs" => :majority_protein_acs,
         "Protein names" => :protein_names, "Gene names" => :gene_names,
         "id" => :protgroup_id])

#=
using VegaLite
p = perseus_tests_df |>
@vlplot(
    :bar,
    transform=[
        {filter="datum.is_signif"},
    ],
    column="timepoint:o",
    y={"sum(is_signif)", axis={title="N regulated", grid=false}},
    x={"change:n", axis={title=""}},
    color={"change:n", scale={range=["firebrick", "dodgerblue"]}},
    spacing=10,
    config={
        view={stroke=:transparent},
        axis={domainWidth=1, labelAngle=0}
    }
)
save(joinpath(plots_path, "test.pdf"), p)

using RCall
R"require(ggplot2)"
R"require(dplyr)"
R"source('~/projects/R/misc/ggplot_ext.R')"
R"p = ggplot(filter($perseus_tests_df, is_signif) %>% group_by(timepoint, change) %>% summarise(n_regulated = n()),
         aes(x = timepoint, y=n_regulated)) +
geom_bar(aes(color=change, fill=change), stat='identity', show.legend=FALSE, position=position_dodge(width=15)) +
geom_text(aes(label=n_regulated, color=change), size=6, vjust=-0.8, show.legend=FALSE, position=position_dodge(width=15)) +
scale_color_manual(values=c('▲' = 'firebrick', '▼' = 'dodgerblue')) +
scale_fill_manual(values=c('▲' = 'firebrick', '▼' = 'dodgerblue')) +
scale_x_continuous(breaks=unique($(perseus_tests_df.timepoint))) +
scale_y_continuous(limits = c(0, 1750)) +
theme_bw_ast()"
R"ggsave(file.path($plots_path, 'cov2_DDA_$(prefix)_stats.pdf'), p, width=6, height=4, device=cairo_pdf)"

R"require(stringr)"
R"require(tidyr)"
R"p = ggplot(filter($obj_contrasts_df, treatment_lhs != treatment_rhs & timepoint_lhs == timepoint_rhs &
                std_type=='median' & is_hit_nomschecks) %>%
             rename(timepoint = timepoint_lhs, change_old=change, change=change_alt) %>%
             group_by(contrast, timepoint, change) %>% summarise(n_regulated = n()),
         aes(x = timepoint, y=n_regulated)) +
geom_bar(aes(color=change, fill=change), stat='identity', show.legend=FALSE, position=position_dodge(width=3)) +
geom_text(aes(label=n_regulated, color=change), size=6, vjust=-0.8, show.legend=FALSE, position=position_dodge(width=3)) +
scale_color_manual(values=c('▲' = 'firebrick', '▼' = 'dodgerblue')) +
scale_fill_manual(values=c('▲' = 'firebrick', '▼' = 'dodgerblue')) +
scale_x_continuous(breaks=unique($(obj_contrasts_df.timepoint_lhs))) +
scale_y_continuous(limits = c(0, 500)) +
theme_bw_ast()"
R"ggsave(file.path($plots_path, 'cov2ts_DIA_phospho_stats.pdf'), p, width=6, height=4, device=cairo_pdf)"
=#

Revise.includet(joinpath(base_scripts_path, "misc_jl", "fasta_reader.jl"))
Revise.includet(joinpath(base_scripts_path, "misc_jl", "protgroup_crossmatch.jl"))

protgroups_cov2ts_phospho_df =
    rename!(unique(select(msglm_rdata[:cov2ts_phospho].input_data["msdata"]["ptmgroups"],
                          [:majority_protein_acs, :PG_ProteinGroups, :gene_names,
                           :PG_Organisms, :protein_names, :is_reverse, :is_contaminant, :is_viral])),
            :PG_Organisms => :organism, :PG_ProteinGroups => :protein_acs)
protgroups_cov2ts_phospho_df[!, :protgroup_id] = 1:nrow(protgroups_cov2ts_phospho_df)

protgroups_dfs = Dict{Symbol, DataFrame}()
for (ds, rdata) in msglm_rdata
    msdata = rdata.input_data["msdata"]
    if haskey(msdata, "protregroups")
        pg_df = rename(msdata["protregroups"], :protregroup_id => :protgroup_id)
    elseif haskey(msdata, "protgroups")
        pg_df = msdata["protgroups"]
    else
        @warn "No protgroups dataframe found for $ds"
        pg_df = nothing
    end
    if pg_df !== nothing
        push!(protgroups_dfs, ds => pg_df)
    end
end
push!(protgroups_dfs, :cov2el_phospho => cov2el_new_protgroups_df)
push!(protgroups_dfs, :cov2ts_phospho => protgroups_cov2ts_phospho_df)
push!(protgroups_dfs, :cov2ts_rnaseq => rnaseq_protgroups_df)
# FIXME? not importing the proteins of cov2ts_phospho

proteins_df = Fasta.read_uniprot(joinpath(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta"))
proteins_df.protein_ac_noiso = Fasta.strip_uniprot_isoform.(proteins_df.protein_ac)
rename!(proteins_df, :genename => :gene_name)

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
pg_matches_df = ProtgroupXMatch.match_protgroups(collect(pairs(protgroups_dfs)), protein_ac_ranks);
pg_matches_df[!, :protgroup_id_united] = 1:nrow(pg_matches_df)
# add protgroup ids of each dataset and add protgroup_id_common to each dataset
for (dsname, pg_df) in pairs(protgroups_dfs)
    rowix_col = Symbol("rowix_", dsname)
    if !hasproperty(pg_matches_df, rowix_col)
        @warn "Dataset $dsname key does not exist, skipping"
        continue
    end
    pg_col = Symbol("protgroup_id_", dsname)
    pg_matches_df[!, pg_col] = missings(nonmissingtype(eltype(pg_df.protgroup_id)), nrow(pg_matches_df))
    pg_df[!, :protgroup_id_united] = missings(Int, nrow(pg_df))
    for (i, rowix) in enumerate(pg_matches_df[!, rowix_col])
        if !ismissing(rowix)
            pg_matches_df[i, pg_col] = pg_df[rowix, :protgroup_id]
            pg_df[rowix, :protgroup_id_united] = pg_matches_df[i, :protgroup_id_united]
        end
    end
    select!(pg_matches_df, Not(rowix_col))
end
pg_matches_long_df = FrameUtils.pivot_longer(pg_matches_df, [:protgroup_id_united, :protein_acs, :majority_protein_acs, :gene_names,
                                             #=:is_contaminant, :is_reverse=#],
                        measure_vars_regex=r"^(?<value>rowix|protgroup_id|pgrank|acrank)_(?<var>[^u].+)$",
                        var_col=:dataset)
pg_matches_long_df.dataset = Symbol.(pg_matches_long_df.dataset)
pg_matches_long_expanded_df = DelimDataUtils.expand_delim_column(
    pg_matches_long_df, key_col=[:dataset, :protgroup_id_united, :protgroup_id],
    list_col=:protein_acs, elem_col=:protein_ac)

obj2protac_df = unique!(dropmissing!(select(pg_matches_long_expanded_df,
                                            [objid_col, :protein_ac])))

nonunique_matches_dfs = [dsname => begin
    pg_col = Symbol("protgroup_id_", dsname)
    combine(groupby(filter(r -> !ismissing(r[pg_col]), pg_matches_df), pg_col)) do df
        return nrow(df) > 1 ? df : df[1:0, :]
    end
end for (dsname, _) in pairs(protgroups_dfs)]

obj_contrasts_df = vcat([begin
    @info "Process $ds results"
    if ds != :cov2ts_phospho
        obj_contrasts_df = rdata.fit["object_contrasts.df"]
        obj_contrasts_df.protgroup_id = obj_contrasts_df.object_id
        obj_contrasts_df = join(obj_contrasts_df,
             select(filter(r -> r.dataset == ds, pg_matches_long_df), [:protgroup_id, :protgroup_id_united]),
             on=:protgroup_id, kind=:left)
    else
        obj_contrasts_df = join(rdata.fit["object_contrasts.df"],
                                select(ptmcollapsed2psitep_best_df,
                                    [:PTM_collapse_key, :PTM_collapse_protein_ac, :PTM_collapse_ref_pos]),
                                on=(:ptmgroup_id => :PTM_collapse_key))
        rename!(obj_contrasts_df, :PTM_collapse_protein_ac => :protein_ac)
        obj_contrasts_df = join(obj_contrasts_df,
             filter(r -> r.dataset == ds, pg_matches_long_expanded_df),
             on=:protein_ac, kind=:left)
    end
    obj_contrasts_df[!, :dataset] .= ds
    obj_contrasts_df[!, :ptm_type] .= datasets[ds].ptm !== nothing ? datasets[ds].ptm : missing
    select(obj_contrasts_df,
        [:dataset, :ptm_type, :contrast, :std_type,
         :object_id, :protgroup_id, :protgroup_id_united,
         :gene_names, :median_log2, :p_value,
         :is_hit, :change, :is_hit_nomschecks, :is_signif])
end for (ds, rdata) in pairs(msglm_rdata)]...)

obj_contrasts_pg_agg_df = combine(groupby(filter(r -> !ismissing(r.protgroup_id_united), obj_contrasts_df),
   [:dataset, :ptm_type, :contrast, :std_type, :protgroup_id_united])) do compXpg_df
   top_signif_ix = findmin(compXpg_df.p_value)[2]
   df = compXpg_df[top_signif_ix:top_signif_ix, :]
   df[!, :nptm_hits] .= ismissing(df.ptm_type[1]) ? 0 : sum(c -> c != ".", compXpg_df.change)
   return df
end

rnaseq_pg_agg_df = combine(groupby(innerjoin(rnaseq_contrasts_df, filter(r -> r.dataset == :cov2ts_rnaseq, pg_matches_long_expanded_df),
                        on=[:dataset, :protein_ac]),
                        [:dataset, :contrast, :protgroup_id_united])) do compXpg_df
   top_signif_ix = findmin(compXpg_df.pvalue)[2]
   df = compXpg_df[top_signif_ix:top_signif_ix, :]
   return df
end

perseus_protgroups_df = unique(select(perseus_analysis_orig_df, [:protein_acs]))
perseus_protgroups_df.ptmprotgroup_id = 1:nrow(perseus_protgroups_df)
perseus_analysis_df = join(perseus_analysis_orig_df, perseus_protgroups_df, on=:protein_acs)
perseus_ptmprotgroup2protein_df = DelimDataUtils.expand_delim_column(perseus_protgroups_df,
    key_col=:ptmprotgroup_id, list_col=:protein_acs, elem_col=:protein_ac)
perseus_ptmprotgroup2protgroup_df = select(
    join(perseus_ptmprotgroup2protein_df,
         filter(r -> r.dataset == :cov2el_phospho, pg_matches_long_expanded_df), on=:protein_ac),
    [:ptmprotgroup_id, :protgroup_id_united]) |> unique!
perseus_analysis_df = join(perseus_analysis_df,
                           perseus_ptmprotgroup2protgroup_df,
                           on=:ptmprotgroup_id, kind=:left)

perseus_pg_agg_df = combine(groupby(filter(r -> !ismissing(r.protgroup_id_united) && !ismissing(r.ptm_poses), perseus_analysis_df),
        [:dataset, :ptm_type, :contrast, :protgroup_id_united])) do compXpg_df
   top_signif_ix = findmax(compXpg_df.mlog10_pvalue)[2]
   df = compXpg_df[top_signif_ix:top_signif_ix, :]
   df[!, :nptm_hits] .= sum(c -> c != ".", compXpg_df.change)
   return df
end
sum(perseus_pg_agg_df.nptm_hits)

contrasts_df = vcat(unique!(select(obj_contrasts_df, [:dataset, :contrast, :change])),
                    unique!(select(perseus_analysis_df, [:dataset, :contrast, :change])),
                    unique!(select(rnaseq_contrasts_df, [:dataset, :contrast, :change])))
contrast_matches = match.(Ref(r"(.+)@(\d+)h_vs_(.+)@(\d+)h"), contrasts_df.contrast)
contrasts_df.dataset = string.(contrasts_df.dataset)
contrasts_df.treatment_lhs = string.(getindex.(contrast_matches, 1))
contrasts_df.timepoint_lhs = parse.(Int, getindex.(contrast_matches, 2))
contrasts_df.treatment_rhs = string.(getindex.(contrast_matches, 3))
contrasts_df.timepoint_rhs = parse.(Int, getindex.(contrast_matches, 4))
contrasts_df.change_alt = getindex.(Ref(Dict("+" => "▲", "-" => "▼", "." => ".")),
                                        contrasts_df.change)

Revise.includet(joinpath(misc_scripts_path, "optcover_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "gmt_reader.jl"));
Revise.includet(joinpath(misc_scripts_path, "omics_collections.jl"));

@info "Loading Human annotations..."
# human mappings from http://download.baderlab.org/EM_Genesets/December_01_2018/Human/UniProt/
genesets_df, genesets_coll = GMT.read(String,
        joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_April_01_2019_UniProt.gmt"),
        id_col = :term_id, src_col = :term_src);

pcomplexes_df, pcomplex_iactors_df, pcomplex_iactor2ac_df =
    OmicsCollections.ppicollection(joinpath(party3rd_data_path, "complexes_20191217.RData"), seqdb=:uniprot);
pcomplexes_df[!, :coll_id] .= "protein_complexes";

# make complexes collections, keep complexes with at least 2 participants
pcomplex_coll = FrameUtils.frame2collection(join(pcomplex_iactors_df, pcomplex_iactor2ac_df,
    on=[:file, :entry_index, :interaction_id, :interactor_id], kind=:inner),
            set_col=:complex_id, obj_col=:protein_ac, min_size=2)
protac_sets = merge(genesets_coll, pcomplex_coll)

terms_df = vcat(rename(genesets_df[!, [:term_src, :term_id, :name, :descr]],
                       :term_src => :coll_id, :name=>:term_name, :descr=>:term_descr),
                #rename(goterm_info_df[[:id, :name, :def]], Dict(:onto => :coll_id, :id=>:term_id, :name=>:term_name, :def=>:term_descr)),
                rename(pcomplexes_df[!, [:coll_id, :complex_id, :interaction_label, :interaction_name]],
                       :complex_id=>:term_id, :interaction_label=>:term_name, :interaction_name=>:term_descr));
protac2term_df = FrameUtils.collection2frame(protac_sets, terms_df,
                                             setid_col=:term_id, objid_col=:protein_ac)

obj2term_df = select!(join(obj2protac_df, protac2term_df, on = :protein_ac, kind=:inner),
                      Not([:protein_ac])) |> unique!
protac_colls = FrameUtils.frame2collections(protac2term_df, obj_col=:protein_ac,
                                            set_col=:term_id, coll_col=:coll_id)
obj_colls = FrameUtils.frame2collections(obj2term_df, obj_col=objid_col,
                                         set_col=:term_id, coll_col=:coll_id)

@info "Preparing hit sets"
ObjectType = eltype(obj2protac_df[!, objid_col])
obj_hit_sets = Dict{Tuple{String, String, String, String}, Set{ObjectType}}()
for hits_df in groupby(filter(r -> coalesce(r.is_hit_nomschecks, false), obj_contrasts_df),
                       [:dataset, :std_type, :contrast, :change])
    obj_hit_sets[(string(hits_df[1, :dataset]), string(hits_df[1, :std_type], "_std"),
                  hits_df[1, :contrast], hits_df[1, :change])] =
        Set(skipmissing(hits_df[!, objid_col]))
end
for hits_df in groupby(filter(r -> coalesce(r.is_signif, false), perseus_pg_agg_df),
                       [:dataset, :contrast, :change])
    obj_hit_sets[(string(hits_df[1, :dataset]), "median_std",
                  hits_df[1, :contrast], hits_df[1, :change])] =
        Set(skipmissing(hits_df[!, objid_col]))
end
for hits_df in groupby(filter(r -> coalesce(r.is_signif, false), rnaseq_pg_agg_df),
                       [:dataset, :contrast, :change])
    obj_hit_sets[(string(hits_df[1, :dataset]), "median_std",
                  hits_df[1, :contrast], hits_df[1, :change])] =
        Set(skipmissing(hits_df[!, objid_col]))
end

# only relevant ones
#sel_std_type = "replicate_std"
sel_std_type = "median_std"
obj_hit_selsets = filter(kv -> (kv[1][2] == sel_std_type) && (
    occursin(r"SARS.+_vs_mock", kv[1][3])) && kv[1][4] != ".",
  obj_hit_sets)

@info "Preparing mosaics..."
observed_protacs = Set(obj2protac_df.protein_ac) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls,
                                      protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);

# remove broad terms larger than 150 elements (too unspecific)
obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200, verbose=true)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

ENV["MKL_NUM_THREADS"] = 1

obj_hit_covers = Dict(begin
    @info "Covering $mosaic_name by hits..."
    mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=Threads.nthreads()-1, MaxRestarts=200),
            true)
    end for (mosaic_name, masked_mosaic) in pairs(obj_hit_mosaics))

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_united_hit_covers_$(proj_info.oesc_ver).jld2")
@save(hit_covers_filename,
      proj_info, datasets, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      #objects_df,
      protgroups_dfs, pg_matches_df,
      #obj_effects_df,
      obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
if !@isdefined(obj_hit_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, datasets, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      #objects_df,
      protgroups_dfs, pg_matches_df,
      #obj_effects_df,
      obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
end

include(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup↦gene_name map..."
obj_id2name = Dict(r.protgroup_id_united => !ismissing(r.gene_names) ? DelimDataUtils.rejoin_unique_substrings([r.gene_names], ";") :
                    !ismissing(r.majority_protein_acs) ? DelimDataUtils.rejoin_unique_substrings([r.majority_protein_acs], ";") :
                    string(r.protgroup_id_united)
                   for r in eachrow(filter(r -> !ismissing(r.protgroup_id_united), pg_matches_df)))

obj_hit_covers_df = innerjoin(OptCoverUtils.covers_report(
    obj_hit_covers, obj_hit_selsets, obj_colls, obj_mosaics, obj_id2name,
    terms_df,
    maskid_col=[:dataset, :std_type, :contrast, :change],
    maskedset_col_prefix="contrast"),
    contrasts_df, on=[:dataset, :contrast, :change])

# don't remove the sets since they are timecourses timepoints
obj_hit_covers_signif_df = combine(groupby(obj_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:dataset, :std_type, :contrast, :change],
                                                   max_term_pvalue=3E-4, max_set_pvalue=1.0, max_entry_pvalue=1.0),
                   Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_united_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_df[obj_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_united_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_signif_df[obj_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_plots.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

using PlotlyJS, TextWrap, ORCA

heatmap_layout_attrs = Dict(
    ("GO_CC", true) => Dict(:margin_l => 200),
    ("GO_CC", false) => Dict(:margin_l => 200),
    ("Reactome", true) => Dict(:margin_l => 220),
)

for (plot_mosaic, cover_coll) in obj_hit_covers
    isempty(cover_coll.results) && continue
    @info "Plotting $plot_mosaic Pareto front"
    paretofront_plot = OptCoverPlots.plot_paretofront(cover_coll.results[1], plot_unfolded=true)
    plot_filename = joinpath(plots_path, "oesc_$(sel_std_type)_std", "paretofront",
                             "$(proj_info.id)_$(plot_mosaic)_X_treatment_$(sel_std_type)_pareto")
    savefig(paretofront_plot.plot, "$plot_filename.svg")
    PlotlyJS.savehtml(paretofront_plot, "$plot_filename.html")
end

stylize_dataset(ds) =
    "<span style=\"font-weight: bold; color: $(datasets[Symbol(ds)].color);\">" * datasets[Symbol(ds)].label * "</span>"

stylize_hit(str) = foldl(replace, [
    "bait_id" => "<span style=\"color: #808080;\">bait:</span>&nbsp;",
    ":orgcode" => "&nbsp;CoV-2&nbsp;<span style=\"color: #808080;\">vs</span>&nbsp;",
    ],
    init = str)

stylize_contrast(str) = foldl(replace, [
    r"(SARS_C[oO]V2+)@(\d+)h_vs_mock@(\d+)h" => s"<span style=\"font-weight: bold; color: black;\">\2</span>h",
    ],
    init = str)

function process_contrast_axis(contrast_df)
    contrast_df,
    stylize_dataset.(contrast_df.dataset) .* ":&nbsp;" .*
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change),
    stylize_dataset.(contrast_df.dataset) .* ":&nbsp;" .*
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change)#stylize_effect.(effect_df.effect)
end

datatype_order = Dict(:rnaseq => 1, :proteome => 2, :ubi => 3, :phospho => 4)

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")hit heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = signif ? obj_hit_covers_signif_df : obj_hit_covers_df
    for outformat in ("html", "pdf", "svg")
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            Symbol(term_coll), elements_label="protein",
            maskedset_axis_title = "contrast",
            maskedset_cols = [:dataset, :contrast, :change, :ncontrast],
            process_maskedset_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 120),
            transpose=false,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",#outformat in ["svg", "pdf"] ? "#000" : "#BBB",
            cell_width=25, cell_height=20, gridwidth=1,
            row_order=contrasts -> begin
                contrast_matches = match.(Ref(r"SARS_C[oO]V2@(\d+)h"), contrasts.contrast)
                contrasts.timepoint = parse.(Int, getindex.(contrast_matches, 1))
                contrasts.datatype_order = [datatype_order[datasets[Symbol(r.dataset)].datatype]
                                            for r in eachrow(contrasts)]
                return sortperm(contrasts, [:datatype_order, :change, :timepoint, :dataset])
            end)
        (coll_heatmap === nothing) && continue
        for (k, v) in [#:width=>800, :height=>400,
                       :margin_r=>100, :margin_t=>20,
                       :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
            coll_heatmap.plot.layout[k] = v
        end
        plot_fname = joinpath(plots_path, "united_hits_oesc_$(sel_std_type)",
                              "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_contrast$(signif ? "_signif" : "")_heatmap.$(outformat)")
        if outformat == "html"
            PlotlyJS.savehtml(coll_heatmap, plot_fname, :embed);
        else
            try
                savefig(coll_heatmap.plot, plot_fname);
            catch e
                @warn e
            end
        end
    end
end
