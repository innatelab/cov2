proj_info = (id = "cov2",
             data_ver = "20200428",
             fit_ver = "20200428",
             oesc_ver = "20200512",
             modelobj = "ptmgroup",
             msfolder = "cov2timecourse_phospho_dia_20200423")

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, CSV, DataFrames, FastaIO
using JLD2
using StatsBase

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));

objid_col = Symbol(string(proj_info.modelobj, "_id"));

input_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(proj_info.msfolder)_$(proj_info.data_ver).RData"), convert=true)
full_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.msfolder)_$(proj_info.data_ver).RData"), convert=true)
fit_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(proj_info.msfolder)_$(proj_info.fit_ver).RData"), convert=true)
@load(joinpath(scratch_path, "$(proj_info.msfolder)_data_$(proj_info.data_ver).jld2"),
      ptmcollapsed2psitep_best_df)

effects_df = copy(input_rdata["effects.df"]);
contrasts_df = unique!(select!(copy(input_rdata["contrastXmetacondition.df"]), [:contrast, :contrast_type]));
objects_df = copy(input_rdata["msdata"][string(proj_info.modelobj, "s")]);
obj_contrasts_df = copy(fit_rdata["object_contrasts.df"]);
obj_effects_df = copy(fit_rdata["object_effects.df"]);
for df in [obj_contrasts_df, obj_effects_df]
    df |> MSGLMUtils.fix_quantile_columns!
    df.ptmgroup_shortid = replace.(df.ptmgroup_id, Ref(r"_M\d$" => ""))
end
# pick the most significant comparison/effect
obj_contrasts_df = by(obj_contrasts_df, [:contrast, :std_type, :ptmgroup_shortid]) do df
    selix= findmin(df.p_value)[2]
    df[selix:selix, :]
end
obj_contrasts_df.object_id = copy(obj_contrasts_df.ptmgroup_shortid)
obj_effects_df = by(obj_effects_df, [:effect, :std_type, :ptmgroup_shortid]) do df
    selix= findmin(df.p_value)[2]
    df[selix:selix, :]
end
obj_effects_df.object_id = copy(obj_effects_df.ptmgroup_shortid)

Revise.includet(joinpath(misc_scripts_path, "optcover_utils.jl"));

using BioSequences, BioAlignments # without it the loading will fail
@load(joinpath(scratch_path, "phosphositeplus_annotations_20200502.jld2"), psitep_annot_dfs)
for df in values(psitep_annot_dfs)
    df.ptm_id = df.protein_ac .* "_" .* uppercase.(df.ptm_AA) .* string.(df.ptm_pos)
end
psitep_annots_df = let kinsub_df = copy(psitep_annot_dfs[:KinaseSubstrate])
    filter!(r -> r.organism=="human", dropmissing!(kinsub_df, :kinase_gene_name))
    kinsub_df[!, :coll_id] .= :KinaseSubstrate
    kinsub_df.term_id = uppercase.(kinsub_df.kinase_gene_name)
    kinsub_df
end

Revise.includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));

perseus_report_df = CSV.read(joinpath(data_path, proj_info.msfolder,
        "COV2_DIA_phospho_0.75probablity_no normalization_psitep_annotated_Perseus_output.txt"),
        header=1, datarow=3, comment="#", delim='\t')
perseus_report_df.ptmgroup_shortid = replace.(perseus_report_df.PTM_collapse_key, Ref(r"_M\d$" => ""))
perseus_annots_collapsed_df = unique!(select(perseus_report_df,
        [:ptmgroup_shortid, :Motifs, :psitep_protein_ac, :psitep_ptm_AA, :psitep_ptm_pos]))
perseus_annots_collapsed_df.ptm_id = [!ismissing(r.psitep_ptm_AA) && !ismissing(r.psitep_protein_ac) ?
    r.psitep_protein_ac * "_" * uppercase(r.psitep_ptm_AA) * string(r.psitep_ptm_pos) : missing
    for r in eachrow(perseus_annots_collapsed_df)]

perseus_obj_annots_df = dropmissing!(DelimDataUtils.expand_delim_column(perseus_annots_collapsed_df,
                        list_col=:Motifs, elem_col=:term_id, key_col=:ptmgroup_shortid),
                        [:term_id, :ptmgroup_shortid])
perseus_ptm_annots_df = dropmissing!(DelimDataUtils.expand_delim_column(perseus_annots_collapsed_df,
                        list_col=:Motifs, elem_col=:term_id, key_col=:ptm_id),
                        [:term_id, :ptm_id])
for df in [perseus_ptm_annots_df, perseus_obj_annots_df]
    df[!, :coll_id] .= :seq_motifs
    df.term_id = replace.(df.term_id, Ref(r"\s+(?:motif|sequence)$" => ""))
end
perseus_colls = FrameUtils.frame2collections(perseus_obj_annots_df, set_col=:term_id, obj_col=:ptmgroup_shortid, coll_col=:coll_id)

#CSV.write(joinpath(data_path, proj_info.msfolder, "Perseus_motifs.txt"),
#          DataFrame(motif = collect(keys(perseus_colls[:seq_motifs]))),
#          delim='\t')
motifXregulators_collapsed_df = CSV.read(joinpath(data_path, proj_info.msfolder, "Perseus_motifsXregulators.txt"))
motifXregulators_df = DelimDataUtils.expand_delim_column(motifXregulators_collapsed_df,
                            list_col=:regulators, elem_col=:regulator, key_col=:motif, delim=",")

function motif2regulators(df::DataFrame)
    dfnew = rename!(select!(join(df, motifXregulators_df,
                         on=[:term_id=>:motif], kind=:left), Not(:term_id)),
            :regulator => :term_id) |> unique!
    dropmissing!(dfnew, [:term_id, hasproperty(dfnew, :ptm_id) ? :ptm_id : :ptmgroup_shortid])
    dfnew[!, :coll_id] .= :regulators
    return dfnew
end
perseus_ptm_regulators_df = motif2regulators(perseus_ptm_annots_df)
perseus_obj_regulators_df = motif2regulators(perseus_obj_annots_df)

# https://linkphinder.insight-centre.org/
linkphinder_df = CSV.read(joinpath(party3rd_data_path, "linkPhinder_data.csv"))
countmap(linkphinder_df.KinaseLabel)
linkphinder_annots_df = filter(r -> r.Score >= 0.75, rename(linkphinder_df, :KinaseLabel => :term_id, :ProteinSubstrate_ID => :protein_ac))
linkphinder_annots_df.ptm_pos = parse.(Int, [r.Site[2:end] for r in eachrow(linkphinder_annots_df)])
linkphinder_annots_df.ptm_AA = getindex.(linkphinder_annots_df.Site, Ref(1:1))
linkphinder_annots_df.ptm_id = [
    r.protein_ac * "_" * uppercase(r.ptm_AA) * string(r.ptm_pos)
    for r in eachrow(linkphinder_annots_df)]
linkphinder_annots_df[!, :coll_id] .= "LinkPhinder"

ptm_annots_df = vcat(
        select(psitep_annots_df, [:coll_id, :term_id, :ptm_id]),
        select(linkphinder_annots_df, [:coll_id, :term_id, :ptm_id]),
        perseus_ptm_annots_df,
        perseus_ptm_regulators_df)
ptm_colls = FrameUtils.frame2collections(ptm_annots_df,
    set_col=:term_id, obj_col=:ptm_id, coll_col=:coll_id)

ptmcollapsed2psitep_best_df.ptmgroup_shortid = replace.(ptmcollapsed2psitep_best_df.PTM_collapse_key, Ref(r"_M\d$" => ""))
ptmcollapsed_short2psitep_best_df = unique!(select(ptmcollapsed2psitep_best_df, Not([:PTM_collapse_key, :pepmodstate_seq])))
ptmcollapsed_short2psitep_best_df.ptm_id = [
    !ismissing(r.psitep_ptm_AA) ? r.protein_ac * "_" * uppercase(r.psitep_ptm_AA) * string(r.psitep_ptm_pos) : missing
    for r in eachrow(ptmcollapsed_short2psitep_best_df)]
obj_annots_df = vcat(
    select(join(psitep_annots_df, select(ptmcollapsed_short2psitep_best_df, Not(:gene_name)),
         on=[:protein_ac, :ptm_pos=>:psitep_ptm_pos, :ptm_id]),
        [:coll_id, :term_id, :ptmgroup_shortid]),
    select(join(linkphinder_annots_df, select(ptmcollapsed_short2psitep_best_df, Not(:gene_name)),
        on=[:protein_ac, :ptm_pos=>:psitep_ptm_pos, :ptm_id]),
        [:coll_id, :term_id, :ptmgroup_shortid]),
    perseus_obj_annots_df,
    perseus_obj_regulators_df)

obj_colls = FrameUtils.frame2collections(obj_annots_df,
    set_col=:term_id, obj_col=:ptmgroup_shortid, coll_col=:coll_id)

terms_df = unique!(select(obj_annots_df, [:coll_id, :term_id]));
terms_df.term_name = copy(terms_df.term_id)
terms_df[!, :term_descr] .= missing

@info "Preparing sets of hits"
ObjectType = eltype(obj_annots_df.ptmgroup_shortid)
obj_hit_sets = Dict{Tuple{String, String, String}, Set{ObjectType}}()
for hit_df in groupby(obj_contrasts_df[coalesce.(obj_contrasts_df.is_hit_nomschecks, false), :], [:std_type, :contrast, :change])
    obj_hit_sets[(string.(hit_df[1, :std_type], "_std"), hit_df[1, :contrast], hit_df[1, :change])] = Set(skipmissing(hit_df.object_id))
end
# only relevant ones
sel_std_type = "median"
obj_hit_selsets = filter(kv -> (kv[1][1] == sel_std_type*"_std") &&
                         occursin(r"SARS.+_vs_mock", kv[1][2]), obj_hit_sets);

@info "Preparing mosaics..."
observed_ptms = Set(skipmissing(join(ptmcollapsed_short2psitep_best_df, ptm_annots_df, on=:ptm_id, kind=:semi).ptm_id)) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, ptm_colls, observed_ptms,
                                    setXset_frac_extra_elms=0.05,
                                    verbose=true);

obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets,
                                          max_sets=2000, min_nmasked=2, max_setsize=2000)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

obj_hit_covers = Dict(begin
    @info "Covering $mosaic_name by hits..."
    mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.02], MaxSteps=3_000_000, WeightDigits=2, NWorkers=Threads.nthreads()-1, MaxRestarts=200),
            true)
    end for (mosaic_name, masked_mosaic) in pairs(obj_hit_mosaics))

using JLD2

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_$(proj_info.oesc_ver).jld2")
@save(hit_covers_filename,
      proj_info, ptm_colls, obj_colls, obj_mosaics,
      terms_df,
      objects_df, obj_effects_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
if !@isdefined(obj_effect_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(effect_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      objects_df, obj_effects_df, obj_effects_weak_df,
      obj_effect_sets, obj_effect_selsets, obj_effect_mosaics, obj_effect_selmosaics, #obj_effect_weak_sets,
      cover_params, obj_effect_covers, obj_effect_selcovers)
end

Revise.includet(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup↦gene_name map..."
obj_id2name = Dict(r.ptmgroup_shortid => r.gene_name * "_" * uppercase(r.psitep_ptm_AA) * string(r.psitep_ptm_pos)
                   for r in eachrow(filter(r -> !ismissing(r.psitep_ptm_AA), ptmcollapsed_short2psitep_best_df)))

obj_hit_covers_df = join(OptCoverUtils.covers_report(
    obj_hit_covers, obj_hit_selsets, obj_colls, obj_mosaics, obj_id2name,
    terms_df,
    maskid_col=[:std_type, :contrast, :change],
    maskedset_col_prefix="contrast"),
    contrasts_df, on=[:contrast], kind=:inner)

obj_hit_covers_signif_df = by(obj_hit_covers_df, :term_collection) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:std_type, :contrast, :change],
                                                   max_term_pvalue=1E-3, max_set_pvalue=1E-2, max_entry_pvalue=1.0),
                   Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_df[obj_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_signif_df[obj_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_plots.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

using PlotlyJS, TextWrap, ORCA

heatmap_layout_attrs = Dict(
    ("GO_CC", true) => Dict(:margin_l => 200),
    ("GO_CC", false) => Dict(:margin_l => 200),
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

stylize_contrast(str) = foldl(replace, [
    r"(SARS_COV2+)@(\d+)h" => s"<span style=\"font-weight: bold\">\2</span>h: SARS-CoV-2",
    r"mock@(\d)h" => "mock",
    "_vs_" => "&nbsp;<span style=\"color: #808080;\">vs</span>&nbsp;",
    ],
    init = str)

function process_contrast_axis(contrast_df)
    contrast_df,
    stylize_contrast.(contrast_df.contrast) .* " " .* OptCoverHeatmap.stylize_change.(contrast_df.change),
    stylize_contrast.(contrast_df.contrast) .* " " .* OptCoverHeatmap.stylize_change.(contrast_df.change)#stylize_effect.(effect_df.effect)
end

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")hit heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = signif ? obj_hit_covers_signif_df : obj_hit_covers_df
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            Symbol(term_coll), elements_label="PTM",
            maskedset_axis_title = "contrast",
            maskedset_cols = [:contrast, :change, :ncontrast],
            process_maskedset_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 160),
            cell_width=25, cell_height=25,
            row_order=contrasts -> begin
                contrast_matches = match.(Ref(r"SARS_COV2@(\d+)h"), contrasts.contrast)
                contrasts.timepoint = parse.(Int, getindex.(contrast_matches, 1))
                return sortperm(contrasts, [:change, :timepoint])
            end,
            transpose=false)
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    plotname = joinpath(plots_path, "$(proj_info.msfolder)", "oesc_hits_$(sel_std_type)",
                        "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_X_contrast$(signif ? "_signif" : "")_heatmap")
    PlotlyJS.savehtml(coll_heatmap, "$(plotname).html", :embed);
    try
        savefig(coll_heatmap.plot, "$(plotname).pdf");
    catch e
        @warn e
    end
end

obj_contrasts_report_df = join(join(
    select(filter(r -> occursin(r"SARS.+_vs_mock", r.contrast), obj_contrasts_df),
           [:contrast, :std_type, :ptmgroup_id, :ptmgroup_shortid, :median_log2, :p_value, :is_signif, :is_hit_nomschecks, :is_hit]),
    select(ptmcollapsed2psitep_best_df, [:ptmgroup_shortid, :gene_name, :protein_ac, :data_ptm_pos, :psitep_ptm_AA, :psitep_ptm_pos, :flanking_15AAs]),
    on=[:ptmgroup_shortid]
), select(psitep_annots_df, [:protein_ac, :ptm_pos, :domain, :kinase_gene_names,
                          :reg_function, :reg_prot_iactions, :reg_other_iactions,
                          :diseases, :diseases_pubmed_ids]),
  on=[:protein_ac, :psitep_ptm_pos => :ptm_pos], kind=:left
)
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_contrasts_$(proj_info.fit_ver).txt"),
          obj_contrasts_report_df, missingstring="", delim='\t')
