proj_info = (id = "cov2",
             data_ver = "20201012",
             fit_ver = "20201012",
             oesc_ver = "20201029",
             modelobj = "ptmn",
             countobj = "ptmgroup",
             msfolder = "snaut_parsars_phospho_20201005")

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, CSV, CodecZlib, DataFrames, FastaIO
using JLD2
using StatsBase

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

includet(joinpath(misc_scripts_path, "frame_utils.jl"));
includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));

objid_col = Symbol(string(proj_info.countobj, "_id"));

input_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(proj_info.msfolder)_$(proj_info.data_ver).RData"), convert=true)
full_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.msfolder)_$(proj_info.data_ver).RData"), convert=true)
fit_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(proj_info.msfolder)_$(proj_info.fit_ver).RData"), convert=true)

ptmns_df = open(GzipDecompressorStream, joinpath(data_path, proj_info.msfolder, "ptm_extractor_$(proj_info.fit_ver)", "ptmns_grouped.txt.gz"), "r") do io
    CSV.read(io, delim='\t')
end
ptm2group_df = select(ptmns_df, [:ptm_id, :ptmgroup_id]) |> unique!

effects_df = copy(input_rdata["effects.df"]);
contrasts_df = copy(input_rdata["contrasts.df"])
contrasts_df.timepoint_lhs = parse.(Int, contrasts_df.timepoint_lhs)
contrasts_df.timepoint_rhs = parse.(Int, contrasts_df.timepoint_rhs)
objects_df = copy(input_rdata["msdata"][string(proj_info.modelobj, "s")]);
obj_contrasts_df = copy(fit_rdata["object_contrasts.df"]);
obj_effects_df = copy(fit_rdata["object_effects.df"]);
for df in [obj_contrasts_df, obj_effects_df]
    df |> MSGLMUtils.fix_quantile_columns!
end
# pick the most significant comparison/effect
obj_contrasts_df.rowix = 1:nrow(obj_contrasts_df)
objcontr_selixs = Vector{Int}()
for df in groupby(obj_contrasts_df, [:contrast, :std_type, objid_col])
    selix = findmin(df.p_value)[2]
    push!(objcontr_selixs, df.rowix[selix])
end
obj_contrasts_df = obj_contrasts_df[objcontr_selixs, :]
obj_contrasts_df.object_id = copy(obj_contrasts_df[!, objid_col])

ptm2gene_df = copy(full_rdata["msdata_full"]["ptm2gene"])
ptm2gene_df.ptm_id = convert(Vector{Int}, ptm2gene_df.ptm_id)
ptm2gene_df.ptm_pos = convert(Vector{Int}, ptm2gene_df.ptm_pos)
ptm2gene_df.ptm_site = ["$(r.ptm_AA_seq)$(r.ptm_pos)" for r in eachrow(ptm2gene_df)]

includet(joinpath(misc_scripts_path, "optcover_utils.jl"));

@load(joinpath(scratch_path, "phosphositeplus_annotations_20200828.jld2"), psitep_annot_dfs, psitep_annots_df, psitepseqs_df)
for df in values(psitep_annot_dfs)
    df.ptm_id = df.protein_ac .* "_" .* uppercase.(df.ptm_AA) .* string.(df.ptm_pos)
end
psitep_annots_df = let kinsub_df = copy(psitep_annot_dfs[:KinaseSubstrate])
    filter!(r -> r.organism=="human", dropmissing!(kinsub_df, :kinase_gene_name))
    kinsub_df[!, :coll_id] .= :KinaseSubstrate
    kinsub_df.term_id = uppercase.(kinsub_df.kinase_gene_name)
    kinsub_df
end

perseus_report_df = CSV.read(joinpath(data_path, proj_info.msfolder,
        "cov2_snaut_parsars_phospho_20201005_contrasts_report_20201012_long_psp.txt"),
        header=1, datarow=3, comment="#", delim='\t')
perseus_motifs_df = dropmissing!(DelimDataUtils.expand_delim_column(select(perseus_report_df, [:protein_ac, :ptm_site, :Motifs]) |> unique!,
                                 list_col=:Motifs, elem_col=:term_id, key_col=[:protein_ac, :ptm_site]),
                                 [:term_id, :protein_ac])
perseus_motifs_df.term_id = replace.(perseus_motifs_df.term_id, Ref(r"\s+(?:motif|sequence)$" => ""))
perseus_motifs_df.ptm_AA = [r.ptm_site[1:1] for r in eachrow(perseus_motifs_df)]
perseus_motifs_df.ptm_pos = [parse(Int, r.ptm_site[2:end]) for r in eachrow(perseus_motifs_df)]
perseus_motifs_df[!, :coll_id] .= :Motifs

motif_colls = FrameUtils.frame2collections(ptmgroup2motif_df, set_col=:term_id, obj_col=:ptmgroup_id, coll_col=:coll_id)

motifXregulators_collapsed_df = CSV.read(joinpath(data_path, proj_info.msfolder, "Perseus_motifsXregulators.txt"))
motifXregulators_df = DelimDataUtils.expand_delim_column(motifXregulators_collapsed_df,
                            list_col=:regulators, elem_col=:regulator, key_col=:motif, delim=",")

function motif2regulators(df::DataFrame; keycol::Union{Symbol, AbstractVector}=:ptmgroup_id)
    dfnew = rename!(select!(join(df, motifXregulators_df,
                         on=[:term_id=>:motif], kind=:left), Not(:term_id)),
            :regulator => :term_id) |> unique!
    dropmissing!(dfnew, [:term_id; keycol])
    dfnew[!, :coll_id] .= :Regulators
    return dfnew
end
perseus_regulators_df = motif2regulators(perseus_motifs_df, keycol=[:protein_ac, :ptm_AA, :ptm_pos])

# https://linkphinder.insight-centre.org/
linkphinder_df = CSV.read(joinpath(party3rd_data_path, "linkPhinder_data.csv"))
countmap(linkphinder_df.KinaseLabel)
linkphinder_annots_df = rename!(filter(r -> r.Score >= 0.75, linkphinder_df), :KinaseLabel => :term_id, :ProteinSubstrate_ID => :protein_ac)
linkphinder_annots_df.ptm_pos = [parse(Int, r.Site[2:end]) for r in eachrow(linkphinder_annots_df)]
linkphinder_annots_df.ptm_AA = getindex.(linkphinder_annots_df.Site, Ref(1:1))
linkphinder_annots_df.coll_id = :LinkPhinder

ptm_annot_cols = [:coll_id, :term_id, :protein_ac, :ptm_AA, :ptm_pos]
ptm_annots_df = vcat(
    [select(df, ptm_annot_cols, copycols=false) for df in
        (psitep_annots_df, perseus_motifs_df,
         perseus_regulators_df, linkphinder_annots_df)]...)

ptmgroup_annots_df = dropmissing!(select!(leftjoin(leftjoin(ptm_annots_df,
            rename!(select(ptm2gene_df, [:ptm_id, :protein_ac, :ptm_AA_seq, :ptm_pos], copycols=false), :ptm_AA_seq => :ptm_AA),
            on=[:protein_ac, :ptm_AA, :ptm_pos]),
            ptm2group_df, on=:ptm_id),
            [:ptmgroup_id, :coll_id, :term_id]), [:ptmgroup_id, :term_id]) |> unique!
countmap(ptmgroup_annots_df.coll_id)

obj_colls = FrameUtils.frame2collections(ptmgroup_annots_df,
    set_col=:term_id, obj_col=:ptmgroup_id, coll_col=:coll_id)

terms_df = unique!(select(ptmgroup_annots_df, [:coll_id, :term_id]));
terms_df.term_name = copy(terms_df.term_id)
terms_df[!, :term_descr] .= missing

@info "Preparing sets of hits"
sel_std_type = "median"
ObjectType = eltype(ptmgroup_annots_df.ptmgroup_id)
obj_hit_sets = Dict{Tuple{String, String}, Set{ObjectType}}()
for hit_df in groupby(filter(r -> coalesce(r.is_hit_composed, false) && (r.std_type == sel_std_type) && !r.is_contaminant, obj_contrasts_df), [:contrast, :change])
    obj_hit_sets[(hit_df[1, :contrast], hit_df[1, :change])] = Set(skipmissing(hit_df.ptmgroup_id))
end
# only relevant ones
sel_contrasts_df = filter(r -> r.contrast_kind == "treatment_vs_treatment" && r.treatment_rhs == "mock"
                          && r.treatment_lhs != "infected", contrasts_df)
obj_hit_selsets = filter(kv -> (kv[1][1] ∈ sel_contrasts_df.contrast), obj_hit_sets)

@info "Preparing mosaics..."
observed_ptms = Set(skipmissing(semijoin(obj_contrasts_df, ptmgroup_annots_df, on=:ptmgroup_id).ptmgroup_id))
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, setXset_frac_extra_elms=0.05, verbose=true);

obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets, max_log10_pvalue=0.0,
                                          max_sets=2000, min_nmasked=1, max_setsize=2000, verbose=true)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.1,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

obj_hit_mosaics_v = collect(pairs(obj_hit_mosaics))
obj_hit_covers_v = similar(obj_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_hit_mosaics_v[i]
    @info "Covering $mosaic_name by hits..."
    obj_hit_covers_v[i] = mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.02], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=1, #Threads.nthreads()-1,
                                    MaxRestarts=200),
            true)
end
obj_hit_covers = Dict(k => v for (k, v) in obj_hit_covers_v)

using JLD2

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_$(proj_info.oesc_ver).jld2")
@save(hit_covers_filename,
      proj_info, obj_colls, obj_mosaics,
      terms_df, ptm_annots_df, ptmgroup_annots_df, contrasts_df,
      ptmns_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
if !@isdefined(obj_effect_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, obj_colls, obj_mosaics,
      terms_df, ptm_annots_df, ptmgroup_annots_df, contrasts_df,
      ptmns_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
end

includet(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup↦gene_name map..."
obj_id2name = Dict(r.ptmgroup_id => replace(r.ptm_label, r"^[^_]+_" => "")
                   for r in eachrow(filter(r -> r.ptmid_is_reference, ptmns_df)))

obj_hit_covers_df = innerjoin(OptCoverUtils.covers_report(
    obj_hit_covers, obj_hit_selsets, obj_colls, obj_mosaics, obj_id2name,
    terms_df,
    maskid_col=[:contrast, :change],
    maskedset_col_prefix="contrast"),
    contrasts_df, on=[:contrast])
filter!(r -> r.term_id != "CHEK2" && r.term_id != "CK1E", obj_hit_covers_df) # redundant with Chk2 and CSK1E

obj_hit_covers_signif_df = combine(groupby(obj_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:contrast, :change],
                                                   max_term_pvalue=1E-2, max_set_pvalue=nothing, min_set_overlap=nothing),
                   Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_$(proj_info.oesc_ver).txt"),
          filter(r -> r.nmasked > 0, obj_hit_covers_df),
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          filter(r -> r.nmasked > 0, obj_hit_covers_signif_df),
          missingstring="", delim='\t');

includet(joinpath(misc_scripts_path, "frame_utils.jl"))
includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

using PlotlyJS, TextWrap

heatmap_layout_attrs = Dict(
    (:Motifs, true) => Dict(:margin_l => 300, :margin_b => 200),
    (:Motifs, false) => Dict(:margin_l => 300, :margin_b => 200),
    (:Regulators, true) => Dict(:margin_l => 150, :margin_b => 200),
    (:Regulators, false) => Dict(:margin_l => 150, :margin_b => 200),
)

stylize_contrast(str) = foldl(replace, [
    r"(SARS_CoV2?)@(\d+)h_vs_mock@(\d+)h" => s"\1:<span style=\"font-weight: bold; color: black;\">\2</span>h",
    "SARS_CoV:" => "<span style=\"font-wieght: bold; color: #811A02;\">SARS</span> ",
    "SARS_CoV2:" => "<span style=\"font-wieght: bold; color: #F4982A;\">CoV2</span> ",
    ],
    init = str)

function process_contrast_axis(contrast_df)
    contrast_df,
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change),
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change)#stylize_effect.(effect_df.effect)
end

heatmaps_path = joinpath(plots_path, "$(proj_info.msfolder)_$(proj_info.fit_ver)", "hits_oesc_$(sel_std_type)_$(proj_info.oesc_ver)")
isdir(heatmaps_path) || mkdir(heatmaps_path)

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")hit heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_hit_covers_signif_df : obj_hit_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end

    for outformat in ("html", "pdf", "svg")
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="PTM",
            mask_axis_title = "contrast",
            mask_cols = [:contrast, :treatment_lhs, :timepoint_lhs, :change, :ncontrast],
            process_mask_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 80),
            margin_b=get(layout_attrs, :margin_b, 100),
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",#outformat in ["svg", "pdf"] ? "#000" : "#BBB",
            zmin=-10, cell_width=18, cell_height=14, gridwidth=1,
            mask_order=contrasts -> begin
                return reverse(sortperm(contrasts, [order(:change, rev=true), :timepoint_lhs, :treatment_lhs]))
            end,
            transpose=true)
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>60,
                   :margin_t=>20,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    plot_fname = joinpath(heatmaps_path,
                          "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_X_contrast$(signif ? "_signif" : "")_heatmap.$(outformat)")
    if outformat == "html"
        PlotlyJS.savehtml(coll_heatmap, plot_fname, :embed);
    else
        try
            savefig(coll_heatmap, plot_fname, width=coll_heatmap.plot.layout[:width], height=coll_heatmap.plot.layout[:height]);
        catch e
            @warn e
        end
    end
    end
end

includet(joinpath(misc_scripts_path, "subheatmap_utils.jl"))
using VegaLite

# triheatmap
triheatmaps_path = joinpath(plots_path, "$(proj_info.msfolder)_$(proj_info.fit_ver)", "hits_oesc_$(sel_std_type)_$(proj_info.oesc_ver)_tri")
isdir(triheatmaps_path) || mkdir(triheatmaps_path)

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "") tri-heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    cover_df = filter(r -> r.term_collection == term_coll, signif ? obj_hit_covers_signif_df : obj_hit_covers_df)
    if nrow(cover_df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end
    scores_mtx, _, term_axis, contrast_axis = OptCoverHeatmap.heatmap_matrices(cover_df,
            elements_label="protein",
            mask_axis_title = "contrast",
            mask_cols = [:contrast, :treatment_lhs, :timepoint_lhs, :treatment_rhs, :timepoint_rhs, :change, :ncontrast],
            process_mask_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            mask_order=contrasts -> begin
                return sortperm(contrasts, [:change, :timepoint_lhs, :treatment_lhs])
        end)
    rename!(contrast_axis, :axis_label => :axisXtri_label, :axis_tip => :axisXtri_tip)
    contrast_axis.timepointXdir_label = string.(contrast_axis.timepoint_lhs) .* "h " .* ifelse.(contrast_axis.change .== "+", "▲", "▼")
    contrast_axis.treatment_lhs = levels!(categorical(contrast_axis.treatment_lhs), ["mock", "SARS_CoV2", "SARS_CoV"])
    contrast_axis.treatment_rhs = FrameUtils.matchcategorical(contrast_axis.treatment_rhs, contrast_axis.treatment_lhs)
    contrast_axis.contrast_label = stylize_contrast.(contrast_axis.contrast)
    
    triheatmap_df, triheatmap_rows_df, triheatmap_cols_df = SubheatmapUtils.subheatmap_frame(scores_mtx, term_axis, contrast_axis, #tips=tips_mtx,
        row_label_col=:term_id, col_label_col=:timepointXdir_label,
        col_cols=[:timepointXdir_label, :timepoint_lhs, :change],
        col_sub_cols = [:treatment_lhs, :contrast, :ncontrast, :timepointXdir_label])
    rename!(triheatmap_df, :row_index=>:col_index, :col_index=>:row_index, :row_label=>:col_label, :col_label=>:row_label)
    triheatmap_plot = SubheatmapUtils.vegalite_subheatmap(triheatmap_df, value_domain=(-10, 0),
                                          yaxis_label="Time, h.p.i", xaxis_label="Term",
                                          subaxis_label="Virus", coloraxis_label="log₁₀(P-value)",
                                          subshapes=reverse(SubheatmapUtils.DiagonalTriangles(scale=0.6)),
                                          xaxis_tick_col=:term_id, yaxis_tick_col=:timepointXdir_label,
                                          labelLimit_l=50, labelLimit_b=get(layout_attrs, :margin_b, 120),
                                          cell_width=12, cell_height=12)

    plot_fname = joinpath(triheatmaps_path,
        "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_contrast$(signif ? "_lowsignif" : "")_triheatmap")
    save(plot_fname * ".html", triheatmap_plot)
    save(plot_fname * ".svg", triheatmap_plot)
    save(plot_fname * ".pdf", triheatmap_plot)
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
