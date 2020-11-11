proj_info = (id = "cov2", ver="20201028")
using Pkg
Pkg.activate(@__DIR__)
using Revise
using DataFrames, CSV, Statistics, StatsBase, RData
using PlotlyBase, PlotlyJS, TextWrap
using Printf: @sprintf

@info "Project '$(proj_info.id)'"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

includet(joinpath(misc_scripts_path, "frame_utils.jl"))
includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

#ura_ipa_df = CSV.read(joinpath(data_path, "20200519_URA_IPA_table.txt"), delim='\t', header=3)

encode_tfea_df = CSV.read(joinpath(analysis_path, "reports", "cov2_rnaseq_tfea_20201028.txt"), missingstrings=["NA"], delim='\t')
encode_tfea_df.log10_FET_pvalue = log10.(max.(encode_tfea_df[!, "FET.p.value"], 1E-20))
rename!(encode_tfea_df, :contrast => :contrastXdir)
contrast_matches = match.(Ref(r"^((SARS_CoV2?)@(\d+)h_vs_(SARS_CoV|mock)@(\d+)h)_(up|down)$"), encode_tfea_df.contrastXdir)
@assert all(!isnothing, contrast_matches)
encode_tfea_df.contrast = String.(getindex.(contrast_matches, 1))
encode_tfea_df.treatment_lhs = String.(getindex.(contrast_matches, 2))
encode_tfea_df.timepoint_lhs = parse.(Int, getindex.(contrast_matches, 3))
encode_tfea_df.treatment_rhs = String.(getindex.(contrast_matches, 4))
encode_tfea_df.timepoint_rhs = parse.(Int, getindex.(contrast_matches, 5))
encode_tfea_df.change = ifelse.(getindex.(contrast_matches, 6) .== "up", "+", "-")
TF_Set_usecount = countmap(encode_tfea_df.Set_name)
encode_tfea_df.Set_n_used = getindex.(Ref(TF_Set_usecount), encode_tfea_df.Set_name)

encode_tfea_mtx, (TF_axis, contrast_axis) =
    FrameUtils.frame2array(encode_tfea_df, [[:TF], [:contrastXdir, :contrast, :treatment_lhs, :timepoint_lhs, :treatment_rhs, :timepoint_rhs, :change]],
                           data_col=:log10_FET_pvalue)

stylize_contrast(str) = foldl(replace, [
    r"(SARS_CoV2?)@(\d+)h_vs_mock@(\d+)h" => s"\1:<span style=\"font-weight: bold; color: black;\">\2</span>h",
    "SARS_CoV:" => "<span style=\"font-wieght: bold; color: #811A02;\">SARS</span> ",
    "SARS_CoV2:" => "<span style=\"font-wieght: bold; color: #F4982A;\">CoV2</span> ",
    ],
    init = str)

encode_tfea_df.tooltip = [begin
    genes = sort!(split(coalesce(r.Overlapping_Genes, ""), ','))
    "TF=<b>$(r.TF)</b> ($(r.Set_name))<br>\n" .*
    "contrast=<b>$(stylize_contrast(r.contrast))</b> $(OptCoverHeatmap.stylize_change(r.change))<br>\n" .*
    "P-value=$(@sprintf("%.4e", r["FET.p.value"]))<br>\n" .*
    "Genes (N=$(length(genes))):<br>\n" *
    replace(TextWrap.wrap(join(genes, ' '), break_long_words=false), '\n' => "<br>\n")
end for r in eachrow(encode_tfea_df)]

encode_tfea_tip_mtx, _ =
    FrameUtils.frame2array(encode_tfea_df, [[:TF], [:contrastXdir]], data_col=:tooltip)

contrast_axis.treatment_lhs = levels!(categorical(contrast_axis.treatment_lhs), ["SARS_CoV2", "SARS_CoV", "mock"])
contrast_axis.treatment_rhs = FrameUtils.matchcategorical(contrast_axis.treatment_rhs, contrast_axis.treatment_lhs)
contrast_axis.contrast_label = stylize_contrast.(contrast_axis.contrast)
contrast_axis.contrastXchange_label = contrast_axis.contrast_label .* "&nbsp;:" .* OptCoverHeatmap.stylize_change(contrast_axis.change)

contrast_perm = sortperm(contrast_axis, [:change, :timepoint_lhs, :treatment_rhs, :treatment_lhs], rev=true)
contrast_axis = contrast_axis[contrast_perm, :]
encode_tfea_mtx = encode_tfea_mtx[:, contrast_perm]
encode_tfea_tip_mtx = encode_tfea_tip_mtx[:, contrast_perm]

#condition_axis = condition_axis[conditions_perm, :]
#encode_tfea_mtx = encode_tfea_mtx[conditions_perm, :]
TF_mask = vec(minimum(coalesce.(encode_tfea_mtx, Inf), dims=2) .<= log10(1E-4))
sel_encode_tfea_mtx = encode_tfea_mtx[TF_mask, :]
sum(TF_mask)
size(encode_tfea_mtx)

duplicate_TFs = Set(TF_axis[TF_mask, :TF][nonunique(TF_axis[TF_mask, [:TF]])])
TF_axis.Set_label = copy(TF_axis.TF)
#TF_axis[TF_axis.TF .∈ Ref(duplicate_TFs), :Set_label] .= TF_axis[TF_axis.TF .∈ Ref(duplicate_TFs), :Set_name]

using Distances, Clustering
TF_dist = pairwise(Euclidean(), coalesce.(sel_encode_tfea_mtx, 0.0), dims=1)
TF_dist[isnan.(TF_dist)] .= 100.0
TF_hclu = hclust(TF_dist, linkage=:ward)
TF_order = axes(encode_tfea_tip_mtx, 1)[TF_mask][TF_hclu.order]

# classical heatmap
encode_heatmap = plot(heatmap(z=coalesce.(encode_tfea_mtx[TF_order, :], 0.0),
                              text=encode_tfea_tip_mtx[TF_order, :],
            hoverinfo="text",
            y=contrast_axis.contrastXchange_label, x=TF_axis.Set_label[TF_order],
            zmin=-10, zmax=0, colorscale="Hot", reversescale=false,
        #zmax=zmax, zmin=zmin, zauto=false, connectgaps=false,
        xtype="array", ytype="array", xgap=1, ygap=1),
        Layout(#paper_bgcolor=paper_bgcolor,
               #plot_bgcolor=plot_bgcolor,
               modebar_bgcolor="#FFF",
               xaxis = attr(tickson="boundaries", gridcolor="#DDD", gridwidth=1, tickangle=45),
               yaxis = attr(tickson="boundaries", gridcolor="#DDD", gridwidth=1),
               width = 1800, height = 300, margin_l = 120, margin_b = 100, margin_t=20))
savefig(encode_heatmap, joinpath(plots_path, "$(proj_info.id)_ENCODE_TFEA_signif_$(proj_info.ver).html"),
        width=encode_heatmap.plot.layout[:width], height=encode_heatmap.plot.layout[:height]);

includet(joinpath(misc_scripts_path, "subheatmap_utils.jl"))
using VegaLite

# triheatmap
contrast_axis.timepointXdir_label = string.(contrast_axis.timepoint_lhs) .* "h " .* ifelse.(contrast_axis.change .== "+", "▲", "▼")
contrast_axis.axisXtri_label = contrast_axis.contrastXchange_label
triheatmap_df, triheatmap_cols_df, triheatmap_rows_df = SubheatmapUtils.subheatmap_frame(encode_tfea_mtx[TF_order, nrow(contrast_axis):-1:1],
                                                                                         TF_axis[TF_order, :], contrast_axis[nrow(contrast_axis):-1:1, :], #tips=tips_mtx,
        col_label_col=:timepointXdir_label, row_label_col=:TF,
        col_sub_cols=[:treatment_lhs, :axisXtri_label],
        col_cols=[:timepointXdir_label, :timepoint_lhs, :change])
rename!(triheatmap_df, :row_index=>:col_index, :col_index=>:row_index, :row_label=>:col_label, :col_label=>:row_label)
triheatmap_df.rowblock_index = min.(fld1.(triheatmap_df.col_index, nrow(triheatmap_cols_df) ÷ 2), 2)
triheatmap_plot = SubheatmapUtils.vegalite_subheatmap(triheatmap_df, value_domain=(-10, 0),
                                      xaxis_label="Transcription Factors", yaxis_label="",
                                      subaxis_label="Virus", coloraxis_label="log₁₀(P-value)",
                                      xaxis_tick_col=:TF, yaxis_tick_col=:timepointXdir_label,
                                      subshapes=SubheatmapUtils.DiagonalTriangles(scale=0.6),
                                      labelLimit_l=120, cell_width=12, cell_height=12)
for ext in ["html", "pdf", "svg"]
    save(joinpath(plots_path, "$(proj_info.id)_ENCODE_TFEA_signif_$(proj_info.ver)_tri.$ext"), triheatmap_plot)
end