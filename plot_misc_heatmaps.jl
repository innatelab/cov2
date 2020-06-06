proj_info = (id = "cov2",)
using Pkg
Pkg.activate(@__DIR__)
using Revise
using DataFrames, CSV, Statistics, RData
using PlotlyBase, PlotlyJS

@info "Project '$(proj_info.id)'"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

ura_ipa_df = CSV.read(joinpath(data_path, "20200519_URA_IPA_table.txt"), delim='\t', header=3)

using ORCA

encode_tfea_df = CSV.read(joinpath(data_path, "20200519_ENCODE_TFEA_all.txt"), delim='\t')
encode_tfea_df.log10_FET_pvalue = log10.(encode_tfea_df.FET_pvalue)
encode_tfea_mtx, (condition_axis, term_axis) =
    FrameUtils.frame2array(encode_tfea_df, [[:timepoint, :direction], [:TF, :Set_name]],
                           data_col=:log10_FET_pvalue)
condition_axis.change = ifelse.(condition_axis.direction .== "up", "+", "-")
condition_axis.condition = string.(condition_axis.timepoint) .* "h " .* condition_axis.direction
condition_axis.condition_label = "<span style=\"font-weight: bold\">" .* string.(condition_axis.timepoint) .*
    "</span>h:&nbsp;" .* OptCoverHeatmap.stylize_change.(condition_axis.change)

conditions_perm = sortperm(condition_axis, [:direction, :timepoint], rev=true)
#condition_axis = condition_axis[conditions_perm, :]
#encode_tfea_mtx = encode_tfea_mtx[conditions_perm, :]
term_mask = vec(minimum(coalesce.(encode_tfea_mtx, Inf), dims=1) .<= log10(0.005))
encode_tfea_mtx = encode_tfea_mtx[:, term_mask]
term_axis = term_axis[term_mask, :]
sum(term_mask)
size(encode_tfea_mtx)

using Distances, Clustering
term_dist = pairwise(Euclidean(), coalesce.(encode_tfea_mtx, 0.0), dims=2)
term_dist[isnan.(term_dist)] .= 100.0
term_hclu = hclust(term_dist, linkage=:ward)

encode_heatmap = plot(heatmap(z=coalesce.(encode_tfea_mtx[conditions_perm, term_hclu.order]', 0.0),
            hoverinfo="text",
            y=condition_axis.condition_label[conditions_perm], x=term_axis.TF[term_hclu.order],
        colorscale="Blackbody", reversescale=true,
        #zmax=zmax, zmin=zmin, zauto=false, connectgaps=false,
        xtype="array", ytype="array", xgap=1, ygap=1),
        Layout(#paper_bgcolor=paper_bgcolor,
               #plot_bgcolor=plot_bgcolor,
               modebar_bgcolor="#FFF",
               xaxis = attr(tickson="boundaries", tickangle=45),
               yaxis = attr(tickson="boundaries"),
               width = 1200, height = 250, margin_l = 60, margin_t=20))
savefig(encode_heatmap.plot, joinpath(plots_path, "$(proj_info.id)_ENCODE_TFEA_all_20200519.pdf"));
