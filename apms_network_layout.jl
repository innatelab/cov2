proj_info = (id = "cov2",
             data_ver = "20200329",
             fit_ver = "20200402",
             modelobj = "protgroup",
             mq_folder = "mq_apms_20200329",
             network_ver = "20200402")
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

Revise.includet(joinpath(misc_scripts_path, "clustering_utils.jl"))

network_rdata = load(joinpath(networks_path, "$(proj_info.id)_4graph_$(proj_info.mq_folder)_$(proj_info.fit_ver).RData"))

objects_orig_df = network_rdata["objects_4graphml.df"]
iactions_orig_df = network_rdata["iactions_ex_4graphml.df"]

nodesize_scale = 0.015
node_type_props = Dict("bait" => (mass=2.0, width=6nodesize_scale, height=6nodesize_scale),
                       "prey" => (mass=2.0, width=3nodesize_scale, height=3nodesize_scale))
objects_df = copy(objects_orig_df)
objects_df.mass = getproperty.(getindex.(Ref(node_type_props), objects_df.exp_role), :mass);
objects_df.width = getproperty.(getindex.(Ref(node_type_props), objects_df.exp_role), :width);
objects_df.height = getproperty.(getindex.(Ref(node_type_props), objects_df.exp_role), :height);
sort!(objects_df, :object_id)

ppi_weight_scales = Dict("ppi_low" => 0.05,
                         "ppi_medium" => 0.05,
                         "ppi_strong" => 0.05,
                         "complex" => 0.05,
                         "experiment" => 1.0)
iactions_df = filter(r -> r.src_object_id != r.dest_object_id, iactions_orig_df)
iactions_df.src_object_id = convert(Vector{Int}, iactions_df.src_object_id)
iactions_df.dest_object_id = convert(Vector{Int}, iactions_df.dest_object_id)
#filter!(r -> coalesce(r.ppi_type, "") ∉ ["enriched_complex_gocc"], iactions_df)
iactions_df.edge_weight =
    get.(Ref(ppi_weight_scales), coalesce.(iactions_df.type, "experiment"), 0.1) .*
    iactions_df.weight
for obj_iactions_df in groupby(iactions_df, :src_object_id)
    obj_id = obj_iactions_df.src_object_id[1]
    node_ix = searchsortedfirst(objects_df.object_id, obj_id)
    @assert objects_df.object_id[node_ix] == obj_id
    if objects_df.exp_role[node_ix] != "bait"
        #@warn "Node $obj_id is not a bait"
        continue
    end
    exp_iactions_mask = ismissing.(obj_iactions_df.ppi_type)
    max_weight = quantile(obj_iactions_df[exp_iactions_mask, :edge_weight], 0.75)
    obj_iactions_df[exp_iactions_mask, :edge_weight] .= clamp.(obj_iactions_df[exp_iactions_mask, :edge_weight] ./ max_weight, 0.1, 1.25)
end

#using StatsBase
#countmap(iactions_df.ppi_type[.!ismissing.(iactions_df.ppi_type)])

using LightGraphs

Revise.includet(joinpath(misc_scripts_path, "forceatlas3_layout.jl"))
FA = ForceAtlas3

gr_apms = SimpleGraph(FA.Graph(iactions_df[ismissing.(iactions_df.ppi_type), :], objects_df,
                       src_node_col=:src_object_id, dest_node_col=:dest_object_id, weight_col=:edge_weight,
                       node_col=:object_id))

#objects_df.object_label[objects_df.object_label .∈ Ref(Set(["4", "9", "12", "28"]))]
#node_dislikes_baits[objects_df.object_label .∈ Ref(Set(["4", "9", "12", "28"])),
#              objects_df.object_label .∈ Ref(Set(["4", "9", "12", "28"]))]

node_dislikes = FA.socioaffinity(gr_apms, p=(0.0, 1.0), q=1.0)
node_dislikes_baits = FA.socioaffinity(gr_apms, p=(1.25, 1.25), q=2.0)
baits_mask = objects_df.exp_role .== "bait"
node_dislikes[baits_mask, baits_mask] .= 4 * node_dislikes_baits[baits_mask, baits_mask]

gr = FA.Graph(iactions_df, objects_df,
              src_node_col=:src_object_id, dest_node_col=:dest_object_id, weight_col=:edge_weight,
              #gravitable_col=:gravitable,
              node_col=:object_id, mass_col=:mass,
              width_col=:width, height_col=:height)
FA.layout!(gr, FA.ForceAtlas3Settings(gr,
            outboundAttractionDistribution=false,
            attractionStrength=10.0, attractionEdgeWeightInfluence=0.5, jitterTolerance=0.1,
            repulsionStrength=0.1.*(1.0 .+ node_dislikes),
            repulsionNodeModel=:Point,
            gravity=1.0, gravityFalloff=1.3, gravityShape=:Rod,
            gravityRodCorners=((0.0, -4.0), (0.0, 4.0)), gravityRodCenterWeight=0.1),
            nsteps=1000, progressbar=true)
FA.layout!(gr, FA.ForceAtlas3Settings(gr,
            outboundAttractionDistribution=false,
            attractionStrength=8.0, attractionEdgeWeightInfluence=0.5, jitterTolerance=0.1,
            repulsionStrength=4.0.*(1.0 .+ node_dislikes),
            repulsionNodeModel=:Circle,
            gravity=1.5, gravityFalloff=1.3, gravityShape=:Rod,
            gravityRodCorners=((0.0, -4.0), (0.0, 4.0)), gravityRodCenterWeight=0.1),
            nsteps=5000, progressbar=true)

objects_df.layout_x = FA.extract_layout(gr)[1] .* 20
objects_df.layout_y = FA.extract_layout(gr)[2] .* 20

includet(joinpath(misc_scripts_path, "graphml_writer.jl"))

objects_df.object_idstr = string.(objects_df.object_id)
iactions_df.src_object_idstr = string.(iactions_df.src_object_id)
iactions_df.dest_object_idstr = string.(iactions_df.dest_object_id)

apms_graph = GraphML.import_graph(objects_df, iactions_df,
                             node_col=:object_idstr,
                             source_col=:src_object_idstr, target_col=:dest_object_idstr,
                             node_attrs = [:object_idstr => "Protein Group ID",
                                           :object_label => "Protein Group Label",
                                           :majority_protein_acs => "Majority ACs",
                                           :gene_names => "Gene Names",
                                           :exp_role => "Experimental Role",
                                           :protein_names => "Protein Names",
                                           :protein_description => "Protein Description",
                                           :protein_class => "Protein Class",
                                           :organism => "Organism",
                                           :seqlen => "Seq Length",
                                           :is_detected => "Detected",
                                           :crispr_plasmid_ids => "CRISPR Plasmid IDs",
                                           :layout_x, :layout_y],
                             edge_attrs = [Symbol("prob_nonpos_min_vs_background") => "P-value (vs Background)",
                                           Symbol("median_log2_vs_background") => "Enrichment (vs Background)",
                                           Symbol("edge_weight") => "Weight",
                                           :type => "type",
                                           :krogan_is_hit => "Krogan hit",
                                           :krogan_MIST => "Krogan MIST score",
                                           :krogan_avg_spec => "Krogan Average SC",
                                           :krogan_fold_change => "Krogan Fold Change",
                                           :virhostnet_confidence => "VirHostNet Confidence",
                                           :virhostnet_references => "VirHostNet PubMeds",
                                           #`Known types` = "known_types",
                                           :iaction_ids => "Known Interaction IDs"])
 #edge.attrs = c( `P-value (vs Background)` = 'p_value_min.vs_background',
 #                 `P-value (WT vs Mock)` = 'p_value.SC35MWT',
 #                 `P-value (delNS1 vs Mock)` = 'p_value.SC35MdelNS1',
 #                 `Enrichment (vs Background)` = 'median_log2.vs_background',
 #                 `Enrichment (WT vs Mock)` = 'median_log2.SC35MWT',
 #                 `Enrichment (delNS1 vs Mock)` = 'median_log2.SC35MdelNS1',
 #                 `Weight` = 'weight',
 #                 `type` = "type" )
                             # verbose=verbose)
open(joinpath(networks_path, "$(proj_info.id)_4graph_$(proj_info.mq_folder)_$(proj_info.fit_ver)_FA3.graphml"), "w") do io
    write(io, apms_graph)
end
