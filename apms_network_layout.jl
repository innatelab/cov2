proj_info = (id = "cov2",
             data_ver = "20200525",
             fit_ver = "20200525",
             msfolder = "mq_apms_20200525",
             prev_network_ver = "20200515",
             network_ver = "20200606")
using Pkg
Pkg.activate(@__DIR__)
using Revise
using DataFrames, CSV, Statistics, StatsBase, RData

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const networks_path = joinpath(analysis_path, "networks")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "clustering_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "delimdata_utils.jl"))

input_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(proj_info.msfolder)_$(proj_info.fit_ver).RData"))
network_rdata = load(joinpath(networks_path, "$(proj_info.id)_4graph_$(proj_info.msfolder)_$(proj_info.fit_ver).RData"))

objects_orig_df = network_rdata["objects_4graphml.df"]
iactions_orig_df = network_rdata["iactions_ex_4graphml.df"]

modelobjs_df = copy(input_rdata["msdata"]["protregroups"])
modelobjs_df.object_id = copy(modelobjs_df.protregroup_id)
modelobjs_df.object_label = copy(modelobjs_df.protregroup_label)
modelobjs_df.chunk = 1:nrow(modelobjs_df)
processed_object_ids = Set(network_rdata["object_contrasts_slim.df"].object_id)
toprocess_objs_df = vcat(select(filter(r -> #(r.object_id ∉ processed_object_ids) ||
        occursin(r"(?:^|;)ITG", coalesce(r.gene_names, "")), modelobjs_df),
       [:chunk, :object_id, :object_label, :gene_names]),
#=toprocess_objs_df = =#
select(semijoin(modelobjs_df, objects_orig_df, on=:object_id),
      [:chunk, :object_id, :object_label, :gene_names]))
unique(toprocess_objs_df.chunk)

nodesize_scale = 0.006
node_type_props = Dict("bait" => (mass=2.0, width=6nodesize_scale, height=6nodesize_scale),
                       "prey" => (mass=2.0, width=3nodesize_scale, height=3nodesize_scale))
objects_df = copy(objects_orig_df)
objects_df.object_label = ifelse.(objects_df.exp_role .== "bait",
                                  replace.(objects_df.object_label, r"(CoV\d*)_" => s"\1\n"),
                                  objects_df.object_label)
objects_df.organism = ifelse.(objects_df.is_bait .& endswith.(objects_df.object_label, "?"),
                              objects_df.organism .* "?", objects_df.organism)
objects_df.mass = getproperty.(getindex.(Ref(node_type_props), objects_df.exp_role), :mass);
objects_df.width = getproperty.(getindex.(Ref(node_type_props), objects_df.exp_role), :width);
objects_df.height = getproperty.(getindex.(Ref(node_type_props), objects_df.exp_role), :height);
sort!(objects_df, :object_id)

ppi_weight_scales = Dict("ppi_low" => 0.05,
                         "ppi_medium" => 0.05,
                         "ppi_strong" => 0.05,
                         "complex" => 0.5,
                         "homology" => 1.0,
                         "experiment" => 1.0)
iactions_df = filter(r -> r.src_object_id != r.dest_object_id, iactions_orig_df)
iactions_df.contrast = iactions_df.bait_full_id .* "_vs_others"
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
    obj_iactions_df[exp_iactions_mask, :edge_weight] .= clamp.(obj_iactions_df[exp_iactions_mask, :edge_weight] ./ max_weight, 0.5, 1.25)
end

#using StatsBase
#countmap(iactions_df.ppi_type[.!ismissing.(iactions_df.ppi_type)])

using LightGraphs

Revise.includet(joinpath(misc_scripts_path, "forceatlas3_layout.jl"))
includet(joinpath(misc_scripts_path, "graphml_writer.jl"))

FA = ForceAtlas3

gr_apms = SimpleGraph(FA.Graph(iactions_df[ismissing.(iactions_df.ppi_type), :], objects_df,
                       src_node_col=:src_object_id, dest_node_col=:dest_object_id, weight_col=:edge_weight,
                       node_col=:object_id))

#objects_df.object_label[objects_df.object_label .∈ Ref(Set(["4", "9", "12", "28"]))]
#node_dislikes_baits[objects_df.object_label .∈ Ref(Set(["4", "9", "12", "28"])),
#              objects_df.object_label .∈ Ref(Set(["4", "9", "12", "28"]))]

node_dislikes = FA.socioaffinity(gr_apms, p=(0.0, 1.0), q=1.0)
node_dislikes_baits = FA.socioaffinity(gr_apms, p=(1.25, 1.25), q=2.0)
for (i, r1) in enumerate(eachrow(objects_df)), (j, r2) in enumerate(eachrow(objects_df))
    if r1.exp_role == "bait" && r2.exp_role == "bait" &&
       replace(r1.gene_label, r"\?$" => "") != replace(r2.gene_label, r"\?$" => "")
       node_dislikes[i, j] = 15 * node_dislikes_baits[i, j]
    end
end

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
            attractionStrength=4.0, attractionEdgeWeightInfluence=0.75, jitterTolerance=0.1,
            repulsionStrength=3 .* (1.0 .+ node_dislikes),
            repulsionNodeModel=:Circle,
            gravity=0.5, gravityFalloff=1.5, gravityShape=:Rod,
            gravityRodCorners=((0.0, -6.0), (0.0, 6.0)), gravityRodCenterWeight=0.1),
            nsteps=5000, progressbar=true)

objects_df.object_label = replace.(objects_df.object_label, Ref(r"\.+$" => ""))
objects_df.layout_x = FA.extract_layout(gr)[1] .* 30
objects_df.layout_y = FA.extract_layout(gr)[2] .* 30

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
                                           :oeproteome_is_hit => "OE Proteome Hit",
                                           :oeproteome_bait_full_ids => "OE Proteome Baits",
                                           :oeproteome_p_value => "OE Proteome: Most Signif. P-value",
                                           :oeproteome_median_log2 => "OE Proteome: Most Signif. Log2(Fold-Change)",
                                           :cov2ts_proteome_timepoints => "CoV-2 Proteome: Timepoints of Significant Changes",
                                           :cov2ts_proteome_p_value => "CoV-2 Proteome: Most Signif. P-value",
                                           :cov2ts_proteome_median_log2 => "CoV-2 Proteome: Most Signif. Median Log2",
                                           :cov2ts_proteome_is_hit => "CoV-2 Proteome: Is Hit",
                                           :cov2ts_phospho_ptms => "CoV-2 Phospho: PTMs and Timepoints of Significant Changes",
                                           :cov2ts_phospho_p_value => "CoV-2 Phospho: Most Signif. P-value",
                                           :cov2ts_phospho_median_log2 => "CoV-2 Phospho: Most Signif. Median Log2",
                                           :cov2ts_phospho_is_hit => "CoV-2 Phospho: Is Hit",
                                           :layout_x, :layout_y],
                             edge_attrs = [Symbol("prob_nonpos") => "P-value (vs Background)",
                                           Symbol("median_log2") => "Enrichment (vs Background)",
                                           Symbol("edge_weight") => "Weight",
                                           :type => "type",
                                           :contrast_carryover => "Carryover test",
                                           :median_log2_carryover => "Carryover Log2(Fold-Change)",
                                           :p_value_carryover => "Carryover P-value",
                                           :contrast_batch => "Batch-specific test",
                                           :median_log2_batch => "Batch-specific Log2(Fold-Change)",
                                           :p_value_batch => "Batch-specific P-value",
                                           :oeproteome_is_hit => "OE Proteome Hit",
                                           :oeproteome_p_value => "OE Proteome P-value",
                                           :oeproteome_median_log2 => "OE Proteome Median Log2",
                                           Symbol(string("is_in_", proj_info.prev_network_ver, "_apms")) => "Is In $(proj_info.prev_network_ver) AP-MS Network",
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
open(joinpath(networks_path, "$(proj_info.id)_4graph_$(proj_info.msfolder)_$(proj_info.fit_ver)_FA3.graphml"), "w") do io
    write(io, apms_graph)
end

cov1_organisms = ["SARS-CoV", "SARS-CoV-GZ02"]
cov2_organisms = ["SARS-CoV-2"]

baits_df = network_rdata["bait_labels.df"]
bait_nodes_df = join(baits_df, objects_df, on=:bait_full_id=>:object_bait_full_id)
bait_nodes_df.is_alt = occursin.(Ref(r"\?$"), bait_nodes_df.bait_full_id)
bait_nodes_df.is_merged = occursin.(Ref("SARS"), bait_nodes_df.organism) .& .! bait_nodes_df.is_alt
bait_nodes_df.bait_homid = copy(bait_nodes_df.bait_id) # stricter "homology" match
bait_nodes_df.is_used = (bait_nodes_df.is_merged .| occursin.(Ref(r"ORF[34]"), String.(bait_nodes_df.bait_full_id))) .&
    .!bait_nodes_df.is_alt
SkipBaitId = -10000
bait_nodes_df[!, :baitmerge_object_id] .= SkipBaitId
bait_nodes_df[bait_nodes_df.is_used, :baitmerge_object_id] .= bait_nodes_df.object_id[bait_nodes_df.is_used]
for hombaits_df in groupby(bait_nodes_df, [:bait_homid])
    if any(hombaits_df.is_merged)
        hombaits_df[hombaits_df.is_merged, :baitmerge_object_id] .=
            minimum(hombaits_df.object_id[hombaits_df.is_merged])
    end
end
objects_genemerge_map_df = combine(groupby(objects_df, :object_label)) do df
    (nrow(df) == 1) && return DataFrame(object_id = Vector{eltype(df.object_id)}(),
                                        baitmerge_object_id = Vector{eltype(df.object_id)}())
    DataFrame(object_id = df.object_id,
              baitmerge_object_id = repeat([minimum(df.object_id)], nrow(df)))
end
select!(objects_genemerge_map_df, Not(:object_label))
objects_baitmerge_map_df = leftjoin(objects_df, vcat(select(bait_nodes_df, [:object_id, :baitmerge_object_id]),
                                                     objects_genemerge_map_df),
                                    on=:object_id)
objects_baitmerge_map_df.old_object_id = copy(objects_baitmerge_map_df.object_id)
objects_baitmerge_map_df.object_label = ifelse.(ismissing.(objects_baitmerge_map_df.baitmerge_object_id),
                                                objects_baitmerge_map_df.object_label,
                                                objects_baitmerge_map_df.gene_label)
objects_baitmerge_map_df.object_id = ifelse.(ismissing.(objects_baitmerge_map_df.baitmerge_object_id),
                                             objects_baitmerge_map_df.object_id,
                                             objects_baitmerge_map_df.baitmerge_object_id)
print(select(filter(r -> r.object_id == SkipBaitId, objects_baitmerge_map_df),
       [:object_label, :object_id, :baitmerge_object_id, :old_object_id]))
filter!(r -> r.object_id != SkipBaitId, objects_baitmerge_map_df)
objects_baitmerge_df = combine(groupby(objects_baitmerge_map_df, :object_id)) do df
    res = df[1:1, :]
    res[!, :organism] .= join(sort(unique(df.organism)), " ")
    return res
end

# expand interactions so that if the interaction is present in SARS/SARS-CoV-2, we also get the
# p-values for the corresponding interaction of the other virus
object_contrasts_df = network_rdata["object_contrasts_slim.df"]
apms_iactions_ex_df = innerjoin(unique!(select!(filter(r -> occursin("experiment", String(r.type)), iactions_df),
                                        [:bait_id, :dest_object_id])),
                                rename!(filter(r -> (r.std_type == "replicate") && occursin("_vs_others", String(r.contrast)),
                                        object_contrasts_df), :object_id => :dest_object_id),
                                on=[:bait_id, :dest_object_id])
unique!(rename!(select(bait_nodes_df, [:bait_full_id, :object_id]), :object_id=>:src_object_id))
apms_iactions_ex_df = leftjoin(apms_iactions_ex_df,
    unique!(rename!(select(bait_nodes_df, [:bait_full_id, :object_id]), :object_id=>:src_object_id)),
    on=:bait_full_id)
iaction_key_cols = [:contrast, :bait_id, :bait_full_id, :src_object_id, :dest_object_id]
select!(apms_iactions_ex_df, union(setdiff(propertynames(apms_iactions_ex_df), propertynames(iactions_df)), iaction_key_cols))
intersect(propertynames(apms_iactions_ex_df), propertynames(iactions_df))
countmap(apms_iactions_ex_df.bait_full_id)

iactions_ex_df = vcat(
    leftjoin(iactions_df, apms_iactions_ex_df, on=iaction_key_cols),
    filter(r -> ismissing(r.type),
           rightjoin(iactions_df, apms_iactions_ex_df, on=iaction_key_cols))
)
iactions_ex_df[!, :type] .= coalesce.(iactions_ex_df.type, "experiment")
countmap(iactions_ex_df.type)

nrow(unique(iactions_ex_df[!, [:src_object_id, :dest_object_id]]))

iactions_baitmerge_pre_df = filter!(r -> r.src_object_id != r.dest_object_id,
        rename!(join(rename!(filter(r -> r.type != "homology", iactions_ex_df),
                         :src_object_id => :old_src_object_id),
                select(filter(r -> r.object_id != SkipBaitId, objects_baitmerge_map_df),
                        [:old_object_id, :object_id, :organism]),
                 on=:old_src_object_id => :old_object_id),
                 :object_id => :src_object_id))
countmap(iactions_baitmerge_pre_df.type)

# remove unattached nodes
filter!(r -> r.object_id ∈ union(Set(filter(r -> r.type == "experiment", iactions_baitmerge_pre_df).src_object_id),
                                 Set(filter(r -> r.type == "experiment" && coalesce(r.is_hit, false), iactions_baitmerge_pre_df).dest_object_id)),
        objects_baitmerge_df)

# merge interactions of SARS/SARS-CoV-2
iactions_baitmerge_df = combine(groupby(filter(r -> r.src_object_id ∈ Set(objects_baitmerge_df.object_id) &&
                                                    r.dest_object_id ∈ Set(objects_baitmerge_df.object_id),
                                               iactions_baitmerge_pre_df),
                                        [:src_object_id, :dest_object_id])) do df
    cov1ix = findfirst(in(["SARS-CoV", "SARS-CoV-GZ02"]), df.organism)
    cov2ix = findfirst(==("SARS-CoV-2"), df.organism)
    min_pval = findmin(coalesce.(df.prob_nonpos, 2.0))
    min_pval_oeprot = findmin(coalesce.(df.oeproteome_p_value, 2.0))
    merged_objix = findfirst(!=(df.src_object_id[1]), df.old_src_object_id)

    DataFrame(merged_src_object_id = isnothing(merged_objix) ? missing : df.old_src_object_id[merged_objix],
              p_value = min_pval[1] == 2.0 ? missing : min_pval[1],
              median_log2 = min_pval[1] == 2.0 ? missing : df.median_log2[min_pval[2]],
              is_hit = any(x -> coalesce(x, false), df.is_hit),
              p_value_SARS_CoV2 = isnothing(cov2ix) ? missing : df.prob_nonpos[cov2ix],
              p_value_SARS_CoV = isnothing(cov1ix) ? missing : df.prob_nonpos[cov1ix],
              median_log2_SARS_CoV2 = isnothing(cov2ix) ? missing : df.median_log2[cov2ix],
              median_log2_SARS_CoV = isnothing(cov1ix) ? missing : df.median_log2[cov1ix],
              is_hit_SARS_CoV2 = isnothing(cov2ix) ? missing : df.is_hit[cov2ix],
              is_hit_SARS_CoV = isnothing(cov1ix) ? missing : df.is_hit[cov1ix],
              type = df.type[1], # FIXME split and recombine?
              ppi_type = df.ppi_type[1],
              iaction_ids = df.iaction_ids[1],
              edge_weight = maximum(w -> coalesce(w, 0.0), df.edge_weight),
              contrast_comparison = isnothing(cov1ix) || isnothing(cov2ix) ? missing :
                df.bait_full_id[cov2ix] * "_vs_" * df.bait_full_id[cov1ix] * "_corrected",
              oeproteome_p_value = min_pval_oeprot[1] == 2.0 ? missing : min_pval_oeprot[1],
              oeproteome_median_log2 = min_pval_oeprot[1] == 2.0 ? missing : df.oeproteome_median_log2[min_pval_oeprot[2]],
              oeproteome_is_hit = min_pval_oeprot[1] == 2.0 ? missing : df.oeproteome_is_hit[min_pval_oeprot[2]],
    )
end
filter!(r -> r.is_hit || r.type != "experiment", iactions_baitmerge_df)

nrow(unique(iactions_baitmerge_pre_df[!, [:src_object_id, :dest_object_id]]))

comparison_cols = [:contrast, :p_value, :median_log2, :is_signif, :is_hit, :is_hit_nomschecks, :change]
iactions_baitmerge_df = leftjoin(iactions_baitmerge_df,
    rename!(select!(filter(r -> r.std_type == "replicate", object_contrasts_df),
                    :object_id, comparison_cols...),
            :object_id=>:dest_object_id,
            [col => Symbol(String(col), "_comparison") for col in comparison_cols]...),
    on=[:contrast_comparison, :dest_object_id])

using CSV
gene_groups_df = CSV.read(joinpath(analysis_path, "networks", "cov2_mq_apms_20200525_hit_oesc_20200604_vigi_20200606_annotations.txt"), delim='\t')
for col in propertynames(gene_groups_df) # FIXME workaround for CSV.Column missing methods
    gene_groups_df[!, col] = convert(Vector, gene_groups_df[!, col])
end
gene_groups_df.annotation_label = replace.(gene_groups_df.annotation_label, "\\n" => "\n")
filter!(r -> r.in_network > 0, gene_groups_df)
gene_groups_df[!, :group_object_id] .= -1000 .- (1:nrow(gene_groups_df))

object2gene_df = semijoin(join(rename!(DelimDataUtils.expand_delim_column(objects_df, list_col=:gene_names, elem_col=:gene_name, key_col=:object_id, delim=";"),
                              :object_id=>:old_object_id),
                      select(objects_baitmerge_map_df, [:object_id, :old_object_id]), on=:old_object_id),
                      objects_baitmerge_df, on=:object_id)
unique!(select!(object2gene_df, [:object_id, :gene_name]))
gene_groups_expanded_df = DelimDataUtils.expand_delim_column(gene_groups_df, list_col=:annotation_genes, elem_col=:term_gene, key_col=:group_object_id, delim=" ")
object_groups_expanded_df = combine(groupby(join(gene_groups_expanded_df, object2gene_df, on=[:term_gene => :gene_name]), [:group_object_id, :object_id])) do df
    df[1:1, :]
end
used_groups_df = semijoin(gene_groups_df, object_groups_expanded_df, on=:group_object_id)

gene_hits_df = CSV.read(joinpath(networks_path, "mq_apms_20200525_20200525_prg_review_hits_vigi_20200607.txt"), delim='\t')
object_hits_df = leftjoin(gene_hits_df, object2gene_df, on=:gene_name)
filter(r -> ismissing(r.object_id), object_hits_df)

objects_baitmerge_df.show_label = objects_baitmerge_df.object_id .∈ Ref(skipmissing(object_hits_df.object_id))

countmap(objects_baitmerge_df.organism)
countmap(iactions_baitmerge_df.type)

iactions_group_matrix_df = join(rename(select(object_groups_expanded_df, :object_id, :group_object_id), :object_id=>:src_object_id),
                                rename(select(object_groups_expanded_df, :object_id, :group_object_id), :object_id=>:dest_object_id),
                                on=:group_object_id)
filter!(r -> r.src_object_id < r.dest_object_id, iactions_group_matrix_df)
iactions_group_matrix_df[!, :edge_weight] .= 8

iactions_baitmerge_layout_df = vcat(select(iactions_baitmerge_df, [:src_object_id, :dest_object_id, :edge_weight]),
                                    select(iactions_group_matrix_df, [:src_object_id, :dest_object_id, :edge_weight]))

bm_nodesize_scale = 0.02
bm_node_type_props = Dict("bait" => (mass=2.0, width=6bm_nodesize_scale, height=6bm_nodesize_scale),
                          "prey" => (mass=2.0, width=1.0bm_nodesize_scale, height=1.0bm_nodesize_scale),
                          "prey hilight" => (mass=0.5, width=2.5bm_nodesize_scale, height=2.5bm_nodesize_scale))
select!(objects_baitmerge_df, Not([:layout_x, :layout_y]))
objects_baitmerge_df.width = getproperty.(getindex.(Ref(bm_node_type_props), objects_baitmerge_df.exp_role), :width);
objects_baitmerge_df.height = getproperty.(getindex.(Ref(bm_node_type_props), objects_baitmerge_df.exp_role), :height);

group_objects_df = repeat(objects_baitmerge_df[1:1, :], nrow(used_groups_df))
for col in propertynames(group_objects_df)
    group_objects_df[!, col] = missings(eltype(group_objects_df[!, col]), nrow(group_objects_df))
end
group_objects_df.object_id = copy(used_groups_df.group_object_id)
group_objects_df.object_label = copy(used_groups_df.annotation_label)
group_objects_df[!, :exp_role] .= "group"
group_objects_df[!, :is_contaminant] .= false
group_objects_df[!, :is_reverse] .= false
group_objects_df[!, :is_viral] .= false
group_objects_df[!, :show_label] .= true

# the graph to define socioaffinity
gr_bm_apms = SimpleGraph(FA.Graph(iactions_baitmerge_layout_df,
                         objects_baitmerge_df,
                         src_node_col=:src_object_id, dest_node_col=:dest_object_id, weight_col=:edge_weight,
                         node_col=:object_id))

node_dislikes = FA.socioaffinity(gr_bm_apms, p=(0.0, 1.0), q=1.0)
node_dislikes_baits = FA.socioaffinity(gr_bm_apms, p=(1.25, 1.25), q=2.0)
for (i, r1) in enumerate(eachrow(objects_baitmerge_df)), (j, r2) in enumerate(eachrow(objects_baitmerge_df))
    r1_label = replace(r1.gene_label, r"\?$" => "")
    r2_label = replace(r2.gene_label, r"\?$" => "")
    if r1.exp_role == "bait" && r2.exp_role == "bait" &&
       (r1_label != r2_label)
       orf3_cluster =
        (r1_label ∈ ["ORF3", "ORF3a", "ORF3b", "ORF4", "ORF4a"] &&
         r2_label ∈ ["ORF3", "ORF3a", "ORF3b", "ORF4", "ORF4a"])
       node_dislikes[i, j] = ifelse(orf3_cluster, 5, 15) * node_dislikes_baits[i, j]
    end
end

bmgr = FA.Graph(iactions_baitmerge_layout_df, objects_baitmerge_df,
              src_node_col=:src_object_id, dest_node_col=:dest_object_id, weight_col=:edge_weight,
              #gravitable_col=:gravitable,
              node_col=:object_id, mass_col=:mass,
              width_col=:width, height_col=:height)
FA.layout!(bmgr, FA.ForceAtlas3Settings(bmgr,
            outboundAttractionDistribution=false,
            attractionStrength=10.0, attractionEdgeWeightInfluence=0.5, jitterTolerance=0.1,
            repulsionStrength=1.0*(1.0 .+ node_dislikes),
            repulsionNodeModel=:Point,
            gravity=1.0, gravityFalloff=1.3, gravityShape=:Point,
            #gravityRodCorners=((0.0, -4.0), (0.0, 4.0)), gravityRodCenterWeight=0.1
            ),
            nsteps=1000, progressbar=true)
FA.layout!(bmgr, FA.ForceAtlas3Settings(bmgr,
            outboundAttractionDistribution=false,
            attractionStrength=5.0, attractionEdgeWeightInfluence=0.75, jitterTolerance=0.1,
            repulsionStrength=2 .* (0.5 .+ node_dislikes),
            repulsionNodeModel=:Circle,
            gravity=1.2, gravityFalloff=1.2, gravityShape=:Point,
            #gravityRodCorners=((0.0, -6.0), (0.0, 6.0)), gravityRodCenterWeight=0.1
            ),
            nsteps=5000, progressbar=true)

objects_baitmerge_df.layout_x = FA.extract_layout(bmgr)[1] .* 20
objects_baitmerge_df.layout_y = FA.extract_layout(bmgr)[2] .* 20
group_objects_df[!, :layout_x] .= NaN
group_objects_df[!, :layout_y] .= NaN
group_objects_df[!, :group_object_id] = missings(Int, nrow(group_objects_df))

# assign object to a single group
objects_baitmerge_and_groups_df = leftjoin(objects_baitmerge_df,
                                           combine(df -> df[1:1, :], groupby(select(object_groups_expanded_df, [:object_id, :group_object_id]), :object_id)),
                                           on=:object_id)
objects_baitmerge_and_groups_df = vcat(objects_baitmerge_and_groups_df,
                                       filter(r -> r.object_id ∈ Set(objects_baitmerge_and_groups_df.group_object_id), group_objects_df))
objects_baitmerge_and_groups_df.exp_role = objects_baitmerge_and_groups_df.exp_role .* ifelse.(objects_baitmerge_and_groups_df.show_label, " hilight", "")

objects_baitmerge_and_groups_df.object_idstr = string.(objects_baitmerge_and_groups_df.object_id)
objects_baitmerge_and_groups_df.group_object_idstr = ifelse.(ismissing.(objects_baitmerge_and_groups_df.group_object_id), missing,
                                                             string.(objects_baitmerge_and_groups_df.group_object_id))
iactions_baitmerge_df.src_object_idstr = string.(iactions_baitmerge_df.src_object_id)
iactions_baitmerge_df.dest_object_idstr = string.(iactions_baitmerge_df.dest_object_id)

bm_apms_graph = GraphML.import_graph(objects_baitmerge_and_groups_df, combine(groupby(iactions_baitmerge_df, [:src_object_idstr, :dest_object_idstr])) do df
    df[1:1, :]
end,
                             node_col=:object_idstr,
                             source_col=:src_object_idstr, target_col=:dest_object_idstr, parent_col=:group_object_idstr,
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
                                           :oeproteome_is_hit => "OE Proteome Hit",
                                           :oeproteome_bait_full_ids => "OE Proteome Baits",
                                           :oeproteome_p_value => "OE Proteome: Most Signif. P-value",
                                           :oeproteome_median_log2 => "OE Proteome: Most Signif. Median Log2",
                                           :cov2ts_proteome_timepoints => "CoV-2 Proteome: Timepoints of Significant Changes",
                                           :cov2ts_proteome_p_value => "CoV-2 Proteome: Most Signif. P-value",
                                           :cov2ts_proteome_median_log2 => "CoV-2 Proteome: Most Signif. Median Log2",
                                           :cov2ts_proteome_is_hit => "CoV-2 Proteome: Is Hit",
                                           :cov2ts_phospho_ptms => "CoV-2 Phospho: PTMs and Timepoints of Significant Changes",
                                           :cov2ts_phospho_p_value => "CoV-2 Phospho: Most Signif. P-value",
                                           :cov2ts_phospho_median_log2 => "CoV-2 Phospho: Most Signif. Median Log2",
                                           :cov2ts_phospho_is_hit => "CoV-2 Phospho: Is Hit",
                                           :layout_x, :layout_y],
                             edge_attrs = [Symbol("edge_weight") => "Weight",
                                           :type => "type",
                                           :oeproteome_is_hit => "OE Proteome Hit",
                                           :oeproteome_p_value => "OE Proteome P-value",
                                           :oeproteome_median_log2 => "OE Proteome Median Log2",

                                           :p_value => "P-value vs background (most significant)",
                                           :median_log2 => "Log2(enrichment vs background) (most significant)",
                                           :is_hit => "Is Hit (vs background)",

                                           :p_value_SARS_CoV2 => "P-value vs background SARS-CoV-2",
                                           :median_log2_SARS_CoV2 => "Log2(enrichment vs background) SARS-CoV-2",
                                           :is_hit_SARS_CoV2 => "Is Hit (vs background) SARS-CoV-2",

                                           :p_value_SARS_CoV => "P-value vs background SARS-CoV",
                                           :median_log2_SARS_CoV => "Log2(enrichment vs background) SARS-CoV",
                                           :is_hit_SARS_CoV => "Is Hit (vs background) SARS-CoV",

                                           :p_value_comparison => "P-value of CoV-2 vs SARS",
                                           :median_log2_comparison => "Log2(enrichment of CoV-2 vs SARS)",
                                           :change_comparison => "Is Hit CoV-2 vs SARS",

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
open(joinpath(networks_path, "$(proj_info.id)_baitsmerged_graph_$(proj_info.msfolder)_$(proj_info.fit_ver)_FA3.graphml"), "w") do io
    write(io, bm_apms_graph)
end

# final stats

baits_merged_df = leftjoin(select(semijoin(objects_baitmerge_df,
         filter(r -> occursin("experiment", r.type), iactions_baitmerge_df),
         on=:object_id=>:src_object_id), [:object_id, :object_label, :organism]),
         rename!(select(bait_nodes_df, [:object_id, :bait_full_id, :organism]),
                 :organism => :bait_organism),
         on=:object_id)
baits_merged_df.is_sars = occursin.("SARS_CoV_", baits_merged_df.bait_full_id)
baits_merged_df.is_cov2 = occursin.("SARS_CoV2_", baits_merged_df.bait_full_id)
baits_merged_df.is_hcov = occursin.("HCoV_", baits_merged_df.bait_full_id)

# baits that require special counting
special_baits_df = filter(r -> occursin(r"^SARS_CoV2?_ORF[38][ab]?$", r.bait_full_id), baits_merged_df)
# SARS/CoV-2 interactions

cov_iactions_df = filter(r -> coalesce(r.is_hit, false) && occursin("experiment", r.type),
       innerjoin(iactions_baitmerge_df, select(filter(r -> r.is_sars || r.is_cov2 || r.is_hcov, baits_merged_df), [:object_id, :bait_full_id, :is_sars, :is_cov2]),
                 on=:src_object_id=>:object_id))
length(unique(cov_iactions_df.dest_object_id))

sars_iactions_df = filter(r -> coalesce(r.is_hit, false) && occursin("experiment", r.type),
       innerjoin(iactions_baitmerge_df, select(filter(r -> r.is_sars || r.is_cov2, baits_merged_df), [:object_id, :bait_full_id, :is_sars, :is_cov2]),
                 on=:src_object_id=>:object_id))
length(unique(sars_iactions_df.dest_object_id))
length(unique(filter(r -> coalesce(r.is_hit, false) && occursin("experiment", r.type), iactions_baitmerge_df).dest_object_id))
countmap(sars_iactions_df.change_comparison)
countmap(coalesce.(sars_iactions_df.change_comparison,
    ifelse.(sars_iactions_df.is_sars, "-",
        ifelse.(sars_iactions_df.is_cov2, "+", "?"))))

hcov_iactions_df = filter(r -> coalesce(r.is_hit, false) && occursin("experiment", r.type),
       innerjoin(iactions_baitmerge_df, select(filter(r -> r.is_hcov, baits_merged_df), [:object_id, :bait_full_id, :is_sars, :is_cov2]),
                 on=:src_object_id=>:object_id))
length(unique(hcov_iactions_df.dest_object_id))

# special cases of non-merged baits
special_baits_df.bait_id = replace.(replace.(special_baits_df.bait_full_id, Ref(r"[ab]$" => "")),
                                    Ref(r"SARS_CoV2?_" => ""))
special_bait_contrasts_df = innerjoin(rename!(semijoin(filter(r -> r.contrast_type == "comparison" &&
                        r.std_type == "replicate" &&
                        occursin(r"^SARS_.+_vs_SARS.+_corrected", r.contrast),
                                            object_contrasts_df),
                                     special_baits_df, on=[:bait_id]), :object_id=>:old_object_id),
                                     objects_baitmerge_map_df, on=:old_object_id)
special_bait_contrasts_df.bait_full_id_rhs = replace.(replace.(special_bait_contrasts_df.contrast, Ref(r".+_vs_" => "")),
                                                      Ref("_corrected" => ""))
countmap(special_bait_contrasts_df.change)

special_bait_iactions_lhs_df = innerjoin(iactions_baitmerge_df, filter(r -> r.bait_organism == "SARS-CoV-2", special_baits_df),
                                     on=:src_object_id => :object_id)
rename!(special_bait_iactions_lhs_df, :is_hit => :is_hit_lhs)
special_bait_iactions_rhs_df = innerjoin(iactions_baitmerge_df, filter(r -> r.bait_organism != "SARS-CoV-2", special_baits_df),
                                        on=:src_object_id => :object_id)
rename!(special_bait_iactions_rhs_df, :bait_full_id => :bait_full_id_rhs, :is_hit => :is_hit_rhs)

special_bait_contrasts_merged_df = combine(groupby(
    leftjoin(leftjoin(special_bait_contrasts_df,
                      select(special_bait_iactions_lhs_df, [:bait_id, :bait_full_id, :dest_object_id, :is_hit_lhs]),
            on = [:bait_id, :bait_full_id, :object_id=>:dest_object_id]),
            select(special_bait_iactions_rhs_df, [:bait_id, :bait_full_id_rhs, :dest_object_id, :is_hit_rhs]),
            on = [:bait_id, :bait_full_id_rhs, :object_id=>:dest_object_id]),
    [:bait_id, :object_id])) do df
    DataFrame(change_comparison = ifelse(any(==("."), df.change[coalesce.(df.is_hit_lhs, false) .| coalesce.(df.is_hit_rhs, false)]),
                                         ".",
                                         ifelse(any(==("+"), df.change[coalesce.(df.is_hit_lhs, false)]), "+", "-")),
              n_iactions = nrow(df),
              contrasts = join(sort(df.contrast), " "),
              is_hit = any(x -> coalesce(x, false), df.is_hit_lhs) || any(x -> coalesce(x, false), df.is_hit_rhs))
end
filter!(r -> r.is_hit, special_bait_contrasts_merged_df)
countmap(special_bait_contrasts_merged_df.change_comparison)
countmap(special_bait_contrasts_merged_df.n_iactions)
countmap(special_bait_contrasts_merged_df.contrasts)

countmap(vcat(antijoin(sars_iactions_df, special_baits_df, on=:bait_full_id).change_comparison,
              special_bait_contrasts_merged_df.change_comparison))

countmap(special_bait_iactions_df.bait_id)
