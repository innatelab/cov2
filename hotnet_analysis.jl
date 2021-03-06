proj_info = (id = "cov2",
             apms_folder = "mq_apms_20200525",
             apms_data_ver = "20200525",
             apms_fit_ver = "20200525",
             oeproteome_folder = "spectronaut_oeproteome_20200527",
             oeproteome_data_ver = "20200527",
             oeproteome_fit_ver = "20200608",
             hotnet_ver = "20201022")
using Pkg
Pkg.activate(joinpath(base_scripts_path, "adhoc", proj_info.id))

using Revise
using Distances, DataFrames, CSV, RData, CodecZlib, JLD2

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

includet(joinpath(misc_scripts_path, "frame_utils.jl"))
includet(joinpath(misc_scripts_path, "delimdata_utils.jl"))
includet(joinpath(misc_scripts_path, "graphml_writer.jl"))
includet(joinpath(misc_scripts_path, "forceatlas3_layout.jl"))
includet(joinpath(misc_scripts_path, "hotnet_utils.jl"))

using Revise, DataFrames, CSV, RData, SimpleWeightedGraphs, LightGraphs,
      Statistics, StatsBase, LinearAlgebra
using HierarchicalHotNet
HHN = HierarchicalHotNet

reactomefi_diedges_df = HotnetUtils.import_reactomefi(joinpath(party3rd_data_path, "FIsInGene_020720_with_annotations.txt"), verbose=true)
reactomefi_genes = levels(reactomefi_diedges_df.gene1)
reactomefi_gene2index = Dict(gene => ix for (ix, gene) in enumerate(reactomefi_genes))

# reverse-direction graph
reactomefi_digraph_revfull, reactomefi_digraph_revfull_vertices =
    HHN.import_digraph(reactomefi_diedges_df, src_col=:gene2, dest_col=:gene1, weight_col=:score)
@assert levels(reactomefi_diedges_df.gene1) == reactomefi_digraph_revfull_vertices
# find large-enough connected (not necessarily strongly) components
reactomefi_digraph_revfull_conncomp = HHN.strongly_connected_components(LightGraphs.weights(reactomefi_digraph_revfull), HHN.EdgeTest{Float64}(threshold=nothing))
reactomefi_digraph_revfull_conncomp_used = filter(comp -> length(comp) > 5, reactomefi_digraph_revfull_conncomp)

reactomefi_genes_mask = eachindex(reactomefi_digraph_revfull_vertices) .∈ Ref(reactomefi_digraph_revfull_conncomp_used.elems)
vertex2gene = findall(reactomefi_genes_mask)
gene2vertex = fill(0, length(reactomefi_genes))
@inbounds for (vertex, gene) in enumerate(vertex2gene)
    gene2vertex[gene] = vertex
end

reactomefi_diedges_df[!, :vertex1] = gene2vertex[CategoricalArrays.level.(reactomefi_diedges_df.gene1)]
reactomefi_diedges_df[!, :vertex2] = gene2vertex[CategoricalArrays.level.(reactomefi_diedges_df.gene2)]
reactomefi_diedges_used_df = filter(r -> r.is_valid && (r.vertex1 > 0) && (r.vertex2 > 0),
                                    reactomefi_diedges_df)
reactomefi_digraph_rev, reactomefi_digraph_vertex_indices =
    HHN.import_digraph(reactomefi_diedges_used_df,
                       src_col=:vertex2, dest_col=:vertex1, weight_col=:score)
@assert reactomefi_digraph_vertex_indices == 1:nv(reactomefi_digraph_rev)
reactomefi_digraph_rev

# load experiments data
apms_fit_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(proj_info.apms_folder)_$(proj_info.apms_fit_ver).RData"))
apms_data_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.apms_folder)_$(proj_info.apms_fit_ver).RData"))
apms_network_rdata = load(joinpath(analysis_path, "networks", "$(proj_info.id)_4graph_$(proj_info.apms_folder)_$(proj_info.apms_fit_ver).RData"))
oeproteome_data_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.oeproteome_folder)_$(proj_info.oeproteome_data_ver).RData"))
oeproteome_fit_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(proj_info.oeproteome_folder)_$(proj_info.oeproteome_fit_ver).RData"))

#--- map OE proteome to reactome network
objects_df = copy(oeproteome_data_rdata["msdata_full"]["protgroups"])
objects_df.is_used = objects_df.q_value .<= 0.001
sel_std_type = "median"
oeproteome_contrasts_df = filter!(r -> r.std_type == sel_std_type && occursin(r"_vs_controls$", r.contrast),
                                  semijoin(oeproteome_fit_rdata["object_contrasts.df"],
                                          filter(r -> r.is_used, objects_df), on=:object_id=>:protgroup_id))
oeproteome_batchcontrasts_df = filter!(r -> r.std_type == sel_std_type &&
                                       occursin(r"_vs_B\d+_others$", r.contrast),
                                       semijoin(oeproteome_fit_rdata["object_contrasts.df"],
                                                filter(r -> r.is_used, objects_df), on=:object_id=>:protgroup_id))
oeproteome_contrasts_df = innerjoin(oeproteome_contrasts_df,
    rename!(select(oeproteome_batchcontrasts_df, [:bait_full_id, :contrast, :std_type, :object_id, :median_log2, :p_value]),
            :contrast => :contrast_batch, :median_log2 => :median_log2_batch, :p_value => :p_value_batch),
            on = [:bait_full_id, :std_type, :object_id])

countmap(oeproteome_contrasts_df.bait_full_id)

# calculate vertex weights
oeproteome_contrasts_df.is_source = [
    (coalesce(r.p_value, 1.0) <= 0.001 && abs(r.median_log2) >= 0.25) &&
    (coalesce(r.p_value_batch, 1.0) <= 0.001 && abs(r.median_log2_batch) >= 0.25)
    for r in eachrow(oeproteome_contrasts_df)]
oeproteome_contrasts_df.vertex_weight = [ifelse(
    (coalesce(r.p_value, 1.0) <= 0.01 && abs(r.median_log2) >= 0.25) &&
    (coalesce(r.p_value_batch, 1.0) <= 0.01 && abs(r.median_log2_batch) >= 0.25),
    (-log10(max(1E-20, r.p_value)))^0.5 * abs(r.median_log2)^0.5, 0.0)
    for r in eachrow(oeproteome_contrasts_df)]
extrema(oeproteome_contrasts_df[oeproteome_contrasts_df.vertex_weight .> 0, :vertex_weight])
quantile(oeproteome_contrasts_df[oeproteome_contrasts_df.vertex_weight .> 0, :vertex_weight], 0.5)
oeproteome_contrasts_stats_df = combine(groupby(oeproteome_contrasts_df, [:bait_full_id, :contrast]),
        :vertex_weight => (w -> quantile(filter(>(0), w), 0.5)) => :vertex_weight_median,
        :vertex_weight => (w -> sum(>(0), w)) => :n_nzweights,
        :is_source => (w -> sum(w)) => :n_sources)

flows_path = joinpath(analysis_path, "networks", "apms_oeproteome_flows_$(proj_info.hotnet_ver)")
isdir(flows_path) || mkdir(flows_path)
CSV.write(joinpath(flows_path, "oeproteome_contrast_stats_$(proj_info.hotnet_ver).txt"), oeproteome_contrasts_stats_df, delim='\t')

oeproteome_obj2gene_df = innerjoin(
    filter(r -> r.is_majority, oeproteome_data_rdata["msdata_full"]["protein2protgroup"]),
    filter!(r -> !ismissing(r.gene_name), select(oeproteome_data_rdata["msdata_full"]["proteins"], Not(:protgroup_id))),
    on=:protein_ac)
rename!(oeproteome_obj2gene_df, :protgroup_id => :object_id)
oeproteome_obj2gene_df[!, :gene_id] = get.(Ref(reactomefi_gene2index), oeproteome_obj2gene_df.gene_name, 0)
oeproteome_gene_effects_df = combine(groupby(innerjoin(select(oeproteome_contrasts_df, Not([:is_viral, :is_contaminant])),
                                     oeproteome_obj2gene_df, on=:object_id),
          [:contrast, :gene_name, :gene_id])) do df
    rowix = findmin(df.p_value)[2]
    return df[rowix:rowix, :]
end;
countmap(oeproteome_gene_effects_df.vertex_weight .> 0)

gene_info_df = combine(groupby(filter(r -> !ismissing(r.gene_name),
                      vcat(select(apms_data_rdata["msdata_full"]["proteins"], [:gene_name, :protein_description]),
                           select(oeproteome_data_rdata["msdata_full"]["proteins"], [:gene_name, :protein_description]))),
                  [:gene_name])) do gene_df
                  gene_df[1:1, :]
              end

bait2oeproteome_genes = Dict{String, Vector{Pair{Int, Float64}}}()
bait2oeproteome_vertices = Dict{String, Vector{Pair{Int, Float64}}}()
for bait_genes_df in groupby(filter(r -> r.is_source && r.gene_id > 0,
                                    oeproteome_gene_effects_df), :bait_full_id)
    baitkey = bait_genes_df.bait_full_id[1]
    gene_weights = bait_genes_df.gene_id .=> bait_genes_df.vertex_weight
    bait2oeproteome_genes[baitkey] = gene_weights
    vertex_weights = map(gene_weights) do (g, w)
        v = gene2vertex[g]
        v2 = v > 0 ? searchsortedfirst(reactomefi_digraph_vertex_indices, v) : 0
        @assert v2 == v "v2=$v2 doesn't match v=$v"
        return (0 < v2 <= length(reactomefi_digraph_vertex_indices)) &&
            (reactomefi_digraph_vertex_indices[v2] == v) ? v2 => w : 0 => NaN
    end
    filter!(vw -> vw[1] > 0, vertex_weights)
    bait2oeproteome_vertices[baitkey] = vertex_weights
end

#--- map APMS data to reactome network
apms_iactions_df = filter(r -> coalesce(r.is_hit, false),
                          apms_network_rdata["iactions_4table.df"])
countmap(apms_iactions_df.bait_full_id)
apms_obj2gene_df = innerjoin(
    filter(r -> r.is_majority, apms_data_rdata["msdata_full"]["protein2protregroup"]),
    filter!(r -> !ismissing(r.gene_name), select(apms_data_rdata["msdata_full"]["proteins"], Not(:protgroup_id))),
    on=:protein_ac)
rename!(apms_obj2gene_df, :protregroup_id => :object_id)
apms_obj2gene_df[!, :gene_id] = get.(Ref(reactomefi_gene2index), apms_obj2gene_df.gene_name, 0)
apms_gene_iactions_df = combine(groupby(innerjoin(select(apms_iactions_df, Not([:is_viral, :is_contaminant, :protein_description])),
                                        apms_obj2gene_df, on=:object_id),
          [:contrast, :gene_name, :gene_id])) do df
    rowix = findmin(df.p_value)[2]
    return df[rowix:rowix, :]
end;

bait2apms_genes = Dict{String, Vector{Int}}()
bait2apms_vertices = Dict{String, Vector{Int}}()
for bait_genes_df in groupby(apms_gene_iactions_df, :bait_full_id)
    baitkey = bait_genes_df.bait_full_id[1]
    gene_ids = sort!(unique!(collect(filter(>(0), bait_genes_df.gene_id))))
    bait2apms_genes[baitkey] = gene_ids
    vertices = unique!(filter!(!=(0), gene2vertex[gene_ids]))
    v2s = map(vertices) do v
        v2 = searchsortedfirst(reactomefi_digraph_vertex_indices, v)
        @assert v2 == v
        return (v2 <= length(reactomefi_digraph_vertex_indices)) &&
            (reactomefi_digraph_vertex_indices[v2] == v) ? v2 : 0
    end
    filter!(>(0), v2s)
    bait2apms_vertices[baitkey] = v2s
end

bait_ids = sort!(collect(keys(bait2oeproteome_vertices)))
bait2vertex_weights = Dict{String, Vector{Float64}}()
for bait_id in bait_ids
    vertex2weight = bait2oeproteome_vertices[bait_id]
    vertex_weights = fill(0.0, nv(reactomefi_digraph_rev))
    vertex_weights[first.(vertex2weight)] .= last.(vertex2weight)
    bait2vertex_weights[bait_id] = vertex_weights
end

randomwalk_params = (restart_prob = 0.4,
                     inedge_weight_min = 1.0,
                     source_weight_min = 1.5,
                     flow_weight_max = 0.005)
reactomefi_adjmtx = Matrix(LightGraphs.weights(reactomefi_digraph_rev));
reactomefi_walkmtx = HHN.random_walk_matrix(reactomefi_digraph_rev, randomwalk_params.restart_prob)

### select optimal restart probability
test_restart_probs = 0.05:0.05:0.95
test_diedge_weight_min = 1.0
reactomefi_test_walkmtxs = [HHN.similarity_matrix(
    HHN.stepmatrix(reactomefi_adjmtx,
                   inedge_weights=bait2vertex_weights["SARS_CoV2_M"] .+ randomwalk_params.inedge_weight_min),
    bait2vertex_weights["SARS_CoV2_M"], restart_probability=p) for p in test_restart_probs]

# calculate average probability to stay in the original direct neighborhood of a vertex
# or travel further
vtx_totprob = vec.(sum.(reactomefi_test_walkmtxs, dims=1))
vtx_totprob_sum = sum.(vtx_totprob)
vtx_noselfprobs = (vec.(sum.(reactomefi_test_walkmtxs, dims=1)) - diag.(reactomefi_test_walkmtxs)) ./ vtx_totprob_sum
vtx_neiprobs = HHN.neighborhood_weights.(reactomefi_test_walkmtxs, Ref(reactomefi_digraph_rev)) ./ vtx_totprob_sum
vtx_extprobs = vtx_noselfprobs .- vtx_neiprobs

restartprob_stats_df = vcat(DataFrame(restart_prob = test_restart_probs,
                                      vertices = "neighbors",
                                      trans_probs = sum.(vtx_neiprobs)),
                            DataFrame(restart_prob = test_restart_probs,
                                      vertices = "exterior",
                                      trans_probs = sum.(vtx_extprobs)))
categorical!(restartprob_stats_df)
using PlotlyJS
restartXneighb_plot = plot(restartprob_stats_df, x=:restart_prob, y=:trans_probs, color=:vertices, group=:vertices)
savefig(restartXneighb_plot,
        joinpath(analysis_path, "networks", "$(proj_info.id)_restart_probability_neighborhood_$(proj_info.hotnet_ver).svg"))

###
sel_bait_ids = bait_ids
bait_stepmtxs = [bait_id => HHN.stepmatrix(reactomefi_adjmtx,
                        inedge_weights=bait2vertex_weights[bait_id] .+ randomwalk_params.inedge_weight_min)
                 for bait_id in sel_bait_ids]
bait2stepmtx = Dict(bait_id => stepmtx for (bait_id, stepmtx) in bait_stepmtxs)

bait_walkmtxs = [bait_id => HHN.similarity_matrix(stepmtx, bait2vertex_weights[bait_id],
                                                  restart_probability=randomwalk_params.restart_prob)
                 for (bait_id, stepmtx) in bait_stepmtxs]
bait2walkmtx = Dict(bait_id => walkmtx for (bait_id, walkmtx) in bait_walkmtxs)
sum(>(0), bait2walkmtx["SARS_CoV_NSP3"] |> vec)

reactomefi_trees = Vector{HHN.SCCTree{Float64}}(undef, length(sel_bait_ids));
Threads.@threads for i in eachindex(sel_bait_ids)
    bait_id = sel_bait_ids[i]
    @assert bait_walkmtxs[i][1] == bait_id
    walkmtx = bait_walkmtxs[i][2]
    @info "Bait $bait_id: $(size(walkmtx)) matrix"
    reactomefi_trees[i] = HHN.scctree(walkmtx, verbose=false, method=:bisect)
end
bait2tree = Dict(bait_id => reactomefi_trees[i] for (i, bait_id) in enumerate(sel_bait_ids))

# bin vertices mapped to pg_ids (i.e. quantified) according to their in/out degree
vertex_bins = HHN.vertexbins(reactomefi_digraph_rev, findall(>(-1), vertex2gene),
                             by=:outXin, method=:tree, nbins=100)
bait2perm_vertex_weights = Dict{String, Matrix{Float16}}()
for (bait_id, vertex_weights) in bait2vertex_weights
    vertex_weights_perms = fill(zero(eltype(valtype(bait2perm_vertex_weights))),
                                length(vertex_weights), 1000)
    vertex_perm = collect(eachindex(vertex_weights))
    for wperm in eachcol(vertex_weights_perms)
        HHN.randpermgroups!(vertex_perm, vertex_bins)
        @inbounds for (i, j) in enumerate(vertex_perm)
            wperm[i] = vertex_weights[j]
        end
    end
    bait2perm_vertex_weights[bait_id] = vertex_weights_perms
end

bait2vertex_walkweights = Dict(bait_id => dropdims(sum(walkmtx, dims=2), dims=2)
                               for (bait_id, walkmtx) in bait2walkmtx)

using Serialization, CodecZstd
open(ZstdCompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_prepare_$(proj_info.hotnet_ver).jlser.zst"), "w") do io
    serialize(io, (vertex2gene, gene2vertex,
              reactomefi_genes, reactomefi_digraph_rev,
              bait2apms_vertices, apms_gene_iactions_df,
              bait_ids, bait2vertex_weights,
              randomwalk_params, bait2walkmtx, bait2tree,
              vertex_bins, bait2perm_vertex_weights))
end

# save per-bait for use on a cluster
permtrees_input_prefix = "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver)"
sel_bait_ids = sort(bait_ids)
open(ZstdCompressorStream, joinpath(scratch_path, "$(permtrees_input_prefix).jlser.zst"), "w") do io
    serialize(io, (vertex2gene, reactomefi_genes, reactomefi_digraph_rev,
                   [bait_id => size(bait2perm_vertex_weights[bait_id], 2) for bait_id in sel_bait_ids],
                   randomwalk_params))
end;

isdir(joinpath(scratch_path, permtrees_input_prefix)) || mkdir(joinpath(scratch_path, permtrees_input_prefix))
for (bait_ix, bait_id) in enumerate(sel_bait_ids)
    open(ZstdCompressorStream, joinpath(scratch_path, permtrees_input_prefix, "bait_$(bait_ix)_perms.jlser.zst"), "w") do io
        diedge_ixs = findall(!=(0), vec(bait2walkmtx[bait_id]))
        serialize(io, (bait_id, bait2vertex_weights[bait_id],
                       bait2vertex_walkweights[bait_id],
                       bait2perm_vertex_weights[bait_id],
                       get(bait2apms_vertices, bait_id, Int[]),
                       diedge_ixs,
                       bait2stepmtx[bait_id][diedge_ixs],
                       bait2walkmtx[bait_id][diedge_ixs],
                       bait2tree[bait_id]))
    end
end;

using Serialization, CodecZstd
vertex2gene, gene2vertex,
reactomefi_genes, reactomefi_digraph_rev,
bait2apms_vertices, apms_gene_iactions_df,
bait_ids, bait2vertex_weights,
randomwalk_params, bait2walkmtx, bait2tree,
vertex_bins, _ = open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_prepare_$(proj_info.hotnet_ver).jlser.zst"), "r") do io
    deserialize(io)
end

using LinearAlgebra, HierarchicalHotNet
HHN = HierarchicalHotNet

# collect chunks of network diffusion / tree statistics
using Base.Filesystem, Serialization, CodecZstd

# collect real data-based statistics (see hotnet_treestats_chunk.jl)
chunk_prefix = "$(proj_info.id)_hotnet_treestats_$(proj_info.hotnet_ver)"
tree_stats_dfs = Vector{DataFrame}()
for (root, dirs, files) in Filesystem.walkdir(joinpath(scratch_path, chunk_prefix))
    if isempty(files)
        @warn "Found no files in $root"
        continue # skip the empty folder
    else
        @info "Found $(length(files)) file(s) in $root, processing..."
    end
    for file in files
        (startswith(file, chunk_prefix) && endswith(file, ".jlser.zst")) || continue
        chunk_proj_info, chunk_bait_ids, chunk_tree_stats_df =
            open(ZstdDecompressorStream, joinpath(root, file), "r") do io
            deserialize(io)
        end
        push!(tree_stats_dfs, chunk_tree_stats_df)
    end
end
tree_stats_df = vcat(tree_stats_dfs...)
tree_stats_dfs = nothing
@save(joinpath(scratch_path, "$(proj_info.id)_hotnet_treestats_$(proj_info.hotnet_ver).jld2"),
      tree_stats_df)

# collect assembled permuted data-based statistics (see hotnet_perm_chunk.jl and hotnet_perm_assemble.jl)
chunk_prefix = "$(proj_info.id)_hotnet_perm_assembled_$(proj_info.hotnet_ver)"
perm_tree_stats_dfs = Vector{DataFrame}()
vertex_stats_dfs = similar(perm_tree_stats_dfs)
diedge_stats_dfs = similar(perm_tree_stats_dfs)
chunks_update = Threads.SpinLock()
for (root, dirs, files) in Filesystem.walkdir(joinpath(scratch_path, chunk_prefix))
    if isempty(files)
        @warn "Found no files in $root"
        continue # skip the empty folder
    else
        @info "Found $(length(files)) file(s) in $root, processing..."
    end
    Threads.@threads for i in eachindex(files)
        file = files[i]
        (startswith(file, chunk_prefix) && endswith(file, ".jlser.zst")) || continue
        chunk_tree_stats_df, chunk_vertex_stats_df, chunk_diedge_stats_df =
            open(ZstdDecompressorStream, joinpath(root, file), "r") do io
            deserialize(io)
        end
        lock(chunks_update)
        push!(perm_tree_stats_dfs, chunk_tree_stats_df)
        push!(vertex_stats_dfs, chunk_vertex_stats_df)
        push!(diedge_stats_dfs, chunk_diedge_stats_df)
        unlock(chunks_update)
    end
end
perm_tree_stats_df = vcat(perm_tree_stats_dfs...)
vertex_stats_df = vcat(vertex_stats_dfs...)
diedge_stats_df = vcat(diedge_stats_dfs...)
perm_tree_stats_dfs = nothing
vertex_stats_dfs = nothing
diedge_stats_dfs = nothing
GC.gc()

vertex_stats_df = leftjoin(vertex_stats_df,
                       rename!(select(oeproteome_gene_effects_df, [:bait_full_id, :gene_id, :p_value, :median_log2]),
                               :p_value => :oeproteome_p_value,
                               :median_log2 => :oeproteome_median_log2,
                               :bait_full_id => :bait_id),
                       on=[:bait_id, :gene_id])
vertex_stats_df2 = leftjoin(vertex_stats_df,
                       rename!(select(apms_gene_iactions_df, [:bait_full_id, :gene_name, :p_value, :median_log2]),
                               :p_value => :apms_p_value,
                               :median_log2 => :apms_median_log2,
                               :bait_full_id => :bait_id),
                       on=[:bait_id, :gene_name])
@assert nrow(vertex_stats_df) == nrow(vertex_stats_df2)
vertex_stats_df = vertex_stats_df2
vertex_stats_df.is_hit = vertex_stats_df.prob_perm_walkweight_greater .<= 0.05
vertex_stats_df = leftjoin(vertex_stats_df, gene_info_df, on=:gene_name)

countmap(tree_stats_df.bait_id)
countmap(perm_tree_stats_df.bait_id)
@save(joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_stats_$(proj_info.hotnet_ver).jld2"),
      vertex_stats_df, diedge_stats_df, tree_stats_df, perm_tree_stats_df)

@load(joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_stats_$(proj_info.hotnet_ver).jld2"),
      vertex_stats_df, diedge_stats_df, tree_stats_df, perm_tree_stats_df)

tree_all_stats_df = vcat(tree_stats_df, perm_tree_stats_df)
tree_all_stats_df.is_permuted = vcat(falses(nrow(tree_stats_df)),
                                     trues(nrow(perm_tree_stats_df)))
threshold_range = (0.0, 0.05)
tree_binstats_df = combine(groupby(tree_all_stats_df, :bait_id)) do bait_treestats_df
    @info "binstats($(bait_treestats_df.bait_id[1]))"
    if nrow(bait_treestats_df) == 0
        @warn "cutstats frame is empty, skipping"
        return DataFrame()
    end
    HHN.bin_treecut_stats(bait_treestats_df,
                          by_cols=[:is_permuted, :tree],
                          threshold_range=threshold_range)
end

# annotate tree_stats with threshold_bin indices
tree_stats_df2 = copy(tree_stats_df, copycols=false)
tree_stats_df2.threshold_bin = missings(Int, nrow(tree_stats_df2))
for df in groupby(tree_binstats_df, :bait_id)
    bait_bins = sort!([unique(collect(skipmissing(df.threshold_binmin))); [maximum(skipmissing(df.threshold_binmax))]])
    bait_stats_df = view(tree_stats_df2, tree_stats_df2.bait_id .== df.bait_id[1], :)
    HHN.add_bins!(bait_stats_df, :threshold, bait_bins)
end
tree_stats_df = tree_stats_df2

tree_perm_aggstats_df = combine(groupby(filter(r -> r.is_permuted, tree_binstats_df), :bait_id)) do perm_binstats_df
    @info "aggregate_treecut_binstats($(perm_binstats_df.bait_id[1]))"
    HHN.aggregate_treecut_binstats(perm_binstats_df, by_cols=[:is_permuted, :threshold_bin])
end
tree_perm_aggstats_wide_df = unstack(filter(r -> !ismissing(r.threshold_binmid), tree_perm_aggstats_df),
                                     [:bait_id, :is_permuted, :threshold_bin, :threshold_binmid], :quantile,
                                     intersect(HHN.treecut_metrics, propertynames(tree_perm_aggstats_df)),
                                     namewidecols=(valcol, qtl, sep) -> Symbol(valcol, sep, HHN.quantile_suffix(qtl)))

cut_threshold_range = (1E-4, 1E-2)
tree_extremes_df = HHN.extreme_treecut_stats(
    filter(r -> !r.is_permuted && (cut_threshold_range[1] <= coalesce(r.threshold, -1.0) <= cut_threshold_range[2]), tree_stats_df2),
    tree_perm_aggstats_df,
    extra_join_cols = [:bait_id])

include(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))

sel_quantile = 0.25
bait2cut_threshold = Dict(begin
    bait_id = df[1, :bait_id]
    valtype = df[1, :value_type]
    thresh_df = nothing
    for (metric, stat) in [(:flow_avginvlen, "max"), (:flow_avghopweight, "max"), (:flow_avgweight, "max"), (:maxcomponent_size, "max")]
        metric_df = filter(r -> (r.metric == metric) && (r.value_type == valtype) &&
                                (r.stat == stat) &&
                                (r.quantile == (stat == "max" ? 1 - sel_quantile : sel_quantile)), df)
        (nrow(metric_df) == 0) && continue
        @assert nrow(metric_df) == 1
        delta = metric_df[1, :delta]
        if !ismissing(delta) &&
           (((stat == "max") && (delta > 0.0)) ||
            ((stat == "min") && (delta < 0.0)))
            thresh_df = metric_df
            break
        end
    end
    if thresh_df !== nothing
        @info "$bait_id $valtype: using $(thresh_df[1, :metric])=$(thresh_df[1, :value]) (delta=$(thresh_df[1, :delta])) for threshold=$(thresh_df[1, :threshold])"
        thresh_df
    else
        @warn "No significant difference between real and permuted results for $bait_id"
    end
    (bait_id, valtype) => thresh_df
end for df in groupby(filter(r -> !ismissing(r.quantile), tree_extremes_df), [:bait_id, :value_type]))
#    if (r.type == "max") && !ismissing(r.flow_avgweight_perm_50))
bait_cut_thresholds_df = reduce(vcat, df for df in values(bait2cut_threshold) if !isnothing(df))
bait2treecut = Dict((baitid, valtype) =>
    !isnothing(threshold_df) ? HHN.cut(bait2tree[baitid], threshold_df[1, :threshold], minsize=1) : nothing
    for ((baitid, valtype), threshold_df) in bait2cut_threshold)

using PlotlyJS, PlotlyBase
using Printf: @sprintf

sel_bait_id = "SARS_CoV2_ORF7a"
sel_metric = :flow_avgweight
sel_metric = :maxcomponent_size
sel_metric = :topn_components_sizesum
sel_metric = :components_signif_sizesum_mw
sel_metric = :compflow_distance
PlotlyUtils.permstats_plot(
    filter(r -> !r.is_permuted && (r.bait_id == sel_bait_id) && !ismissing(r[sel_metric]),
           tree_binstats_df),
    filter(r -> r.bait_id == sel_bait_id, tree_perm_aggstats_wide_df),
    filter(r -> r.bait_id == sel_bait_id && r.quantile == (r.stat == "max" ? 1.0 - sel_quantile : sel_quantile) &&
                r.metric == sel_metric && !ismissing(r.value),
           tree_extremes_df),
    threshold_range=(0.00, 0.01),
    yaxis_metric=sel_metric)

includet(joinpath(misc_scripts_path, "graphml_writer.jl"))

function flowgraphml(bait_id::AbstractString, threshold::Number; kwargs...)
    orig_diedges = DataFrames.rename(reactomefi_diedges_used_df,
                                     :vertex1 => :source, :vertex2 => :target,
                                     :score => :weight, :direction => :interaction_type)
    HotnetUtils.fix_reactomefi_iactiontype!(orig_diedges.interaction_type)
    HotnetUtils.fix_reactomefi_iactiontype!(orig_diedges.diedge_type)

    vertices_df = DataFrame(vertex = 1:nv(reactomefi_digraph_rev),
                            gene = ifelse.(vertex2gene .> 0, vertex2gene, missing),
                            gene_name = ifelse.(vertex2gene .> 0, getindex.(Ref(reactomefi_genes), vertex2gene), missing))
    vertices_df = leftjoin(vertices_df, select!(filter(r -> r.bait_id == bait_id, vertex_stats_df),
                                [:vertex, :oeproteome_median_log2, :oeproteome_p_value,
                                 :apms_median_log2, :apms_p_value]), on=:vertex)
    vertices_df.oeproteome_median_log2_signif = ifelse.(coalesce.(vertices_df.oeproteome_p_value, false) .<= 0.01,
                                                        vertices_df.oeproteome_median_log2, missing)
    vertices_df.apms_median_log2_signif = ifelse.(coalesce.(vertices_df.apms_p_value, false) .<= 0.01,
                                                  vertices_df.apms_median_log2, missing)

    return HotnetUtils.flowgraphml(
            bait2tree[bait_id], threshold;
            walkmatrix=bait2walkmtx[bait_id], stepmatrix=bait2stepmtx[bait_id],
            vertices_labels=[v2g > 0 ? reactomefi_genes[v2g] : missing for v2g in vertex2gene],
            vertices_info=vertices_df,
            vertices_weights=bait2vertex_weights[bait_id],
            vertices_stats=@isdefined(vertex_stats_df) ? select!(filter(r -> r.bait_id == bait_id, vertex_stats_df),
                    Not([:gene_name, :oeproteome_median_log2, :oeproteome_p_value, :apms_median_log2, :apms_p_value])) : nothing,
            diedges_stats=@isdefined(diedge_stats_df) ?
                select!(rename!(filter(r -> r.bait_id == bait_id, diedge_stats_df),
                        :src=>:source, :dest=>:target), Not([:weight, :walkweight]))
                : nothing,
            source_threshold=randomwalk_params.source_weight_min,
            sinks=get(bait2apms_vertices, bait_id, valtype(bait2apms_vertices)()),
            flow_metrics = filter(r -> r.bait_id == bait_id, tree_stats_df),
            step_sinks=nothing,
            orig_diedges=orig_diedges,
            extra_node_attrs=[:apms_p_value, :apms_median_log2, :apms_median_log2_signif,
                              :oeproteome_p_value, :oeproteome_median_log2, :oeproteome_median_log2_signif],
            kwargs...)
end

Revise.includet(joinpath(misc_scripts_path, "forceatlas3_layout.jl"))
FA = ForceAtlas3

flows_path = joinpath(analysis_path, "networks", "apms_oeproteome_flows_$(proj_info.hotnet_ver)")
isdir(flows_path) || mkdir(flows_path)
CSV.write(joinpath(flows_path, "apms_oeproteome_flows_$(proj_info.hotnet_ver)_cut_thresholds.txt"), bait_cut_thresholds_df, delim='\t')

sel_bait_ids = unique(first.(collect(keys(bait2cut_threshold))))
sel_bait_ids = filter(r -> r.nflows == 0, innerjoin(unique!(select(bait_cut_thresholds_df, [:bait_id, :value_type, :threshold])),
                                                    tree_stats_df, on=[:bait_id, :threshold])).bait_id |> unique
flowgraph_todo = vec([(bait_id, valtype) for bait_id in sel_bait_ids, valtype in unique(last.(collect(keys(bait2cut_threshold))))])
Threads.@threads for i in 1:length(flowgraph_todo)
    bait_id, valtype = flowgraph_todo[i]
    @info "Flowgraph #$i: bait=$bait_id cut=$valtype"
    cut_threshold_df = get(bait2cut_threshold, (bait_id, valtype), nothing)
    if isnothing(cut_threshold_df)
        @warn "  No threshold found, skipping"
        continue
    end
    cut_threshold = cut_threshold_df.threshold[1]
    flow_graph = flowgraphml(bait_id, cut_threshold,
        component_groups=false,
        flowpaths = :flowattr,
        step_threshold=cut_threshold, maxsteps=2,
        layout_cut_threshold=cut_threshold * 1.25,
        pvalue_mw_max = 0.001, verbose=true)
    isempty(flow_graph.nodes) || HotnetUtils.layout_flowgraph!(flow_graph, scale=80, progressbar=false)
    valtype_flows_path = joinpath(flows_path, valtype)
    isdir(valtype_flows_path) || mkdir(valtype_flows_path)
    open(joinpath(valtype_flows_path, "$(proj_info.id)_$(bait_id)_flow_$(proj_info.hotnet_ver)_$(valtype).graphml"), "w") do io
        write(io, flow_graph)
    end
end
HotnetUtils.layout_flowgraph!(test_flowgraph, scale=80);
open(joinpath(flows_path, "test.graphml"), "w") do io
    write(io, test_flowgraph)
end

include(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))

sel_bait_ids = ["SARS_CoV2_ORF7a"]
metrics_plot_path = joinpath(flows_path, "metrics")
isdir(metrics_plot_path) || mkdir(metrics_plot_path)
for bait_id in sel_bait_ids, (metric, stat) in [(:flow_avginvlen, "max"), #=(:flow_avghopweight, "max"),
                                                (:flow_avgweight, "max"),=# (:maxcomponent_size, "max")]
    @info "Plotting $metric stats for $bait_id"
    realstats_df = filter(r -> !r.is_permuted && (r.bait_id == bait_id) && !ismissing(r[metric]), tree_stats_df)
    aggstats_df = filter(r -> r.bait_id == bait_id, tree_perm_aggstats_wide_df)
    isempty(aggstats_df) && continue
    extremes_df = filter(r -> (r.bait_id == bait_id) && (r.stat == stat) && (coalesce(r.quantile, NaN) == (stat == "max" ? 1.0 - sel_quantile : sel_quantile)) &&
                    (r.metric == metric) && !ismissing(r.value) && (r.value_type == "opt"), tree_extremes_df)
    bait_stats_plot = PlotlyUtils.permstats_plot(
        realstats_df, aggstats_df, extremes_df,
        threshold_range=(max(0.0, cut_threshold_range[1]*0.75), cut_threshold_range[2]*1.25),
        yaxis_metric=metric)
    #bait_stats_plot.plot.layout[:width] = "300px"
    #bait_stats_plot.plot.layout[:height] = "200px"
    #bait_stats_plot.plot.layout[:font_size] = "20"
    try
        savefig(bait_stats_plot.plot,
            #joinpath(metrics_plot_path, "$(proj_info.id)_$(bait_id)_$(metric)_$(proj_info.hotnet_ver).pdf"))
            joinpath(metrics_plot_path, "$(bait_id)_$(metric).pdf"))
    catch e
        if e isa InterruptException
            rethrow(e)
        else
            @warn e
        end
    end
end

sel_bait_id = "SARS_CoV2_ORF7a"
sel_metric_x = :threshold
sel_metric_y = :flow_avginvlen
sel_metric_z = :flow_avghopweight
max_threshold = 0.01
df = filter(r -> r.bait_id == sel_bait_id && r.threshold <= max_threshold, tree_stats_df)
metric3d_plot = plot([
    scatter3d(combine(groupby(sort!(filter(r -> r.bait_id == sel_bait_id && r.threshold <= max_threshold, perm_tree_stats_df), [:tree, :threshold]),
                              [:tree])) do tree_gdf
                tree_df = push!(allowmissing(tree_gdf), tree_gdf[1, :])
                tree_df[nrow(tree_df), :] .= missing
                tree_df[nrow(tree_df), :tree] = tree_df[1, :tree]
                return tree_df
            end, connectgaps=false,
            x=sel_metric_x, y=sel_metric_y, z=sel_metric_z, mode="lines", line_color="gray", opacity=0.25, line_size=0.5),
    scatter3d(filter(r -> r.bait_id == sel_bait_id && r.threshold <= max_threshold, tree_stats_df),
              x=sel_metric_x, y=sel_metric_y, z=sel_metric_z, marker_color="firebrick", marker_size=2.0)],
    Layout(title=attr(text=sel_bait_id),
           scene=attr(xaxis=attr(title=attr(text=string(sel_metric_x))),
                      yaxis=attr(title=attr(text=string(sel_metric_y))),
                      zaxis=attr(title=attr(text=string(sel_metric_z))))))
savehtml(metric3d_plot, joinpath(metrics_plot_path, "metrics3d_$(sel_bait_id)_2.html"))
