# job_info = (project = "cov2",
#            name = "cov2_permtrees_stats",
#            hotnet_ver = "20200327",
#            id = 0,
#            chunk = 2,
#            ntrees_perchunk = 2)
job_info = (project = ARGS[1],
            name = ARGS[2],
            hotnet_ver = ARGS[3],
            id = parse(Int, ARGS[4]),
            chunk = parse(Int, ARGS[5]),
            ntrees_perchunk = 10)

@info "Tree stats for $(job_info.project) (job '$(job_info.name)' id=$(job_info.id))" *
      " chunk #$(job_info.chunk)"

proj_info = (id = job_info.project,
             hotnet_ver = job_info.hotnet_ver)

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(base_scratch_path, proj_info.id)
const plots_path = joinpath(analysis_path, "plots")

using JLD2, DataFrames, Serialization, CodecZstd
using LinearAlgebra, HierarchicalHotNet
HHN = HierarchicalHotNet

@load(joinpath(scratch_path, "$(proj_info.id)_hotnet_prepare_$(proj_info.hotnet_ver).jld2"),
      bait_ids, bait2vertex_weights,
      bait2apms_vertices,
      reactomefi_walkmtx, bait2reactomefi_tree,
      bait2vertex_weights_perms)
@load(joinpath(scratch_path, "$(proj_info.id)_hotnet_permtrees_$(proj_info.hotnet_ver).jld2"),
      bait2reactomefi_perm_trees)
@load(joinpath(scratch_path, "$(proj_info.id)_hotnet_vertex_stats_$(proj_info.hotnet_ver).jld2"),
      vertex_stats_df)

ntrees_pergroup = size.(getindex.(Ref(bait2vertex_weights_perms), bait_ids), 2)
pushfirst!(ntrees_pergroup, length(bait_ids)) # non-permuted trees
nchunks_pergroup = fld1.(ntrees_pergroup, job_info.ntrees_perchunk)
group2lastchunk = cumsum(nchunks_pergroup)
chunkgroup = searchsortedfirst(group2lastchunk, job_info.chunk)
@assert (chunkgroup <= length(group2lastchunk)) "Chunk #$(job_info.chunk) is outside of the valid chunks range"
localchunk = job_info.chunk - (chunkgroup > 1 ? group2lastchunk[chunkgroup-1] : 0)
# indices of chunk trees to process (within bait_id perm_trees vector)
chunk_trees = ((localchunk-1)*job_info.ntrees_perchunk + 1):min(
               localchunk*job_info.ntrees_perchunk,
               ntrees_pergroup[chunkgroup])

treestats_dfs = sizehint!(Vector{DataFrame}(), length(chunk_trees))
W = eltype(reactomefi_walkmtx)
tree_mtx = similar(reactomefi_walkmtx)

if chunkgroup > 1
bait_id = bait_ids[chunkgroup-1]
apms_vertices = bait2apms_vertices[bait_id]
vertex_weights_perms = bait2vertex_weights_perms[bait_id]
empty!(bait2vertex_weights_perms) # release unneeded
permtrees = bait2reactomefi_perm_trees[bait_id]
empty!(bait2reactomefi_perm_trees) # release unneeded
vertex_bait_stats_df = filter(r -> r.bait_full_id == bait_id, vertex_stats_df)
@info "Processing $bait_id permuted trees $(chunk_trees)"
for i in chunk_trees
    @info "Treecut statistics for $bait_id permuted tree #$i"
    vertex_weights = convert(Vector{W}, view(vertex_weights_perms, :, i))
    treestats_df = HHN.treecut_stats(permtrees[i], vertex_bait_stats_df,
                                     mannwhitney_tests=true,
                                     walkmatrix=mul!(tree_mtx, reactomefi_walkmtx, Diagonal(vertex_weights)),
                                     sources=findall(>(0), vertex_weights),
                                     sinks=bait2apms_vertices[bait_id],
                                     nflows_ratio=0.99,
                                     pools=nothing)
    treestats_df[!, :tree] .= i
    treestats_df[!, :bait_full_id] .= bait_id
    treestats_df[!, :is_permuted] .= true
    push!(treestats_dfs, treestats_df)
end
else
@info "Processing non-permuted trees"
empty!(bait2vertex_weights_perms) # release unneeded
empty!(bait2reactomefi_perm_trees) # release unneeded
for i in chunk_trees
    bait_id = bait_ids[i]
    @info "Treecut statistics for non-permuted tree of $bait_id"
    vertex_weights = convert(Vector{W}, bait2vertex_weights[bait_id])
    treestats_df = HHN.treecut_stats(bait2reactomefi_tree[bait_id],
                                     filter(r -> r.bait_full_id == bait_id, vertex_stats_df),
                                     mannwhitney_tests=true,
                                     walkmatrix=mul!(tree_mtx, reactomefi_walkmtx, Diagonal(vertex_weights)),
                                     sources=findall(>(0), vertex_weights),
                                     sinks=bait2apms_vertices[bait_id],
                                     nflows_ratio=0.999,
                                     pools=nothing)
    treestats_df[!, :tree] .= 0
    treestats_df[!, :bait_full_id] .= bait_id
    treestats_df[!, :is_permuted] .= false
    push!(treestats_dfs, treestats_df)
end
bait_id = nothing
end
chunkstats_df = reduce(vcat, treestats_dfs)
chunk_prefix = "$(proj_info.id)_hotnet_perm_treecut_stats_$(proj_info.hotnet_ver)"
open(ZstdCompressorStream, joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(job_info.chunk).jlser.zst"), "w") do io
    serialize(io, (job_info, bait_id, chunk_trees, chunkstats_df))
end
@info "Permuted trees treecut statistics chunk #$(job_info.chunk) complete"
