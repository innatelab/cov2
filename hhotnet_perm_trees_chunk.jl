# job_info = (project = "cov2",
#             name = "cov2_permtrees",
#             hotnet_ver = "20200415",
#             id = 0,
#             chunk = 23,
#             ntrees_perchunk = 2)
job_info = (project = ARGS[1],
            name = ARGS[2],
            hotnet_ver = ARGS[3],
            id = parse(Int, ARGS[4]),
            chunk = parse(Int, ARGS[5]),
            ntrees_perchunk = 20)
@info "Permuted tree for $(job_info.project) (job '$(job_info.name)' id=$(job_info.id))" *
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
      reactomefi_walkmtx, bait2vertex_weights_perms)

nbaittrees = cumsum(size.(getindex.(Ref(bait2vertex_weights_perms), bait_ids), 2))
chunk_1st_tree = (job_info.chunk-1) * job_info.ntrees_perchunk + 1
bait_ix = searchsortedfirst(nbaittrees, chunk_1st_tree)
@assert (bait_ix <= length(bait_ids)) "Chunk #$(job_info.chunk) outside of permuted trees range"
bait_id = bait_ids[bait_ix]
apms_vertices = bait2apms_vertices[bait_id]
vertex_weights_perms = bait2vertex_weights_perms[bait_id]
bait_idprev_lastree = bait_ix > 1 ? nbaittrees[bait_ix-1] : 0
# indices of chunk trees to process (within bait_id perm_trees vector)
chunk_trees = (chunk_1st_tree - bait_idprev_lastree):min(
            chunk_1st_tree - bait_idprev_lastree - 1 + job_info.ntrees_perchunk,
            size(vertex_weights_perms, 2))
@info "Processing $bait_id permuted trees $(chunk_trees)"
permtrees = sizehint!(Vector{HHN.SCCTree{Float64}}(), length(chunk_trees))
W = HHN.weighttype(eltype(permtrees))
permtree_mtx = similar(reactomefi_walkmtx)
seedling = HHN.SCCSeedling(permtree_mtx)
for i in chunk_trees
    @info "Growing permuted tree $bait_id-#$i"
    vertex_weights = convert(Vector{W}, view(vertex_weights_perms, :, i))
    mul!(permtree_mtx, reactomefi_walkmtx, Diagonal(vertex_weights))
    push!(permtrees, HHN.scctree(permtree_mtx, seedling=seedling,
                                 verbose=false, method=:bisect))
end

chunk_prefix = "$(proj_info.id)_hotnet_permtrees_$(proj_info.hotnet_ver)"
open(ZstdCompressorStream, joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(job_info.chunk).jlser.zst"), "w") do io
    serialize(io, (job_info, bait_id, chunk_trees, permtrees))
end
@info "Permuted trees chunk #$(job_info.chunk) complete"
