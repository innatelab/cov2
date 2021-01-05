#=
job_info = (project = "cov2",
            name = "cov2_permtrees",
            hotnet_ver = "20200901",
            id = 0,
            chunk = 23,
            ntrees_perchunk = 20)
=#
job_info = (project = ARGS[1],
            name = ARGS[2],
            hotnet_ver = ARGS[3],
            id = parse(Int, ARGS[4]),
            chunk = parse(Int, ARGS[5]),
            ntrees_perchunk = parse(Int, ARGS[6]))
@info "Permuted tree for $(job_info.project) (job '$(job_info.name)' id=$(job_info.id))" *
      " chunk #$(job_info.chunk)"
flush(stdout); flush(stderr)

proj_info = (id = job_info.project,
             hotnet_ver = job_info.hotnet_ver)

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(base_scratch_path, proj_info.id)
const plots_path = joinpath(analysis_path, "plots")

using JLD2, DataFrames, CSV, Serialization, CodecZstd
using LinearAlgebra, LightGraphs, HierarchicalHotNet
HHN = HierarchicalHotNet

_, _, reactomefi_digraph_rev, bait2nperms, randomwalk_params =
open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver).jlser.zst"), "r") do io
    deserialize(io)
end;

bait_ids = first.(bait2nperms)
nbaittrees = cumsum(last.(bait2nperms))
chunk_1st_tree = (job_info.chunk-1) * job_info.ntrees_perchunk + 1
bait_ix = searchsortedfirst(nbaittrees, chunk_1st_tree)
@assert (bait_ix <= length(bait_ids)) "Chunk #$(job_info.chunk) outside of permuted trees range"

chunk_prefix = "$(proj_info.id)_hotnet_perm_$(proj_info.hotnet_ver)"
if isfile(joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(bait_ix)_$(job_info.chunk).jlser.zst"))
    @info "Chunk $(chunk_prefix)-#$(job_info.chunk) already exists, skipping"
    exit()
end

bait_prev_lastree = bait_ix > 1 ? nbaittrees[bait_ix-1] : 0
# indices of chunk trees to process (within bait_id perm_trees vector)
chunk_treeixs = (chunk_1st_tree - bait_prev_lastree):min(
            chunk_1st_tree - bait_prev_lastree - 1 + job_info.ntrees_perchunk,
            bait2nperms[bait_ix][2])
@info "Processing bait #$(bait_ix) ($(bait_ids[bait_ix])) permuted trees $(chunk_treeixs)"
flush(stdout); flush(stderr)

reactomefi_mtx = Matrix(LightGraphs.weights(reactomefi_digraph_rev));

bait_id, _, _, perm_vertex_weights, sink_ixs, diedge_ixs, _, _, _ =
open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver)",
                                      "bait_$(bait_ix)_perms.jlser.zst"), "r") do io
    deserialize(io)
end;
@assert bait_id == bait_ids[bait_ix] "Bait #$(bait_ix) is $(bait_ids[bait_ix]), got data for $(bait_id)"

vertex_walkweights = Matrix{Float64}(undef, (size(perm_vertex_weights, 1), length(chunk_treeixs)))
diedge_stepweights = Matrix{Float64}(undef, (length(diedge_ixs), length(chunk_treeixs)))
diedge_walkweights = Matrix{Float64}(undef, (length(diedge_ixs), length(chunk_treeixs)))
trees = Vector{HHN.SCCTree{Float64}}(undef, length(chunk_treeixs))
treestats_dfs = Vector{DataFrame}()
W = HHN.weighttype(eltype(trees))
#seedling = HHN.SCCSeedling(reactomefi_mtx)

for (i, permix) in enumerate(chunk_treeixs)
    @info "Growing permuted tree $bait_id-#$permix ($i of $(length(chunk_treeixs)))..."
    flush(stdout); flush(stderr)
    vertex_weights = convert(Vector{W}, view(perm_vertex_weights, :, permix))
    stepmtx = HHN.stepmatrix(reactomefi_mtx,
                             inedge_weights=vertex_weights .+ randomwalk_params.inedge_weight_min)
    walkmtx = HHN.similarity_matrix(stepmtx, vertex_weights,
        restart_probability=randomwalk_params.restart_prob)
    diedge_stepweights[:, i] .= vec(stepmtx)[diedge_ixs]
    diedge_walkweights[:, i] .= vec(walkmtx)[diedge_ixs]
    vertex_walkweights[:, i] .= dropdims(sum(walkmtx, dims=2), dims=2)
    trees[i] = HHN.scctree(walkmtx, #seedling=seedling,
                           verbose=false, method=:bisect)
    @info "Treecut statistics for $bait_id-#$permix"
    flush(stdout); flush(stderr)
    local treestats_df = HHN.treecut_stats(trees[i], walkmatrix=walkmtx, maxweight=randomwalk_params.flow_weight_max,
                                    sources=findall(>=(randomwalk_params.source_weight_min), vertex_weights),
                                    sinks=sink_ixs,
                                    nflows_ratio=0.99, pools=nothing)
    treestats_df[!, :tree] .= permix
    treestats_df[!, :is_permuted] .= true
    push!(treestats_dfs, treestats_df)
end
treestats_df = vcat(treestats_dfs...)

chunk_prefix = "$(proj_info.id)_hotnet_perm_$(proj_info.hotnet_ver)"
isdir(joinpath(scratch_path, chunk_prefix)) || mkdir(joinpath(scratch_path, chunk_prefix))
open(ZstdCompressorStream, joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(bait_ix)_$(job_info.chunk).jlser.zst"), "w") do io
    serialize(io, (job_info, bait_id, chunk_treeixs, vertex_walkweights,
                   diedge_stepweights, diedge_walkweights, trees, treestats_df))
end
@info "Permuted trees chunk #$(job_info.chunk) complete"
