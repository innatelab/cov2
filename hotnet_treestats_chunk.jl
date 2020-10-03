#=
job_info = (project = "cov2",
            name = "cov2_treestats",
            hotnet_ver = "20200901",
            id = 0,
            chunk = 23,
            nbaits_perchunk = 2)
=#
job_info = (project = ARGS[1],
            name = ARGS[2],
            hotnet_ver = ARGS[3],
            id = parse(Int, ARGS[4]),
            chunk = parse(Int, ARGS[5]),
            nbaits_perchunk = parse(Int, ARGS[6]))
@info "Treestats for $(job_info.project) (job '$(job_info.name)' id=$(job_info.id))" *
      " chunk #$(job_info.chunk)"

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
bait_ix1 = (job_info.chunk-1) * job_info.nbaits_perchunk + 1
@assert (bait_ix1 <= length(bait_ids)) "Chunk #$(job_info.chunk): outside of valid range"
bait_ixs = bait_ix1:min(bait_ix1 + job_info.nbaits_perchunk - 1, length(bait_ids))

reactomefi_mtx = Matrix(LightGraphs.weights(reactomefi_digraph_rev));

treestats_dfs = Vector{DataFrame}()
for (i, bait_ix) in enumerate(bait_ixs)
    @info "Processing bait #$(bait_ix) ($(bait_ids[bait_ix])), $i of $(length(bait_ixs)))..."
    bait_id, vertex_weights, vertex_walkweights, _, sink_ixs, _, _, _, tree =
        open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver)",
                                        "bait_$(bait_ix)_perms.jlser.zst"), "r") do io
        deserialize(io)
    end
    @assert bait_id == bait_ids[bait_ix] "Bait #$(bait_ix) is $(bait_ids[bait_ix]), got data for $(bait_id)"
    # recreate walkmtx since it's too expensive to store it in perm_input
    stepmtx = HHN.stepmatrix(reactomefi_mtx, inedge_weights=vertex_weights .+ randomwalk_params.inedge_weight_min)
    walkmtx = HHN.similarity_matrix(stepmtx, vertex_weights, restart_probability=randomwalk_params.restart_prob)
    @info "Treecut statistics for bait #$(bait_ix) ($(bait_id), $i of $(length(bait_ixs)))..."
    local treestats_df = HHN.treecut_stats(tree, walkmatrix=walkmtx, maxweight=randomwalk_params.flow_weight_max,
                                    sources=findall(>=(randomwalk_params.source_weight_min), vertex_weights),
                                    sinks=sink_ixs,
                                    nflows_ratio=0.999, pools=nothing)
    treestats_df[!, :tree] .= 0
    treestats_df[!, :is_permuted] .= false
    treestats_df[!, :bait_id] .= bait_id
    push!(treestats_dfs, treestats_df)
end
treestats_df = vcat(treestats_dfs...)

chunk_prefix = "$(proj_info.id)_hotnet_treestats_$(proj_info.hotnet_ver)"
isdir(joinpath(scratch_path, chunk_prefix)) || mkdir(joinpath(scratch_path, chunk_prefix))
open(ZstdCompressorStream, joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(job_info.chunk).jlser.zst"), "w") do io
    serialize(io, (job_info, bait_ids[bait_ixs], treestats_df))
end
@info "Treestats chunk #$(job_info.chunk) complete"
