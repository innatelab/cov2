proj_info = (id = "cov2",
             data_ver = "20200907",
             msfolder = "snaut_parsars_ptm_20200907",
             ptm_locprob_min = 0.75,
             )
using Pkg
Pkg.activate(@__DIR__)
using Revise
using StatsBase, DataFrames, CSV, FastaIO, JLD2, CodecZlib, Base.Filesystem

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl")
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "spectronaut_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "fasta_reader.jl"));
Revise.includet(joinpath(misc_scripts_path, "ptm_extractor.jl"));
Revise.includet(joinpath(misc_scripts_path, "phosphositeplus_reader.jl"));

# read direct Spectronaut output
pmsreports = Dict(name => begin
    report_df, colgroups = SpectronautUtils.read_report(joinpath(data_path, proj_info.msfolder, filename),
            import_data=[:protgroup, :peptide, :pepmodstate, :quantity, :metrics], delim=',')
    (data = report_df, colgroups = colgroups)
end for (name, filename) in [(:phospho, "phospho/COVID_phospho_EG level.csv"),
                             (:ubi, "ubi/Ubi_COVID_EG_level_Report.csv")])

# fix contaminant names (strip "CON__<SRC>:")
for report in values(pmsreports), col in [:protgroup_sn_id, :majority_protein_acs, :majority_protein_acs_alt, :peptide_protein_acs]
    df = report.data
    mask = .!ismissing.(df[!, col])
    df[mask, col] .= replace.(df[mask, col], Ref(r"CON__(?:[^:;]+):" => ""))
end

samples_df = CSV.read(joinpath(data_path, proj_info.msfolder, "msruns_info.txt"))
samples_df.treatment = levels!(categorical(samples_df.treatment), ["mock", "SARS_CoV2", "SARS_CoV"])

rawfiles_df = vcat([begin
    df = select!(SpectronautUtils.metrics_colinfo(last.(report.colgroups[:quantity])), [:rawfile_ix, :rawfile]) |> unique!
    df[!, :dataset] .= dsname
    df.msrun_ix = [begin
            m = match(r"_(\d+)(?:_\d+)?\.raw$", fname)
            parse(Int, m[1])
        end for fname in df.rawfile]
    df
end for (dsname, report) in pmsreports]...)
msruns_df = sort!(innerjoin(rawfiles_df, samples_df, on=:msrun_ix),
                  [:dataset, :treatment, :timepoint, :replicate])

msfasta_path = joinpath(data_path, proj_info.msfolder, "fasta used")
proteins_df = let
    human_df = Fasta.read_uniprot(joinpath(msfasta_path, "Human_July2019_with_isoforms_only_swissprot.fasta"))
    human_df[!, :is_contaminant] .= false
    human_df[!, :is_viral] .= false

    sars_df = Fasta.read_uniprot(joinpath(msfasta_path, "SARS_CoV_20200817.fasta"))
    sars_df.genename = "SARS_CoV_" .* sars_df.genename
    sars_df[!, :is_contaminant] .= false
    sars_df[!, :is_viral] .= true

    sars2_df = Fasta.read_uniprot(joinpath(msfasta_path, "SARS_CoV2_20200908.fasta"))
    sars2_df.genename = "SARS_CoV2_" .* sars2_df.genename
    sars2_df[!, :is_contaminant] .= false
    sars2_df[!, :is_viral] .= true

    contaminants_df = Fasta.read_contaminants(joinpath(msfasta_path, "contaminants_20200405.fasta"))
    contaminants_df[!, :is_contaminant] .= true
    contaminants_df[!, :is_viral] .= false
    contaminants_df[!, :protein_code] .= missing
    contaminants_df[!, :genename] .= missing
    contaminants_df[!, :protein_existence] .= missing
    contaminants_df[!, :seq_version] .= missing

    res = vcat(human_df, sars_df, sars2_df, contaminants_df)
    res.orig_genename = res.genename
    res.genename = coalesce.(res.orig_genename, res.protein_ac)
    res.protein_ac_isoform = [(isnothing(m) ? 1 : parse(Int, m[1])) for m in match.(Ref(r"-(\d+)(?:#.+)?$"), res.protein_ac)]
    res
end
protein_seqs = Dict(r.protein_ac => r.seq for r in eachrow(proteins_df))

# find all modified_peptide -to- sequence matches
protgroups_df = vcat([begin
    df = unique(report.data[!, last.(report.colgroups[:protgroup])])
    df[!, :dataset] .= dsname
    df
end for (dsname, report) in pairs(pmsreports)]...)
peptide2protgroups_df = vcat([begin
    df = unique(report.data[!, [last.(report.colgroups[:peptide]); ["protgroup_id", "majority_protein_acs"]]])
    df[!, :dataset] .= dsname
    df
end for (dsname, report) in pairs(pmsreports)]...)
peptide2protgroups_df.peptide_id = FrameUtils.indexunique(peptide2protgroups_df.peptide_seq)
# global peptides
peptides_df = sort!(unique(select(peptide2protgroups_df, [:peptide_id, :peptide_seq], copycols=false)), :peptide_id)
# global pepmods
pepmodstates_df = vcat([begin
    df = innerjoin(unique(report.data[!, [last.(report.colgroups[:pepmodstate]); ["peptide_seq"]]]),
                   peptides_df, on=:peptide_seq)
    select!(df, Not([:peptide_seq, :pepmodstate_id]))
end for (dsname, report) in pairs(pmsreports)]...) |> unique!
pepmodstates_df.pepmodstate_id = FrameUtils.indexunique(pepmodstates_df.pepmodstate_seq)
sort!(pepmodstates_df, :pepmodstate_id)
# replace pepmodstate_id and peptide_id with the global ones
pmsreports = Dict(dsname => begin
    (data = innerjoin(innerjoin(
            select!(report.data, Not([:peptide_id, :pepmodstate_id])),
            select(peptides_df, [:peptide_seq, :peptide_id], copycols=false), on=:peptide_seq),
            select(pepmodstates_df, [:pepmodstate_seq, :pepmodstate_id], copycols=false), on=:pepmodstate_seq),
     colgroups = report.colgroups)
end for (dsname, report) in pairs(pmsreports))
protein2protgroup_df = DelimDataUtils.expand_delim_column(protgroups_df, list_col=:majority_protein_acs,
                                                          elem_col=:protein_ac, key_col=[:dataset, :protgroup_id])

# read reference sequences of PhosphoSitePlus
psitepseqs_df = Fasta.read_phosphositeplus(joinpath(party3rd_data_path, "PhosphoSitePlus", "Phosphosite_PTM_seq_noheader.fasta"))
psitepseqs_df.is_ref_isoform = .!occursin.(Ref(r"-\d+$"), psitepseqs_df.protein_ac)

peptide2protein_df = PTMExtractor.peptide_matches(peptide2protgroups_df)

ptm2pms_df = innerjoin(unique!(PTMExtractor.extract_ptms(pepmodstates_df, objid_col=:pepmodstate_id)),
                       select(pepmodstates_df, [:pepmodstate_id, :peptide_id], copycols=false), on=:pepmodstate_id) |> unique!
select!(ptm2pms_df, Not(:peptide_seq))
sort!(ptm2pms_df, [:pepmodstate_id, :ptm_offset])
ptm2pms2protein_df = innerjoin(ptm2pms_df, peptide2protein_df, on=:peptide_id)
ptm2pms2protein_df.ptm_pos = ptm2pms2protein_df.peptide_pos .+ ptm2pms2protein_df.ptm_offset
ptm2protein_df = leftjoin(unique!(select(ptm2pms2protein_df, [:ptm_type, :ptm_AAs, :ptm_AA_seq, :ptm_pos, :protein_ac])),
                          select(proteins_df, [:protein_ac, :genename, :organism, :is_viral, :is_contaminant], copycols=false),
                          on=:protein_ac)
#=
ptm2protgroup_df = PTMExtractor.group_aaobjs(ptm2protein_df, innerjoin(proteins_df, protein2protgroup_df, on=:protein_ac),
                                             seqgroup_col=[:protgroup_id, :is_viral, :is_contaminant],
                                             seqid_col=[:protein_ac, :is_viral, :is_contaminant],
                                             seqrank_col=:protein_ac_isoform,
                                             obj_prefix=:ptm_, objid_col=[:ptm_type, :ptm_pos, :ptm_AA_seq],
                                             force_refseqs=true, verbose=true)
=#
ptm2gene_df = PTMExtractor.group_aaobjs(ptm2protein_df, proteins_df,
                                        seqgroup_col=[:genename, :is_viral, :is_contaminant],
                                        seqid_col=:protein_ac, seqrank_col=:protein_ac_isoform,
                                        obj_prefix=:ptm_, objid_col=[:ptm_type, :ptm_pos, :ptm_AA_seq],
                                        force_refseqs=true, verbose=true)
ptm2gene_df.flanking_15AAs = PTMExtractor.flanking_sequence.(getindex.(Ref(protein_seqs), ptm2gene_df.protein_ac), ptm2gene_df.ptm_pos, flanklen=15)

filter(r -> coalesce(r.genename, "") == "SYNPO", ptm2gene_df)
filter(r -> r.ptm_is_reference && coalesce(r.is_viral, false), ptm2gene_df) |> print
countmap(filter(r -> r.ptm_is_reference, ptm2gene_df).ptm_type)

ptmn2pms_df = select!(innerjoin(select(ptm2gene_df, [:ptm_id, :ptm_label, :protein_ac, :ptm_type, :ptm_AA_seq, :ptm_pos]),
                                ptm2pms2protein_df, on=[:protein_ac, :ptm_type, :ptm_AA_seq, :ptm_pos]),
                      [:ptm_id, :ptm_label, :ptm_type, :nptms, :pepmodstate_id]) |> unique!
ptmns_df = unique!(select(ptmn2pms_df, [:ptm_id, :ptm_label, :ptm_type, :nptms]))
sort!(ptmns_df, [:ptm_id, :nptms])
ptmns_df.ptmn_id = 1:nrow(ptmns_df)
ptmns_df.ptmn_label = categorical(ptmns_df.ptm_label .* "_M" .* string.(ptmns_df.nptms))
ptmn2pms_df = innerjoin(ptmn2pms_df, select(ptmns_df, Not(:ptm_label), copycols=false), on=[:ptm_id, :ptm_type, :nptms])
select!(ptmn2pms_df, Not([:ptm_type, :nptms]))
sort!(ptmn2pms_df, [:ptmn_id, :pepmodstate_id])

print(countmap(collect(values(countmap(ptmn2pms_df.pepmodstate_id)))))
print(countmap(collect(values(countmap(filter(r -> r.nptms == 1, ptmn2pms_df).ptmn_id)))))

print(filter(p -> p[2] >= 5, pairs(countmap(ptmn2pms_df.ptmn_id))))
print(innerjoin(filter(r -> r.ptm_id == "Phospho_MARCKS_S170", ptmn2pms_df), pepmodstates_df, on=:pepmodstate_id))

ptm2psitep_df = PTMExtractor.map_aapos(ptm2protein_df, proteins_df, psitepseqs_df,
                                       destmap_prefix=:psitep_, obj_prefix=:ptm_, objid_col=[:ptm_type, :ptm_pos, :ptm_AA_seq], verbose=true)
filter!(r -> !ismissing(r.destseq_ix), ptm2psitep_df)
ptm2psitep_df.genename = [!ismissing(r.destseq_ix) && (coalesce(psitepseqs_df.genename[r.destseq_ix], "-") != "-") ? psitepseqs_df.genename[r.destseq_ix] :
                          !ismissing(r.srcseq_ix) ? ptm2protein_df.genename[r.srcseq_ix] :
                          missing
                          for r in eachrow(ptm2psitep_df)]
ptm2psitep_df.pos_match = coalesce.(ptm2psitep_df.ptm_pos .== ptm2psitep_df.psitep_ptm_pos, false)
ptm2psitep_df.AA_match = coalesce.(ptm2psitep_df.ptm_AA_seq .== uppercase.(coalesce.(ptm2psitep_df.psitep_ptm_AA, '?')), false)
ptm2psitep_df.protein_ac_isoform = [(isnothing(m) ? 1 : parse(Int, m[1])) for m in match.(Ref(r"-(\d+)(?:#.+)?$"), ptm2psitep_df.protein_ac)]
sort!(ptm2psitep_df, [:genename, :protein_ac_isoform, order(:pos_match, rev=true), order(:AA_match, rev=true)])

pms_intensities_df = vcat([begin
    df = SpectronautUtils.pivot_longer(report.data, [:pepmodstate_id, :pepmod_seq])
    df[!, :dataset] .= dsname
    df
end for (dsname, report) in pmsreports]...)
rename!(pms_intensities_df, :EG_normfactor => :normfactor, :EG_intensity => :intensity_norm, :EG_qvalue => :qvalue)
pms_intensities_df.intensity = pms_intensities_df.intensity_norm ./ pms_intensities_df.normfactor

ptm_locprobs_df = PTMExtractor.extract_ptm_locprobs(pms_intensities_df, pepmodseq_col=:pepmod_seq, msrun_col=[:dataset, :rawfile_ix])
ptm_locprobs_df.ptm_locprob_match = ptm_locprobs_df.ptm_locprob .>= proj_info.ptm_locprob_min
ptm_locprobs_df = innerjoin(innerjoin(ptm_locprobs_df, select(pepmodstates_df, [:pepmodstate_id, :peptide_id], copycols=false), on=:pepmodstate_id),
                            peptide2protein_df, on=:peptide_id)
ptm_locprobs_df.ptm_pos = ptm_locprobs_df.ptm_offset .+ ptm_locprobs_df.peptide_pos
countmap(collect(zip(ptm_locprobs_df.dataset, ptm_locprobs_df.ptm_locprob_match)))
ptm_locprobs_df

# strip rows/cols that are no longer required
filter!(r -> !ismissing(r.intensity), select!(pms_intensities_df, Not([:pepmod_seq, :EG_locprob_seq])))

# save the files
output_path = joinpath(data_path, proj_info.msfolder, "ptm_extractor")
isdir(output_path) || mkdir(output_path)
for (fname, df) in ["rawfiles_info.txt" => msruns_df,
                    "proteins.txt.gz" => proteins_df,
                    "ptm_locprobs.txt.gz" => ptm_locprobs_df,
                    "pms_intensities.txt.gz" => pms_intensities_df,
                    "protgroups.txt.gz" => protgroups_df,
                    "peptide_to_protgroup.txt.gz" => peptide2protgroups_df,
                    "peptides.txt.gz" => peptides_df,
                    "pepmodstates.txt.gz" => pepmodstates_df,
                    "protein_to_protgroup.txt.gz" => protein2protgroup_df,
                    "peptide_to_protein.txt.gz" => peptide2protein_df,
                    "ptm_to_protein.txt.gz" => ptm2protein_df,
                    "ptm_to_gene.txt.gz" => ptm2gene_df,
                    "ptmns.txt.gz" => ptmns_df,
                    "ptmn_to_pepmodstate.txt.gz" => ptmn2pms_df]
    @info "Saving $fname..."
    if endswith(fname, ".gz")
        open(GzipCompressorStream, joinpath(output_path, fname), "w") do io
            CSV.write(io, df, delim='\t')
        end
    else
        CSV.write(joinpath(output_path, fname), df, delim='\t')
    end
end

psitep_annot_dfs = PhosphoSitePlus.read_annotations(joinpath(party3rd_data_path, "PhosphoSitePlus"), verbose=true)
psitep_annots_df = PhosphoSitePlus.combine_annotations(psitep_annot_dfs)

@save(joinpath(scratch_path, "phosphositeplus_annotations_20200502.jld2"), psitep_annot_dfs, psitep_annots_df)
