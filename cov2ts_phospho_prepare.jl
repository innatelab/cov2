proj_info = (id = "cov2",
             data_ver = "20200428",
             msfolder = "cov2timecourse_phospho_dia_20200423",
             )
using Pkg
Pkg.activate(@__DIR__)
using Revise
using StatsBase, DataFrames, CSV, FastaIO, JLD2

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl")
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "spectronaut_utils.jl"));

# read collapsed PTM report that we use
ptmcollapsed_df = CSV.read(joinpath(data_path, proj_info.msfolder, "COV2_DIA_phospho_0.75probablity_no normalization.txt"))
rename!(col -> replace(string(col), r"^\w\:\s" => ""), ptmcollapsed_df)
PTM_collapse_key_matches = match.(Ref(r"(\w+)_\w(\d+)_M\d"), ptmcollapsed_df.PTM_collapse_key)
ptmcollapsed_df.PTM_collapse_gene_name = string.(getindex.(PTM_collapse_key_matches, 1))
ptmcollapsed_df.PTM_collapse_ref_pos = parse.(Int, getindex.(PTM_collapse_key_matches, 2))
#ptmcollapsed_df.PTM_collapse_ref_protein_ac = getindex.(match.(Ref(r"(?:^|;)(\w+(?:-\d+)?)(?:;|$)"),
#                                                  ptmcollapsed_df[!, Symbol("PG.UniProtIds")]), 1)

# read direct Spectronaut output
pmsreport_df, colgroups = SpectronautUtils.read_peptides_report(joinpath(data_path, proj_info.msfolder,
        "20200422_071729_200421 COVID phosphoenrichment NO normalization _Report.txt"))

# find all modified_peptide -to- sequence matches
pepmodstates_df = unique!(pmsreport_df[:, colgroups[:precursor]])
pepmodstates_df.pepmodstate_id = 1:nrow(pepmodstates_df)
pepmodstates_df.peptide = replace.(pepmodstates_df.pepmodstate_seq, Ref(r"^_|_\.\d$|\[[^]]+\]" => ""))
length(unique(pepmodstates_df.pepmodstate_seq))
pmsreport_df = join(pmsreport_df, select(pepmodstates_df, [:pepmodstate_id, :pepmodstate_seq, :peptide]), kind=:left, on=:pepmodstate_seq)

pms2pepmatches_df = DelimDataUtils.expand_delim_column(pepmodstates_df, list_col=:majority_protein_acs, elem_col=:protein_ac, key_col=:pepmodstate_id)
pms2pepmatches_df.pms2pepmatches_ix = 1:nrow(pms2pepmatches_df)
pms2pepmatches_df.peptide_poses = DelimDataUtils.expand_delim_column(pepmodstates_df, list_col=:peptide_poses, elem_col=:peptide_pos, key_col=:pepmodstate_id).peptide_pos

pms2pepmatch_df = DelimDataUtils.expand_delim_column(pms2pepmatches_df, list_col=:peptide_poses, elem_col=:peptide_pos,
                                                     key_col=:pms2pepmatches_ix, delim=",")
pms2pepmatch_df.peptide_pos = parse.(Int, pms2pepmatch_df.peptide_pos)
pms2pepmatch_df.protein_ac = pms2pepmatches_df[pms2pepmatch_df.pms2pepmatches_ix, :protein_ac]
pms2pepmatch_df.pepmodstate_id = pms2pepmatches_df[pms2pepmatch_df.pms2pepmatches_ix, :pepmodstate_id]
pms2pepmatch_df.peptide = pepmodstates_df[pms2pepmatch_df.pepmodstate_id, :peptide]

pepmatches_df = unique!(select(pms2pepmatch_df, [:protein_ac, :peptide, :peptide_pos]))

# extract inidividual PTMs and their offsets from the modified sequence
function append_ptms!(out::AbstractDataFrame,
                      id::Integer,
                      seq::AbstractString)
    curoffset = 1 # offset in real AA seq
    for m in eachmatch(r"\[(\w+)\s\((\w+)\)\]", seq)
        push!(out, (pepmodstate_id=id, ptm_type=m[1],
                    ptm_AAs=m[2], data_ptm_AA=seq[m.offset-1],
                    ptm_offset=m.offset-1 - curoffset - 1))
        curoffset += length(m.match)
    end
    return out
end

function append_ptms_from_locprob_seq!(out::AbstractDataFrame,
            pmsreport_ix::Integer,
            pepmodstate_id::Integer,
            seq::AbstractString)
    curoffset = 1 # offset in real AA seq
    for m in eachmatch(r"\[(\w+)\s\((\w+)\):\s(\d+(?:\.\d+)?)%\]", seq)
        push!(out, (pmsreport_ix=pmsreport_ix, pepmodstate_id=pepmodstate_id,
                    peptide=replace(replace(seq, r"^_|_\.\d$|\[[^]]+\]" => ""), "_" => ""),
                    ptm_type=m[1],
                    ptm_AAs=m[2], data_ptm_AA=seq[m.offset-1],
                    ptm_offset=m.offset-1 - curoffset - 1,
                    ptm_locprob=parse(Float64, m[3])))
        curoffset += length(m.match)
    end
    return out
end
ptmoffsets_with_locprob_df = DataFrame()
for (i, r) in enumerate(eachrow(pmsreport_df))
    append_ptms_from_locprob_seq!(ptmoffsets_with_locprob_df, i, r.pepmodstate_id, r.pepmodstate_locprob_seq)
end
categorical!(ptmoffsets_with_locprob_df, [:ptm_type, :ptm_AAs])
ptmoffsets_df = unique!(select(ptmoffsets_with_locprob_df, Not([:pmsreport_ix, :ptm_locprob, :pepmodstate_id])))
countmap(ptmoffsets_df.data_ptm_AA)

ptm2pepmatch_df = join(ptmoffsets_df, pepmatches_df, on=:peptide, kind=:inner)
ptm2pepmatch_df.data_ptm_pos = ptm2pepmatch_df.peptide_pos .+ ptm2pepmatch_df.ptm_offset

# expand protein acs and peptide_start columns (note: the position in PTM_collapse_key may not match the position in given sequence)
ptmcollapsed_expanded_df = join(join(DelimDataUtils.expand_delim_column(ptmcollapsed_df, list_col=Symbol("PG.UniProtIds"),
                                                                   elem_col=:PTM_collapse_protein_ac, key_col=:PTM_collapse_key),
                                DelimDataUtils.expand_delim_column(ptmcollapsed_df, list_col=Symbol("PEP.PeptidePosition"),
                                                                   elem_col=:peptide_pos, key_col=:PTM_collapse_key),
                                on=:PTM_collapse_key),
                                DelimDataUtils.expand_delim_column(ptmcollapsed_df, list_col=Symbol("EG.PrecursorId"),
                                                                   elem_col=:pepmodstate_seq, key_col=:PTM_collapse_key),
                                on=:PTM_collapse_key)
# handle multiple peptide matches per protein
ptmcollapsed_expanded_df.tmp_ix = 1:nrow(ptmcollapsed_expanded_df)
ptmcollapsed_expanded_df = select!(join(select(ptmcollapsed_expanded_df, Not(:peptide_pos)),
                                DelimDataUtils.expand_delim_column(ptmcollapsed_expanded_df, list_col=:peptide_pos, elem_col=:peptide_pos, key_col=:tmp_ix, delim=","),
                                on=:tmp_ix), Not(:tmp_ix))
ptmcollapsed_expanded_df.peptide_pos = parse.(Int, ptmcollapsed_expanded_df.peptide_pos)
ptmcollapsed_expanded_df.peptide = replace.(ptmcollapsed_expanded_df.pepmodstate_seq, Ref(r"^_|_\.\d$|\[[^]]+\]" => ""))
PTM_collapse_key_expaned_matches = match.(Ref(r"(\w+)_\w(\d+)_M\d"), ptmcollapsed_expanded_df.PTM_collapse_key)
ptmcollapsed_expanded_df.PTM_collapse_gene_name = string.(getindex.(PTM_collapse_key_expaned_matches, 1))
ptmcollapsed_expanded_df.PTM_collapse_ref_pos = parse.(Int, getindex.(PTM_collapse_key_expaned_matches, 2))

# control examples
#filter(r -> r.PTM_collapse_key == "EGFR_S946_M1", ptmcollapsed_expanded_df)
#filter(r -> r.protein_ac == "E9PFD7" && r.data_ptm_pos == 946, ptm2pms_df)

# merge with ptm2pms_df to eliminate false expansion
ptmcollapsed_expanded_df = join(ptmcollapsed_expanded_df, ptm2pepmatch_df,
                                on=[:PTM_collapse_protein_ac => :protein_ac, :peptide, :peptide_pos,
                                :PTM_collapse_ref_pos => :data_ptm_pos],
                                kind=:inner)
#filter(r -> r.PTM_collapse_key == "AAAS_T397_M2", ptmcollapsed_expanded_df)
# check if some PTMs in the report was not matched to exact protein sequence
setdiff(ptmcollapsed_df.PTM_collapse_key, Set(ptmcollapsed_expanded_df.PTM_collapse_key))

# read fasta sequences used for quanting
dataseqs = FastaIO.readfasta(joinpath(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta"))
dataseqs_header = match.(Ref(r"(?:sp|tr)\|([^|]+)\|(\w+)\s"), first.(dataseqs))
dataseqs_df = DataFrame(dataseq_ix = eachindex(dataseqs),
                        dataseq_len = length.(last.(dataseqs)),
                        protein_ac = string.(getindex.(dataseqs_header, 1)),
                        protein_name = string.(getindex.(dataseqs_header, 2)))
dataseqs_df.has_ptms = dataseqs_df.protein_ac .∈ Ref(Set(pms2pepmatch_df.protein_ac))
#dataseqs_df.in_psitep = dataseqs_df.protein_ac .∈ Ref(Set(psitepseqs_df.protein_ac))
#sum(data_seq_df.has_ptms .& data_seq_df.in_psp)
#sum(data_seq_df.has_ptms .& .!data_seq_df.in_psp)

# read reference sequences of PhosphoSitePlus
psitepseqs = FastaIO.readfasta(joinpath(party3rd_data_path, "PhosphoSitePlus", "Phosphosite_PTM_seq_noheader.fasta"))
psitepseqs_header_matches = match.(Ref(r"GN:([^|]+)\|([^|]+)\|([^|]+)\|([^|]*)$"), first.(psitepseqs))
psitepseqs_df = DataFrame(
    psitepseq_ix = eachindex(psitepseqs),
    psitepseq_len = length.(last.(psitepseqs)),
    gene_name = string.(getindex.(psitepseqs_header_matches, 1)),
    protein_name = string.(getindex.(psitepseqs_header_matches, 2)),
    organism = string.(getindex.(psitepseqs_header_matches, 3)),
    protein_ac = string.(getindex.(psitepseqs_header_matches, 4)),
)
psitepseqs_df.is_ref_isoform = .!occursin.(Ref(r"-\d+$"), psitepseqs_df.protein_ac)

# align MS sequences to PhosphoSitePlus sequences
using BioSequences, BioAlignments
dataseq2psitepseq_df = join(select(filter(r -> r.has_ptms, dataseqs_df), [:dataseq_ix, :dataseq_len, :protein_ac]),
                            select(psitepseqs_df, [:psitepseq_ix, :protein_ac, :gene_name, :is_ref_isoform]),
                            on=:protein_ac, kind=:inner)
dataseq2psitepseq_df.data2psitepagn_ix = eachindex(dataseq2psitepseq_df.protein_ac)
data2psitep_agns = [pairalign(GlobalAlignment(),
                              LongAminoAcidSeq(dataseqs[r.dataseq_ix][2]),
                              LongAminoAcidSeq(psitepseqs[r.psitepseq_ix][2]),
                              AffineGapScoreModel(BLOSUM90, gap_open=-10, gap_extend=-1))
                    for r in eachrow(dataseq2psitepseq_df)]

ptm2pepmatch_df = join(ptm2pepmatch_df, select(dataseqs_df, [:dataseq_ix, :protein_ac]), on=:protein_ac, kind=:left)
ptm2psitepseq_df = join(ptm2pepmatch_df, dataseq2psitepseq_df, on=[:protein_ac, :dataseq_ix], kind=:left)
rename!(ptm2psitepseq_df, :peptide_pos => :aligned_peptide_pos)
ptm2psitepseq_df.data_ptm_pos = ptm2psitepseq_df.aligned_peptide_pos .+ ptm2psitepseq_df.ptm_offset
ptm2psitepseq_df.psitep_ptm_pos = [ismissing(r.data2psitepagn_ix) || r.data_ptm_pos > r.dataseq_len ? missing :
                                   seq2ref(alignment(data2psitep_agns[r.data2psitepagn_ix]), r.data_ptm_pos)[1]
                                   for r in eachrow(ptm2psitepseq_df)]
ptm2psitepseq_df.psitep_ptm_match = [ismissing(r.data2psitepagn_ix) || r.data_ptm_pos > r.dataseq_len ? missing :
                                     seq2ref(alignment(data2psitep_agns[r.data2psitepagn_ix]), r.data_ptm_pos)[2]
                                     for r in eachrow(ptm2psitepseq_df)]
ptm2psitepseq_df.psitep_ptm_AA = [ismissing(r.psitep_ptm_pos) ? missing : psitepseqs[r.psitepseq_ix][2][r.psitep_ptm_pos]
                                  for r in eachrow(ptm2psitepseq_df)]
countmap(ptm2psitepseq_df.psitep_ptm_match)
filter(r -> !ismissing(r.psitep_ptm_pos) && r.data_ptm_pos != r.psitep_ptm_pos, ptm2psitepseq_df)

ptmcollapsed2psitep_df = join(ptmcollapsed_expanded_df, ptm2psitepseq_df,
                              on=[:ptm_type, :ptm_AAs, :data_ptm_AA, :peptide, :ptm_offset], kind=:left)
ptmcollapsed2psitep_df.ptm_pos_dist = abs.(ptmcollapsed2psitep_df.peptide_pos .- ptmcollapsed2psitep_df.aligned_peptide_pos)
ptmcollapsed2psitep_df.has_match = .!ismissing.(ptmcollapsed2psitep_df.psitep_ptm_pos)
ptmcollapsed2psitep_df.protein_match_rank =
    ifelse.(coalesce.(ptmcollapsed2psitep_df.PTM_collapse_gene_name .== ptmcollapsed2psitep_df.gene_name, false), 4, 0) .+
    ifelse.(coalesce.(ptmcollapsed2psitep_df.PTM_collapse_protein_ac .== ptmcollapsed2psitep_df.protein_ac, false), 2, 0) .+
    ifelse.(coalesce.(ptmcollapsed2psitep_df.is_ref_isoform, false), 1, 0)
countmap(ptmcollapsed2psitep_df.protein_match_rank)

df = filter(r -> r.PTM_collapse_key == "AAAS_S371_M1", ptmcollapsed2psitep_df)[:, [:PTM_collapse_key, :ptm_type, :ptm_AAs, :data_ptm_AA, :gene_name, :protein_ac, :psitep_ptm_pos, :protein_match_rank]]
filter(r -> occursin(r"EGFR_", r.PTM_collapse_key), ptmcollapsed_df)[:, [:PTM_collapse_key, Symbol("EG.PrecursorId"), Symbol("PG.UniProtIds")]]

ptmcollapsed2psitep_best_df = by(ptmcollapsed2psitep_df, :PTM_collapse_key) do matches_df
    # sort by best match
    match_order = sortperm(matches_df, (order(:has_match, rev=true), order(:protein_match_rank, rev=true), :ptm_pos_dist))
    return matches_df[match_order[1]:match_order[1], :]
end

function flanking_sequence(seq::AbstractString, pos::Integer; flanklen::Integer=15)
    start_flank = pos - flanklen
    res = repeat('_', length(start_flank:0))
    res *= uppercase(seq[max(start_flank, 1):pos-1])
    res *= uppercase(seq[pos])
    end_flank = pos + flanklen
    res *= uppercase(seq[pos+1:min(end_flank, length(seq))])
    res *= repeat('_', length(length(seq)+1:end_flank))
    return res
end

countmap(length.(skipmissing(ptmcollapsed2psitep_best_df.flanking_15AAs)))

ptmcollapsed2psitep_best_df.flanking_15AAs = [
    ismissing(r.dataseq_ix) || r.data_ptm_pos > length(dataseqs[r.dataseq_ix][2]) ? missing :
        flanking_sequence(dataseqs[r.dataseq_ix][2], r.data_ptm_pos, flanklen=15)
    for r in eachrow(ptmcollapsed2psitep_best_df)]

ptmcollapsed2psitep_matches_df = select(ptmcollapsed2psitep_best_df, [:PTM_collapse_key, :ptm_type, :ptm_AAs, :data_ptm_AA, :gene_name, :protein_ac, :psitep_ptm_AA, :psitep_ptm_pos, :flanking_15AAs])
countmap((ptmcollapsed2psitep_best_df.PTM_collapse_protein_ac .!= ptmcollapsed2psitep_best_df.protein_ac) .&
         (ptmcollapsed2psitep_best_df.PTM_collapse_ref_pos .!= ptmcollapsed2psitep_best_df.psitep_ptm_pos))

ptmcollapsed_fixed_df = join(ptmcollapsed_df,
                             rename!(select(ptmcollapsed2psitep_best_df, [:PTM_collapse_key, :gene_name, :protein_ac, :psitep_ptm_AA, :data_ptm_pos, :psitep_ptm_pos, :flanking_15AAs]),
                                     :gene_name => :psitep_gene_name, :protein_ac => :psitep_protein_ac),
                             on=:PTM_collapse_key)
CSV.write(joinpath(data_path, proj_info.msfolder, "COV2_DIA_phospho_0.75probablity_no normalization_psitep.txt"),
          ptmcollapsed_fixed_df, delim='\t')

CSV.write(joinpath(data_path, proj_info.msfolder, "COV2_DIA_phospho_0.75probablity_no normalization_psitep_nodata.txt"),
          ptmcollapsed2psitep_best_df, delim='\t')

@save(joinpath(scratch_path, "$(proj_info.msfolder)_data_$(proj_info.data_ver).jld2"),
      ptmcollapsed2psitep_best_df, ptmcollapsed_fixed_df)

psitep_annot_files = [:KinaseSubstrate => "Kinase_Substrate_Dataset",
                      :DiseaseAssoc => "Disease-associated_sites",
                      :Phosphosites => "Phosphorylation_site_dataset",
                      :Regsites => "Regulatory_sites"]
psitep_colmap = Dict(
    :GENE => :gene_name, :SUB_GENE => :gene_name,
    :PROTEIN => :protein_name, :KINASE => :kinase_protein_name, :SUBSTRATE => :protein_name,
    :ACC_ID => :protein_ac, :KIN_ACC_ID => :kinase_protein_ac, :SUB_ACC_ID => :protein_ac,
    :GENE_ID => :entrez_id, :SUB_GENE_ID=> :entrez_id,
    :ORGANISM => :organism, :KIN_ORGANISM => :kinase_organism, :SUB_ORGANISM => :organism,
    Symbol("SITE_+/-7_AA") => :flanking_AAs,
    :DOMAIN => :domain,
    :DISEASE => :disease
)
psitep_custom_colmap = Dict(
    :KinaseSubstrate => Dict(
        :GENE => :kinase_gene_name,
    ),
    :Regsites => Dict(
        :ON_FUNCTION => :reg_function,
        :ON_PROT_INTERACT => :reg_prot_iactions,
        :ON_OTHER_INTERACT => :reg_other_iactions,
        :PMIDs => :reg_pubmed_ids,
    ),
    :DiseaseAssoc => Dict(
        :PMIDs => :disease_pubmed_ids,
    )
)

psitep_annot_dfs = Dict(begin
    df = CSV.read(joinpath(party3rd_data_path, "PhosphoSitePlus", filename), header=4)
    custom_colmap = get(psitep_custom_colmap, dsname, nothing)
    for col in names(df)
        newcol = nothing
        isnothing(custom_colmap) || (newcol = get(custom_colmap, col, nothing))
        isnothing(newcol) && (newcol = get(psitep_colmap, col, nothing))
        isnothing(newcol) || rename!(df, col => newcol)
    end
    if hasproperty(df, :MOD_RSD)
        ptm_matches = match.(Ref(r"(\w)(\d+)-(\w+)"), df.MOD_RSD)
        df.ptm_AA = string.(getindex.(ptm_matches, 1))
        df.ptm_pos = parse.(Int, getindex.(ptm_matches, 2))
        df.ptm_code = string.(getindex.(ptm_matches, 3))
    elseif hasproperty(df, :SUB_MOD_RSD)
        ptm_matches = match.(Ref(r"(\w)(\d+)"), df.SUB_MOD_RSD)
        df.ptm_AA = string.(getindex.(ptm_matches, 1))
        df.ptm_pos = parse.(Int, getindex.(ptm_matches, 2))
    end
    dsname => df
end for (dsname, filename) in psitep_annot_files)

collapsevals(vals; delim="; ") =
    all(ismissing, vals) ? missing : join(skipmissing(vals), delim)

#psitep_annots_df = join(rename!(select(ptmcollapsed2psitep_best_df, [:PTM_collapse_key, :gene_name, :protein_ac, :data_ptm_pos, :psitep_ptm_AA, :psitep_ptm_pos, :flanking_15AAs]),
#                             :psitep_ptm_AA=>:ptm_AA, :psitep_ptm_pos=>:ptm_pos),

psitep_annots_df = select(psitep_annot_dfs[:Phosphosites], [:protein_ac, :ptm_pos, :domain])
psitep_annots_df = join(psitep_annots_df,
                     by(psitep_annot_dfs[:KinaseSubstrate], [:protein_ac, :ptm_pos],
                        kinase_gene_names = :kinase_gene_name => collapsevals),
                     on=[:protein_ac, :ptm_pos], kind=:outer)
psitep_annots_df = join(psitep_annots_df,
                     by(psitep_annot_dfs[:Regsites], [:protein_ac, :ptm_pos],
                        reg_function = :reg_function => collapsevals,
                        reg_prot_iactions = :reg_prot_iactions => collapsevals,
                        reg_other_iactions= :reg_other_iactions => collapsevals,
                        reg_pubmed_ids = :reg_pubmed_ids => collapsevals),
                     on=[:protein_ac, :ptm_pos], kind=:outer)
psitep_annots_df = join(psitep_annots_df,
                     by(psitep_annot_dfs[:DiseaseAssoc], [:protein_ac, :ptm_pos],
                        diseases = :disease => collapsevals,
                        diseases_pubmed_ids = :disease_pubmed_ids => collapsevals),
                     on=[:protein_ac, :ptm_pos], kind=:outer)

@save(joinpath(scratch_path, "phosphositeplus_annotations_20200502.jld2"), psitep_annot_dfs, psitep_annots_df)

CSV.write(joinpath(data_path, proj_info.msfolder, "COV2_DIA_phospho_psitep_annotations.txt"),
          ptm_annots_df, delim='\t')

ptmcollapsed_annot_df = join(ptmcollapsed_fixed_df,
                             select(ptm_annots_df, Not([:gene_name, :protein_ac, :ptm_AA, :ptm_pos, :data_ptm_pos, :flanking_15AAs])),
                             on=:PTM_collapse_key)
CSV.write(joinpath(data_path, proj_info.msfolder, "COV2_DIA_phospho_0.75probablity_no normalization_psitep_annotated.txt"),
          ptmcollapsed_annot_df, delim='\t')

filter(r -> !ismissing(r.psitep_ptm_pos), ptmcollapsed_annot_df)[nonunique(filter(r -> !ismissing(r.psitep_ptm_pos), ptmcollapsed_annot_df), [:PTM_collapse_key, :psitep_protein_ac, :psitep_ptm_pos]),
    [:PTM_collapse_key, :psitep_protein_ac, :psitep_ptm_pos]]
filter(r -> r.PTM_collapse_key == "MAP2K6_S207_M1", ptmcollapsed_annot_df)[:, [:reg_function, :reg_prot_iactions, :reg_other_iactions, :reg_pubmed_ids]]
[:PTM_collapse_key, :psitep_protein_ac, :psitep_ptm_pos]]
