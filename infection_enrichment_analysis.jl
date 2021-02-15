proj_info = (id = "cov2",
             oesc_ver = "20201024")

datasets = Dict(
    :fp => (folder="snaut_parsars_fp_20200829",
            fit_ver="20200830",
            analysis=:msglm, ptm=nothing,
            datatype=:proteome, entity=:protein,
            label="Proteome",
            color="#226430"),
    :phospho => (folder="snaut_parsars_phospho_20201005",
                 fit_ver="20201012",
                 analysis=:msglm, ptm=:phospho,
                 datatype=:phospho, entity=:ptm,
                 label="Phosphoproteome",
                 color="#5e268f"),
    :ubi => (folder="snaut_parsars_ptm_20200907",
             fit_ver="20201012",
             analysis=:msglm, ptm=:ubiquitin, prefix=:ubi,
             datatype=:ubi, entity=:ptm,
             label="Ubiquitinome",
             color="#bf1c2c"),
    :rnaseq => (folder="parsars_rnaseq_20201020",
                fit_ver="20201020",
                analysis=:limma, ptm=nothing,
                datatype=:rnaseq, entity=:gene,
                label="Transcriptome",
                color="#383c9b"),
)

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, CSV, DataFrames, FastaIO
using JLD2
using StatsBase

@info "Project '$(proj_info.id)' analysis version=$(proj_info.oesc_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

includet(joinpath(misc_scripts_path, "frame_utils.jl"));
includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));
includet(joinpath(misc_scripts_path, "fasta_reader.jl"));
includet(joinpath(misc_scripts_path, "ms_import.jl"));
includet(joinpath(misc_scripts_path, "protgroup_assembly.jl"));
includet(joinpath(misc_scripts_path, "protgroup_crossmatch.jl"));

# load data and fits
msglm_rdata = Dict(begin
    @info "Loading $dsname analysis (fit_ver=$(dsinfo.fit_ver))..."
    dsname => (input_data = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(dsinfo.folder)_$(dsinfo.fit_ver).RData"), convert=true),
               full_data = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(dsinfo.folder)_$(dsinfo.fit_ver).RData"), convert=true),
               fit = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(dsinfo.folder)_$(dsinfo.fit_ver).RData"), convert=true))
end for (dsname, dsinfo) in pairs(datasets) if dsinfo.analysis == :msglm);

# fix some data differences
for (dsname, rdata) in pairs(msglm_rdata)
    @info "Fixing $dsname"
    msdata_full = rdata.full_data["msdata_full"]

    proteins_df = msdata_full["proteins"]
    proteins_df.organism = [ismissing(org) ? org : replace(org, r"\s+OX=\d+$" => "") for org in proteins_df.organism]
    
    hasproperty(proteins_df, :genename) && rename!(proteins_df, :genename => :gene_name)
    if hasproperty(proteins_df, :protein_name) && hasproperty(proteins_df, :protein_description)
        rename!(proteins_df, :protein_name => :protein_code, :protein_description => :protein_name)
    end

    peptides_df = msdata_full["peptides"]
    hasproperty(peptides_df, :peptide_seq) || rename!(peptides_df, :seq => :peptide_seq)

    pepmodstates_df = msdata_full["pepmodstates"]
    if haskey(msdata_full, "pepmods")
        pepmods_df = msdata_full["pepmods"]
        missing_pms_cols = intersect(setdiff([:peptide_id, :nselptms], propertynames(pepmodstates_df)), propertynames(pepmods_df))
        if !isempty(missing_pms_cols)
            pepmodstates_df = innerjoin(pepmodstates_df, pepmods_df[!, [[:pepmod_id]; missing_pms_cols]], on=:pepmod_id)
            @assert nrow(pepmodstates_df) == nrow(msdata_full["pepmodstates"])
            msdata_full["pepmodstates"] = pepmodstates_df
        end
    end

    pms_intensities_df = msdata_full["pepmodstate_intensities"]
    hasproperty(pms_intensities_df, :psm_qvalue) || rename!(pms_intensities_df, :qvalue => :psm_qvalue)
    if !hasproperty(pms_intensities_df, :peptide_id)
        pms_intensities_df = innerjoin(pms_intensities_df, pepmodstates_df[!, [:pepmodstate_id, :peptide_id]], on=:pepmodstate_id)
        @assert nrow(pms_intensities_df) == nrow(msdata_full["pepmodstate_intensities"])
        msdata_full["pepmodstate_intensities"] = pms_intensities_df
    end
end

peptide2protein_df = combine(groupby(reduce(vcat, begin
    @info "Processing peptides of $dsname"
    isptm = !isnothing(datasets[dsname].ptm)
    msdata = dsdata.full_data["msdata_full"]
    pep2prot_df = copy(msdata["peptide2protein"], copycols=false)
    if !hasproperty(pep2prot_df, :peptide_seq)
        pep2prot_df = innerjoin(msdata["peptide2protein"], msdata["peptides"],
                                on=:peptide_id)
    end
    if isptm # only use modified peptides for PTM datasets
        pepmodstates_df = filter(r -> r.nselptms > 0, msdata["pepmodstates"])
        # remove phosphosites from ubi data
        (dsname == :ubi) && filter!(r -> occursin("[GlyGly", r.pepmod_seq), pepmodstates_df)
        select!(pepmodstates_df, [:pepmodstate_id, :peptide_id])
    else
        pepmodstates_df = select(msdata["pepmodstates"], [:pepmodstate_id, :peptide_id])
    end
    pms_intensities_df = copy(msdata["pepmodstate_intensities"], copycols=false)
    if !hasproperty(pms_intensities_df, :peptide_seq)
        pms_intensities_df = innerjoin(semijoin(pms_intensities_df, pepmodstates_df, on=:pepmodstate_id),
                                       select(msdata["peptides"], [:peptide_id, :peptide_seq]), on=:peptide_id)
    end
    pepstats_df = combine(groupby(pms_intensities_df, :peptide_seq),
                          :psm_qvalue => (x -> minimum(skipmissing(x))) => :psm_qvalue)
    # set the rank of non-PTM peptides to -1, so that proteome-only observed peptides are ignored (was: still used to define
    # protein groups, but they would not split the protein group derived from PTM datasets)
    pepstats_df.peptide_rank = ifelse.(coalesce.(pepstats_df.psm_qvalue, 1.0) .<= 1E-3, ifelse(isptm, 1, 2), -1)
    pep2prot_df = innerjoin(pep2prot_df, pepstats_df, on=:peptide_seq)
    pep2prot_df.peptide_seq = replace.(pep2prot_df.peptide_seq, Ref('_' => ""))
    select!(pep2prot_df, [:protein_ac, :peptide_seq, :peptide_rank])
    unique!(pep2prot_df)
    pep2prot_df[!, :dataset] .= dsname
    @info "  $(nrow(pep2prot_df)) associations of $(length(unique(pep2prot_df.peptide_seq))) peptide(s) to $(length(unique(pep2prot_df.protein_ac))) protein(s)"
    pep2prot_df
end for (dsname, dsdata) in pairs(msglm_rdata)), [:protein_ac, :peptide_seq]),
    :peptide_rank => (x -> any(>(0), x) ? minimum(filter(>(0), x)) : -1) => :peptide_rank,
    :dataset => (x -> join(sort(x), ' ')) => :datasets)

peptide2proteins = Dict(df.peptide_seq[1] => (Set{String}(df.protein_ac), df.peptide_rank[1])
                        for df in groupby(peptide2protein_df, :peptide_seq))

proteins_df = reduce(vcat, begin
    select(rdata.full_data["msdata_full"]["proteins"], [:protein_ac, :gene_name, :protein_code, :protein_name, :protein_existence, :src_db, :is_contaminant, :is_viral, :organism], copycols=false)
end for (dsname, rdata) in msglm_rdata) |> unique!
proteins_df.protein_ac_noiso = Fasta.strip_uniprot_isoform.(proteins_df.protein_ac)
for noiso_group_df in groupby(proteins_df, :protein_ac_noiso) # extrapolate best PE across isoforms
    pe = noiso_group_df.protein_existence
    if any(ismissing, pe) && any(!ismissing, pe)
        noiso_group_df.protein_existence .= coalesce.(pe, minimum(skipmissing(pe)))
    elseif any(ismissing, pe) && any(noiso_group_df.is_viral)
        noiso_group_df.protein_existence .= 1
    end
end
# calculate AC ranks - the more canonical proptein is, the better
protein_ac_ranks = Dict(r.protein_ac => ProtgroupXMatch.rank_uniprot(r) for r in eachrow(proteins_df))

ptm_protgroups = ProtgroupAssembly.assemble_protgroups(peptide2proteins, verbose=true,
                                                       nspec_peptides=2, rspec_peptides=0.25)
ptm_protgroup2protein_df = reduce(vcat, [
    DataFrame(protgroup_id = i-1,
              protein_ac = collect(prg.major_prots))
    for (i, prg) in enumerate(ptm_protgroups)])
ptm_protgroups_df = ProtgroupAssembly.dataframe(ptm_protgroups, protein_ranks=protein_ac_ranks,
                            protgroup_col=:protgroup_id, protein_col=:protein_ac,
                            proteins_info=proteins_df)
rename!(ptm_protgroups_df, [:gene_name => :gene_names, :protein_name => :protein_names])

using CodecZlib
ptmn_dfs = Dict(dsname => begin
    open(GzipDecompressorStream, joinpath(data_path, dsinfo.folder, "ptm_extractor_$(dsinfo.fit_ver)", "ptmns_grouped.txt.gz"), "r") do io
        CSV.read(io, delim='\t')
    end
end for (dsname, dsinfo) in datasets if !isnothing(dsinfo.ptm))

ptmgroup2protgroup_dfs = Dict(dsname => begin
    msdata_full = msglm_rdata[dsname].full_data["msdata_full"]
    obj_contrasts_df = msglm_rdata[dsname].fit["object_contrasts.df"]
    ptm2ptmgroup_df = semijoin(select(ptmn_dfs[dsname], [:ptmgroup_id, :ptm_id], copycols=false),
                               select(obj_contrasts_df, [:ptmgroup_id], copycols=false), on=:ptmgroup_id)
    ptmgroup2proteinac_df = innerjoin(ptm2ptmgroup_df, msdata_full["ptm2protein"], on=:ptm_id)
    select!(innerjoin(ptmgroup2proteinac_df, ptm_protgroup2protein_df, on=:protein_ac),
            [:ptmgroup_id, :protgroup_id]) |> unique!
end for (dsname, dsinfo) in datasets if !isnothing(dsinfo.ptm))

for (dsname, df) in ptmgroup2protgroup_dfs
    dsinfo = datasets[dsname]
    open(GzipCompressorStream, joinpath(data_path, dsinfo.folder, "ptm2protgroup_$(dsinfo.fit_ver).txt.gz"), "w") do io
        CSV.write(io, df, delim='\t')
    end
end

fp_protgroups_df = rename!(copy(msglm_rdata[:fp].full_data["msdata_full"]["protregroups"]),
                           :protregroup_id => :protgroup_id)

rnaseq_rdata = load(joinpath(scratch_path, "$(proj_info.id)_$(datasets[:rnaseq].folder)_$(datasets[:rnaseq].fit_ver).RData"), convert=true)
rnaseq_contrasts_df = copy(rnaseq_rdata["object_contrasts.df"], copycols=false)
                            
rnaseq_genename2proteinac_df = rename!(copy(rnaseq_rdata["genename2proteinac.df"], copycols=false),
                                       :protein_ac => :protein_ac_noiso, :GeneName => :gene_name)
rnaseq_genename2proteinac_ex_df = semijoin(innerjoin(
        filter(r -> r.genename_source == "NCBI_Symbol" || !r.has_genesymbol, rnaseq_genename2proteinac_df),
        select(proteins_df, [:protein_ac, :protein_ac_noiso], copycols=false), on=:protein_ac_noiso),
        select(rnaseq_contrasts_df, :gene_name), on=:gene_name)
rnaseq_protgroups_df = combine(groupby(filter!(r -> !ismissing(r.protein_ac),
                unique!(select(rnaseq_genename2proteinac_ex_df, [:gene_name, :protein_ac]))), :gene_name)) do df
    res = DataFrame(protein_acs = join(sort!(unique(df.protein_ac)), ";"),
                    gene_names = df.gene_name[1:1])
    res.majority_protein_acs = copy(res.protein_acs)
    res[!, :is_contaminant] .= false
    res[!, :is_reverse] .= false
    res.organism = missings(String, nrow(res))
    return res
end
rnaseq_protgroups_df.protgroup_id = 1:nrow(rnaseq_protgroups_df)

protgroups_dfs = Dict(
    :rnaseq => rnaseq_protgroups_df,
    :ptm => ptm_protgroups_df,
    :fp => fp_protgroups_df,
)

pg_matches_df = ProtgroupXMatch.match_protgroups(collect(pairs(protgroups_dfs)), protein_ac_ranks);
objid_col = :protgroup_id_united;
pg_matches_df.protgroup_id_united = 1:nrow(pg_matches_df)

# add protgroup ids of each dataset and add protgroup_id_common to each dataset
for (dsname, pg_df) in pairs(protgroups_dfs)
    rowix_col = Symbol("rowix_", dsname)
    if !hasproperty(pg_matches_df, rowix_col)
        @warn "Dataset $dsname key does not exist, skipping"
        continue
    end
    pg_col = Symbol("protgroup_id_", dsname)
    pg_matches_df[!, pg_col] = missings(nonmissingtype(eltype(pg_df.protgroup_id)), nrow(pg_matches_df))
    pg_df[!, :protgroup_id_united] = missings(Int, nrow(pg_df))
    for (i, rowix) in enumerate(pg_matches_df[!, rowix_col])
        if !ismissing(rowix)
            pg_matches_df[i, pg_col] = pg_df[rowix, :protgroup_id]
            pg_df[rowix, :protgroup_id_united] = pg_matches_df[i, :protgroup_id_united]
        end
    end
    select!(pg_matches_df, Not(rowix_col))
end
pg_matches_long_df = FrameUtils.pivot_longer(pg_matches_df, [:protgroup_id_united, :protein_acs, :majority_protein_acs, :gene_names,
                                             #=:is_contaminant, :is_reverse=#],
                        measure_vars_regex=r"^(?<value>rowix|protgroup_id|pgrank|acrank)_(?<var>[^u].+)$",
                        var_col=:dataset)
pg_matches_long_df.dataset = Symbol.(pg_matches_long_df.dataset)
pg_matches_long_expanded_df = DelimDataUtils.expand_delim_column(
    pg_matches_long_df, key_col=[:dataset, :protgroup_id_united, :protgroup_id],
    list_col=:protein_acs, elem_col=:protein_ac)
length(union(
    unique(filter(r -> occursin(r"SARS_.+_vs_mock", r.contrast) && r.std_type == "median" && r.change in ["+", "-"],
              innerjoin(msglm_rdata[:cov2ts_proteome].fit["object_contrasts.df"],
                        select(pg_matches_df, [:protgroup_id_united, :protgroup_id_cov2ts_proteome]),
                        on=[:protgroup_id => :protgroup_id_cov2ts_proteome])).protgroup_id_united),
    unique(filter(r -> occursin(r"SARS_.+_vs_mock", r.contrast) && r.std_type == "median" && r.change in ["+", "-"],
              innerjoin(msglm_rdata[:cov2el_proteome].fit["object_contrasts.df"],
                        select(pg_matches_df, [:protgroup_id_united, :protgroup_id_cov2el_proteome]),
                        on=[:protregroup_id => :protgroup_id_cov2el_proteome])).protgroup_id_united)
    ))


obj2protac_df = unique!(dropmissing!(select(pg_matches_long_expanded_df,
                                            [objid_col, :protein_ac])))

nonunique_matches_dfs = [dsname => begin
    pg_col = Symbol("protgroup_id_", dsname)
    combine(groupby(filter(r -> !ismissing(r[pg_col]), pg_matches_df), pg_col)) do df
        return nrow(df) > 1 ? df : df[1:0, :]
    end
end for (dsname, _) in pairs(protgroups_dfs)]

contrast_cols = [:contrast, :timepoint_lhs, :timepoint_rhs, :treatment_lhs, :treatment_rhs]

obj_contrasts_msglm_dfs = Dict(ds => begin
    @info "Processing $ds results..."
    orig_objcontrasts_df = rdata.fit["object_contrasts.df"]
    sel_cols = [:std_type; contrast_cols;
        [:object_id,
        :median_log2, :p_value,
        :is_hit_composed, :is_hit, :change, :is_hit_nomschecks, :is_signif]]
    for col in [:ptm_type, :ptmgroup_id, :ptmn_id, :protregroup_id]
        hasproperty(orig_objcontrasts_df, col) && push!(sel_cols, col)
    end
    objcontrasts_df = select(orig_objcontrasts_df, sel_cols)
    objcontrasts_df[!, :dataset] .= ds
    objcontrasts_df.change = ifelse.(objcontrasts_df.is_hit_composed,
            ifelse.(objcontrasts_df.median_log2 .> 0, "+", "-"), ".")
    @show nrow(objcontrasts_df)
    isptm = !isnothing(datasets[ds].ptm)
    if isptm
        @info "  Associating PTM groups to protein groups..."
        objcontrasts_df = innerjoin(objcontrasts_df,
                                    ptmgroup2protgroup_dfs[ds], on=:ptmgroup_id)
        @assert hasproperty(objcontrasts_df, :ptm_type)
        #objcontrasts_df.ptm_type = datasets[ds].ptm
    else
        if hasproperty(objcontrasts_df, :protregroup_id) &&
          !hasproperty(objcontrasts_df, :protgroup_id)
            rename!(objcontrasts_df, :protregroup_id => :protgroup_id)
        end
        objcontrasts_df.ptmgroup_id = missings(Int, nrow(objcontrasts_df))
        objcontrasts_df.ptmn_id = missings(Int, nrow(objcontrasts_df))
        objcontrasts_df.ptm_type = missings(String, nrow(objcontrasts_df))
    end
    @show nrow(objcontrasts_df)
    @info "  Associating dataset-specific protein groups to united protein groups..."
    objcontrasts_df = leftjoin(objcontrasts_df,
                               select(protgroups_dfs[isptm ? :ptm : ds], [:protgroup_id, :protgroup_id_united], copycols=false),
                               on=:protgroup_id)
    @show nrow(objcontrasts_df)
    objcontrasts_df
end for (ds, rdata) in pairs(msglm_rdata))

sel_std_type = "median"
# in FP one protgroup could belong to a single +/-/. group
obj_contrasts_fp_agg_df = combine(groupby(filter(r -> (r.std_type == sel_std_type) && !ismissing(r.protgroup_id_united), obj_contrasts_msglm_dfs[:fp]),
   [[:std_type, :dataset, :ptm_type]; contrast_cols; :protgroup_id_united])
) do df
    min_ix = findmin(df.p_value)[2]
    return df[min_ix:min_ix, :]
end

obj_contrasts_ptm_agg_df = combine(groupby(filter!(r -> (r.std_type == sel_std_type) && !ismissing(r.protgroup_id_united),
                                                   vcat(obj_contrasts_msglm_dfs[:ubi], obj_contrasts_msglm_dfs[:phospho])),
   [:std_type, :dataset, :ptm_type,
    :contrast, :timepoint_lhs, :timepoint_rhs, :treatment_lhs, :treatment_rhs,
    :protgroup_id_united, :change])
) do df
    return DataFrame(is_hit_composed = any(df.is_hit_composed),
                     nrows=nrow(df),
                     nobjs = length(unique(df.object_id)),
                     nobj_hits = length(unique(df.object_id[coalesce.(df.is_hit_composed, false)])))
end

obj_contrasts_rnaseq_agg_df = combine(groupby(innerjoin(rnaseq_contrasts_df, protgroups_dfs[:rnaseq],
                                             on=[:gene_name]),
                        [contrast_cols; :protgroup_id_united])) do df
   min_ix = findmin(df.p_value)[2]
   return df[min_ix:min_ix, :]
end
obj_contrasts_rnaseq_agg_df[!, :dataset] .= :rnaseq

united_cols = [:dataset; contrast_cols; :protgroup_id_united; :is_hit_composed; :change]
obj_contrasts_united_df = vcat(select(obj_contrasts_fp_agg_df, united_cols, copycols=false),
                               select(obj_contrasts_ptm_agg_df, united_cols, copycols=false),
                               select(obj_contrasts_rnaseq_agg_df, united_cols, copycols=false))
filter!(r -> r.treatment_rhs == "mock" && r.treatment_lhs != "infected" && r.timepoint_lhs == r.timepoint_rhs,
        obj_contrasts_united_df)
countmap(collect(zip(obj_contrasts_united_df.is_hit_composed, obj_contrasts_united_df.change)))

contrasts_df = unique!(select(obj_contrasts_united_df, [:dataset; contrast_cols; :change]))
contrasts_df.dataset = String.(contrasts_df.dataset)
contrasts_df.timepoint_lhs = parse.(Int, contrasts_df.timepoint_lhs)
contrasts_df.timepoint_rhs = parse.(Int, contrasts_df.timepoint_rhs)
contrasts_df.change_alt = getindex.(Ref(Dict("+" => "▲", "-" => "▼", "." => ".")),
                                        contrasts_df.change)

includet(joinpath(misc_scripts_path, "optcover_utils.jl"));
includet(joinpath(misc_scripts_path, "gmt_reader.jl"));
includet(joinpath(misc_scripts_path, "omics_collections.jl"));

@info "Loading Human annotations..."
# human mappings from http://download.baderlab.org/EM_Genesets/December_01_2018/Human/UniProt/
genesets_df, genesets_coll = GMT.read(String,
        joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_October_01_2020_UniProt.gmt"),
        id_col = :term_id, src_col = :term_src);

pcomplexes_df, pcomplex_iactors_df, pcomplex_iactor2ac_df =
    OmicsCollections.ppicollection(joinpath(party3rd_data_path, "complexes_20191217.RData"), seqdb=:uniprot);
pcomplexes_df[!, :coll_id] .= "protein_complexes";

# make complexes collections, keep complexes with at least 2 participants
pcomplex_coll = FrameUtils.frame2collection(innerjoin(pcomplex_iactors_df, pcomplex_iactor2ac_df,
    on=[:file, :entry_index, :interaction_id, :interactor_id]),
            set_col=:complex_id, obj_col=:protein_ac, min_size=2)
protac_sets = merge(genesets_coll, pcomplex_coll)

terms_df = vcat(rename(genesets_df[!, [:term_src, :term_id, :name, :descr]],
                       :term_src => :coll_id, :name=>:term_name, :descr=>:term_descr),
                #rename(goterm_info_df[[:id, :name, :def]], :onto => :coll_id, :id=>:term_id, :name=>:term_name, :def=>:term_descr),
                rename(pcomplexes_df[!, [:coll_id, :complex_id, :interaction_label, :interaction_name]],
                       :complex_id=>:term_id, :interaction_label=>:term_name, :interaction_name=>:term_descr));
protac2term_df = FrameUtils.collection2frame(protac_sets, terms_df,
                                             setid_col=:term_id, objid_col=:protein_ac)

obj2term_df = select!(innerjoin(obj2protac_df, protac2term_df, on = :protein_ac),
                      Not([:protein_ac])) |> unique!
protac_colls = FrameUtils.frame2collections(protac2term_df, obj_col=:protein_ac,
                                            set_col=:term_id, coll_col=:coll_id)
obj_colls = FrameUtils.frame2collections(obj2term_df, obj_col=objid_col,
                                         set_col=:term_id, coll_col=:coll_id)

@info "Preparing hit sets"
ObjectType = eltype(obj2protac_df[!, objid_col])
obj_hit_sets = Dict{Tuple{String, String, String}, Set{ObjectType}}()
for hits_df in groupby(filter(r -> coalesce(r.is_hit_composed, false), obj_contrasts_united_df),
                       [:dataset, :contrast, :change])
    obj_hit_sets[(string(hits_df[1, :dataset]),
                  hits_df[1, :contrast], hits_df[1, :change])] =
        Set(skipmissing(hits_df[!, objid_col]))
end

# only relevant ones
obj_hit_selsets = obj_hit_sets

@info "Preparing mosaics..."
observed_protacs = Set(obj2protac_df.protein_ac) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls,
                                      protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);

# remove broad terms larger than 150 elements (too unspecific)
obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200, verbose=true)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

ENV["MKL_NUM_THREADS"] = 1

obj_hit_mosaics_v = collect(pairs(obj_hit_mosaics))
obj_hit_covers_v = similar(obj_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_hit_mosaics_v[i]
    @info "Covering $mosaic_name by hits..."
    obj_hit_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=1,#Threads.nthreads()-1,
                                    MaxRestarts=200),
            true)
end
obj_hit_covers = Dict(k => v for (k, v) in obj_hit_covers_v)

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_united_hit_covers_$(proj_info.oesc_ver).jld2")
@save(hit_covers_filename,
      proj_info, datasets, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      #objects_df,
      protgroups_dfs, pg_matches_df,
      #obj_effects_df,
      obj_contrasts_united_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
if !@isdefined(obj_hit_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, datasets, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      #objects_df,
      protgroups_dfs, pg_matches_df,
      #obj_effects_df,
      obj_contrasts_united_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
end

includet(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup↦gene_name map..."
obj_id2name = Dict(r.protgroup_id_united => !ismissing(r.gene_names) ? DelimDataUtils.rejoin_unique_substrings([r.gene_names], ";") :
                    !ismissing(r.majority_protein_acs) ? DelimDataUtils.rejoin_unique_substrings([r.majority_protein_acs], ";") :
                    string("PGU_", r.protgroup_id_united)
                   for r in eachrow(filter(r -> !ismissing(r.protgroup_id_united), pg_matches_df)))

obj_hit_covers_df = innerjoin(
    OptCoverUtils.covers_report(
    obj_hit_covers, obj_hit_selsets, obj_colls, obj_mosaics, obj_id2name,
    terms_df,
    maskid_col=[:dataset, :contrast, :change],
    maskedset_col_prefix="contrast"),
    contrasts_df, on=[:dataset, :contrast, :change])

# don't remove the sets since they are timecourses timepoints
obj_hit_covers_signif_df = combine(groupby(obj_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:dataset, :contrast, :change],
                        max_term_pvalue=1E-4, max_set_pvalue=nothing, min_set_overlap=nothing),
                    Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_united_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_df[obj_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_united_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_signif_df[obj_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

function xlsx_compatible(df::AbstractDataFrame)
    res = copy(df, copycols=false)
    for (coln, col) in zip(names(res), eachcol(res))
        if nonmissingtype(eltype(col)) <: Symbol || nonmissingtype(eltype(col)) <: CategoricalValue
            res[!, coln] = String.(col)
        end
    end
    return res
end

obj_hit_covers_4xl_df = xlsx_compatible(obj_hit_covers_df[obj_hit_covers_df.nmasked .> 0, :])
using XLSX
XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_united_hit_covers_$(proj_info.oesc_ver).xlsx"),
                collect(DataFrames.eachcol(obj_hit_covers_4xl_df)),
                DataFrames.names(obj_hit_covers_4xl_df), overwrite=true)

obj_hit_covers_signif_4xl_df = xlsx_compatible(obj_hit_covers_signif_df[obj_hit_covers_signif_df.nmasked .> 0, :])
XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_united_hit_covers_signif_$(proj_info.oesc_ver).xlsx"),
                collect(DataFrames.eachcol(obj_hit_covers_signif_4xl_df)),
                DataFrames.names(obj_hit_covers_signif_4xl_df), overwrite=true)

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_plots.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

using PlotlyJS, TextWrap

heatmap_layout_attrs = Dict(
    ("GO_CC", true) => Dict(:margin_l => 200),
    ("GO_CC", false) => Dict(:margin_l => 200),
    ("Reactome", true) => Dict(:margin_l => 350),
)

#=
for (plot_mosaic, cover_coll) in obj_hit_covers
    isempty(cover_coll.results) && continue
    @info "Plotting $plot_mosaic Pareto front"
    paretofront_plot = OptCoverPlots.plot_paretofront(cover_coll.results[1], plot_unfolded=true)
    plot_filename = joinpath(plots_path, "oesc_$(sel_std_type)_std", "paretofront",
                             "$(proj_info.id)_$(plot_mosaic)_X_treatment_$(sel_std_type)_pareto")
    savefig(paretofront_plot.plot, "$plot_filename.svg")
    PlotlyJS.savehtml(paretofront_plot, "$plot_filename.html")
end
=#

stylize_dataset(ds) =
    "<span style=\"font-weight: bold; color: $(datasets[Symbol(ds)].color);\">" * datasets[Symbol(ds)].label * "</span>"

stylize_contrast(str) = foldl(replace, [
    r"(SARS_CoV2?)@(\d+)h_vs_mock@(\d+)h" => s"\1:<span style=\"font-weight: bold; color: black;\">\2</span>h",
    "SARS_CoV:" => "<span style=\"font-wieght: bold; color: #811A02;\">SARS</span> ",
    "SARS_CoV2:" => "<span style=\"font-wieght: bold; color: #F4982A;\">CoV2</span> ",
    ],
    init = str)

function process_contrast_axis(contrast_df)
    contrast_df,
    stylize_dataset.(contrast_df.dataset) .* ":&nbsp;" .*
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change),
    stylize_dataset.(contrast_df.dataset) .* ":&nbsp;" .*
        stylize_contrast.(contrast_df.contrast) .*
        " " .* OptCoverHeatmap.stylize_change.(contrast_df.change)#stylize_effect.(effect_df.effect)
end

datatype_order = Dict(:rnaseq => 1, :proteome => 2, :ubi => 3, :phospho => 4)

heatmaps_path = joinpath(plots_path, "united_hits_oesc_$(sel_std_type)_$(proj_info.oesc_ver)")
isdir(heatmaps_path) || mkdir(heatmaps_path)

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")hit heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_hit_covers_signif_df : obj_hit_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end

    for outformat in ("html", "pdf", "svg")
        coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            mask_axis_title = "contrast",
            mask_cols = [:dataset, :contrast, :treatment_lhs, :timepoint_lhs, :change, :ncontrast],
            process_mask_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 150),
            transpose=false,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",#outformat in ["svg", "pdf"] ? "#000" : "#BBB",
            cell_width=25, cell_height=20, gridwidth=1,
            mask_order=contrasts -> begin
                contrasts.datatype_order = [datatype_order[datasets[Symbol(r.dataset)].datatype]
                                            for r in eachrow(contrasts)]
                return sortperm(contrasts, [:datatype_order, :change, :timepoint_lhs, :treatment_lhs, :dataset])
            end)
        (coll_heatmap === nothing) && continue
        for (k, v) in [#:width=>800, :height=>400,
                       :margin_r=>130, :margin_t=>20,
                       :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
            coll_heatmap.plot.layout[k] = v
        end
        plot_fname = joinpath(heatmaps_path,
                              "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_contrast$(signif ? "_signif" : "")_heatmap.$(outformat)")
        if outformat == "html"
            PlotlyJS.savehtml(coll_heatmap, plot_fname, :embed);
        else
            savefig(coll_heatmap, plot_fname, width=coll_heatmap.plot.layout[:width], height=coll_heatmap.plot.layout[:height]);
        end
    end
end

using VegaLite
p = perseus_tests_df |>
@vlplot(
    :bar,
    transform=[
        {filter="datum.is_signif"},
    ],
    column="timepoint:o",
    y={"sum(is_signif)", axis={title="N regulated", grid=false}},
    x={"change:n", axis={title=""}},
    color={"change:n", scale={range=["firebrick", "dodgerblue"]}},
    spacing=10,
    config={
        view={stroke=:transparent},
        axis={domainWidth=1, labelAngle=0}
    }
)
save(joinpath(plots_path, "test.pdf"), p)

triheatmaps_path = joinpath(plots_path, "united_hits_oesc_$(sel_std_type)_$(proj_info.oesc_ver)_tri")
isdir(triheatmaps_path) || mkdir(triheatmaps_path)

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")tri-heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    cover_df = filter(r -> r.term_collection == term_coll, signif ? obj_hit_covers_signif_df : obj_hit_covers_df)
    if nrow(cover_df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end
    scores_mtx, tips_mtx, term_axis, contrast_axis = OptCoverHeatmap.heatmap_matrices(cover_df,
            elements_label="protein",
            mask_axis_title = "contrast",
            mask_cols = [:dataset, :contrast, :treatment_lhs, :timepoint_lhs, :change, :change_alt, :ncontrast],
            process_mask_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            mask_order=contrasts -> begin
                contrasts.datatype_order = [datatype_order[datasets[Symbol(r.dataset)].datatype]
                                            for r in eachrow(contrasts)]
                return sortperm(contrasts, [:datatype_order, :change, :timepoint_lhs, :treatment_lhs, :dataset])
        end)
    rename!(contrast_axis, :axis_label => :axisXtri_label, :axis_tip => :axisXtri_tip)
    contrast_axis.axis_label = [datasets[Symbol(r.dataset)].label * ": " * string(r.timepoint_lhs) * "h " * r.change_alt for r in eachrow(contrast_axis)]

    triheatmap_df, triheatmap_rows_df, triheatmap_cols_df = subheatmap_frame(scores_mtx, term_axis, contrast_axis, #tips=tips_mtx,
        col_sub_cols = [:treatment_lhs, :contrast, :ncontrast, :axisXtri_label, :axisXtri_tip])
    triheatmap_plot = vegalite_subheatmap(triheatmap_df,
                                          xaxis_label="Time, h.p.i", yaxis_label="Term",
                                          subaxis_label="Virus", coloraxis_label="log₁₀(P-value)",
                                          sublabels=["SARS_CoV2", "SARS_CoV"],
                                          labelLimit_l=get(layout_attrs, :margin_l, 400))

    plot_fname = joinpath(triheatmaps_path,
        "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_contrast$(signif ? "_signif" : "")_triheatmap")
    save(plot_fname * ".html", triheatmap_plot)
    save(plot_fname * ".svg", triheatmap_plot)
    save(plot_fname * ".pdf", triheatmap_plot)
end

### global stats

infection_stats_df = combine(groupby(vcat([begin
    df = filter(r -> r.std_type == "median" && r.timepoint_lhs == r.timepoint_rhs &&
                     r.treatment_rhs == "mock" && r.treatment_lhs != "infected",
                msglm_rdata[dsname].fit["object_contrasts.df"])
    df[!, :dataset] .= dsname
    df[!, :entity] .= datasets[dsname].entity
    df.entity_id = datasets[dsname].entity == "PTM" ? df.ptmgroup_id : df.object_id
    rename!(df, :timepoint_lhs=>:timepoint, :treatment_lhs=>:treatment)
    select!(df, [:entity, :dataset, :treatment, :timepoint, :composed_hit_type, :change, :is_viral, :is_contaminant, :entity_id])
end for dsname in [:ubi, :phospho, :fp]]..., begin
    df = filter(r -> r.timepoint_lhs == r.timepoint_rhs &&
                     r.treatment_rhs == "mock" && r.treatment_lhs != "infected",
                rnaseq_rdata["object_contrasts.df"])
    df[!, :dataset] .= :rnaseq
    df[!, :entity] .= datasets[:rnaseq].entity
    df.entity_id = df.gene_name
    df[!, :is_viral] .= false
    df[!, :is_contaminant] .= false
    rename!(df, :timepoint_lhs=>:timepoint, :treatment_lhs=>:treatment)
    select!(df, [:entity, :dataset, :treatment, :timepoint, :composed_hit_type, :change, :is_viral, :is_contaminant, :entity_id])
end),
    [:entity, :dataset, :treatment, :timepoint, :composed_hit_type, :change, :is_viral, :is_contaminant]),
    :entity_id => (x -> length(unique(x))) => :n_regulated)
infection_stats_df.composed_hit_type = levels!(categorical(infection_stats_df.composed_hit_type),
                                               ["none", "shared", "SARS_CoV2", "SARS_CoV"])
infection_stats_df.timepoint = parse.(Int, infection_stats_df.timepoint)
infection_stats_df.dataset_label = [datasets[r.dataset].label for r in eachrow(infection_stats_df)]
infection_stats_df[!, :n_regulated_stacked] .= 0
sort!(infection_stats_df, [:entity, :dataset, :treatment, :timepoint, :is_viral, :is_contaminant, :change, :composed_hit_type])
for df in groupby(infection_stats_df, [:entity, :dataset, :treatment, :timepoint, :change, :is_viral, :is_contaminant])
    cumsum!(view(df, :, :n_regulated_stacked), df.n_regulated)
end
infection_stats_df.n_regulated_stacked_prev = infection_stats_df.n_regulated_stacked .- infection_stats_df.n_regulated
for col in ["n_regulated", "n_regulated_stacked", "n_regulated_stacked_prev"]
    infection_stats_df[!, col*"_vis"] = [ifelse(ch == "+", n, -n)
            for (ch, n) in zip(infection_stats_df.change, infection_stats_df[!, col])]
end
infection_stats_df.composed_hit_type_str = get.(infection_stats_df.composed_hit_type)

infection_stats_wide_df = FrameUtils.unstack(infection_stats_df, [:treatment, :timepoint, :composed_hit_type, :is_contaminant],
                            :dataset, [:n_regulated, :n_regulated_vis, :n_regulated_stacked_vis, :n_regulated_stacked_prev_vis])

select(infection_stats_df, Not(:composed_hit_type), copycols=false) |>
@vlplot(
    transform=[
        {filter="(datum.change != '.') & !datum.is_contaminant & !datum.is_viral"},
    ],
    vconcat=[{
        transform=[{
            filter="datum.dataset=='fp'"
        }],
        y={"n_regulated_stacked_prev_vis:q",
            axis={title="", grid=false}},
        y2={"n_regulated_stacked_vis:q", axis={grid=false}},
        x={"treatment:n", axis={title="", labelAngle=45},
            scale={domain=["SARS_CoV2", "SARS_CoV"]}},
        color={"composed_hit_type_str:o", scale={
            domain=["none", "shared", "SARS_CoV2", "SARS_CoV"],
            range=["lightgray", "gray", "#F4982A", "#811A02"]},
            axis={title="Regulation"}},
        config={
            view={stroke=:transparent},
            axis={domainWidth=1, labelAngle=0},
            bar={binSpacing=0},
        },
        layer=[{
            column={"timepoint:o", align="each"},
            mark={:bar},
        }, {
            column={"timepoint:o", align="each"},
            mark={
                :text,
            },
            text={"n_regulated:q"},
        }],
    },{
        transform=[{
            filter="datum.dataset=='rnaseq'"
        }],
        column={"timepoint:o", align="each"},
        row={"dataset_label:n", axis={title=""}},
        y={"n_regulated_stacked_prev_vis:q",
            axis={title="", grid=false}},
        y2={"n_regulated_stacked_vis:q", axis={grid=false}},
        x={"treatment:n", axis={title="", labelAngle=45},
            scale={domain=["SARS_CoV2", "SARS_CoV"]}},
        color={"composed_hit_type_str:o", scale={
            domain=["none", "shared", "SARS_CoV2", "SARS_CoV"],
            range=["lightgray", "gray", "#F4982A", "#811A02"]},
            axis={title="Regulation"}},
        config={
            view={stroke=:transparent},
            axis={domainWidth=1, labelAngle=0},
            bar={binSpacing=0},
        },
        layer=[{
            mark={:bar},
        }, {
            mark={
                :text,
            },
            text={"n_regulated:q"},
        }],
    }],
    spacing=2,
    width=100,
    height=200
)

mutate(dataset_label = sapply(dataset, function(ds) datasets[[ds]]$label),
             dataset_type = factor(sapply(dataset, function(ds) datasets[[ds]]$type),
                                   levels = datatype_order),
             protocol = factor(sapply(dataset, function(ds) datasets[[ds]]$protocol)),
             n_regulated_vis = n_regulated * if_else(change == "+", 1, -1),
             timepoint_num = parse_integer(timepoint),
             timepoint = factor(timepoint_num, ordered = TRUE)) %>%
    dplyr::arrange(dataset_type, timepoint, treatment, dataset) %>%
    dplyr::mutate(data_slice = str_c(dataset_label, "@", timepoint, "h"),
                  data_slice = factor(data_slice, levels=unique(data_slice)))



plots_path
using RCall
R"require(ggplot2)"
R"require(dplyr)"
R"require(stringr)"
R"require(tidyr)"
R"require(readr)"
R"source('~/projects/R/misc/ggplot_ext.R')"

length(unique(rnaseq_contrasts_df.GeneName))
length(unique(filter(r -> r.timepoint_lhs == r.timepoint_rhs && r.change in ["+", "-"], rnaseq_contrasts_df).GeneName))

countmap(perseus_analysis_orig_df.comparison)

nrow(unique(select(filter(r -> r.dataset == :cov2el_phospho,
                          perseus_analysis_orig_df),
                     [:lead_protein_acs, :ptm_poses])))
nrow(unique(select(filter(r -> r.dataset == :cov2el_phospho &&
                          occursin(r"norm_SARS_.+_ph_norm_mock_.+", r.comparison) &&
                          r.change in ["+" , "-"],
                          perseus_analysis_orig_df),
                     [:lead_protein_acs, :ptm_poses])))

prefix = "cov2el_ubi"
R"stats.df <- filter($perseus_analysis_orig_df, change %in% c('+', '-') & dataset == $prefix) %>%
             group_by(contrast, timepoint, change) %>%
             summarise(n_regulated = if_else(change[[1]] == '+', n(), -n())) %>% ungroup()"
R"stats.df <- filter($perseus_analysis_orig_df, change %in% c('+', '-') & dataset == $prefix) %>%
             group_by(contrast, timepoint) %>%
             summarise(n_regulated = n_distinct(lead_protein_acs)) %>% ungroup()"

prefix = "cov2ts_phospho"
obj_contrasts_df = msglm_rdata[Symbol(prefix)].fit["object_contrasts.df"]
R"stats.df <- filter($obj_contrasts_df, change %in% c('+', '-') & std_type == 'median' &
                     str_detect(contrast, 'SARS_.+_vs_mock')) %>%
              extract(contrast, c('timepoint_lhs', 'timepoint_rhs'), 'SARS_COV2@(\\d+)h_vs_mock@(\\d+)h', remove=FALSE) %>%
              mutate(timepoint = parse_integer(timepoint_lhs)) %>%
              group_by(contrast, contrast, timepoint, change) %>%
              summarise(n_regulated = if_else(change[[1]] == '+', n_distinct(object_id), -n_distinct(object_id))) %>%
              ungroup()"
R"stats.df <- filter($obj_contrasts_df, change %in% c('+', '-') & std_type == 'median' &
                     str_detect(contrast, 'SARS_.+_vs_mock')) %>%
              extract(contrast, c('timepoint_lhs', 'timepoint_rhs'), 'SARS_COV2@(\\d+)h_vs_mock@(\\d+)h', remove=FALSE) %>%
              mutate(timepoint = parse_integer(timepoint_lhs)) %>%
              group_by(contrast, contrast, timepoint) %>%
              summarise(n_regulated = n_distinct(majority_protein_acs)) %>%
              ungroup()"

prefix = "cov2ts_rnaseq"
R"stats.df <- filter($rnaseq_contrasts_df, change %in% c('+', '-') &
             treatment_lhs != treatment_rhs & timepoint_lhs == timepoint_rhs) %>%
             mutate(timepoint = factor(timepoint_lhs)) %>%
             group_by(contrast, timepoint, change) %>%
             summarise(n_regulated = if_else(change[[1]] == '+', n(), -n())) %>%
             ungroup()"

R"stats.df <- filter($obj_contrasts_df, treatment_lhs != treatment_rhs & timepoint_lhs == timepoint_rhs &
                std_type=='median' & change %in% c('+', '-')) %>%
             rename(timepoint = timepoint_lhs) %>%
             group_by(contrast, timepoint, change) %>% summarise(n_regulated = n())"

R"print(stats.df)"
R"p = ggplot(stats.df,
         aes(x = factor(timepoint), y=n_regulated)) +
geom_bar(aes(), fill='lightgray', color='darkgray', stat='identity', position='identity', show.legend=FALSE) +
geom_text(data=filter(stats.df, n_regulated>0), aes(label=n_regulated),
          size=6, vjust=-0.8, show.legend=FALSE) +
geom_text(data=filter(stats.df, n_regulated<0), aes(label=-(n_regulated)),
          size=6, vjust=1.5, show.legend=FALSE) +
scale_x_discrete('Timepoint') +# breaks=unique(stats.df[['timepoint']])) +
scale_y_continuous('N regulated sites', limits=c(-300, 500)) +
theme_bw_ast(base_size=20) +
theme(panel.grid = element_blank(), axis.line = element_blank())"
R"ggsave(file.path($plots_path, str_c('cov2_', $prefix, '_stats.pdf')), p, width=6, height=4, device=cairo_pdf)"

R"p = ggplot(stats.df,
         aes(x = factor(timepoint), y=n_regulated)) +
geom_bar(aes(), fill='lightgray', color='darkgray', stat='identity', position='identity', show.legend=FALSE) +
geom_text(data=filter(stats.df, n_regulated>0), aes(label=n_regulated),
          size=6, vjust=-0.8, show.legend=FALSE) +
scale_x_discrete('Timepoint') +# breaks=unique(stats.df[['timepoint']])) +
scale_y_continuous('N regulated proteins', limits=c(0, 500)) +
theme_bw_ast(base_size=20) +
theme(panel.grid = element_blank(), axis.line = element_blank())"
R"ggsave(file.path($plots_path, str_c('cov2_', $prefix, '_protein_stats.pdf')), p, width=6, height=4, device=cairo_pdf)"

R"require(stringr)"
R"require(tidyr)"
R"p = ggplot(filter($obj_contrasts_df, treatment_lhs != treatment_rhs & timepoint_lhs == timepoint_rhs &
                std_type=='median' & is_hit_nomschecks) %>%
             rename(timepoint = timepoint_lhs, change_old=change, change=change_alt) %>%
             group_by(contrast, timepoint, change) %>% summarise(n_regulated = n()),
         aes(x = timepoint, y=n_regulated)) +
geom_bar(aes(color=change, fill=change), stat='identity', show.legend=FALSE, position=position_dodge(width=3)) +
geom_text(aes(label=n_regulated, color=change), size=6, vjust=-0.8, show.legend=FALSE, position=position_dodge(width=3)) +
scale_color_manual(values=c('▲' = 'firebrick', '▼' = 'dodgerblue')) +
scale_fill_manual(values=c('▲' = 'firebrick', '▼' = 'dodgerblue')) +
scale_x_continuous(breaks=unique($(obj_contrasts_df.timepoint_lhs))) +
scale_y_continuous(limits = c(0, 500)) +
theme_bw_ast()"
R"ggsave(file.path($plots_path, 'cov2_', $prefix, '_DIA_phospho_stats.pdf'), p, width=6, height=4, device=cairo_pdf)"

