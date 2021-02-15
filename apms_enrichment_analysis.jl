proj_info = (id = "cov2",
             data_ver = "20200525",
             fit_ver = "20200525",
             oesc_ver = "20200604",
             modelobj = "protregroup",
             msfolder = "mq_apms_20200525")
using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, DataFrames
using StatsBase

@info "Project '$(proj_info.id)' dataset $(proj_info.msfolder) version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "msglm_utils.jl"));

objid_col = Symbol(string(proj_info.modelobj, "_id"));

input_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(proj_info.msfolder)_$(proj_info.data_ver).RData"), convert=true)
full_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.msfolder)_$(proj_info.data_ver).RData"), convert=true)
fit_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(proj_info.msfolder)_$(proj_info.fit_ver).RData"), convert=true)
network_rdata = load(joinpath(analysis_path, "networks", "$(proj_info.id)_4graph_$(proj_info.msfolder)_$(proj_info.fit_ver).RData"), convert=true)
contrasts_df = unique!(select!(copy(input_rdata["contrastXmetacondition.df"]), [:contrast, :contrast_type]));
objects_df = copy(input_rdata["msdata"][string(proj_info.modelobj, "s")]) |> MSGLMUtils.fix_object_id!;
protacs_df = copy(full_rdata["msdata_full"]["proteins"]);
obj2protac_df = select!(
         filter(r -> r.is_majority, full_rdata["msdata_full"][string("protein2", proj_info.modelobj)]),
         [objid_col, :protein_ac]) |> unique! |> MSGLMUtils.fix_object_id!; # is_majority?
#obj2protac_df = obj2protac_df[obj2protac_df.is_majority, [:object_id, :protein_ac]];
obj_effects_df = filter(r -> r.var == "obj_effect", fit_rdata["object_effects.df"]);
obj_contrasts_df = copy(fit_rdata["object_contrasts.df"]);
for df in [obj_effects_df, obj_contrasts_df]
    df |> MSGLMUtils.fix_quantile_columns! |> MSGLMUtils.fix_object_id!
end
obj_iactions_df = copy(network_rdata["iactions_4table.df"]);

Revise.includet(joinpath(misc_scripts_path, "gmt_reader.jl"));
Revise.includet(joinpath(misc_scripts_path, "optcover_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "omics_collections.jl"));

@info "Loading Human annotations..."
# human mappings from http://download.baderlab.org/EM_Genesets/December_01_2018/Human/UniProt/
genesets_df, genesets_coll = GMT.read(String,
        joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_April_01_2019_UniProt.gmt"),
        id_col = :term_id, src_col = :term_src);

pcomplexes_df, pcomplex_iactors_df, pcomplex_iactor2ac_df =
    OmicsCollections.ppicollection(joinpath(party3rd_data_path, "complexes_20191217.RData"), seqdb=:uniprot);
pcomplexes_df[!, :coll_id] .= "protein_complexes";
# make complexes collections, keep complexes with at least 2 participants
pcomplex_coll = FrameUtils.frame2collection(innerjoin(pcomplex_iactors_df, pcomplex_iactor2ac_df,
            on=[:file, :entry_index, :interaction_id, :interactor_id]),
            set_col=:complex_id, obj_col=:protein_ac, min_size=2)
protac_sets = merge!(genesets_coll, pcomplex_coll)

terms_df = vcat(rename(genesets_df[!, [:term_src, :term_id, :name, :descr]],
                       :term_src => :coll_id, :name=>:term_name, :descr=>:term_descr),
                #rename(goterm_info_df[[:id, :name, :def]], Dict(:onto => :coll_id, :id=>:term_id, :name=>:term_name, :def=>:term_descr)),
                rename(pcomplexes_df[!, [:coll_id, :complex_id, :interaction_label, :interaction_name]],
                       :complex_id=>:term_id, :interaction_label=>:term_name, :interaction_name=>:term_descr));
protac2term_df = FrameUtils.collection2frame(protac_sets, terms_df,
                                             setid_col=:term_id, objid_col=:protein_ac)

# link protein group IDs to annots and create protgroup collections
obj2term_df = select!(innerjoin(obj2protac_df, protac2term_df, on = :protein_ac),
                      Not([:protein_ac])) |> unique!
protac_colls = FrameUtils.frame2collections(protac2term_df, obj_col=:protein_ac,
                                            set_col=:term_id, coll_col=:coll_id)
obj_colls = FrameUtils.frame2collections(obj2term_df, obj_col=objid_col,
                                         set_col=:term_id, coll_col=:coll_id)

@info "Preparing hit sets"
ObjectType = eltype(obj2protac_df.object_id)
obj_hit_sets = Dict{String, Set{ObjectType}}()
#=
for contrast_df in groupby(filter(r -> coalesce(r.is_hit, false), obj_contrasts_df), [:std_type, :contrast, :change])
    obj_hit_sets[(string.(contrast_df[1, :std_type], "_std"), contrast_df[1, :contrast], contrast_df[1, :change])] = Set(skipmissing(contrast_df.object_id))
end
# leave only true interactors in bait-to-bait comparison
for ((std_type, contrast, change), objs) in obj_contrast_sets
    conditions = match(r"(.+)_vs_(.+)", contrast)
    if !in(conditions[2], ["controls", "others"])
        cond1_objs = obj_hit_sets[(std_type, string(conditions[1], "_vs_controls"), "+")]
        cond2_objs = obj_hit_sets[(std_type, string(conditions[2], "_vs_controls"), "+")]
        new_objs = intersect(objs, union(cond1_objs, cond2_objs))
        @info "$contrast: was $(length(objs)), new $(length(new_objs))"
        obj_hit_sets[(std_type, contrast, change)] = new_objs
    end
end
=#
for bait_iactors_df in groupby(filter(r -> coalesce(r.is_hit, false), obj_iactions_df), [:contrast, :bait_full_id])
    obj_hit_sets[bait_iactors_df[1, :contrast]] = Set(skipmissing(bait_iactors_df.object_id))
end
# only relevant ones
obj_hit_selsets = copy(obj_hit_sets)

@info "Preparing mosaics..."
observed_protacs = Set(obj2protac_df.protein_ac) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);

obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by contrasts..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets,
                                          max_sets=2000, min_nmasked=2, max_setsize=150)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

obj_hit_covers = Dict(begin
    @info "Covering $mosaic_name by contrasts..."
    mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=Threads.nthreads()-1, MaxRestarts=200),
            true)
    end for (mosaic_name, masked_mosaic) in pairs(obj_hit_mosaics))

using JLD2

@info "Saving data and analysis results"
contrast_covers_filename = joinpath(scratch_path,
                                    "$(proj_info.id)_$(proj_info.msfolder)_contrasts_$(proj_info.oesc_ver)_covers.jld2")
@save(contrast_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      objects_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics, #obj_contrast_weak_sets,
      cover_params, obj_hit_covers)
if !@isdefined(obj_hit_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
contrast_covers_jld = jldopen(contrast_covers_filename)
@load(contrast_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      objects_df, obj_contrasts_df,
      obj_contrast_sets, obj_contrast_selsets, obj_contrast_mosaics, #obj_contrast_weak_sets,
      cover_params, obj_contrast_covers)
end

Revise.includet(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup↦gene_name map..."
obj_id2name = Dict(r.object_id => r[Symbol(proj_info.modelobj, "_label")]
                   for r in eachrow(objects_df))

obj_hit_covers_df = innerjoin(OptCoverUtils.covers_report(
    obj_hit_covers, obj_hit_selsets, obj_colls, obj_mosaics, obj_id2name,
    terms_df,
    maskid_col=[:contrast],
    maskedset_col_prefix="hit"),
    contrasts_df, on=[:contrast])

obj_hit_covers_signif_df = combine(groupby(filter(r -> !occursin("?", r.contrast), obj_hit_covers_df), :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:contrast],
                                                   max_term_pvalue=5E-4, max_set_pvalue=1E-1, max_entry_pvalue=1.0),
                   Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_oesc_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_df[obj_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.msfolder)_hit_oesc_signif_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_signif_df[obj_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_plots.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

using PlotlyJS, TextWrap, ORCA

heatmap_layout_attrs = Dict(
    ("GO_CC", true) => Dict(:margin_l => 200),
    ("GO_CC", false) => Dict(:margin_l => 200),
)

for (plot_mosaic, cover_coll) in obj_hit_covers
    isempty(cover_coll.results) && continue
    @info "Plotting $plot_mosaic Pareto front"
    paretofront_plot = OptCoverPlots.plot_paretofront(cover_coll.results[1], plot_unfolded=true)
    plot_filename = joinpath(plots_path, "oesc_$(sel_std_type)_std", "paretofront",
                             "$(proj_info.id)_$(plot_mosaic)_X_treatment_$(sel_std_type)_pareto")
    savefig(paretofront_plot.plot, "$plot_filename.svg")
    PlotlyJS.savehtml(paretofront_plot, "$plot_filename.html")
end

stylize_contrast(str) = foldl(replace, [
    r"_vs_(?:controls|others)" => "",#"&nbsp;<span style=\"color: #808080;\">vs</span>&nbsp;",
    r"SARS_CoV2_(\w+)" => s"\1<sub><span style=\"color: orange;\">CoV2</span></sub>",
    r"SARS_CoV_(\w+)" => s"\1<sub><span style=\"color: sienna;\">CoV</span></sub>",
    r"HCoV_(\w+)" => s"\1<sub><span style=\"color: seagreen;\">HCoV</span></sub>",
    #"controls" => "<span style=\"color: #808080;\">controls</span>",
    #"others" => "<span style=\"color: #808080;\">others</span>",
    ],
    init = str)

function process_contrast_axis(contrast_df)
    contrast_df,
    string.(stylize_contrast.(contrast_df.contrast)),
    string.(stylize_contrast.(contrast_df.contrast))
end

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = signif ? obj_hit_covers_signif_df : obj_hit_covers_df
    for outformat in ["svg", "pdf", "html"]# ("html", "pdf", "svg")
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            Symbol(term_coll), elements_label="proteins",
            maskedset_axis_title = "Comparison",
            maskedset_cols = [:contrast, :nhit],
            process_maskedset_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 100),
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",#outformat in ["svg", "pdf"] ? "#000" : "#BBB",
            cell_width=25, cell_height=14, gridwidth=1,
            transpose=false)
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>100, :margin_t=>20,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    plot_fname = joinpath(plots_path, "$(proj_info.msfolder)_$(proj_info.fit_ver)", "oesc_hits",
                        "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_X_hits$(signif ? "_signif" : "")_heatmap.$(outformat)")
    if outformat == "html"
        PlotlyJS.savehtml(coll_heatmap, plot_fname, :embed);
    else
        try
            savefig(coll_heatmap.plot, plot_fname);
        catch e
            @warn e
        end
    end
    end
end
