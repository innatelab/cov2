proj_info = (id = "cov2",
             data_ver = "20200429",
             fit_ver = "20200429",
             oesc_ver = "20200429",
             modelobj = "protgroup",
             ms_folder = "cov2timecourse_dia_20200423")
using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, DataFrames
using StatsBase

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"));
Revise.includet(joinpath(misc_scripts_path, "msglm_utils.jl"));

objid_col = Symbol(string(proj_info.modelobj, "_id"));

input_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_data_$(proj_info.ms_folder)_$(proj_info.data_ver).RData"), convert=true)
full_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msdata_full_$(proj_info.ms_folder)_$(proj_info.data_ver).RData"), convert=true)
fit_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(proj_info.ms_folder)_$(proj_info.fit_ver).RData"), convert=true)
effects_df = copy(input_rdata["effects.df"]);
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

contrast_matches = match.(Ref(r"(.+)@(\d+)h_vs_(.+)@(\d+)h"), obj_contrasts_df.contrast)
obj_contrasts_df.treatment_lhs = string.(getindex.(contrast_matches, 1))
obj_contrasts_df.timepoint_lhs = parse.(Int, getindex.(contrast_matches, 2))
obj_contrasts_df.treatment_rhs = string.(getindex.(contrast_matches, 3))
obj_contrasts_df.timepoint_rhs = parse.(Int, getindex.(contrast_matches, 4))
obj_contrasts_df.change_alt = getindex.(Ref(Dict("+" => "▲", "-" => "▼", "." => ".")),
                                        obj_contrasts_df.change)

comparisons_df = vcat(hcat(DataFrame(comparison_type = fill("effect", nrow(effects_df))),
                           rename(select!(effects_df, [:effect, :effect_label]),
                                  :effect=>:comparison, :effect_label=>:comparison_label)),
                      hcat(DataFrame(comparison_type = fill("contrast", nrow(contrasts_df))),
                           rename(select!(contrasts_df, [:contrast]),
                                  :contrast=>:comparison),
                           rename(select!(contrasts_df, [:contrast]),
                                  :contrast=>:comparison_label)))

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
pcomplex_coll = FrameUtils.frame2collection(join(pcomplex_iactors_df, pcomplex_iactor2ac_df,
            on=[:file, :entry_index, :interaction_id, :interactor_id], kind=:inner),
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
obj2term_df = select!(join(obj2protac_df, protac2term_df, on = :protein_ac, kind=:inner),
                      Not([:protein_ac])) |> unique!
protac_colls = FrameUtils.frame2collections(protac2term_df, obj_col=:protein_ac,
                                            set_col=:term_id, coll_col=:coll_id)
obj_colls = FrameUtils.frame2collections(obj2term_df, obj_col=objid_col,
                                         set_col=:term_id, coll_col=:coll_id)

@info "Preparing hit sets"
ObjectType = eltype(obj2protac_df.object_id)
obj_hit_sets = Dict{Tuple{String, String, String#=, String=#}, Set{ObjectType}}()
for hits_df in groupby(filter(r -> coalesce(r.is_hit_nomschecks, false), obj_effects_df), [:std_type, :effect#=, :change=#])
    obj_hit_sets[("effect", string.(hits_df[1, :std_type], "_std"), hits_df[1, :effect]#=, hits_df[1, :change]=#)] =
        Set(skipmissing(hits_df.object_id))
end
for hits_df in groupby(filter(r -> coalesce(r.is_hit_nomschecks, false), obj_contrasts_df), [:std_type, :contrast#=, :change=#])
    obj_hit_sets[("contrast", string.(hits_df[1, :std_type], "_std"), hits_df[1, :contrast]#=, hits_df[1, :change]=#)] =
        Set(skipmissing(hits_df.object_id))
end
# only relevant ones
#sel_std_type = "replicate_std"
sel_std_type = "median_std"
obj_hit_selsets = filter(kv -> (kv[1][2] == sel_std_type) && (
    ((kv[1][1] == "contrast") && occursin(r"SARS.+_vs_mock", kv[1][3])) ||
    false && ((kv[1][1] == "effect") && occursin(":treatment", kv[1][3]))),
  obj_hit_sets);

@info "Preparing mosaics..."
observed_protacs = Set(obj2protac_df.protein_ac) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);

obj_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_hit_selsets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.5,
                           uncovered_factor=0.0, covered_factor=0.0)#, covered_factor=0.002)

obj_hit_covers = Dict(begin
    @info "Covering $mosaic_name by bits..."
    mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=Threads.nthreads()-1, MaxRestarts=200),
            true)
    end for (mosaic_name, masked_mosaic) in pairs(obj_hit_mosaics))

using JLD2

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_hit_$(proj_info.oesc_ver)_covers.jld2")
@save(hit_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      objects_df, obj_effects_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics,
      cover_params, obj_hit_covers)
if !@isdefined(obj_hit_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      objects_df, obj_effects_df, obj_contrasts_df,
      obj_hit_sets, obj_hit_selsets, obj_hit_mosaics, obj_hit_selmosaics,
      cover_params, obj_hit_covers)
end

Revise.includet(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing protgroup↦gene_name map..."
obj_id2name = Dict(r.object_id => r[Symbol(proj_info.modelobj, "_label")]
                   for r in eachrow(objects_df))

obj_hit_covers_df = join(OptCoverUtils.covers_report(
   obj_hit_covers, obj_hit_sets, obj_colls, obj_mosaics, obj_id2name,
   terms_df,
   maskid_col=[:comparison_type, :std_type, :comparison#=, :change=#],
   maskedset_col_prefix="hit"),
   comparisons_df, on=[:comparison_type, :comparison], kind=:inner)

obj_hit_covers_signif_df = by(obj_hit_covers_df, :term_collection) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:comparison_type, :comparison_label, :comparison#=, :change=#],
                                                   max_term_pvalue=1E-3, max_set_pvalue=1E-2, max_entry_pvalue=1.0),
                   Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.ms_folder)_hit_oesc_$(sel_std_type)_std_$(proj_info.oesc_ver).txt"),
          obj_hit_covers_df[obj_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_$(proj_info.ms_folder)_hit_oesc_$(sel_std_type)_std_signif_$(proj_info.oesc_ver).txt"),
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
    r"(SARS_COV2+)@(\d+)h" => s"<span style=\"font-weight: bold\">\2</span>h: SARS-CoV-2",
    r"mock@(\d)h" => "mock",
    "_vs_" => "&nbsp;<span style=\"color: #808080;\">vs</span>&nbsp;",
    ],
    init = str)

function process_comparison_axis(comparison_df)
    comparison_df,
    stylize_contrast.(comparison_df.comparison_label),
    stylize_contrast.(comparison_df.comparison_label)#stylize_effect.(effect_df.effect)
end

for term_coll in unique(obj_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "") heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = signif ? obj_hit_covers_signif_df : obj_hit_covers_df
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            Symbol(term_coll), elements_label="proteins",
            maskedset_axis_title = "Comparison",
            maskedset_cols = [:comparison_type, :comparison_label, :comparison#=, :change=#, :nhit],
            process_maskedset_axis=process_comparison_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 160),
            cell_width=25, cell_height=25,
            transpose=false,
            row_order=contrasts -> begin
                contrast_matches = match.(Ref(r"SARS_COV2@(\d+)h"), contrasts.comparison)
                contrasts.timepoint = parse.(Int, getindex.(contrast_matches, 1))
                return sortperm(contrasts, [:timepoint])
            end)
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    plotname = joinpath(plots_path, "$(proj_info.ms_folder)", "oesc_hits_$(sel_std_type)",
                        "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_X_hits$(signif ? "_signif" : "")_heatmap")
    PlotlyJS.savehtml(coll_heatmap, "$(plotname).html", :embed);
    try
        savefig(coll_heatmap.plot, "$(plotname).pdf");
    catch e
        @warn e
    end
end
