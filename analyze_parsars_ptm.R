source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)

object_contrasts_trunc.df <- tibble(
    contrast_kind = c("treatment_vs_treatment"),
    #p_value_threshold = 1E-2,
    #median_log2_threshold = c(1.0),
    median_log2_max = c(3)
)

object_contrasts_4show.df <- object_contrasts.df %>%
    dplyr::inner_join(object_contrasts_trunc.df) %>%
    dplyr::filter(treatment_lhs != "infected") %>%
    dplyr::inner_join(dplyr::select(modelobjs_df, object_id, ptmn_label_no_ptm_type, ptm_id, protein_code)) %>%
    dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median_log2) >= median_log2_threshold),
                  is_hit_nomschecks = is_signif, is_hit=is_hit_nomschecks & is_msvalid_object & !is_contaminant, # & !is_decoy # FIXME
                  p_value_compressed = 10^(-mlog10_pvalue_compress(-log10(p_value))),
                  #object_label = if_else(is_viral & !is.na(obj_bait_full_id), obj_bait_full_id, object_label),
                  show_label = coalesce(is_hit_nomschecks, FALSE),
                  mean_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, mean_log2)),
                  median_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, median_log2)),
                  truncation = scatter_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit | !is_signif),
                  truncation_type = point_truncation_type(truncation, is_signif)) %>%
    dplyr::group_by(contrast, std_type, ptm_type) %>%
    dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 400L, n()), is_hit, show_label),
                  orgcode = str_remove(str_remove(protein_code, ".+_"), "\\.\\.\\.$"),
                  orgcode = case_when(is_contaminant ~ "contaminant",
                                      is.na(orgcode) ~ "HUMAN",
                                      TRUE ~ orgcode)) %>%
    dplyr::ungroup()

orgcode_palette <- c(HUMAN="black", contaminant="gray", SARS2 = "goldenrod", CVHSA = "brown")

group_by(object_contrasts_4show.df, ptm_type, std_type, contrast) %>% do({
    sel_object_contrast.df <- .
    sel_std_type <- sel_object_contrast.df$std_type[[1]]
    sel_ptm_type <- sel_object_contrast.df$ptm_type[[1]]
    sel_contrast <- sel_object_contrast.df$contrast[[1]]
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df,
                                                   dplyr::select(sel_object_contrast.df, contrast_kind))
    message("Plotting ", sel_contrast, " ptm_type=", sel_ptm_type, " std_type=", sel_std_type)
    nlabels <- nrow(dplyr::filter(sel_object_contrast.df, is_signif & show_label))
    
    p <- ggplot(sel_object_contrast.df,
                aes(x=median_log2_trunc, y=p_value_compressed, shape=truncation, size=truncation_type, color=orgcode)) +
        geom_hline(data=sel_object_contrast_thresholds.df,
                   aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
        #geom_hline(data=sel_object_contrast_thresholds.df,
        #           aes(yintercept = p_value_max), linetype=1, color="darkgray") +
        geom_vline(data=sel_object_contrast_thresholds.df,
                   aes(xintercept = contrast_offset_log2), linetype=1, color="darkgray") +
        geom_vline(data=sel_object_contrast_thresholds.df,
                   aes(xintercept = contrast_offset_log2 + median_log2_threshold), linetype=2, color="darkgray") +
        geom_vline(data=sel_object_contrast_thresholds.df,
                   aes(xintercept = contrast_offset_log2 - median_log2_threshold), linetype=2, color="darkgray") +
        geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_signif),
                        alpha=0.1, size=0.5, color="darkgray") +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & !is_hit), shape=1) +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & is_hit)) +
        geom_text_repel(data=dplyr::filter(sel_object_contrast.df, is_signif & show_label),
                        aes(label = ptmn_label_no_ptm_type),
                        size=if_else(nlabels > 20, 2.5, 3.5),
                        force=if_else(nlabels > 20, 0.25, 1.0),
                        label.padding=if_else(nlabels > 20, 0.1, 0.25),
                        show.legend = FALSE, segment.color = "gray") +
        scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
        #scale_fill_gradient(low="gray75", high="black") +
        #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
        scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
        scale_size_manual(values=point_truncation_size_palette, guide="none") +
        scale_color_manual(values=orgcode_palette, guide="none") +
        #facet_grid(p_value_range ~ contrast, scales = "free_y") +
        ggtitle(sel_contrast, subtitle=str_c("ptm_type=", sel_ptm_type, " std_type=", sel_std_type)) +
        theme_bw_ast()
    plot_path <- file.path(analysis_path, 'plots', str_c(msfolder,'_', fit_version),
                           str_c("volcanos_contrasts_", sel_std_type, modelobj_suffix))
    if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

    ggsave(filename = file.path(plot_path,
                                str_c(project_id, '_', fit_version, '_volcano_', sel_ptm_type, "_",
                                      str_replace_all(sel_contrast, ":|@", "_"), '.pdf')),
           plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
    tibble()
})

object_iactions_4show.df <- filter(object_contrasts.df, contrast_type=="comparison" &
                                   contrast_kind == "treatment_vs_treatment" &
                                   treatment_lhs %in% c("SARS_CoV2", "SARS_CoV") &
                                   treatment_rhs %in% c("SARS_CoV2", "SARS_CoV")) %>%
    dplyr::select(std_type, contrast_kind, contrast, timepoint_lhs, timepoint_rhs,
                  contrast_median_log2 = median_log2, contrast_p_value = p_value,
                  is_signif, is_hit_nomschecks, is_hit, object_id, object_label) %>%
    dplyr::left_join(dplyr::select(modelobjs_df, ptm_type, object_id, protein_ac, protein_code, ptmn_label_no_ptm_type, is_contaminant, is_viral)) %>%
    #dplyr::filter(!is.na(ptm_type)) %>%
    dplyr::mutate(condition_lhs = str_remove(contrast, "_vs_.+"),
                  condition_rhs = str_remove(contrast, ".+_vs_"),
                  orgcode = str_remove(str_remove(protein_code, ".+_"), "\\.\\.\\.$"),
                  orgcode = case_when(is_contaminant ~ "contaminant",
                                      is.na(orgcode) ~ "HUMAN",
                                      TRUE ~ orgcode),
                  is_hit_nomschecks = if_else(!is.na(is_hit_nomschecks), is_hit_nomschecks, is_signif & is_viral),
                  is_hit = if_else(!is.na(is_hit), is_hit, is_hit_nomschecks)) %>%
    dplyr::mutate(contrast_lhs = str_c(condition_lhs, "_vs_mock@", timepoint_lhs, "h"),
                  contrast_rhs = str_c(condition_rhs, "_vs_mock@", timepoint_rhs, "h")) %>%
    dplyr::left_join(dplyr::select(object_contrasts_trunc.df, contrast_kind, contrast_median_log2_max = median_log2_max)) %>%
    dplyr::left_join(dplyr::filter(object_contrasts.df, contrast_type=="comparison") %>%
                         dplyr::select(std_type, contrast_lhs = contrast, object_id, contrast_lhs_p_value = p_value, contrast_lhs_median_log2 = median_log2,
                                       is_signif_lhs = is_signif, is_hit_nomschecks_lhs = is_hit_nomschecks, is_hit_lhs = is_hit)) %>%
    dplyr::left_join(dplyr::filter(object_contrasts.df, contrast_type=="comparison") %>%
                         dplyr::select(std_type, contrast_rhs = contrast, object_id, contrast_rhs_p_value = p_value, contrast_rhs_median_log2 = median_log2,
                                       is_signif_rhs = is_signif, is_hit_nomschecks_rhs = is_hit_nomschecks, is_hit_rhs = is_hit)) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_lhs=condition, lhs_median_log2=median_log2 + global_pepmodstate_labu_shift/log(2))) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_rhs=condition, rhs_median_log2=median_log2 + global_pepmodstate_labu_shift/log(2))) %>%
    dplyr::mutate(is_foreground = coalesce(is_hit_nomschecks_lhs, FALSE) | coalesce(is_hit_nomschecks_rhs, FALSE),
                  is_hilite = is_foreground & is_hit_nomschecks,
                  contrast_lhs_median_log2_trunc = if_else(is.na(contrast_median_log2_max), contrast_lhs_median_log2,
                                                           pmax(pmin(contrast_median_log2_max, contrast_lhs_median_log2), -contrast_median_log2_max)),
                  contrast_rhs_median_log2_trunc = if_else(is.na(contrast_median_log2_max), contrast_rhs_median_log2,
                                                           pmax(pmin(contrast_median_log2_max, contrast_rhs_median_log2), -contrast_median_log2_max)),
                  show_label = is_foreground,# & (is_hit_nomschecks | is_viral),
                  truncation = scatter_truncation(contrast_lhs_median_log2, contrast_lhs_median_log2_trunc,
                                                  contrast_rhs_median_log2, contrast_rhs_median_log2_trunc,
                                                  is_hit),
                  truncation_type = point_truncation_type(truncation, is_foreground))

object_iactions_4show.df %>% #filter(bait_id %in% c("NSP2", "ORF8")) %>%
group_by(contrast, ptm_type, std_type) %>% do({
    sel_object_contrast.df <- dplyr::filter(., is_hilite | between(percent_rank(contrast_lhs_median_log2), 0.001, 0.999))
    sel_ptm_type <- sel_object_contrast.df$ptm_type[[1]]
    sel_std_type <- sel_object_contrast.df$std_type[[1]]
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
    message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]],
            " ptm_type=", sel_ptm_type, " std_type=", sel_std_type,
            " (", sum(sel_object_contrast.df$show_label), " label(s))")
    manylabels <- sum(sel_object_contrast.df$show_label) > 300
    p <- ggplot(sel_object_contrast.df,
                aes(x = contrast_lhs_median_log2_trunc, y = contrast_rhs_median_log2_trunc, color=orgcode,
                    shape=truncation, size=truncation_type)) +
        geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_foreground),
                        alpha=0.1, size=0.5, color="darkgray", shape=16L) +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]], color="firebrick", linetype="dashed") +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] - sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                    color="firebrick", linetype="dotted", size=0.5) +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] + sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                    color="firebrick", linetype="dotted", size=0.5) +
        geom_text_repel(data=dplyr::filter(sel_object_contrast.df, show_label),
                        aes(label = ptmn_label_no_ptm_type),
                        size=ifelse(manylabels, 2.5, 3.5),
                        show.legend = FALSE, segment.color = "gray") +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_foreground)) +
        scale_color_manual(values=orgcode_palette, na.value="black", guide="none") +
        scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
        scale_size_manual(values=if_else(manylabels, 0.5, 1.0) * point_truncation_size_palette, guide="none", ) +
        xlab(str_c("log2(fold-change) ", sel_object_contrast.df$contrast_lhs)) +
        ylab(str_c("log2(fold-change) ", sel_object_contrast.df$contrast_rhs)) +
        coord_fixed() +
        theme_bw_ast()
    plot_path <- file.path(analysis_path, 'plots', str_c(msfolder,'_', fit_version),
                           str_c("scatter_contrasts_", sel_std_type, modelobj_suffix))
    if (!dir.exists(plot_path)) dir.create(plot_path)
    ggsave(filename = file.path(plot_path,
                                str_c(project_id, '_', fit_version, '_scatter_contrasts_',
                                      sel_ptm_type, "_", sel_std_type, "_", 
                                      str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'), '.pdf')),
           plot = p, width=16, height=16, device=cairo_pdf, family="Arial")
    tibble()
})

sel_std_type <- "replicate"
sel_objects.df <- dplyr::filter(modelobjs_df, object_label == "GlyGly_AHNAK_K3609_M1")
sel_objects.df <- dplyr::semi_join(modelobjs_df, object_contrasts.df)
sel_objects.df <- dplyr::semi_join(modelobjs_df,
    dplyr::distinct(
    bind_rows(
        #dplyr::select(dplyr::filter(object_effects.df, #std_type==sel_std_type & 
        #                            is_hit_nomschecks & effect_type == "treatmentXtimepoint"), object_id),
        dplyr::select(dplyr::filter(object_contrasts.df, #std_type==sel_std_type & 
                                    is_hit_nomschecks & contrast_kind == "treatment_vs_treatment"), object_id)
    )))

treatment_palette <- c(mock="gray", SARS_CoV2 = "goldenrod", SARS_CoV = "brown")

dplyr::left_join(sel_objects.df, dplyr::select(msdata_full$ptmn_stats, ptmn_id, n_pepmodstates)) %>%
dplyr::left_join(dplyr::select(msdata$proteins, protein_ac, protein_description=protein_name)) %>%
group_by(object_id) %>% do({
    sel_obj.df <- .
    sel_ptm_type <- sel_obj.df$ptm_type
    message("Plotting ", sel_obj.df$object_label, " time course")
    sel_var <- if_else(sel_std_type == "median", "iaction_labu", "iaction_labu_replCI")
    sel_obj_iactions.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id),
                                             dplyr::filter(fit_stats$iactions, var == sel_var)) %>%
        dplyr::inner_join(conditions.df) %>%
        dplyr::mutate_at(vars(mean, ends_with("%")),
                         list(~exp(. + global_labu_shift)))
    sel_obj_msdata.df <- sel_obj.df %>%
        dplyr::inner_join(msdata_full$ptmn2pepmodstate) %>%
        dplyr::left_join(dplyr::select(dplyr::filter(fit_stats$subobjects, var == "suo_shift"),
                                       ptmn_id, pepmodstate_id, glm_subobject_ix, suo_shift = `50%`)) %>%
        dplyr::inner_join(msdata_full$pepmodstate_intensities) %>%
        dplyr::inner_join(msdata_full$ptmn_locprobs) %>%
        dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
        dplyr::mutate(intensity_norm_orig = intensity_norm,
                      intensity_norm_orig_scaled = intensity_norm_orig*exp(-suo_shift),
                      intensity_norm = intensity*exp(-total_msrun_shift),
                      intensity_norm_scaled = intensity_norm*exp(-suo_shift)) %>%
        dplyr::inner_join(msdata_full$msruns) %>%
        dplyr::mutate(locprob_valid = coalesce(ptm_locprob, 0) >= data_info$locprob_min,
                      qvalue_valid = coalesce(qvalue, 1) <= data_info$qvalue_max,
                      ms_status = case_when(locprob_valid & qvalue_valid ~ "valid",
                                            !locprob_valid & qvalue_valid ~ "bad PTM loc",
                                            locprob_valid & !qvalue_valid ~ "bad ident",
                                            TRUE ~ 'bad ident&loc'))
    #print(sel_obj_msdata.df)
if (nrow(sel_obj_iactions.df) > 0) {
p <-
ggplot(data=sel_obj_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
    geom_ribbon(aes(x = timepoint_num, ymin = `2.5%`, ymax=`97.5%`),
                 alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
    geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                 alpha=0.5, stat = "identity", size=0.5) +
    geom_path(aes(x = timepoint_num, y = `50%`), alpha=0.5, size=1, stat="identity") +
    geom_point(data=sel_obj_msdata.df,
               aes(y = intensity_norm_scaled, shape=ms_status),
               position = position_jitter(width = 0.75, height = 0), size=0.5, alpha=0.5, show.legend=FALSE) +
    geom_point(data=sel_obj_msdata.df,
               aes(y = intensity_norm_orig_scaled, shape=ms_status),
               position = position_jitter(width = 0.75, height = 0), size=1.5) +
    theme_bw_ast(base_family = "", base_size = 8) +
    scale_x_continuous(breaks=unique(msdata$msruns$timepoint_num)) +
    scale_color_manual(values=treatment_palette) +
    scale_shape_manual(values=c("valid"=19, "bad PTM loc"=8, "bad ident"=1, "bad ident&loc"=4)) +
    scale_fill_manual(values=treatment_palette) +
    scale_y_log10() +
    ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
            subtitle=str_c(sel_obj.df$protein_description, " npepmodstates=", sel_obj.df$n_pepmodstates)) +
    facet_wrap( ~ object_label, scales = "free")
    plot_path <- file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                           str_c("timecourse_", sel_ptm_type, "_", sel_std_type,
                                 modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
    if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
    print(plot_path)
    ggsave(p, file = file.path(plot_path, str_c(project_id, "_", msfolder, '_', fit_version, "_",
                               str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], ".pdf")),
       width=8, height=6, device = cairo_pdf)
}
    tibble()
})

viral_ptm2ptm.df <- read_tsv(file.path(data_path, msfolder, str_c("ptm_extractor_", data_version),
                                       str_c("viral_ptm2ptm_", data_version, ".txt"))) %>%
    left_join(dplyr::select(msdata_full$proteins, protein_ac, is_viral, is_contaminant, genename)) %>%
    left_join(dplyr::rename_all(dplyr::select(msdata_full$proteins, protein_ac, is_viral, is_contaminant, genename), ~ str_c("hom_", .))) %>%
    mutate(protein_ac_pair = if_else(protein_ac < hom_protein_ac,
                                str_c(protein_ac, "-", hom_protein_ac),
                                str_c(hom_protein_ac, "-", protein_ac)),
           genename_short = str_remove(genename, "SARS_CoV\\d?_"),
           hom_genename_short = str_remove(hom_genename, "SARS_CoV\\d?_"),
           genename_pair = if_else(is.na(hom_genename_short) | (genename_short == hom_genename_short),
                                   genename_short,
                                   if_else(genename_short < hom_genename_short,
                                           str_c(genename_short, "_", hom_genename_short),
                                           str_c(hom_genename_short, "_", genename_short))),
           ptm_correct = (ptm_type == "GlyGly" & ptm_AA_seq %in% "K") |
               (ptm_type == "Phospho" & ptm_AA_seq %in% c("S","T","Y")),
           hom_ptm_correct = (ptm_type == "GlyGly" & hom_ptm_AA %in% "K") |
               (ptm_type == "Phospho" & hom_ptm_AA %in% c("S","T","Y")))

ptm_qvalue_max <- 1E-3

viral_ptmn_aligned.df <- select(viral_ptm2ptm.df, -agn_match_ratio) %>% distinct() %>%
    dplyr::left_join(distinct(select(msdata_full$ptm2gene, ptm_id, protein_ac, is_viral, is_contaminant, ptm_pos, ptm_AA_seq))) %>%
    dplyr::left_join(distinct(select(msdata_full$ptm2gene, hom_ptm_id=ptm_id, hom_protein_ac=protein_ac, hom_ptm_pos=ptm_pos, hom_ptm_AA=ptm_AA_seq))) %>%
    dplyr::left_join(select(msdata_full$ptmns, ptmn_id, ptm_id)) %>%
    #dplyr::left_join(select(msdata_full$ptmns, hom_ptmn_id=ptmn_id, hom_ptm_id=ptm_id)) %>%
    dplyr::left_join(dplyr::select(msdata_full$ptmn_stats, ptmn_id, pms_qvalue_min)) %>%
    dplyr::left_join(dplyr::select(msdata_full$ptmn_stats, hom_ptm_id=ptm_id, hom_pms_qvalue_min=pms_qvalue_min) %>%
                     dplyr::group_by(hom_ptm_id) %>% dplyr::summarise(hom_pms_qvalue_min=min(hom_pms_qvalue_min, na.rm=TRUE), .groups="drop")) %>%
    dplyr::filter(pmin(coalesce(pms_qvalue_min, 1), coalesce(hom_pms_qvalue_min, 1)) <= ptm_qvalue_max) %>%
    dplyr::mutate(ptm_status = case_when(is_observed ~ "observed",
                                         ptm_correct ~ "potential",
                                         TRUE ~ "N/A"),
                  organism_short = case_when(organism == "Severe acute respiratory syndrome coronavirus 2 OX=2697049" ~ "SARS_CoV2",
                            organism == "Human SARS coronavirus OX=694009" ~ "SARS_CoV",
                            TRUE ~ organism))
View(filter(viral_ptm2ptm.df, genename_short == "ORF9b" & ptm_pos == 81))

viral_ptmn_stats.df <- dplyr::full_join(viral_ptmn_aligned.df, msdata_full$msruns, by=character()) %>%
    dplyr::left_join(dplyr::select(msdata_full$ptmns, ptmn_id, ptm_id)) %>%
    dplyr::left_join(msdata_full$ptmn2pepmodstate) %>%
    dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
    dplyr::group_by(ptmn_id, treatment) %>%
    dplyr::summarise(ptm_treatment_qvalue_min = min(qvalue, na.rm=TRUE),
                     .groups="drop")

viral_ptmn_intensity.df <- dplyr::full_join(viral_ptmn_aligned.df, conditions.df, by=character()) %>%
    dplyr::left_join(viral_ptmn_stats.df) %>%
    dplyr::left_join(dplyr::select(filter(fit_stats$iactions, var == "iaction_labu"), ptmn_id, nselptms, condition, median_log2)) %>%
    dplyr::mutate(ptm_id = if_else(is.na(ptm_id), -hom_ptm_id, ptm_id)) %>% # assign some ptm_id for unobserved ptm
    dplyr::filter(treatment == organism_short) %>%
    # show only the highest ptmn_id per ptm_id
    dplyr::group_by(ptm_id, condition) %>% dplyr::filter(row_number(-coalesce(100 + median_log2, 0)) == 1L) %>% dplyr::ungroup() %>%
    dplyr::mutate(timepoint_label = factor(str_c(timepoint, "h p.i."), levels=str_c(sort(unique(conditions.df$timepoint)), "h p.i.")),
                  shown_ptm_pos = if_else(is.na(agn_ptm_pos), ptm_pos, agn_ptm_pos),
                  shown_median_log2 = coalesce(global_labu_shift/log(2) + median_log2, 1),
                  shown_median_log2 = if_else(organism_short == "SARS_CoV2", shown_median_log2, -shown_median_log2),
                  ptm_label = str_c(ptm_AA_seq, ptm_pos, if_else(coalesce(nselptms, 1) > 1, str_c(" (M", nselptms,")"), "")),
                  ptm_status = case_when(ptm_status == "observed" & ptm_treatment_qvalue_min > ptm_qvalue_max ~ "observed (high q-value)",
                                         TRUE ~ ptm_status))

View(filter(viral_ptmn_intensity.df, genename_short == "S" & ptm_pos == 478))

ptm_palette <- c("Phospho" = "royalblue", "GlyGly" = "black", "N/")
ptm_status_palette <- c("observed"=19L, "observed (high q-value)"=10L, "potential" = 21L, "N/A" = 4L)

group_by(viral_ptmn_intensity.df, genename_pair, agn_len) %>%
group_walk(~{
    sel_iactions.df <- .x
    viral_gene <- .y$genename_pair[[1]]
    genenames_df <- semi_join(viral_ptm2ptm.df, .y) %>% dplyr::select(genename, hom_genename) %>% dplyr::distinct()
    genenames <- c(genenames_df$genename, genenames_df$hom_genename)
    genenames <- sort(unique(genenames[!is.na(genenames)]), decreasing = TRUE)
    message("Plotting PTMs of ", viral_gene)
    p <- ggplot(sel_iactions.df, aes(x = agn_ptm_pos, y = shown_median_log2, color = organism_short)) +
        geom_hline(yintercept = 0, color = "black", size = 1.5) +
        geom_segment(aes(xend = agn_ptm_pos, y = 0, yend = shown_median_log2)) +
        geom_point(aes(shape=ptm_status), fill="white") +
        scale_color_manual("Virus", values = treatment_palette) +
        new_scale_color() +
        geom_text_repel(data=filter(sel_iactions.df, treatment == "SARS_CoV2" & ptm_status != "N/A"),
                        aes(label = ptm_label, color=ptm_type), min.segment.length=0.2, segment.alpha=0.5, segment.size=0.25, nudge_y = 1) +
        geom_text_repel(data=filter(sel_iactions.df, treatment == "SARS_CoV" & ptm_status != "N/A"),
                        aes(label = ptm_label, color=ptm_type), min.segment.length=0.2, segment.alpha=0.5, segment.size=0.25, nudge_y = -1) +
        scale_color_manual("PTM Type", values = ptm_palette) +
        scale_shape_manual("PTM Status", values = ptm_status_palette) +
        scale_x_continuous("Aligned Position", limits=c(0L, .y$agn_len[[1]])) +
        facet_grid(timepoint_label ~ .) +
        theme_bw_ast(base_family = "") +
        ggtitle(str_c("PTMs on ", str_c(genenames, collapse="/"), " proteins"))
    plot_path <- file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version), "viral_ptm")
    if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
    ggsave(p, file = file.path(plot_path, str_c(project_id, "_", msfolder, '_', fit_version, "_viral_",
                                                viral_gene, #"_qvalue_1E-3", 
                                                ".pdf")),
           width=14, height=16, device = cairo_pdf)
})
