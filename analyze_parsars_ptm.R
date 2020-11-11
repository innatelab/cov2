source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)

object_contrasts_trunc.df <- tibble(
    contrast_kind = c("treatment_vs_treatment"),
    #p_value_threshold = 1E-2,
    #median_log2_threshold = c(1.0),
    median_log2_max = c(3)
)

object_contrasts_4show.df <- object_contrasts.df %>%
    dplyr::inner_join(object_contrasts_trunc.df) %>%
    dplyr::filter(treatment_lhs != "infected") %>%
    dplyr::inner_join(dplyr::select(dplyr::filter(modelobjs_df, ptmid_is_reference), object_id, ptm_label_no_ptm_type, ptm_id, protein_code)) %>%
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

orgcode_palette <- c(HUMAN="black", contaminant="gray", SARS2 = "#F4982A", CVHSA = "#811A02")
composed_hit_type_palette <- c(shared="black", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02", none="gray")

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
                aes(x=median_log2_trunc, y=p_value_compressed, shape=truncation, size=truncation_type, color=composed_hit_type)) +
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
                        aes(label = ptm_label_no_ptm_type),
                        size=if_else(nlabels > 20, 2.5, 3.5),
                        force=if_else(nlabels > 20, 0.25, 1.0),
                        label.padding=if_else(nlabels > 20, 0.1, 0.25),
                        show.legend = FALSE, segment.color = "gray") +
        scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
        #scale_fill_gradient(low="gray75", high="black") +
        #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
        scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
        scale_size_manual(values=point_truncation_size_palette, guide="none") +
        scale_color_manual(values=composed_hit_type_palette, na.value="black") +
        #scale_color_manual(values=orgcode_palette, guide="none") +
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
                  is_signif, is_hit_nomschecks, is_hit, is_hit_composed, composed_hit_type,
                  object_id, object_label) %>%
    dplyr::inner_join(dplyr::select(dplyr::filter(modelobjs_df, ptmid_is_reference),
                                    ptm_type, object_id, protein_ac, protein_code, ptm_label_no_ptm_type, is_contaminant, is_viral)) %>%
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
                                       is_signif_lhs = is_signif, is_hit_nomschecks_lhs = is_hit_nomschecks, is_hit_lhs = is_hit,
                                       is_hit_composed_lhs = is_hit_composed)) %>%
    dplyr::left_join(dplyr::filter(object_contrasts.df, contrast_type=="comparison") %>%
                         dplyr::select(std_type, contrast_rhs = contrast, object_id, contrast_rhs_p_value = p_value, contrast_rhs_median_log2 = median_log2,
                                       is_signif_rhs = is_signif, is_hit_nomschecks_rhs = is_hit_nomschecks, is_hit_rhs = is_hit,
                                       is_hit_composed_rhs = is_hit_composed)) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_lhs=str_replace(condition, "_(\\d+)h", "@\\1h"),
                                      lhs_median_log2=median_log2 + global_pepmodstate_labu_shift/log(2))) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_rhs=str_replace(condition, "_(\\d+)h", "@\\1h"),
                                      rhs_median_log2=median_log2 + global_pepmodstate_labu_shift/log(2))) %>%
    dplyr::mutate(is_foreground = coalesce(is_hit_composed_lhs, FALSE) | coalesce(is_hit_composed_rhs, FALSE),
                  # coalesce(is_hit_nomschecks_lhs, FALSE) | coalesce(is_hit_nomschecks_rhs, FALSE),
                  is_hilite = is_foreground,# & is_hit_nomschecks,
                  contrast_lhs_median_log2_trunc = if_else(is.na(contrast_median_log2_max), contrast_lhs_median_log2,
                                                           pmax(pmin(contrast_median_log2_max, contrast_lhs_median_log2), -contrast_median_log2_max)),
                  contrast_rhs_median_log2_trunc = if_else(is.na(contrast_median_log2_max), contrast_rhs_median_log2,
                                                           pmax(pmin(contrast_median_log2_max, contrast_rhs_median_log2), -contrast_median_log2_max)),
                  show_label = is_foreground,# & (is_hit_nomschecks | is_viral),
                  truncation = scatter_truncation(contrast_lhs_median_log2, contrast_lhs_median_log2_trunc,
                                                  contrast_rhs_median_log2, contrast_rhs_median_log2_trunc,
                                                  composed_hit_type != "none"),
                  truncation_type = point_truncation_type(truncation, is_foreground))

scatter_type <- "contrast"
object_iactions_4show_kde.df <- object_iactions_4show.df %>% filter(truncation == "none") %>%
group_by(contrast, ptm_type, std_type) %>% group_modify(~{
    kde2d4plot(.x, if (scatter_type == "contrast") { "contrast_lhs_median_log2_trunc" } else { "lhs_median_log2" },
               if (scatter_type == "contrast") { "contrast_rhs_median_log2_trunc" } else { "rhs_median_log2" }, n = 400)$density_df
}) %>% ungroup()

composed_hit_type_palette <- c(shared="black", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02", none="gray")

show_labels <- TRUE
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
    sel_kde.df <- semi_join(object_iactions_4show_kde.df,
                            dplyr::select(sel_object_contrast.df, contrast, std_type, ptm_type)[1, ]) %>%
                  dplyr::filter(bin2d_density >= 0.01)
    p <- if (scatter_type == "contrast"){
        ggplot(sel_object_contrast.df,
               aes(x = contrast_lhs_median_log2_trunc, y = contrast_rhs_median_log2_trunc))
    } else {
        ggplot(sel_object_contrast.df,
               aes(x = lhs_median_log2, y = rhs_median_log2))
    }
    p <- p +
        #geom_raster(data=sel_kde.df, aes(fill=bin2d_density), color=NA) +
        stat_contour_filled(data=sel_kde.df, aes(z=bin2d_density, fill=after_stat(level_mid)), bins=20) +
        stat_contour(data=sel_kde.df, aes(z=bin2d_density, color=after_stat(level))) +
        scale_color_gradient("density", low="gray75", high="black", trans=power_trans(0.25), guide=FALSE) +
        scale_fill_gradient("density", low="gray95", high="slategray4", trans=power_trans(0.25)) +
        geom_vline(xintercept=0, size=1, color="dodgerblue4") +
        geom_hline(yintercept=0, size=1, color="dodgerblue4") +
        new_scale_color() +
        #geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_foreground),
        #                alpha=0.2, size=0.5, color="darkgray", shape=16L) +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]],
                    color="dodgerblue4", linetype="dashed") +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] - sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                    color="dodgerblue4", linetype="dotted", size=0.5) +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] + sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                    color="dodgerblue4", linetype="dotted", size=0.5)
    if (show_labels) {
        p <- p +
        geom_text_repel(data=dplyr::filter(sel_object_contrast.df, show_label),
                        aes(label = ptm_label_no_ptm_type),
                        size=ifelse(manylabels, 2.5, 3.5),
                        show.legend = FALSE, segment.color = "gray")
    }
    p <- p +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_foreground),
                   aes(color=composed_hit_type, shape=truncation, size=truncation_type)) +
        scale_color_manual(values=composed_hit_type_palette, na.value="black") +
        #scale_color_manual(values=orgcode_palette, na.value="black", guide="none") +
        scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
        scale_size_manual(values=if_else(manylabels, 0.5, 1.0) * point_truncation_size_palette, guide="none", ) +
        xlab(str_c("log2(fold-change) ", sel_object_contrast.df$contrast_lhs)) +
        ylab(str_c("log2(fold-change) ", sel_object_contrast.df$contrast_rhs)) +
        coord_fixed() +
        theme_bw_ast()
    plot_path <- file.path(analysis_path, 'plots', str_c(msfolder,'_', fit_version),
                           str_c("scatter_", scatter_type, "s_", sel_std_type, modelobj_suffix))
    if (!dir.exists(plot_path)) dir.create(plot_path)
    ggsave(filename = file.path(plot_path,
                                str_c(project_id, '_', fit_version, "_scatter_", scatter_type, "s_",
                                      sel_ptm_type, "_", sel_std_type, "_",
                                      str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'),
                                      if_else(show_labels, "", "_nolabels"), '.pdf')),
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

treatment_palette <- c(mock="gray", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02")

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
    sel_intensity_range.df <- dplyr::group_by(sel_obj_iactions.df, object_id) %>% dplyr::summarise(
                                               intensity_min = 0.5*quantile(`25%`, 0.05, na.rm=TRUE),
                                               intensity_max = 1.5*quantile(`75%`, 0.95, na.rm=TRUE),
                                               .groups="drop")
    sel_obj_msdata.df <- sel_obj.df %>%
        dplyr::inner_join(sel_intensity_range.df) %>%
        dplyr::inner_join(msdata_full$ptmn2pepmodstate) %>%
        dplyr::left_join(dplyr::select(dplyr::filter(fit_stats$subobjects, var == "suo_shift"),
                                       ptmn_id, pepmodstate_id, glm_subobject_ix, suo_shift = `50%`)) %>%
        dplyr::inner_join(msdata_full$pepmodstate_intensities) %>%
        dplyr::inner_join(msdata_full$ptmn_locprobs) %>%
        dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
        dplyr::mutate(intensity_norm_orig = intensity_norm,
                      intensity_norm_orig_scaled = intensity_norm_orig*exp(-suo_shift),
                      intensity_norm = intensity*exp(-total_msrun_shift),
                      intensity_norm_scaled = intensity_norm*exp(-suo_shift),
                      intensity_norm_orig_scaled_trunc = pmin(pmax(intensity_norm_orig_scaled,
                                                                   intensity_min), intensity_max),
                      intensity_norm_scaled_trunc = pmin(pmax(intensity_norm_scaled,
                                                                   intensity_min), intensity_max),
        ) %>%
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
               aes(y = intensity_norm_scaled_trunc, shape=ms_status),
               position = position_jitter(width = 0.75, height = 0), size=0.5, alpha=0.5, show.legend=FALSE) +
    geom_point(data=sel_obj_msdata.df,
               aes(y = intensity_norm_orig_scaled_trunc, shape=ms_status),
               position = position_jitter(width = 0.75, height = 0), size=1.5) +
    theme_bw_ast(base_family = "", base_size = 8) +
    scale_x_continuous(breaks=unique(msdata$msruns$timepoint_num)) +
    scale_color_manual(values=treatment_palette) +
    scale_shape_manual(values=c("valid"=19, "bad PTM loc"=8, "bad ident"=1, "bad ident&loc"=4)) +
    scale_fill_manual(values=treatment_palette) +
    scale_y_log10() +
    ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
            subtitle=str_c(sel_obj.df$protein_description, " (npepmodstates=", sel_obj.df$n_pepmodstates, ")")) +
    facet_wrap( ~ object_label, scales = "free")
    plot_path <- file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                           str_c("timecourse_", sel_ptm_type, "_", sel_std_type,
                                 modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
    if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
    ggsave(p, file = file.path(plot_path, str_c(project_id, "_", msfolder, '_', fit_version, "_",
                               str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], ".pdf")),
       width=8, height=6, device = cairo_pdf)
}
    tibble()
})

sel_objects.df <- inner_join(msdata_full$proteins, msdata_full$protein2protgroup) %>%
    dplyr::filter(is_contaminant) %>%
    dplyr::mutate(contaminant_type = case_when(str_detect(fasta_header, "Bos") ~ "Cow",
                                               str_detect(fasta_header, "Sus") ~ "Pig",
                                               TRUE ~ "Other")) %>%
    dplyr::select(protgroup_id, contaminant_type) %>% dplyr::distinct()

sel_data.df <- dplyr::inner_join(sel_objects.df, msdata_full$peptide2protgroup) %>%
    dplyr::inner_join(msdata_full$pepmodstates) %>%
    dplyr::filter(nselptms == 0) %>%
    dplyr::inner_join(dplyr::select(msdata_full$pepmodstate_intensities, msrun, pepmodstate_id, intensity_norm)) %>%
    dplyr::inner_join(dplyr::arrange(msdata_full$msruns, dataset, timepoint, treatment, replicate) %>%
                      dplyr::mutate(sample = str_c(condition, "_", replicate),
                                    sample = factor(sample, levels=as.character(sample) %>% unique))) %>%
    dplyr::mutate(intensity_norm_log2 = log2(intensity_norm)) %>%
    dplyr::group_by(dataset, pepmodstate_id) %>%
    dplyr::mutate(intensity_mock0_log2_avg = median(intensity_norm_log2[condition == "mock_0h"], na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(intensity_norm_log2_rel = intensity_norm_log2 - intensity_mock0_log2_avg)
    
contaminant_palette = c("Other" = "lightgray", "Cow" = "firebrick", "Pig" = "goldenrod")
contaminants_plot = ggplot(sel_data.df) +
    geom_boxplot(aes(y = intensity_norm_log2_rel, x = sample, fill=contaminant_type, color=contaminant_type),
                 alpha=0.5, position=position_dodge()) +
    scale_y_continuous(limits=c(-2.5, 2.5)) +
    scale_fill_manual(values = contaminant_palette) +
    scale_color_manual(values = contaminant_palette) +
    theme_bw_ast(base_family = "") +
    facet_grid(dataset ~ .) +
    theme(axis.text.x = element_text(hjust=0, angle=-45))
contaminants_plot
ggsave(contaminants_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                           str_c(project_id, "_", msfolder, '_', fit_version, "_contaminants.pdf")),
       width=16, height=12, device = cairo_pdf)

sel_std_type <- "median"
sel_objects.df <- dplyr::filter(modelobjs_df, object_id == 26850L)# str_detect(object_label, "Phospho_CON__P12763_S325_M2"))
sel_objects.df <- modelobjs_df

treatment_palette <- c(mock="gray", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02")
fp_treatment_palette <- c(mock="lightgray", SARS_CoV2 = "gold", SARS_CoV = "salmon")

msdata_full$pepmodstates <- dplyr::mutate(msdata_full$pepmodstates,
                                          pepmodstate_seq = str_c(pepmod_seq, ".", charge))

require(multidplyr)
plot_cluster <- new_cluster(16)
cluster_library(plot_cluster, c("dplyr", "ggplot2", "stringr", "Cairo", "readr", "rlang"))
cluster_copy(plot_cluster, c("total_msrun_shifts.df", "conditions.df", "fit_stats", "msdata_full", "theme_bw_ast",
                             "sel_std_type", "treatment_palette", "fp_treatment_palette",
                             "global_labu_shift", "project_id", "data_info",
                             "analysis_path", "msfolder", "fit_version", "modelobj_suffix", "fp.env", "ptm2protregroup.df"))

dplyr::left_join(sel_objects.df, dplyr::select(msdata_full$ptmn_stats, ptmn_id, n_pepmodstates)) %>%
dplyr::left_join(dplyr::select(msdata_full$proteins, protein_ac, protein_description=protein_name)) %>%
group_by(object_id) %>% partition(plot_cluster) %>%
do({
    sel_obj.df <- .
    sel_ptm_type <- sel_obj.df$ptm_type
    plot_path <- file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                           str_c("timecourse_", sel_ptm_type, "_", sel_std_type,
                                 modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
    message("Plotting ", sel_obj.df$object_label, " time course")
    sel_var <- if_else(sel_std_type == "median", "iaction_labu", "iaction_labu_replCI")
    sel_obj_iactions.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id),
                                             dplyr::filter(fit_stats$iactions, var == sel_var)) %>%
        dplyr::inner_join(conditions.df) %>%
        dplyr::left_join(dplyr::select(dplyr::filter(fit_stats$subobjects, var == "suo_shift"),
                                       ptmn_id, pepmodstate_id, glm_subobject_ix, suo_shift = `50%`)) %>%
        dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, pepmodstate_seq)) %>%
        dplyr::mutate_at(vars(mean, ends_with("%")),
                         list(~exp(. + suo_shift + global_labu_shift)))
    sel_intensity_range.df <- dplyr::group_by(sel_obj_iactions.df, object_id, pepmodstate_seq) %>% dplyr::summarise(
        intensity_min = 0.25*quantile(`25%`, 0.05, na.rm=TRUE),
        intensity_max = 1.75*quantile(`75%`, 0.95, na.rm=TRUE),
        .groups="drop")
    sel_obj_msdata.df <- sel_obj.df %>%
        dplyr::inner_join(sel_intensity_range.df) %>%
        dplyr::inner_join(msdata_full$ptmn2pepmodstate) %>%
        dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, pepmodstate_seq)) %>%
        dplyr::inner_join(msdata_full$pepmodstate_intensities) %>%
        dplyr::inner_join(msdata_full$ptmn_locprobs) %>%
        dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
        dplyr::mutate(#intensity_norm_orig = intensity_norm,
                      intensity_norm = intensity*exp(-total_msrun_shift),
                      #intensity_norm_orig_trunc = pmin(pmax(intensity_norm_orig, intensity_min), intensity_max),
                      intensity_norm_trunc = pmin(pmax(intensity_norm, intensity_min), intensity_max),
                      intensity_used = intensity_norm_trunc#,
                      #intensity_alt = intensity_norm_orig_trunc
        ) %>%
        dplyr::inner_join(msdata_full$msruns) %>%
        dplyr::mutate(locprob_valid = coalesce(ptm_locprob, 0) >= data_info$locprob_min,
                      qvalue_valid = coalesce(psm_qvalue, 1) <= data_info$qvalue_max,
                      #qvalue_valid = coalesce(psm_pvalue, 1) <= data_info$pvalue_max,
                      ms_status = case_when(locprob_valid & qvalue_valid ~ "valid",
                                            !locprob_valid & qvalue_valid ~ "bad PTM loc",
                                            locprob_valid & !qvalue_valid ~ "bad ident",
                                            TRUE ~ 'bad ident&loc'))
    if (nrow(sel_obj_iactions.df) > 0) {
        p <- ggplot(data=sel_obj_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
            geom_ribbon(aes(x = timepoint_num, ymin = `2.5%`, ymax=`97.5%`),
                        alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
            geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                        alpha=0.5, stat = "identity", size=0.5) +
            geom_path(aes(x = timepoint_num, y = `50%`), alpha=0.5, size=1, stat="identity") +
            #geom_point(data=sel_obj_msdata.df,
            #           aes(y = intensity_alt, shape=ms_status),
            #           position = position_jitter(width = 0.75, height = 0, seed=12323), size=0.5, alpha=0.5, show.legend=FALSE) +
            geom_point(data=sel_obj_msdata.df,
                       aes(y = intensity_used, shape=ms_status),
                       position = position_jitter(width = 0.75, height = 0, seed=12323), size=1.5) +
            geom_text(data=sel_obj_msdata.df,
                      aes(y = intensity_used, label=replicate),
                      position = position_jitter(width = 0.75, height = 0, seed=12323), size=1, color="lightgray", show.legend=FALSE) +
            theme_bw_ast(base_family = "", base_size = 8) +
            scale_x_continuous(breaks=unique(msdata_full$msruns$timepoint_num)) +
            scale_color_manual(values=treatment_palette) +
            scale_shape_manual(values=c("valid"=19, "bad PTM loc"=8, "bad ident"=1, "bad ident&loc"=4)) +
            scale_fill_manual(values=treatment_palette) +
            scale_y_log10("log10(PTM Intensity)") +
            facet_wrap( ~ object_label + pepmodstate_seq, ncol=1, scales = "free")
        h <- 2+3*n_distinct(sel_obj_iactions.df$pepmodstate_id)
        if (exists("fp.env")) {
            sel_fp_iactions.df <- semi_join(fp.env$fit_stats$iactions, semi_join(ptm2protregroup.df, sel_obj.df)) %>%
                filter(var == if_else(sel_std_type == "replicate", "iaction_labu","iaction_labu_replCI")) %>%
                dplyr::inner_join(conditions.df)
            if (nrow(sel_fp_iactions.df) > 0L) {
            fp_plot <- ggplot(data=sel_fp_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
                geom_ribbon(aes(x = timepoint_num, ymin = `2.5%`, ymax=`97.5%`),
                            alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
                geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                            alpha=0.5, stat = "identity", size=0.5) +
                geom_path(aes(x = timepoint_num, y = `50%`), alpha=0.5, size=1, stat="identity") +
                theme_bw_ast(base_family = "", base_size = 8) +
                scale_x_continuous(breaks=unique(msdata_full$msruns$timepoint_num)) +
                scale_color_manual(values=fp_treatment_palette) +
                scale_fill_manual(values=fp_treatment_palette) +
                scale_y_log10("log10(Protein Intensity)")
            p <- ggpubr::ggarrange(p, fp_plot, ncol=1, heights=c(0.6, 0.3))
            h <- h + 2L
            }
        }
        p <- p + ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
                         subtitle=str_c(sel_obj.df$protein_description, " (npepmodstates=", sel_obj.df$n_pepmodstates, ")"))
        if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
        ggsave(p, file = file.path(plot_path, str_c(project_id, "_", msfolder, '_', fit_version, "_",
                                                    str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], "_long.pdf")),
               width=8, height=h, limitsize=FALSE, device = cairo_pdf)
    }
    tibble()
})

object_effects.df <- dplyr::mutate(object_effects.df,
                                   timepointXtreatment_effect_type = case_when(as.logical(replace_na(infected, FALSE)) ~ "common",
                                                                               replace_na(strain == "SARS_CoV2", FALSE) ~ "specific",
                                                                               TRUE ~ NA_character_))
object_effects4corr.df <- dplyr::filter(object_effects.df, !is.na(timepointXtreatment_effect_type) & (std_type == "replicate") & !is.na(timepoint)) %>%
    pivot_wider(c(ptmn_id, timepoint), names_from = timepointXtreatment_effect_type, values_from = c("median_log2")) %>%
    dplyr::filter(between(common, -4, 4) & between(specific, -0.5, 0.5)) %>%
    # & !(between(common, -0.1, 0.1) & between(specific, -0.1, 0.1))) %>%
    dplyr::mutate(common_bin = cut(common, 80)) %>%
    dplyr::group_by(common_bin) %>%
    dplyr::slice_sample(n = 500) %>%
    dplyr::ungroup()

effect_scatter_plot <- ggplot(object_effects4corr.df) +
    geom_density2d_filled(aes(x = common, y = specific, alpha=after_stat(level_mid)), fill="midnightblue", bins=100, contour_var = "density") +
    #scale_fill_viridis_c(option="cividis", trans = power_trans(0.5), limits=c(0.04, 3), na.value = "#00000000") +
    scale_alpha_continuous(trans = power_trans(0.5), limits=c(0.04, 2), na.value = 0) +
    facet_wrap(~ timepoint) +
    theme_bw_ast()
effect_scatter_plot
ggsave(effect_scatter_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                             str_c(project_id, "_", msfolder, '_', fit_version, "_effects_scatter.pdf")),
       width=16, height=16, device = cairo_pdf)

pepmodstate_intensities_treatwide.df <- pivot_wider(dplyr::inner_join(msdata_full$pepmodstate_intensities, msdata_full$msruns) %>%
                                                   dplyr::filter(timepoint_num > 0) %>%
                                                   dplyr::mutate(intensity_norm_log2 = log2(intensity_norm)),
            c(pepmodstate_id, timepoint, replicate), values_from=c(intensity_norm_log2), names_from="treatment", names_sep=".") %>%
    dplyr::inner_join(dplyr::filter(msdata_full$pepmodstates, between(nselptms, 1, 2)))
pepmodstate_intensities_treatwide.df <- mutate(pepmodstate_intensities_treatwide.df,
                                          CoV2_vs_mock = SARS_CoV2 - mock,
                                          SARS_vs_mock = SARS_CoV - mock,
                                          infected_vs_mock = 0.5*(SARS_vs_mock + CoV2_vs_mock),
                                          CoV2_vs_SARS = 0.5*(CoV2_vs_mock - SARS_vs_mock))
pepmodstate_intensities_treatwide4corr.df <- pepmodstate_intensities_treatwide.df %>%
    dplyr::filter(between(infected_vs_mock, -2, 2) & between(CoV2_vs_SARS, -1, 1)) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::mutate(common_bin = as.integer(cut(infected_vs_mock, 80))) %>%
    dplyr::group_by(replicate, timepoint, nselptms, common_bin) %>%
    dplyr::slice_sample(n = 1000) %>%
    dplyr::ungroup()
intensity_scatter_plot <- ggplot(pepmodstate_intensities_treatwide4corr.df) +
    geom_density2d_filled(aes(x = infected_vs_mock, y = CoV2_vs_SARS, alpha=after_stat(level_mid)), fill="midnightblue", bins=60, contour_var = "density") +
    #scale_fill_viridis_c(option="cividis", trans = power_trans(1.0), limits=c(0.04, 2), na.value = "#00000000") +
    scale_alpha_continuous(trans = power_trans(0.5), limits=c(0.04, 2), na.value = 0) +
    facet_grid(replicate + nselptms ~ timepoint) +
    theme_bw_ast()
intensity_scatter_plot
ggsave(intensity_scatter_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                             str_c(project_id, "_", msfolder, '_', fit_version, "_intensity_scatter_treatment.pdf")),
       width=20, height=40, device = cairo_pdf)

pepmodstate_intensities_timewide.df <- pivot_wider(dplyr::inner_join(msdata_full$pepmodstate_intensities, msdata_full$msruns) %>%
                                                   dplyr::filter(timepoint_num > 0) %>%
                                                   dplyr::mutate(intensity_log2 = log2(intensity_norm)),
                                               c(pepmodstate_id, treatment, replicate), values_from=c(intensity_log2),
                                               names_from="timepoint", names_sep=".", names_glue = "{.value}.{timepoint}h") %>%
    dplyr::inner_join(dplyr::filter(msdata_full$pepmodstates, between(nselptms, 1, 2)))
pepmodstate_intensities_timewide.df <- mutate(pepmodstate_intensities_timewide.df,
                                          intensity_log2.12_vs_6 = intensity_log2.12h - intensity_log2.6h,
                                          intensity_log2.24_vs_12 = intensity_log2.24h - intensity_log2.12h)
pepmodstate_intensities_timewide4corr.df <- pepmodstate_intensities_timewide.df %>%
    dplyr::filter(between(intensity_log2.12_vs_6, -2, 2) & between(intensity_log2.24_vs_12, -2, 2)) %>%
    dplyr::group_by(treatment) %>%
    dplyr::mutate(common_bin = as.integer(cut(intensity_log2.24_vs_12, 80))) %>%
    dplyr::group_by(replicate, treatment, nselptms, common_bin) %>%
    dplyr::slice_sample(n = 500) %>%
    dplyr::ungroup()
intensity_scatter_plot <- ggplot(pepmodstate_intensities_timewide4corr.df) +
    geom_density2d_filled(aes(x = intensity_log2.24_vs_12, y = intensity_log2.12_vs_6, alpha=after_stat(level_mid)), fill="midnightblue", bins=100, contour_var = "density") +
    #scale_fill_viridis_c(option="cividis", trans = power_trans(1.0), limits=c(0.04, 2), na.value = "#00000000") +
    scale_alpha_continuous(trans = power_trans(0.5), limits=c(0.04, 1), na.value = 0) +
    facet_grid(replicate + nselptms ~ treatment) +
    theme_bw_ast()
intensity_scatter_plot
ggsave(intensity_scatter_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                                str_c(project_id, "_", msfolder, '_', fit_version, "_intensity_scatter_time.pdf")),
       width=16, height=40, device = cairo_pdf)

pepmodstate_intensities_wide.df <- pivot_wider(dplyr::inner_join(msdata_full$pepmodstate_intensities, msdata_full$msruns) %>%
                                               dplyr::inner_join(total_msrun_shifts.df) %>%
                                               dplyr::filter(timepoint_num > 0) %>%
                                               dplyr::mutate(intensity_log2 = log2(intensity_norm) - total_msrun_shift/log(2)),
                                               c(pepmodstate_id, condition), values_from=c(intensity_log2),
                                               names_from="replicate", names_sep=".", names_glue = "{.value}.{replicate}") %>%
    dplyr::inner_join(dplyr::filter(msdata_full$pepmodstates, between(nselptms, 1, 1)))
pepmodstate_intensities_wide4corr.df <- pepmodstate_intensities_wide.df %>%
    #dplyr::filter(between(intensity_log2.1, -2, 2) & between(intensity_log2.2, -2, 2)) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(common_bin = as.integer(cut(intensity_log2.1, 80))) %>%
    dplyr::group_by(condition, nselptms, common_bin) %>%
    dplyr::slice_sample(n = 100) %>%
    dplyr::ungroup()
intensity_scatter_plot <- ggplot(pepmodstate_intensities_wide4corr.df) +
    geom_abline(intercept = 0, slope = 1, color="firebrick") +
    geom_density2d_filled(aes(x = intensity_log2.1, y = intensity_log2.4, alpha=after_stat(level_mid)), fill="midnightblue", bins=25, contour_var = "density") +
    #scale_fill_viridis_c(option="cividis", trans = power_trans(1.0), limits=c(0.04, 2), na.value = "#00000000") +
    #scale_alpha_continuous(trans = power_trans(0.5), limits=c(0.005, 0.2), na.value = 0) +
    scale_alpha_continuous(trans = power_trans(0.5), limits=c(0.002, 0.05), na.value = 0) +
    scale_x_continuous(limits=c(12,24)) + scale_y_continuous(limits=c(12,24)) +
    facet_wrap(~ condition, scales = "free", ncol=4) +
    theme_bw_ast()
intensity_scatter_plot
ggsave(intensity_scatter_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                                str_c(project_id, "_", msfolder, '_', fit_version, "_intensity_scatter_14.pdf")),
       width=20, height=16, device = cairo_pdf)

pepmodstate_intensities.df <- msdata_full$pepmodstate_intensities %>%
    dplyr::inner_join(total_msrun_shifts.df) %>%
    dplyr::mutate(intensity_norm = intensity_norm*exp(-total_msrun_shift),
                  intensity_norm_log2 = log2(intensity_norm)) %>%
    inner_join(msdata_full$msruns) %>%
    dplyr::group_by(dataset, pepmodstate_id) %>%
    dplyr::mutate(intensity_median = median(intensity_norm, na.rm=TRUE),
                  intensity_mock0h = median(intensity_norm[condition == "mock_0h"], na.rm=TRUE),
                  .groups="drop") %>%
    dplyr::group_by(dataset, pepmodstate_id, timepoint) %>%
    dplyr::mutate(intensity_mock = median(intensity_norm[treatment == "mock"], na.rm=TRUE),
                  intensity_infected = median(intensity_norm[treatment != "mock"], na.rm=TRUE),
                  .groups="drop") %>%
    dplyr::group_by(dataset, pepmodstate_id, condition) %>%
    dplyr::mutate(intensity_cond = median(intensity_norm, na.rm=TRUE),
                  .groups="drop") %>%
    dplyr::mutate(intensity_cond_log2 = log2(intensity_cond),
                  intensity_infected_log2 = log2(intensity_infected),
                  intensity_mock_log2 = log2(intensity_mock),
                  intensity_mock0h_log2 = log2(intensity_mock0h),
                  intensity_median_log2 = log2(intensity_median),
                  intensity_cond_vs_mock_log2 = intensity_cond_log2 - intensity_mock_log2,
                  intensity_cond_vs_mock0h_log2 = intensity_cond_log2 - intensity_mock0h_log2,
                  intensity_infected_vs_mock_log2 = intensity_infected_log2 - intensity_mock_log2,
                  intensity_infected_vs_mock0h_log2 = intensity_infected_log2 - intensity_mock0h_log2,
                  intensity_vs_mock_log2 = intensity_norm_log2 - intensity_mock_log2,
                  intensity_vs_mock0h_log2 = intensity_norm_log2 - intensity_mock0h_log2)
pepmodstate_intensities.df <- dplyr::mutate(pepmodstate_intensities.df,
                                            infected = treatment != "mock") %>%
    dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, nselptms))

pepmodstate_intensities_pairs.df <- dplyr::inner_join(pepmodstate_intensities.df, pepmodstate_intensities.df,
                                                      by = c("dataset", "timepoint", "infected", "pepmodstate_id", "nselptms")) %>%
    dplyr::filter(between(nselptms, 1, 3))
group_by(dplyr::filter(pepmodstate_intensities_pairs.df, (msrun.x != msrun.y) & (timepoint != 0)), dataset, timepoint, infected) %>% do({
    pms_intensity_scatter.df <- .
    sel_ds <- pms_intensity_scatter.df$dataset[[1]]
    sel_timepoint <- pms_intensity_scatter.df$timepoint[[1]]
    sel_infected <- pms_intensity_scatter.df$infected[[1]]
    message("Plotting dataset=", sel_ds, " timepoint=", sel_timepoint, "h infected=", sel_infected)
    intensity_scatter_plot <- ggplot(pms_intensity_scatter.df) +
        geom_abline(intercept = 0, slope = 1, color="firebrick") +
        geom_point_rast(aes(x = intensity_vs_mock_log2.x, y = intensity_vs_mock_log2.y, color=factor(nselptms)),
                        fill="midnightblue", size=1, alpha=0.25) +
        scale_color_manual(values=c("1" = "black", "2" = "midnightblue", "3" = "sienna4")) +
        #geom_density2d_filled(aes(x = intensity_vs_mock_log2.2, y = intensity_vs_mock_log2.4, alpha=after_stat(level_mid)),
        #                      fill="midnightblue", bins=100, contour_var = "density") +
        #scale_fill_viridis_c(option="cividis", trans = power_trans(1.0), limits=c(0.04, 2), na.value = "#00000000") +
        #scale_alpha_continuous(trans = power_trans(1.0), limits=c(0.05, 1.0), na.value = 0) +
        scale_x_continuous(limits=c(-3,3)) + scale_y_continuous(limits=c(-3,3)) +
        facet_grid(msrun.x ~ msrun.y) +
        theme_bw_ast()
    ggsave(intensity_scatter_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                                    str_c(project_id, "_", msfolder, '_', fit_version, "_intensity_vs_mock_scatter_",
                                                          sel_ds, "_", sel_timepoint, "h_", if_else(sel_infected, "infected", "mock"), "_raster.pdf")),
           width=20, height=20, device = cairo_pdf)
    tibble()
})

pepmodstate_intensities_wide.df <- pepmodstate_intensities.df %>%
    dplyr::filter(intensity_vs_mock_log2 != 0) %>%
    pivot_wider(c(pepmodstate_id, condition, treatment, timepoint),
                values_from=c(intensity_norm_log2, intensity_vs_mock_log2, intensity_vs_mock0h_log2),
                names_from="replicate", names_sep=".")

pepmodstate_intensities_wide4corr.df <- pepmodstate_intensities_wide.df %>%
    dplyr::inner_join(msdata_full$pepmodstates) %>%
    dplyr::filter(nselptms == 1 & between(intensity_vs_mock_log2.1, -3, 3) & between(intensity_vs_mock_log2.1, -3, 3)) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(common_bin = as.integer(cut(intensity_vs_mock_log2.1, 80))) %>%
    dplyr::group_by(condition, common_bin) %>%
    dplyr::slice_sample(n = 200) %>%
    dplyr::ungroup()

intensity_scatter_plot <- ggplot(dplyr::filter(pepmodstate_intensities_wide4corr.df, condition != "mock_0h")) +
    geom_abline(intercept = 0, slope = 1, color="firebrick") +
    geom_point_rast(aes(x = intensity_vs_mock_log2.3, y = intensity_vs_mock_log2.4),
                    fill="midnightblue", size=1, alpha=0.25) +
    #geom_density2d_filled(aes(x = intensity_vs_mock_log2.2, y = intensity_vs_mock_log2.4, alpha=after_stat(level_mid)),
    #                      fill="midnightblue", bins=100, contour_var = "density") +
    #scale_fill_viridis_c(option="cividis", trans = power_trans(1.0), limits=c(0.04, 2), na.value = "#00000000") +
    scale_alpha_continuous(trans = power_trans(1.0), limits=c(0.05, 1.0), na.value = 0) +
    scale_x_continuous(limits=c(-3,3)) + scale_y_continuous(limits=c(-3,3)) +
    facet_wrap(~ condition, scales = "free", ncol=4) +
    theme_bw_ast()
intensity_scatter_plot
ggsave(intensity_scatter_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                                str_c(project_id, "_", msfolder, '_', fit_version, "_intensity_vs_mock_scatter_34_raster.pdf")),
       width=16, height=12, device = cairo_pdf)

intensity_scatter_plot <- ggplot(pepmodstate_intensities_wide4corr.df %>% dplyr::filter(treatment == "mock" & timepoint != 0)) +
    geom_abline(intercept = 0, slope = 1, color="firebrick") +
    geom_point_rast(aes(x = intensity_vs_med_log2.2, y = intensity_vs_med_log2.4),
                    fill="midnightblue", size=1, alpha=0.25) +
    #geom_density2d_filled(aes(x = intensity_vs_mock0h_log2.1, y = intensity_vs_mock0h_log2.3, alpha=after_stat(level_mid)),
    #                      fill="midnightblue", bins=100, contour_var = "density") +
    #scale_fill_viridis_c(option="cividis", trans = power_trans(1.0), limits=c(0.04, 2), na.value = "#00000000") +
    scale_alpha_continuous(trans = power_trans(1.0)) +#, limits=c(0.05, 1.0), na.value = 0) +
    scale_x_continuous(limits=c(-3,3)) + scale_y_continuous(limits=c(-3,3)) +
    facet_wrap(~ condition, scales = "free", ncol=2) +
    theme_bw_ast()
intensity_scatter_plot
ggsave(intensity_scatter_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                                str_c(project_id, "_", msfolder, '_', fit_version, "_intensity_vs_mock0h_scatter_13_raster.pdf")),
       width=12, height=12, device = cairo_pdf)

require(broom)

pepmodstate_intensities4glm.df <- pepmodstate_intensities.df %>%
    dplyr::semi_join(dplyr::filter(msdata_full$pepmodstates, nselptms == 1)) %>%
    dplyr::filter(between(intensity_vs_mock_log2, -2, 2) &
                  between(intensity_infected_vs_mock_log2, -2, 2)) %>%
    dplyr::filter(treatment != "mock") %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(intensity_infected_vs_mock_bin = as.integer(cut(intensity_infected_vs_mock_log2, 50))) %>%
    dplyr::group_by(condition, intensity_infected_vs_mock_bin) %>%
    dplyr::slice_sample(n = 500) %>%
    dplyr::ungroup()

infection_scale.glm <- glm(intensity_vs_mock_log2 ~ 0 + intensity_infected_vs_mock_log2:msrun,
                           data=pepmodstate_intensities4glm.df)

infection_scales.df <- broom::tidy(infection_scale.glm) %>%
    tidyr::extract(term, c("msrun"), ":msrun(.+)$", remove=FALSE) %>%
    dplyr::inner_join(msdata_full$msruns) %>%
    dplyr::mutate(msrun_scale = estimate,
                  msrun_scale_log2 = log2(estimate)) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::mutate(msrun_scale_log2_median = median(msrun_scale_log2)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(msrun_scale_norm = 2^(msrun_scale_log2 - msrun_scale_log2_median))

infection_scales_plot <- ggplot(infection_scales.df %>% dplyr::mutate(msrun_short = str_remove(msrun, "phospho_SARS_"))) +
    geom_bar(aes(x = msrun_short, y = msrun_scale_norm, fill=treatment), stat="identity") +
    scale_fill_manual(values=treatment_palette) +
    scale_y_log10() +
    facet_wrap(~ timepoint, scales="free_x", ncol=2) +
    theme_bw_ast() +
    theme(axis.text.x = element_text(angle = -45, hjust=0))
ggsave(infection_scales_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                                str_c(project_id, "_", msfolder, '_', fit_version, "_infection_scales.pdf")),
       width=8, height=8, device = cairo_pdf)

pepmodstate_intensities4mock_glm.df <- pepmodstate_intensities.df %>%
    dplyr::semi_join(dplyr::filter(msdata_full$pepmodstates, nselptms == 1)) %>%
    dplyr::filter(between(intensity_vs_mock0h_log2, -2, 2) &
                      between(intensity_cond_vs_mock0h_log2, -2, 2)) %>%
    dplyr::filter(treatment == "mock") %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(intensity_infected_vs_mock_bin = as.integer(cut(intensity_cond_vs_mock0h_log2, 50))) %>%
    dplyr::group_by(condition, intensity_infected_vs_mock_bin) %>%
    dplyr::slice_sample(n = 500) %>%
    dplyr::ungroup()

mock_scale.glm <- glm(intensity_vs_mock0h_log2 ~ 0 + intensity_cond_vs_mock0h_log2:msrun,
                      data=pepmodstate_intensities4mock_glm.df)

mock_scales.df <- broom::tidy(mock_scale.glm) %>%
    tidyr::extract(term, c("msrun"), ":msrun(.+)$", remove=FALSE) %>%
    dplyr::inner_join(msdata_full$msruns) %>%
    dplyr::mutate(msrun_scale = estimate,
                  msrun_scale_log2 = log2(estimate)) %>%
    dplyr::group_by(timepoint) %>%
    dplyr::mutate(msrun_scale_log2_median = median(msrun_scale_log2)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(msrun_scale_norm = 2^(msrun_scale_log2 - msrun_scale_log2_median))

mock_scales_plot <- ggplot(mock_scales.df %>% dplyr::mutate(msrun_short = str_remove(msrun, "phospho_")) %>% dplyr::filter(condition != "mock_0h")) +
    geom_bar(aes(x = msrun_short, y = msrun_scale_norm, fill=treatment), stat="identity") +
    scale_fill_manual(values=treatment_palette) +
    scale_y_log10() +
    facet_wrap(~ timepoint, scales="free_x", ncol=2) +
    theme_bw_ast() +
    theme(axis.text.x = element_text(angle = -45, hjust=0))
mock_scales_plot
ggsave(mock_scales_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                          str_c(project_id, "_", msfolder, '_', fit_version, "_mock_scales.pdf")),
       width=8, height=8, device = cairo_pdf)

msrunXeffect.df <- inner_join(conditionXeffect.df, msdata$msruns) %>%
    dplyr::left_join(bind_rows(dplyr::select(mock_scales.df, msrun, msrun_scale=msrun_scale_norm),
                               dplyr::select(infection_scales.df, msrun, msrun_scale=msrun_scale_norm))) %>%
    dplyr::mutate(scaled_mult = mult * msrun_scale)
msrunXeffect_wide.df <- pivot_wider(bind_rows(msrunXeffect.df, tibble(msrun = setdiff(msdata_full$msruns$msrun, msrunXeffect.df$msrun),
                                                                      effect="tmp_intercept", scaled_mult = 0)),
                                    msrun, values_from = scaled_mult, names_from = effect, values_fill = 0) %>%
    dplyr::select(-tmp_intercept)
msrunXeffect.mtx <- as.matrix(dplyr::select(msrunXeffect_wide.df, -msrun))
rownames(msrunXeffect.mtx) <- msrunXeffect_wide.df$msrun
msrunXeffect.mtx <- msrunXeffect.mtx[msdata_full$msruns$msrun, ]

plots_path <- file.path(analysis_path, "plots", str_c(msfolder, "_", fit_version))
if (!dir.exists(plots_path)) dir.create(plots_path)

pheatmap(msrunXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plots_path, paste0(project_id, "_msrun_exp_design_", msfolder, "_", fit_version, ".pdf")),
         width = 8, height = 12)

prev.env <- new_environment()
load("/pool/analysis/astukalov/cov2/scratch/cov2_msglm_fit_snaut_parsars_phospho_20201005_20201008.RData", envir = prev.env)
c2c.df <- left_join(object_contrasts.df, prev.env$object_contrasts.df,
                    by=c("object_label", "std_type", "contrast", "ptm_type", "ptm_AA_seq", "ptm_pos", "nselptms"),
                    suffix=c("", ".prev"))
i2i.df <- left_join(fit_stats$iactions, prev.env$fit_stats$iactions,
                    by=c("object_label", "condition", "ptm_type", "ptm_AA_seq", "ptm_pos", "nselptms", "var"),
                    suffix=c("", ".prev")) %>% dplyr::filter(var == "iaction_labu") %>%
    inner_join(conditions.df)
e2e.df <- left_join(fit_stats$object_effects, prev.env$fit_stats$object_effects,
                    by=c("object_label", "effect", "ptm_type", "ptm_AA_seq", "ptm_pos", "nselptms", "var"),
                    suffix=c("", ".prev")) %>% dplyr::filter(var == "obj_effect") %>%
    inner_join(effects.df)

hits.df <- filter(c2c.df, contrast=="SARS_CoV2@6h_vs_mock@6h" & std_type=="median" & is_hit.prev & median_log2.prev < -0.25)
View(dplyr::select(hits.df, object_label, object_label, contrast, is_hit, is_hit.prev,
                   median_log2, median_log2.prev, sd_log2, sd_log2.prev, p_value, p_value.prev))
hit_iactions.df <- semi_join(i2i.df, dplyr::select(hits.df, object_label)) %>% filter(timepoint_num <= 6)
View(dplyr::select(hit_iactions.df, object_label, object_label, timepoint, treatment, condition,
                   median_log2, median_log2.prev, sd_log2, sd_log2.prev))
hit_effs.df <- semi_join(e2e.df, dplyr::select(hits.df, object_label))
View(dplyr::select(hit_effs.df, object_label, object_label, effect, effect_label,
                   median_log2, median_log2.prev, sd_log2, sd_log2.prev, p_value, p_value.prev))
