source(file.path(project_scripts_path, 'cov2_plots_common.R'))

narrow_bait_ids <- c()

object_contrasts_4show.df <- object_contrasts.df %>% #filter(std_type == "replicate") %>%
    #dplyr::left_join(dplyr::select(bait_checks.df, object_id, obj_bait_full_id = bait_full_id, obj_organism = bait_organism)) %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
                  show_label = coalesce(is_viral, FALSE) | coalesce(is_hit_nomschecks, FALSE),
                  median_log2_max = case_when(bait_full_id %in% narrow_bait_ids ~ 8,
                                              TRUE ~ 18),
                  mean_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, mean_log2)),
                  median_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, median_log2)),
                  truncation = scatter_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit | !is_signif),
                  truncation_type = point_truncation_type(truncation, is_signif)) %>%
    dplyr::group_by(std_type, contrast) %>%
    dplyr::mutate(show_label = if_else(rep.int(sum(show_label, na.rm=TRUE) >= 400L, n()),
                                       rep.int(FALSE, n()), show_label)) %>%
    dplyr::ungroup()

object_contrasts_4show.df %>%
group_by(contrast, std_type) %>% do({
    sel_object_contrast.df <- .
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df,
                                                   dplyr::select(sel_object_contrast.df, std_type, contrast))
    message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]], " std_type=", sel_object_contrast.df$std_type[[1]])
    nlabels <- nrow(dplyr::filter(sel_object_contrast.df, is_signif & show_label))

    p <- ggplot(sel_object_contrast.df,
                aes(x=median_log2_trunc, y=p_value_compressed, color=is_viral, shape=truncation, size=truncation_type)) +
        geom_hline(data=sel_object_contrast_thresholds.df,
                   aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
        #geom_hline(data=sel_object_contrast_thresholds.df,
        #           aes(yintercept = p_value_max), linetype=1, color="darkgray") +
        geom_vline(data=sel_object_contrast_thresholds.df,
                   aes(xintercept = 0), linetype=1, color="darkgray") +
        geom_vline(data=sel_object_contrast_thresholds.df,
                   aes(xintercept = median_log2_threshold), linetype=2, color="darkgray") +
        geom_vline(data=sel_object_contrast_thresholds.df,
                   aes(xintercept = -median_log2_threshold), linetype=2, color="darkgray") +
        geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_signif),
                        alpha=0.1, size=0.5, color="darkgray") +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & !is_hit)) +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & is_hit)) +
        geom_text_repel(data=dplyr::filter(sel_object_contrast.df, is_signif & show_label),
                  aes(label = object_label),
                  vjust=-0.5,
                  size=if_else(nlabels > 20, 2.5, 3.5),
                  force=if_else(nlabels > 20, 0.25, 1.0),
                  label.padding=if_else(nlabels > 20, 0.1, 0.25),
                  show.legend = FALSE, segment.color = "gray") +
        scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
        scale_color_manual(values=c("TRUE" = "red", "FALSE" = "black"), na.value="black") +
        #scale_fill_gradient(low="gray75", high="black") +
        #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
        scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
        scale_size_manual(values=point_truncation_size_palette, guide="none") +
        #facet_grid(p_value_range ~ contrast, scales = "free_y") +
        ggtitle(sel_object_contrast_thresholds.df$contrast[[1]],
                subtitle=str_c("std_type=", sel_object_contrast_thresholds.df$std_type[[1]])) +
        theme_bw_ast()
    ggsave(filename = file.path(analysis_path, 'plots', str_c(ms_folder,'_', fit_version),
                                str_c("volcanos_contrasts_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                paste0(project_id, '_', fit_version, '_volcano_',
                                       str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'), '.pdf')),
           plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
    tibble()
})

object_iactions_4show.df <- filter(object_contrasts.df, contrast_type=="comparison" &
                                    std_type=="replicate" & str_detect(contrast, "others|controls", negate = TRUE)) %>%
    dplyr::left_join(dplyr::select(bait_checks_protgroup.df, bait_id, bait_full_id)) %>%
    dplyr::select(std_type, bait_id, contrast, contrast_median_log2 = median_log2, contrast_p_value = p_value,
                  is_signif, is_hit_nomschecks, is_hit, object_id, object_label) %>%
    dplyr::left_join(dplyr::select(modelobjs_df, object_id, organism, majority_protein_acs, is_viral, is_object_used = is_used)) %>%
    dplyr::filter(is_object_used) %>%
    dplyr::mutate(organism = case_when(is.na(organism) & str_detect(majority_protein_acs, "SARS_CoV2") ~ "SARS-CoV-2",
                                       is.na(organism) & str_detect(majority_protein_acs, "SARS_CoV") ~ "SARS-CoV",
                                       TRUE ~ organism),
                  is_viral = case_when(!is.na(is_viral) ~ is_viral,
                                       TRUE ~ str_detect(organism, "CoV")),
                  is_hit_nomschecks = if_else(!is.na(is_hit_nomschecks), is_hit_nomschecks, is_signif & is_viral),
                  is_hit = if_else(!is.na(is_hit), is_hit, is_hit_nomschecks),
                  contrast_median_log2_max = case_when(bait_id %in% c("NSP2", "ORF8") ~ 5,
                                                       TRUE ~ NA_real_)) %>%
    dplyr::inner_join(dplyr::select(dplyr::filter(contrastXcondition.df, weight>0), contrast, condition_lhs=condition, bait_full_id_lhs=bait_full_id)) %>%
    dplyr::inner_join(dplyr::select(dplyr::filter(contrastXcondition.df, weight<0), contrast, condition_rhs=condition, bait_full_id_rhs=bait_full_id)) %>%
    dplyr::mutate(contrast_lhs = str_c(bait_full_id_lhs, "_vs_controls"),
                  contrast_rhs = str_c(bait_full_id_rhs, "_vs_controls")) %>%
    dplyr::left_join(dplyr::filter(object_contrasts.df, contrast_type=="comparison" & std_type == "replicate") %>%
                     dplyr::select(contrast_lhs = contrast, object_id, contrast_lhs_p_value = p_value, contrast_lhs_median_log2 = median_log2,
                                   is_signif_lhs = is_signif, is_hit_nomschecks_lhs = is_hit_nomschecks, is_hit_lhs = is_hit)) %>%
    dplyr::left_join(dplyr::filter(object_contrasts.df, contrast_type=="comparison" & std_type == "replicate") %>%
                     dplyr::select(contrast_rhs = contrast, object_id, contrast_rhs_p_value = p_value, contrast_rhs_median_log2 = median_log2,
                                   is_signif_rhs = is_signif, is_hit_nomschecks_rhs = is_hit_nomschecks, is_hit_rhs = is_hit)) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_lhs=condition, lhs_median_log2=median_log2 + global_protgroup_labu_shift/log(2))) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_rhs=condition, rhs_median_log2=median_log2 + global_protgroup_labu_shift/log(2))) %>%
    dplyr::mutate(is_foreground = coalesce(is_hit_nomschecks_lhs, FALSE) | coalesce(is_hit_nomschecks_rhs, FALSE) |
                    (bait_full_id_rhs == object_label) | (bait_full_id_lhs == object_label) |
                    (bait_full_id_rhs == majority_protein_acs) | (bait_full_id_lhs == majority_protein_acs),
                  is_hilite = is_foreground & is_hit_nomschecks,
                  contrast_lhs_median_log2_trunc = if_else(is.na(contrast_median_log2_max), contrast_lhs_median_log2,
                                                           pmax(pmin(contrast_median_log2_max, contrast_lhs_median_log2), -contrast_median_log2_max)),
                  contrast_rhs_median_log2_trunc = if_else(is.na(contrast_median_log2_max), contrast_rhs_median_log2,
                                                           pmax(pmin(contrast_median_log2_max, contrast_rhs_median_log2), -contrast_median_log2_max)),
                  show_label = is_foreground,# & (is_hit_nomschecks | is_viral),
                  truncation = scatter_truncation(contrast_lhs_median_log2, contrast_lhs_median_log2_trunc,
                                                  contrast_rhs_median_log2, contrast_rhs_median_log2_trunc,
                                                  is_hilite | !is_foreground),
                  truncation_type = point_truncation_type(truncation, is_foreground))

object_iactions_4show.df %>%
group_by(contrast, std_type) %>% do({
    sel_object_contrast.df <- dplyr::filter(., is_hilite |
                                            (lhs_median_log2 >= quantile(lhs_median_log2, na.rm = TRUE, 0.01) &
                                             rhs_median_log2 >= quantile(rhs_median_log2, na.rm = TRUE, 0.01)))
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
    message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]], " std_type=", sel_object_contrast.df$std_type[[1]])

    p <- ggplot(sel_object_contrast.df,
       aes(x = lhs_median_log2, y = rhs_median_log2, color=is_viral)) +
    geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset[[1]]/log(2), color="firebrick", linetype="dashed") +
    geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset[[1]]/log(2) - sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                color="firebrick", linetype="dotted", linewidth=0.5) +
    geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset[[1]]/log(2) + sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                color="firebrick", linetype="dotted", linewidth=0.5) +
    geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_hilite),
                    alpha=0.1, size=0.5, color="darkgray", shape=16L) +
    geom_point(data=dplyr::filter(sel_object_contrast.df, is_hilite),
               aes(shape=is_hilite, size=show_label)) +
    geom_text_repel(data=dplyr::filter(sel_object_contrast.df, show_label),
                    aes(label = object_label),
                    size=3.5,
                    show.legend = FALSE, segment.color = "gray") +
    scale_color_manual(values=organism_type_palette, na.value="black") +
    scale_shape_manual(values=c("TRUE" = 16L, "FALSE" = 1L)) +
    scale_size_manual(values=c("TRUE" = 2, "FALSE" = 1)) +
    coord_fixed() +
    theme_bw_ast()
    ggsave(filename = file.path(analysis_path, 'plots', str_c(ms_folder,'_', fit_version),
                                str_c("scatter_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                paste0(project_id, '_', fit_version, '_scatter_',
                                       str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'), '.pdf')),
           plot = p, width=10, height=10, device=cairo_pdf, family="Arial")
    tibble()
})

viral_prefix = c("SARS-CoV" = "SARS_CoV_", "SARS-CoV-2" = "SARS_CoV2_")

object_iactions_4show.df %>% #filter(bait_id %in% c("NSP2", "ORF8")) %>%
  group_by(contrast, std_type) %>% do({
    sel_object_contrast.df <- dplyr::filter(., is_hilite | between(percent_rank(contrast_lhs_median_log2), 0.001, 0.999)) %>%
      mutate(organism_type = organism_type_map[organism],
             organism_type = if_else(!is_viral | (str_c(viral_prefix[organism], object_label) == bait_full_id_lhs) |
                                       (majority_protein_acs == bait_full_id_lhs)
                                     | (str_c(viral_prefix[organism], object_label) == bait_full_id_rhs) |
                                       (majority_protein_acs == bait_full_id_rhs),
                                     organism_type, "control"))
    print(select(filter(sel_object_contrast.df, is_viral & is_foreground | object_label == "NSP2"),
                 bait_full_id_rhs, bait_full_id_lhs, is_viral, contrast, object_label, organism, majority_protein_acs, organism_type))
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
    message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]], " std_type=", sel_object_contrast.df$std_type[[1]],
            " (", sum(sel_object_contrast.df$show_label), " label(s))")
    manylabels <- sum(sel_object_contrast.df$show_label) > 50
    p <- ggplot(sel_object_contrast.df,
                aes(x = contrast_lhs_median_log2_trunc, y = contrast_rhs_median_log2_trunc, color=organism_type,
                    shape=truncation, size=truncation_type)) +
      geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_foreground),
                      alpha=0.1, size=0.5, color="darkgray", shape=16L) +
      geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]], color="firebrick", linetype="dashed") +
      geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] - sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                  color="firebrick", linetype="dotted", size=0.5) +
      geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] + sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                  color="firebrick", linetype="dotted", size=0.5) +
      geom_text_repel(data=dplyr::filter(sel_object_contrast.df, show_label),
                      aes(label = object_label),
                      size=ifelse(manylabels, 2.5, 3.5),
                      show.legend = FALSE, segment.color = "gray") +
      geom_point(data=dplyr::filter(sel_object_contrast.df, is_foreground)) +
      scale_color_manual(values=organism_type_palette, na.value="black", guide="none") +
      scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
      scale_size_manual(values=if_else(manylabels, 0.5, 1.0) * point_truncation_size_palette, guide="none", ) +
      xlab(str_c("log2(fold-change) ", sel_object_contrast.df$bait_full_id_lhs, " VS background")) +
      ylab(str_c("log2(fold-change) ", sel_object_contrast.df$bait_full_id_rhs, " VS background")) +
      coord_fixed() +
      theme_bw_ast()
    ggsave(filename = file.path(analysis_path, 'plots', str_c(ms_folder,'_', fit_version),
                                str_c("scatter_contrasts_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                paste0(project_id, '_', fit_version, '_scatter_contrasts_',
                                       str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'), '.pdf')),
           plot = p, width=5, height=5, device=cairo_pdf, family="Arial")
    tibble()
  })

object_effects_4show.df <- fit_stats$object_batch_effects %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
                  p_value_threshold = 0.01,
                  median_log2_threshold = case_when(batch_effect == "batch5" ~ 0.75,
                                                    TRUE ~ 0.4),
                  is_signif = p_value <= p_value_threshold & abs(median_log2) >= median_log2_threshold,
                  is_hit = is_signif,
                  show_label = coalesce(is_viral, FALSE) | is_signif,
                  median_log2_max = 5.0,
                  mean_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, mean_log2)),
                  median_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, median_log2)),
                  truncation = volcano_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit, is_signif),
                  truncation_type = volcano_truncation_type(truncation, is_signif)) %>%
    dplyr::group_by(batch_effect) %>%
    dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 1000L, n()), rep.int(FALSE, n()), show_label)) %>%
    dplyr::ungroup()

object_effects_4show.df %>%
    group_by(batch_effect) %>% do({
        sel_object_effect.df <- .
        thresholds_df = select(sel_object_effect.df,
                               batch_effect,
                               p_value_threshold,
                               median_log2_threshold) %>% distinct()
        message("Plotting ", thresholds_df$batch_effect)
        nlabels <- nrow(dplyr::filter(sel_object_effect.df, is_signif & show_label))

        p <- ggplot(sel_object_effect.df,
                    aes(x=median_log2_trunc, y=p_value_compressed, color=is_viral, shape=truncation)) +
            geom_hline(data=thresholds_df,
                       aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
            #geom_hline(data=sel_object_effect_thresholds.df,
            #           aes(yintercept = p_value_max), linetype=1, color="darkgray") +
            geom_vline(data=thresholds_df,
                       aes(xintercept = 0), linetype=1, color="darkgray") +
            geom_vline(data=thresholds_df,
                       aes(xintercept = median_log2_threshold), linetype=2, color="darkgray") +
            geom_vline(data=thresholds_df,
                       aes(xintercept = -median_log2_threshold), linetype=2, color="darkgray") +
            geom_point_rast(data=dplyr::filter(sel_object_effect.df, !is_signif),
                            alpha=0.1, size=0.5, color="darkgray") +
            geom_point(data=dplyr::filter(sel_object_effect.df, is_signif & !is_hit)) +
            geom_point(data=dplyr::filter(sel_object_effect.df, is_signif & is_hit)) +
            geom_text_repel(data=dplyr::filter(sel_object_effect.df, is_signif & show_label),
                            aes(label = object_label),
                            vjust=-0.5,
                            size=if_else(nlabels > 20, 2.5, 3.5),
                            force=if_else(nlabels > 20, 0.25, 1.0),
                            label.padding=if_else(nlabels > 20, 0.1, 0.25),
                            show.legend = FALSE, segment.color = "gray") +
            scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
            scale_color_manual(values=c("TRUE" = "red", "FALSE" = "black"), na.value="black") +
            #scale_fill_gradient(low="gray75", high="black") +
            #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
            scale_shape_manual(values=volcano_truncation_shape_palette, guide="none") +
            scale_size_manual(values=volcano_truncation_size_palette, guide="none") +
            #facet_grid(p_value_range ~ effect, scales = "free_y") +
            ggtitle(thresholds_df$batch_effect[[1]]) +
            theme_bw_ast()
        ggsave(filename = file.path(analysis_path, 'plots', str_c(ms_folder,'_', fit_version),
                                    str_c("volcanos_batch_effects", modelobj_suffix),
                                    paste0(project_id, '_', fit_version, '_volcano_',
                                           thresholds_df$batch_effect, '.pdf')),
               plot = p, width=17, height=22, device=cairo_pdf, family="Arial")
        tibble()
    })

sel_protgroups.df <- msdata$protgroups#dplyr::filter(msdata$protgroups, str_detect(gene_names, "TAOK\\d+"))
sel_protgroups.df <- filter(msdata$protgroups, is_viral)
sel_pg_msdata.df <- sel_protgroups.df %>%
    dplyr::inner_join(msdata_full$protgroup_intensities) %>%
    dplyr::inner_join(dplyr::select(msdata$msruns, -organism)) %>%
    dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
    dplyr::arrange(bait_kind, batch, condition, msrun, protgroup_id) %>%
    dplyr::mutate(batch = str_c("B", batch),
                  intensity_norm = intensity * exp(-total_msrun_shift)) %>%
    dplyr::mutate(bait_full_id = factor(bait_full_id, levels=unique(bait_full_id)))

sel_pg_iactions.df <- dplyr::inner_join(fit_stats$iactions_obsCI, dplyr::select(sel_protgroups.df, protgroup_id, protgroup_label)) %>%
    dplyr::inner_join(inner_join(conditions.df, msdata$msruns %>%
                      mutate(batch = if_else(bait_kind == "control", "control", str_c("B", batch))) %>%
                      dplyr::select(condition, batch) %>% distinct() %>%
                          filter(!(batch == "B2" & condition == "FPMS_SARS_CoV_NSP3_macroD")))) %>%
    dplyr::mutate_at(vars(mean, ends_with("%")),
                     list(~exp(. + global_protgroup_labu_shift))) %>%
    dplyr::mutate(bait_full_id = factor(bait_full_id, levels=levels(sel_pg_msdata.df$bait_full_id)))

batch_palette <- c("control" = "gray", "B1" = "deepskyblue", "B2" = "darkorange", "B3" = "darkseagreen", "B4" = "gold", "B5" = "firebrick")

group_by(sel_pg_msdata.df, protgroup_id) %>% do({
    shown_pg_msdata.df <- .
    shown_pg_iactions.df <- semi_join(sel_pg_iactions.df, distinct(select(shown_pg_msdata.df, protgroup_id)))
    gene_name <- str_remove(shown_pg_iactions.df$protgroup_label[[1]], "\\.\\.\\.$")
    message("Plotting ", gene_name)
    obj_id <- shown_pg_iactions.df$protgroup_id[[1]]
    shown_pg.df <- semi_join(modelobjs_df, shown_pg_iactions.df)
    if (nrow(shown_pg.df)) {
        prot_desc <- shown_pg.df$protein_descriptions[[1]]
        is_viral <- coalesce(shown_pg.df$is_viral[[1]], FALSE)
        if (is_viral) {
            bait_df <- filter(bait_checks_protgroup.df, protgroup_id == obj_id)
            if (nrow(bait_df) > 0) {
                gene_name <- bait_df$bait_full_id[[1]]
            }
        }
    } else {
        prot_desc <- "N/A"
        is_viral <- FALSE
    }
    p <- ggplot(shown_pg_msdata.df, aes(color=batch, fill=batch, x=bait_full_id)) +
    geom_boxplot(data=shown_pg_iactions.df,
                 aes(ymin = `2.5%`, ymax=`97.5%`, lower=`25%`, upper=`75%`, middle=`50%`),
                 alpha=0.5, position = position_dodge(), size=0.2, stat = "identity") +
    geom_point(data=shown_pg_msdata.df,
               aes(x=bait_full_id, y=intensity_norm), position = position_jitter(width = 0.5, height = 0), size=0.3) +
    ggtitle(str_c(gene_name,  " (pg_id=", obj_id,
            ", ac=", shown_pg_iactions.df$majority_protein_acs[[1]], ") model fit"),
            subtitle = prot_desc) +
    theme_bw_ast(base_family = "", base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values=batch_palette) +
    scale_fill_manual(values=batch_palette) +
    scale_y_log10() +
    facet_grid(protgroup_label + protgroup_id ~ ., scales = "free")
ggsave(p, file = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version), "protgroups",
                        paste0(project_id, "_", ms_folder, '_', fit_version, "_", gene_name, "_", obj_id, ".pdf")),
       width=10, height=6, device = cairo_pdf)
     tibble()
})
