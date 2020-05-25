source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)

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
                  truncation = volcano_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit, is_signif),
                  truncation_type = volcano_truncation_type(truncation, is_signif)) %>%
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
        scale_shape_manual(values=volcano_truncation_shape_palette, guide="none") +
        scale_size_manual(values=volcano_truncation_size_palette, guide="none") +
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
    dplyr::select(std_type, contrast, contrast_median_log2 = median_log2, contrast_p_value = p_value,
                  is_signif, is_hit_nomschecks, is_hit, is_viral, object_id, object_label) %>%
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
                                      object_id, condition_lhs=condition, lhs_median_log2=median_log2 + global_protgroup_labu_shift)) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_rhs=condition, rhs_median_log2=median_log2 + global_protgroup_labu_shift)) %>%
    dplyr::mutate(is_hilite = coalesce(is_hit_nomschecks_lhs, FALSE) | coalesce(is_hit_nomschecks_rhs, FALSE),
                  show_label = is_hilite & (coalesce(is_hit_nomschecks, FALSE) | coalesce(is_viral, FALSE))) %>%
    dplyr::filter(is_hilite | (lhs_median_log2 >= quantile(lhs_median_log2, na.rm = TRUE, 0.01) &
                      rhs_median_log2 >= quantile(rhs_median_log2, na.rm = TRUE, 0.01)))

object_iactions_4show.df %>%
group_by(contrast, std_type) %>% do({
    sel_object_contrast.df <- .
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
    message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]], " std_type=", sel_object_contrast.df$std_type[[1]])

    p <- ggplot(sel_object_contrast.df,
       aes(x = lhs_median_log2, y = rhs_median_log2, color=is_viral)) +
    geom_abline(slope=1, intercept=0, color="firebrick", linetype="dashed") +
    geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                color="firebrick", linetype="dotted", linewidth=0.5) +
    geom_abline(slope=1, intercept=-sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                color="firebrick", linetype="dotted", linewidth=0.5) +
    geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_hilite),
                    alpha=0.1, size=0.5, color="darkgray", shape=16L) +
    geom_point(data=dplyr::filter(sel_object_contrast.df, is_hilite),
               aes(shape=show_label, size=show_label)) +
    geom_text_repel(data=dplyr::filter(sel_object_contrast.df, show_label),
                    aes(label = object_label),
                    nudge_y=-0.25,
                    size=3.5,
                    show.legend = FALSE, segment.color = "gray") +
    xlab(str_c(sel_object_contrast.df$condition_lhs[[1]], " median(log2 fold-change)")) +
    ylab(str_c(sel_object_contrast.df$condition_rhs[[1]], " median(log2 fold-change)")) +
    scale_color_manual(values=c("TRUE" = "red", "FALSE" = "black"), na.value="black") +
    scale_shape_manual(values=c("TRUE" = 16L, "FALSE" = 1L)) +
    scale_size_manual(values=c("TRUE" = 2, "FALSE" = 1)) +
    coord_fixed() +
    theme_bw_ast()
    ggsave(filename = file.path(analysis_path, 'plots', str_c(ms_folder,'_', fit_version),
                                str_c("scatter_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                paste0(project_id, '_', fit_version, '_scatter_',
                                       str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'), '.pdf')),
           plot = p, width=18, height=18, device=cairo_pdf, family="Arial")
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

sel_protgroups.df <- dplyr::filter(msdata$protgroups, str_detect(gene_names, "TAOK\\d+"))
sel_pg_msdata.df <- sel_protgroups.df %>%
    dplyr::inner_join(msdata$protgroup_intensities) %>%
    dplyr::inner_join(msdata$msruns)
sel_pg_iactions.df <- dplyr::inner_join(fit_stats$iactions_obsCI, dplyr::select(sel_protgroups.df, protgroup_id, protgroup_label)) %>%
    dplyr::inner_join(conditions.df) %>%
    dplyr::mutate_at(vars(mean, ends_with("%")),
                     list(~exp(. + global_protgroup_labu_shift)))

p <-
ggplot(sel_pg_msdata.df) +
    geom_boxplot(data=sel_pg_iactions.df,
                 aes(x = bait_full_id, ymin = `2.5%`, ymax=`97.5%`, lower=`25%`, upper=`75%`, middle=`50%`),
                 alpha=0.5, color="darkgray", fill="lightgray", position = position_dodge(), stat = "identity") +
    geom_point(data=sel_pg_msdata.df,
               aes(x=bait_full_id, y=intensity), position = position_jitter(width = 0.75, height = 0), size=1) +
    theme_bw_ast(base_family = "", base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    #scale_size_manual(values=c("By MS/MS" = 1, "By matching" = 0.2)) +
    #scale_color_manual(values=c("proteome" = "gray", "AP" = "black")) +
    #scale_shape_manual(values=c("proteome" = 1, "AP" = 16)) +
    #scale_color_manual(values=treatment_palette) +
    scale_y_log10() +
    facet_grid(. ~ protgroup_label + protgroup_id, scales = "free")
ggsave(p, file = file.path(analysis_path, "plots", ms_folder, "protgroups",
                        paste0(project_id, "_", ms_folder, '_', fit_version, "_", sel_protgroups.df$protgroup_label[[1]], "_", sel_protgroups.df$protgroup_id[[1]], ".pdf")),
       width=8, height=6, device = cairo_pdf)

