source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)

mlog10_pvalue_compress <- function(x, threshold = 10) {
    if (x < threshold) {
        return (x)
    } else if (is.finite(x)) {
        t <- x - threshold
        return (threshold + sqrt(t))
    } else {
        return (2.5*threshold)
    }
}

object_contrasts_4show.df <- object_contrasts.df %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
                  object_annotation = is_viral,
                  show_label = is_viral | is_hit_nomschecks,
                  truncation = case_when(#p_value < p_value_capped ~ "p_value",
                                         median_log2 > median_log2_trunc ~ "median_right",
                                         median_log2 < median_log2_trunc ~ "median_left",
                                         is_signif & !is_hit ~ "significant nonhit",
                                         is_hit ~ "hit",
                                         TRUE ~ "none"))
    
object_contrasts_4show.df %>%
group_by(contrast, std_type) %>% do({
    sel_object_contrast.df <- .
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
    message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]], " std_type=", sel_object_contrast.df$std_type[[1]])
    
    p <- ggplot(sel_object_contrast.df,
                aes(x=median_log2_trunc, y=p_value_compressed, color=is_viral, shape=truncation)) +
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
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & !is_hit), aes(size=truncation=="none"), shape=1L) +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_signif & is_hit), aes(size=truncation=="none"), shape=16L) +
        geom_text_repel(data=dplyr::filter(sel_object_contrast.df, is_signif & show_label),
                  aes(label = object_label),
                  vjust=-0.5, size=3, show.legend = FALSE, segment.color = "gray") +
        scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
        scale_color_manual(values=c("TRUE" = "red", "FALSE" = "black"), na.value="black") +
        scale_shape_manual(values=c("p_value"=-9650, "median_left"=-9664, "median_right"=-9654,
                                    "significant nonhit"=1L, "hit"=16L, "none"=16L), guide="none") +
        #scale_fill_gradient(low="gray75", high="black") +
        scale_size_manual(values=c("TRUE"=1, "FALSE"=2.5), guide="none") +
        #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
        #facet_grid(p_value_range ~ contrast, scales = "free_y") +
        ggtitle(sel_object_contrast_thresholds.df$contrast[[1]],
                subtitle=str_c("std_type=", sel_object_contrast_thresholds.df$std_type[[1]])) +
        theme_bw_ast()
    ggsave(filename = file.path(analysis_path, 'plots', ms_folder,
                                str_c("volcanos_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                paste0(project_id, '_', fit_version, '_volcano_',
                                       sel_object_contrast_thresholds.df$contrast[[1]], '.pdf')),
           plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
    tibble()
})

object_effects_4show.df <- object_effects.df %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
                  object_annotation = is_viral,
                  show_label = is_viral | is_hit_nomschecks,
                  truncation = case_when(#p_value < p_value_capped ~ "p_value",
                      median_log2 > median_log2_trunc ~ "median_right",
                      median_log2 < median_log2_trunc ~ "median_left",
                      is_signif & !is_hit ~ "significant nonhit",
                      is_hit ~ "hit",
                      TRUE ~ "none"))

object_effects_4show.df %>%
    group_by(effect, std_type) %>% do({
        sel_object_effect.df <- .
        sel_object_effect_thresholds.df <- semi_join(object_effects_thresholds.df, sel_object_effect.df)
        message("Plotting ", sel_object_effect_thresholds.df$effect[[1]], " std_type=", sel_object_effect.df$std_type[[1]])
        
        p <- ggplot(sel_object_effect.df,
                    aes(x=median_log2_trunc, y=p_value_compressed, color=is_viral, shape=truncation)) +
            geom_hline(data=sel_object_effect_thresholds.df,
                       aes(yintercept = p_value_threshold), linetype=2, color="darkgray") +
            #geom_hline(data=sel_object_effect_thresholds.df,
            #           aes(yintercept = p_value_max), linetype=1, color="darkgray") +
            geom_vline(data=sel_object_effect_thresholds.df,
                       aes(xintercept = 0), linetype=1, color="darkgray") +
            geom_vline(data=sel_object_effect_thresholds.df,
                       aes(xintercept = median_log2_threshold), linetype=2, color="darkgray") +
            geom_vline(data=sel_object_effect_thresholds.df,
                       aes(xintercept = -median_log2_threshold), linetype=2, color="darkgray") +
            geom_point_rast(data=dplyr::filter(sel_object_effect.df, !is_signif),
                            alpha=0.1, size=0.5, color="darkgray") +
            geom_point(data=dplyr::filter(sel_object_effect.df, is_signif & !is_hit), aes(size=truncation=="none"), shape=1L) +
            geom_point(data=dplyr::filter(sel_object_effect.df, is_signif & is_hit), aes(size=truncation=="none"), shape=16L) +
            geom_text_repel(data=dplyr::filter(sel_object_effect.df, is_signif & show_label),
                            aes(label = object_label),
                            vjust=-0.5, size=3, show.legend = FALSE, segment.color = "gray") +
            scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
            scale_color_manual(values=c("TRUE" = "red", "FALSE" = "black"), na.value="black") +
            scale_shape_manual(values=c("p_value"=-9650, "median_left"=-9664, "median_right"=-9654,
                                        "significant nonhit"=1L, "hit"=16L, "none"=16L), guide="none") +
            #scale_fill_gradient(low="gray75", high="black") +
            scale_size_manual(values=c("TRUE"=1, "FALSE"=2.5), guide="none") +
            #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
            #facet_grid(p_value_range ~ effect, scales = "free_y") +
            ggtitle(sel_object_effect_thresholds.df$effect[[1]],
                    subtitle=str_c("std_type=", sel_object_effect_thresholds.df$std_type[[1]])) +
            theme_bw_ast()
        ggsave(filename = file.path(analysis_path, 'plots', ms_folder,
                                    str_c("volcanos_effects_", sel_object_effect_thresholds.df$std_type[[1]], modelobj_suffix),
                                    paste0(project_id, '_', fit_version, '_volcano_',
                                           sel_object_effect_thresholds.df$effect[[1]], '.pdf')),
               plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
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

