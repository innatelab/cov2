require(ggrastr)
require(ggrepel)

mlog10_pvalue_compress <- function(x, threshold = 10) {
    if (x < threshold) {
        return (x)
    } else {
        t <- x - threshold
        return (threshold + sqrt(t))
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
    ggsave(filename = file.path(analysis_path, 'plots', mq_folder,
                                str_c("volcanos_", sel_object_contrast_thresholds.df$std_type[[1]]),
                                paste0(project_id, '_', fit_version, '_volcano_',
                                       sel_object_contrast_thresholds.df$contrast[[1]], '.pdf')),
           plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
    tibble()
})
