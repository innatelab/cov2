source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)

object_contrasts_4show.df <- object_contrasts.df %>% #filter(std_type == "replicate") %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  #object_label = if_else(is_viral & !is.na(obj_bait_full_id), obj_bait_full_id, object_label),
                  show_label = is_viral | is_hit_nomschecks,
                  truncation = volcano_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit, is_signif),
                  truncation_type = volcano_truncation_type(truncation, is_signif)) %>%
    dplyr::group_by(std_type, contrast) %>%
    dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 300L, n()), rep.int(FALSE, n()), show_label)) %>%
    dplyr::ungroup()

object_contrasts_4show.df %>%
group_by(contrast, std_type) %>% do({
    sel_object_contrast.df <- .
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
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
    ggsave(filename = file.path(analysis_path, 'plots', ms_folder,
                                str_c("volcanos_contrasts_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                paste0(project_id, '_', fit_version, '_volcano_',
                                       sel_object_contrast_thresholds.df$contrast[[1]], '.pdf')),
           plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
    tibble()
})

object_effects_4show.df <- object_effects.df %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
                  show_label = is_viral | is_hit_nomschecks,
                  truncation = volcano_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit, is_signif),
                  truncation_type = volcano_truncation_type(truncation, is_signif)) %>%
    dplyr::group_by(std_type, effect) %>%
    dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 300L, n()), rep.int(FALSE, n()), show_label)) %>%
    dplyr::ungroup()

object_effects_4show.df %>%
    group_by(effect, std_type) %>% do({
        sel_object_effect.df <- .
        sel_object_effect_thresholds.df <- semi_join(object_effects_thresholds.df, sel_object_effect.df)
        message("Plotting ", sel_object_effect_thresholds.df$effect[[1]], " std_type=", sel_object_effect.df$std_type[[1]])
        
        p <- ggplot(sel_object_effect.df,
                    aes(x=median_log2_trunc, y=p_value_compressed, color=is_viral, shape=truncation, size=truncation_type)) +
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
            geom_point(data=dplyr::filter(sel_object_effect.df, is_signif & !is_hit)) +
            geom_point(data=dplyr::filter(sel_object_effect.df, is_signif & is_hit)) +
            geom_text_repel(data=dplyr::filter(sel_object_effect.df, is_signif & show_label),
                            aes(label = object_label),
                            vjust=-0.5, size=3, show.legend = FALSE, segment.color = "gray") +
            scale_y_continuous(trans=mlog10_trans(), limits=c(1.0, NA)) +
            scale_color_manual(values=c("TRUE" = "red", "FALSE" = "black"), na.value="black") +
            #scale_fill_gradient(low="gray75", high="black") +
            #scale_alpha_manual(values=c("TRUE"=1.0, "FALSE"=0.5)) +
            scale_shape_manual(values=volcano_truncation_shape_palette, guide="none") +
            scale_size_manual(values=volcano_truncation_size_palette, guide="none") +
            #facet_grid(p_value_range ~ effect, scales = "free_y") +
            ggtitle(sel_object_effect_thresholds.df$effect[[1]],
                    subtitle=str_c("std_type=", sel_object_effect_thresholds.df$std_type[[1]])) +
            theme_bw_ast()
        ggsave(filename = file.path(analysis_path, 'plots', ms_folder,
                                    str_c("volcanos_effects_", sel_object_effect_thresholds.df$std_type[[1]], modelobj_suffix),
                                    paste0(project_id, '_', fit_version, '_volcano_',
                                           str_replace(sel_object_effect_thresholds.df$effect[[1]], ':', '_'), '.pdf')),
               plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
        tibble()
    })

sel_protgroups.df <- dplyr::filter(msdata$protgroups, str_detect(gene_names, "DDX59"))
sel_protgroups.df <- dplyr::semi_join(msdata$protgroups,
    dplyr::distinct(
    bind_rows(
        dplyr::select(dplyr::filter(object_effects.df, std_type=="replicate" & 
                                    is_hit_nomschecks & effect_type == "treatmentXtimepoint"), protgroup_id),
        dplyr::select(dplyr::filter(object_contrasts.df, std_type=="replicate" & 
                                    is_hit_nomschecks & str_detect(contrast, "SARS_COV2@\\d+h_vs_mock@\\d+h")), protgroup_id)
    )))

group_by(sel_protgroups.df, protgroup_id) %>% do({
    sel_protgroup.df <- .
    message("Plotting ", sel_protgroup.df$protgroup_label, " timecourse")
    sel_pg_msdata.df <- sel_protgroup.df %>%
        dplyr::inner_join(msdata$protgroup_intensities) %>%
        dplyr::inner_join(msdata$msruns)
    sel_pg_iactions.df <- dplyr::inner_join(dplyr::select(sel_protgroup.df, protgroup_id, protgroup_label),
                                            dplyr::filter(fit_stats$iactions, var=="iaction_labu_replCI")) %>%
        dplyr::inner_join(conditions.df) %>%
        dplyr::mutate_at(vars(mean, ends_with("%")),
                         list(~exp(. + global_protgroup_labu_shift)))
p <-
ggplot(data=sel_pg_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
    geom_ribbon(aes(x = timepoint_num, ymin = `2.5%`, ymax=`97.5%`),
                 alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
    geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                 alpha=0.5, stat = "identity", size=0.5) +
    geom_path(aes(x = timepoint_num, y = `50%`), alpha=0.5, size=1, stat="identity") +
    geom_point(data=sel_pg_msdata.df,
               aes(y = intensity_norm), position = position_jitter(width = 0.75, height = 0), size=1) +
    theme_bw_ast(base_family = "", base_size = 8) +
    scale_x_continuous(breaks=unique(msdata$msruns$timepoint_num)) +
    scale_color_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
    scale_fill_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
    scale_y_log10() +
    ggtitle(str_c(sel_protgroup.df$protgroup_label, " timecourse"),
            subtitle=sel_protgroup.df$protein_descriptions) +
    facet_wrap( ~ protgroup_label, scales = "free")
ggsave(p, file = file.path(analysis_path, "plots", ms_folder, "protgroups",
                        paste0(project_id, "_", ms_folder, '_', fit_version, "_",
                               sel_protgroup.df$protgroup_label[[1]], "_", sel_protgroup.df$protgroup_id[[1]], ".pdf")),
       width=8, height=6, device = cairo_pdf)
    tibble()
})
