source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)

pgmodel_fit_version <- "20200503"
pgmodel.env <- new.env(parent=base_env())
load(file=file.path(scratch_path, paste0(project_id, '_msglm_data_', ms_folder, '_', pgmodel_fit_version, '.RData')),
     envir=pgmodel.env)

object_contrasts_4show.df <- object_contrasts.df %>% #filter(std_type == "replicate") %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  #object_label = if_else(is_viral & !is.na(obj_bait_full_id), obj_bait_full_id, object_label),
                  show_label = coalesce(is_viral, FALSE) | coalesce(is_hit_nomschecks, FALSE),
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
    fname <- file.path(analysis_path, 'plots', str_c(ms_folder, "_", fit_version),
                       str_c("volcanos_contrasts_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                       paste0(project_id, '_', fit_version, '_volcano_',
                              sel_object_contrast_thresholds.df$contrast[[1]], '.pdf'))
    message(fname)
    ggsave(filename = fname,
           plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
    tibble()
})

object_effects_4show.df <- object_effects.df %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
                  show_label = coalesce(is_viral, FALSE) | coalesce(is_hit_nomschecks, FALSE),
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

sel_std_type <- "median"
sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(gene_names, "MAPK14"))
sel_objects.df <- dplyr::semi_join(modelobjs_df,
    dplyr::distinct(
    bind_rows(
        dplyr::select(dplyr::filter(object_effects.df, std_type==sel_std_type & 
                                    is_hit_nomschecks & effect_type == "treatmentXtimepoint"), object_id),
        dplyr::select(dplyr::filter(object_contrasts.df, std_type==sel_std_type & 
                                    is_hit_nomschecks & str_detect(contrast, "SARS_COV2@\\d+h_vs_mock@\\d+h")), object_id)
    )))

group_by(sel_objects.df, object_id) %>% do({
    sel_obj.df <- .
    message("Plotting ", sel_obj.df$object_label, " time course")
    sel_var <- if_else(sel_std_type == "median", "iaction_labu", "iaction_labu_replCI")
    sel_obj_iactions.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id, object_label),
                                             dplyr::filter(fit_stats$iactions, var == sel_var)) %>%
        dplyr::inner_join(conditions.df) %>%
        dplyr::mutate_at(vars(mean, ends_with("%")),
                         list(~exp(. + global_labu_shift)))
    if (modelobj == "protgroup") {
    sel_obj_msdata.df <- sel_obj.df %>%
        dplyr::inner_join(msdata$protgroup_intensities) %>%
        dplyr::inner_join(msdata$msruns)
    } else if (modelobj == "protregroup") {
        sel_obj_msdata.df <- sel_obj.df %>%
            dplyr::inner_join(msdata_full$protein2protregroup) %>%
            dplyr::filter(is_majority) %>%
            dplyr::select(-is_majority) %>%
            dplyr::inner_join(msdata_full$protein2protgroup) %>%
            dplyr::select(object_id, protgroup_id, protregroup_id) %>% dplyr::distinct() %>%
            dplyr::inner_join(msdata_full$protgroup_intensities) %>%
            dplyr::inner_join(dplyr::select(pgmodel.env$total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
            dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) %>%
            dplyr::inner_join(msdata_full$msruns)
        obj_shifts.df = dplyr::inner_join(sel_obj_msdata.df, sel_obj_iactions.df) %>%
            dplyr::group_by(object_id) %>%
            dplyr::summarise(obj_shift = median(log2(intensity_norm) - log2(`50%`), na.rm=TRUE)) %>%
        dplyr::ungroup()
        sel_obj_msdata.df <- dplyr::left_join(sel_obj_msdata.df, obj_shifts.df) %>%
            dplyr::mutate(intensity_norm = intensity_norm * 2^(-obj_shift))
    }
p <-
ggplot(data=sel_obj_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
    geom_ribbon(aes(x = timepoint_num, ymin = `2.5%`, ymax=`97.5%`),
                 alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.5) +
    geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                 alpha=0.5, stat = "identity", size=0.5) +
    geom_path(aes(x = timepoint_num, y = `50%`), alpha=0.5, size=1, stat="identity") +
    geom_point(data=sel_obj_msdata.df,
               aes(y = intensity_norm), position = position_jitter(width = 0.75, height = 0), size=1) +
    theme_bw_ast(base_family = "", base_size = 8) +
    scale_x_continuous(breaks=unique(msdata$msruns$timepoint_num)) +
    scale_color_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
    scale_fill_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
    scale_y_log10() +
    ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
            subtitle=sel_obj.df$protein_descriptions) +
    facet_wrap( ~ object_label, scales = "free")
ggsave(p, file = file.path(analysis_path, "plots", str_c(ms_folder, '_', fit_version), str_c(modelobj, "s"),
                        paste0(project_id, "_", ms_folder, '_', fit_version, "_",
                               sel_obj.df$object_label[[1]], "_", sel_obj.df$object_id[[1]], "_median.pdf")),
       width=8, height=6, device = cairo_pdf)
    tibble()
})

sel_objects.df <- dplyr::semi_join(modelobjs_df,
                                   dplyr::filter(msdata_full$proteins, str_detect(gene_name, "(^|;)(STAU2)($|;)")) %>%
                                       dplyr::inner_join(msdata_full[[str_c("protein2", modelobj)]]) %>%
                                       dplyr::select(!!modelobj_idcol))
sel_objects.df <- dplyr::filter(modelobjs_df, is_viral)
sel_objects.df <- modelobjs_df # all!!!

sel_pepmods.df <- dplyr::inner_join(sel_objects.df, msdata_full[[str_c(modelobj, "2pepmod")]]) %>%
    dplyr::filter(is_specific) %>%
    dplyr::inner_join(select(msdata_full$pepmods, pepmod_id, peptide_id, mod_seq=pepmod_seq, unmod_seq=peptide_seq)) %>%
    dplyr::select(object_id, pepmod_id, majority_protein_acs, protac_label, gene_label, gene_names,
                  mod_seq, unmod_seq)
if (exists("fit_stats") && has_name(fit_stats, "subobjects")) {
    sel_pepmods.df <- dplyr::left_join(sel_pepmods.df,
                                            dplyr::select(filter(fit_stats$subobjects, var=="suo_shift"),
                                                          pepmod_id, pms_median_log2 = median_log2)) %>%
        dplyr::arrange(object_id, desc(coalesce(pms_median_log2, min(pms_median_log2, na.rm=TRUE))), unmod_seq, mod_seq)
} else {
    sel_pepmods.df <- dplyr::arrange(sel_pepmods.df, object_id, unmod_seq, mod_seq) %>%
        mutate(pms_median_log2 = NA_real_)
}
sel_pepmods.df <- dplyr::mutate(sel_pepmods.df,
                                pepmod_ext = paste0(mod_seq, " (", pepmod_id, ")", if_else(is.na(pms_median_log2), "", "*"))
                                %>% factor(., levels=.))

sel_pepmod_intens.df$msrun
sel_pepmod_intens.df <- dplyr::inner_join(dplyr::select(sel_pepmods.df, object_id, protac_label, gene_label, gene_names,
                                                        pepmod_id, pepmod_ext, pms_median_log2) %>%
                                              dplyr::distinct(),
                                          msdata_full$pepmod_intensities) %>%
    dplyr::inner_join(distinct(select(msdata$msruns, msrun, replicate, treatment, timepoint, condition)) %>%
                          dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
                          dplyr::arrange(treatment, timepoint, replicate) %>%
                          dplyr::mutate(msrun.2 = factor(msrun, levels=msdata$msruns$msrun))) %>%
    dplyr::mutate(intensity_norm = exp(-total_msrun_shift)*intensity,
                  ident_type = "DIA",
                  msrun = msrun.2,
                  msrun.2 = NULL) %>%
    dplyr::group_by(pepmod_id) %>%
    dplyr::mutate(intensity_norm_trunc = pmin(intensity_norm, quantile(intensity_norm, 0.95, na.rm=TRUE))) %>%
    dplyr::ungroup()

group_by(sel_pepmod_intens.df, object_id) %>% do({
    shown_pepmod_intens.df <- mutate(., !!modelobj_idcol := object_id)#sel_pepmod_intens.df
    gene_name <- str_remove(shown_pepmod_intens.df$gene_label[[1]], "\\.\\.\\.$")
    obj_id <- shown_pepmod_intens.df$object_id[[1]]
    is_viral <- semi_join(modelobjs_df, shown_pepmod_intens.df, by="object_id")$is_viral[[1]]
    message("Plotting ", gene_name)
    p <- ggplot(shown_pepmod_intens.df) +
        geom_tile(aes(x=msrun, y=pepmod_ext,
                      fill=intensity_norm_trunc#, color=ident_type
                  ), size=0.5, width=0.95, height=0.95) +
        theme_bw_ast(base_family = "", base_size = 10) +
        theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0)) +
        facet_grid(object_id + gene_label ~ ., scales = "free_y", space="free_y") +
        guides(color=guide_legend("ident_type", override.aes = list(fill=NA, size=2))) +
        ggtitle(str_c(gene_name,  " (pg_id=", obj_id,
                      ", ac=", shown_pepmod_intens.df$protac_label[[1]], ") peptide map"),
                subtitle = semi_join(modelobjs_df, dplyr::select(shown_pepmod_intens.df, !!modelobj_idcol), by=modelobj_idcol)$protein_description[[1]]) +
        scale_x_discrete(drop=FALSE) +
        scale_fill_distiller(na.value="#00000000", type="div", palette = "Spectral", trans="log2") +
        scale_color_manual(na.value="#00000000",
                           values=c("MULTI-MSMS"="black", "MULTI-MATCH-MSMS"="khaki",
                                    "MSMS"="cornflowerblue", "MULTI-SECPEP"="firebrick",
                                    "MULTI-MATCH"="gray", "DIA"="black"))
    fname <- file.path(analysis_path, "plots", str_c(ms_folder, "_", data_version), str_c("peptide_heatmaps", if_else(is_viral, "/viral", "")),
                       paste0(project_id, "_", ms_folder, '_', data_version, "_pepmod_heatmap_", gene_name, "_", obj_id, ".pdf"))
    #message(fname)
    ggsave(filename = fname,
           plot = p, width=16, height=3 + n_distinct(shown_pepmod_intens.df$object_id) + min(20, 0.1*n_distinct(shown_pepmod_intens.df$pepmod_id)),
           device=cairo_pdf, family="Arial")
    tibble()
})
