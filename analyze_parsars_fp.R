source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)

object_contrasts_trunc.df <- tibble(
    contrast_kind = c("treatment_vs_treatment"),
    #p_value_threshold = 1E-2,
    #median_log2_threshold = c(1.0),
    median_log2_max = c(1.5)
)

object_contrasts_4show.df <- object_contrasts.df %>%
    dplyr::inner_join(object_contrasts_trunc.df) %>%
    dplyr::filter(treatment_lhs != "infected") %>%
    dplyr::inner_join(dplyr::select(modelobjs_df, object_id, is_msvalid_object, protein_label)) %>%
    dplyr::mutate(is_signif = (p_value <= p_value_threshold) & (abs(median_log2) >= median_log2_threshold),
                  is_hit_nomschecks = is_signif, is_hit=is_hit_nomschecks & is_msvalid_object & !is_reverse & !is_contaminant,
                  p_value_compressed = 10^(-mlog10_pvalue_compress(-log10(p_value))),
                  #object_label = if_else(is_viral & !is.na(obj_bait_full_id), obj_bait_full_id, object_label),
                  show_label = coalesce(is_hit_nomschecks, FALSE),
                  mean_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, mean_log2)),
                  median_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, median_log2)),
                  truncation = scatter_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit | !is_signif),
                  truncation_type = point_truncation_type(truncation, is_signif)) %>%
    dplyr::group_by(contrast, std_type) %>%
    dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 400L, n()), is_hit, show_label),
                  orgcode = str_remove(str_remove(protein_label, ".+_"), "\\.\\.\\.$"),
                  orgcode = case_when(is_contaminant ~ "contaminant",
                                      is.na(orgcode) ~ "HUMAN",
                                      TRUE ~ orgcode)) %>%
    dplyr::ungroup()

orgcode_palette <- c(HUMAN="black", contaminant="gray", SARS2 = "#F4982A", CVHSA = "#811A02")

group_by(object_contrasts_4show.df, std_type, contrast) %>% do({
    sel_object_contrast.df <- .
    sel_std_type <- sel_object_contrast.df$std_type[[1]]
    sel_contrast <- sel_object_contrast.df$contrast[[1]]
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df,
                                                   dplyr::select(sel_object_contrast.df, contrast_kind))
    message("Plotting ", sel_contrast, " std_type"=sel_std_type)
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
                        aes(label = object_label),
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
        ggtitle(sel_contrast, subtitle=str_c("std_type=", sel_std_type)) +
        theme_bw_ast()
    plot_path <- file.path(analysis_path, 'plots', str_c(msfolder,'_', fit_version),
                           str_c("volcanos_contrasts_", sel_std_type, modelobj_suffix))
    if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)

    ggsave(filename = file.path(plot_path,
                                str_c(project_id, '_', fit_version, '_volcano_',
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
    dplyr::left_join(dplyr::select(modelobjs_df, object_id, majority_protein_acs, protein_label, is_contaminant, is_viral)) %>%
    dplyr::mutate(condition_lhs = str_remove(contrast, "_vs_.+"),
                  condition_rhs = str_remove(contrast, ".+_vs_"),
                  orgcode = str_remove(str_remove(protein_label, ".+_"), "\\.\\.\\.$"),
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
                                      object_id, condition_lhs=condition, lhs_median_log2=median_log2 + global_protgroup_labu_shift/log(2))) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, str_detect(var, "iaction_labu(?:_replCI)?")),
                                      std_type = if_else(var=="iaction_labu", "median", "replicate"),
                                      object_id, condition_rhs=condition, rhs_median_log2=median_log2 + global_protgroup_labu_shift/log(2))) %>%
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

object_iactions_4show_kde.df <- object_iactions_4show.df %>% #filter(bait_id %in% c("NSP2", "ORF8")) %>%
group_by(contrast, std_type) %>% group_modify(~{
    kde2d4plot(.x, "contrast_lhs_median_log2_trunc", "contrast_rhs_median_log2_trunc", n = 400)$density_df
}) %>% ungroup()

show_labels <- FALSE
object_iactions_4show.df %>% #filter(bait_id %in% c("NSP2", "ORF8")) %>%
group_by(contrast, std_type) %>% do({
    sel_object_contrast.df <- dplyr::filter(., is_hilite | between(percent_rank(contrast_lhs_median_log2), 0.001, 0.999))
    sel_std_type <- sel_object_contrast.df$std_type[[1]]
    sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
    message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]], " std_type=", sel_std_type,
            " (", sum(sel_object_contrast.df$show_label), " label(s))")
    manylabels <- sum(sel_object_contrast.df$show_label) > 300
    sel_kde.df <- semi_join(object_iactions_4show_kde.df,
                            dplyr::select(sel_object_contrast.df, contrast, std_type)[1, ]) %>%
                  dplyr::filter(bin2d_density > 0.001)
    p <- ggplot(sel_object_contrast.df,
                aes(x = contrast_lhs_median_log2_trunc, y = contrast_rhs_median_log2_trunc)) +
        #geom_raster(data=sel_kde.df, aes(fill=bin2d_density), color=NA) +
        stat_contour_filled(data=sel_kde.df, aes(z=bin2d_density, fill=after_stat(level_mid))) +
        stat_contour(data=sel_kde.df, aes(z=bin2d_density, color=after_stat(level))) +
        scale_color_gradient(low="gray75", high="black", trans=power_trans(0.25), guide=FALSE) +
        scale_fill_gradient("density", low="gray95", high="slategray4", trans=power_trans(0.25)) +
        geom_vline(xintercept=0, size=1, color="dodgerblue4") +
        geom_hline(yintercept=0, size=1, color="dodgerblue4") +
        new_scale_color() +
        #geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_foreground),
        #                alpha=0.1, size=0.5, color="darkgray", shape=16L) +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]], color="firebrick", linetype="dashed") +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] - sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                    color="firebrick", linetype="dotted", size=0.5) +
        geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] + sel_object_contrast_thresholds.df$median_log2_threshold[[1]],
                    color="firebrick", linetype="dotted", size=0.5)
    if (show_labels) {
        p <- p +
        geom_text_repel(data=dplyr::filter(sel_object_contrast.df, show_label),
                        aes(label = object_label),
                        size=ifelse(manylabels, 2.5, 3.5),
                        show.legend = FALSE, segment.color = "gray")
    }
    p <- p +
        geom_point(data=dplyr::filter(sel_object_contrast.df, is_foreground),
                   aes(color=orgcode, shape=truncation, size=truncation_type)) +
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
                                str_c(project_id, '_', fit_version, '_scatter_contrasts_', sel_std_type, "_",
                                      str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'),
                                      if_else(show_labels, "", "_nolabels"), '.pdf')),
           plot = p, width=16, height=16, device=cairo_pdf, family="Arial")
    tibble()
})

sel_std_type <- "replicate"
sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(gene_names, "STAT1"))
sel_objects.df <- dplyr::semi_join(modelobjs_df, object_contrasts.df)
sel_objects.df <- dplyr::semi_join(modelobjs_df,
    dplyr::distinct(
    bind_rows(
        dplyr::select(dplyr::filter(object_effects.df, #std_type==sel_std_type & 
                                    is_hit_nomschecks & effect_type == "treatmentXtimepoint"), object_id),
        dplyr::select(dplyr::filter(object_contrasts.df, #std_type==sel_std_type & 
                                    is_hit_nomschecks & contrast_kind == "treatment_vs_treatment"), object_id)
    )))

treatment_palette <- c(mock="gray", SARS_CoV2 = "goldenrod", SARS_CoV = "brown")

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
            dplyr::group_by(object_id) %>%
            dplyr::filter(is_majority == any(is_majority)) %>%
            dplyr::ungroup() %>%
            dplyr::select(object_id, protregroup_id, protein_ac) %>%
            dplyr::inner_join(msdata_full$protein2protgroup) %>%
            dplyr::group_by(object_id, protgroup_id) %>%
            dplyr::filter(is_majority == any(is_majority)) %>%
            dplyr::ungroup() %>%
            dplyr::select(object_id, protgroup_id, protregroup_id) %>% dplyr::distinct() %>%
            dplyr::inner_join(msdata_full$protgroup_intensities) %>%
            dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
            dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) %>%
            dplyr::inner_join(msdata_full$msruns)
        #print(sel_obj_msdata.df)
        obj_shifts.df = dplyr::inner_join(sel_obj_msdata.df, sel_obj_iactions.df) %>%
            dplyr::group_by(object_id) %>%
            dplyr::summarise(obj_shift = median(log2(intensity_norm) - log2(`50%`), na.rm=TRUE)) %>%
        dplyr::ungroup()
        #print(obj_shifts.df)
        sel_obj_msdata.df <- dplyr::left_join(sel_obj_msdata.df, obj_shifts.df) %>%
            dplyr::mutate(intensity_norm = intensity_norm * 2^(-obj_shift))
    }
if (nrow(sel_obj_iactions.df) > 0) {
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
    scale_color_manual(values=treatment_palette) +
    scale_fill_manual(values=treatment_palette) +
    scale_y_log10() +
    ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
            subtitle=sel_obj.df$protein_descriptions) +
    facet_wrap( ~ object_label, scales = "free")
ggsave(p, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                           str_c("timecourse_", sel_std_type, modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")),
                           str_c(project_id, "_", msfolder, '_', fit_version, "_",
                                 sel_obj.df$object_label[[1]], "_", sel_obj.df$object_id[[1]], ".pdf")),
       width=8, height=6, device = cairo_pdf)
}
    tibble()
})

sel_objects.df <- dplyr::semi_join(modelobjs_df,
                                   dplyr::filter(msdata_full$proteins, str_detect(gene_name, "(^|;)(CUEDC2|LSG1|CHMP1B|BORCS6|C4A|HBA1|MAP4|MPDU1|HLA|PCLAF|PTPRF|SMIM26|TBC1D10A|EDEM3)($|;)")) %>%
                                       dplyr::inner_join(msdata_full[[str_c("protein2", modelobj)]]) %>%
                                       dplyr::select(!!modelobj_idcol))
sel_objects.df <- dplyr::filter(modelobjs_df, is_viral)
sel_objects.df <- modelobjs_df # all!!!

sel_pepmodstates.df <- dplyr::inner_join(sel_objects.df, msdata_full[[str_c(modelobj, "2pepmod")]]) %>%
    dplyr::filter(is_specific) %>%
    dplyr::inner_join(select(msdata_full$pepmodstates, pepmodstate_id, pepmod_id, charge)) %>%
    dplyr::inner_join(select(msdata_full$pepmods, pepmod_id, peptide_id, pepmod_seq, peptide_seq)) %>%
    dplyr::select(object_id, pepmodstate_id, pepmod_id, majority_protein_acs, object_label=protregroup_label, protac_label, gene_label, gene_names,
                  pepmod_seq, peptide_seq, charge)
if (exists("fit_stats") && has_name(fit_stats, "subobjects")) {
    sel_pepmodstates.df <- dplyr::left_join(sel_pepmodstates.df,
                                            dplyr::select(filter(fit_stats$subobjects, var=="suo_shift"),
                                                          pepmodstate_id, pms_median_log2 = median_log2)) %>%
        dplyr::arrange(object_id, desc(coalesce(pms_median_log2, min(pms_median_log2, na.rm=TRUE))), pepmod_seq, peptide_seq)
} else {
    sel_pepmodstates.df <- dplyr::arrange(sel_pepmodstates.df, object_id, pepmod_seq, peptide_seq) %>%
        mutate(pms_median_log2 = NA_real_)
}
sel_pepmodstates.df <- dplyr::mutate(sel_pepmodstates.df,
                                pepmod_ext = paste0(pepmod_seq, " (", pepmodstate_id, ").", charge,
                                                    if_else(is.na(pms_median_log2), "", "*"))
                                %>% factor(., levels=.))

sel_pepmod_intens.df <- tidyr::expand(msdata_full$pepmodstate_intensities, pepmodstate_id, msrun) %>%
    dplyr::inner_join(dplyr::select(sel_pepmodstates.df, object_id, object_label, protac_label, gene_names,
                                                        pepmodstate_id, pepmod_ext, pms_median_log2) %>%
                                              dplyr::distinct()) %>%
    dplyr::left_join(dplyr::select(msdata_full$pepmodstate_intensities, pepmodstate_id, msrun, intensity, intensity_norm, qvalue)) %>%
    dplyr::inner_join(distinct(select(msdata_full$msruns, sample, msrun, replicate, treatment, timepoint, condition)) %>%
                          dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
                          dplyr::arrange(treatment, timepoint, replicate) %>%
                          dplyr::mutate(msrun.2 = factor(msrun, levels=msdata$msruns$msrun))) %>%
    dplyr::mutate(msrun = msrun.2,
                  msrun.2 = NULL) %>%
    dplyr::group_by(pepmodstate_id) %>%
    dplyr::mutate(intensity_norm_trunc = pmin(intensity_norm, quantile(intensity_norm, 0.95, na.rm=TRUE))) %>%
    dplyr::ungroup()

group_by(sel_pepmod_intens.df, object_id) %>% do({
    shown_pepmod_intens.df <- mutate(., !!modelobj_idcol := object_id)#sel_pepmod_intens.df
    obj_label <- str_remove(shown_pepmod_intens.df$object_label[[1]], "\\.\\.\\.$")
    obj_id <- shown_pepmod_intens.df$object_id[[1]]
    is_viral <- semi_join(modelobjs_df, shown_pepmod_intens.df, by="object_id")$is_viral[[1]]
    message("Plotting ", obj_label)
    p <- ggplot(shown_pepmod_intens.df) +
        geom_tile(aes(x=msrun, y=pepmod_ext, fill=intensity_norm_trunc, color=qvalue), size=0.5, width=0.95, height=0.95) +
        theme_bw_ast(base_family = "", base_size = 10) +
        theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0)) +
        facet_grid(object_id + object_label ~ ., scales = "free_y", space="free_y") +
        guides(color=guide_legend("EG.qvalue", override.aes = list(fill=NA, size=2))) +
        ggtitle(str_c(obj_label, " (pg_id=", obj_id,
                      ", ac=", shown_pepmod_intens.df$protac_label[[1]], ") peptide map"),
                subtitle = semi_join(modelobjs_df, dplyr::select(shown_pepmod_intens.df, !!modelobj_idcol), by=modelobj_idcol)$protein_description[[1]]) +
        scale_x_discrete(drop=FALSE) +
        scale_fill_distiller(na.value="#00000000", type="div", palette = "Spectral", trans="log2") +
        scale_color_distiller(na.value="#00000000", type="div", palette = "Spectral", trans="mlog10")
    fname <- file.path(analysis_path, "plots", str_c(msfolder, "_", data_version), str_c("pepmod_heatmaps", if_else(is_viral, "/viral", "")),
                       paste0(project_id, "_", msfolder, '_', data_version, "_pepmod_heatmap_", obj_label, "_", obj_id, ".pdf"))
    #message(fname)
    ggsave(filename = fname,
           plot = p, width=8, height=3 + n_distinct(shown_pepmod_intens.df$object_id) + min(20, 0.1*n_distinct(shown_pepmod_intens.df$pepmodstate_id)),
           device=cairo_pdf, family="Arial")
    tibble()
})

sel_std_type <- "replicate"
treatment_palette <- c(mock="gray", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02")

dplyr::filter(modelobjs_df, is_viral) %>% do({
    sel_obj.df <- .
    sel_var <- if_else(sel_std_type == "median", "iaction_labu", "iaction_labu_replCI")
    sel_obj_iactions.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id, object_label, gene_label, protein_label),
                                             dplyr::filter(fit_stats$iactions, var == sel_var)) %>%
        dplyr::inner_join(conditions.df) %>%
        dplyr::mutate_at(vars(mean, ends_with("%")),
                         list(~exp(. + global_labu_shift))) %>%
        dplyr::filter((treatment == "SARS_CoV" & str_detect(protein_label, "_CVHSA$")) |
                      (treatment == "SARS_CoV2" & str_detect(protein_label, "_SARS2$")))
    if (modelobj == "protgroup") {
        sel_obj_msdata.df <- sel_obj.df %>%
            dplyr::inner_join(msdata$protgroup_intensities) %>%
            dplyr::inner_join(msdata$msruns)
    } else if (modelobj == "protregroup") {
        sel_obj_msdata.df <- sel_obj.df %>%
            dplyr::inner_join(msdata_full$protein2protregroup) %>%
            dplyr::group_by(object_id) %>%
            dplyr::filter(is_majority == any(is_majority)) %>%
            dplyr::ungroup() %>%
            dplyr::select(object_id, gene_label, protein_label, protregroup_id, protein_ac) %>%
            dplyr::inner_join(msdata_full$protein2protgroup) %>%
            dplyr::group_by(object_id, protgroup_id) %>%
            dplyr::filter(is_majority == any(is_majority)) %>%
            dplyr::ungroup() %>%
            dplyr::select(object_id, gene_label, protein_label, protgroup_id, protregroup_id) %>% dplyr::distinct() %>%
            dplyr::inner_join(msdata_full$protgroup_intensities) %>%
            dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
            dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) %>%
            dplyr::inner_join(msdata_full$msruns)
        #print(sel_obj_msdata.df)
        obj_shifts.df = dplyr::inner_join(sel_obj_msdata.df, sel_obj_iactions.df) %>%
            dplyr::group_by(object_id) %>%
            dplyr::summarise(obj_shift = median(log2(intensity_norm) - log2(`50%`), na.rm=TRUE)) %>%
            dplyr::ungroup()
        #print(obj_shifts.df)
        sel_obj_msdata.df <- dplyr::left_join(sel_obj_msdata.df, obj_shifts.df) %>%
            dplyr::mutate(intensity_norm = intensity_norm * 2^(-obj_shift)) %>%
            dplyr::filter((treatment == "SARS_CoV" & str_detect(protein_label, "_CVHSA$")) |
                          (treatment == "SARS_CoV2" & str_detect(protein_label, "_SARS2$")))
    }
    if (nrow(sel_obj_iactions.df) > 0) {
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
            scale_color_manual(values=treatment_palette) +
            scale_fill_manual(values=treatment_palette) +
            scale_y_log10() +
            ggtitle("Viral proteins timecourse") +
            facet_wrap( ~ gene_label, ncol=4, scales = "free")
        ggsave(p, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                   str_c("timecourse_", sel_std_type, modelobj_suffix), "viral",
                                   str_c(project_id, "_", msfolder, '_', fit_version, "_viral_comparison.pdf")),
               width=12, height=16, device = cairo_pdf)
    }
    tibble()
})
