source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)

bait_checks.df <- get(str_c("bait_checks_", modelobj, ".df"))
bait_checks.df$object_id <- bait_checks.df[[modelobj_idcol]]

narrow_bait_ids <- c()

object_contrasts_4show.df <- object_contrasts.df %>%
    dplyr::left_join(dplyr::select(bait_checks.df, object_id, obj_bait_full_id = bait_full_id, obj_organism = bait_organism)) %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
                  show_label = is_viral | is_hit_nomschecks,
                  median_log2_max = case_when(bait_full_id %in% narrow_bait_ids ~ 8,
                                              TRUE ~ 18),
                  mean_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, mean_log2)),
                  median_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, median_log2)),
                  truncation = volcano_truncation(median_log2, median_log2_trunc, p_value, p_value, is_hit, is_signif),
                  truncation_type = volcano_truncation_type(truncation, is_signif)) %>%
    dplyr::group_by(std_type, contrast) %>%
    dplyr::mutate(show_label = if_else(rep.int(sum(show_label) >= 300L, n()), rep.int(FALSE, n()), show_label)) %>%
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
                                       sel_object_contrast_thresholds.df$contrast[[1]], '.pdf')),
           plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
    tibble()
})

object_iactions_4show.df <- filter(object_contrasts.df, contrast_type=="comparison" & std_type=="replicate") %>%
    dplyr::select(std_type, contrast, contrast_median_log2 = median_log2, contrast_p_value = p_value,
                  is_signif, is_hit_nomschecks, is_hit, is_viral, object_id, object_label) %>%
    dplyr::inner_join(dplyr::select(dplyr::filter(contrastXcondition.df, weight>0), contrast, condition_lhs=condition, bait_full_id_lhs=bait_full_id)) %>%
    dplyr::inner_join(dplyr::select(dplyr::filter(contrastXcondition.df, weight<0), contrast, condition_rhs=condition, bait_full_id_rhs=bait_full_id)) %>%
    dplyr::mutate(contrast_lhs = str_c(bait_full_id_lhs, "_vs_others"),
                  contrast_rhs = str_c(bait_full_id_rhs, "_vs_others")) %>%
    dplyr::left_join(dplyr::filter(object_contrasts.df, contrast_type=="filter" & std_type == "replicate") %>%
                     dplyr::select(contrast_lhs = contrast, object_id, contrast_lhs_p_value = p_value, contrast_lhs_median_log2 = median_log2,
                                   is_signif_lhs = is_signif, is_hit_nomschecks_lhs = is_hit_nomschecks, is_hit_lhs = is_hit)) %>%
    dplyr::left_join(dplyr::filter(object_contrasts.df, contrast_type=="filter" & std_type == "replicate") %>%
                     dplyr::select(contrast_rhs = contrast, object_id, contrast_rhs_p_value = p_value, contrast_rhs_median_log2 = median_log2,
                                   is_signif_rhs = is_signif, is_hit_nomschecks_rhs = is_hit_nomschecks, is_hit_rhs = is_hit)) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, var=="iaction_labu_replCI"),
                                      object_id, condition_lhs=condition, lhs_median_log2=median_log2 + global_pepmodstate_labu_shift)) %>%
    dplyr::left_join(dplyr::transmute(dplyr::filter(fit_stats$iactions, var=="iaction_labu_replCI"),
                                      object_id, condition_rhs=condition, rhs_median_log2=median_log2 + global_pepmodstate_labu_shift)) %>%
    dplyr::mutate(is_hilite = is_hit_nomschecks_lhs | is_hit_nomschecks_rhs,
                  show_label = (is_hit_nomschecks_lhs | is_hit_nomschecks_rhs) & (is_hit_nomschecks | is_viral)) %>%
    dplyr::filter(is_hilite | lhs_median_log2 >= quantile(lhs_median_log2, na.rm = TRUE, 0.01) &
                      rhs_median_log2 >= quantile(rhs_median_log2, na.rm = TRUE, 0.01))

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
                    vjust=-0.5,
                    size=3.5,
                    show.legend = FALSE, segment.color = "gray") +
    scale_color_manual(values=c("TRUE" = "red", "FALSE" = "black"), na.value="black") +
    scale_shape_manual(values=c("TRUE" = 16L, "FALSE" = 1L)) +
    scale_size_manual(values=c("TRUE" = 2, "FALSE" = 1)) +
    coord_fixed() +
    theme_bw_ast()
    ggsave(filename = file.path(analysis_path, 'plots', str_c(ms_folder,'_', fit_version),
                                str_c("scatter_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                paste0(project_id, '_', fit_version, '_scatter_',
                                       sel_object_contrast_thresholds.df$contrast[[1]], '.pdf')),
           plot = p, width=10, height=10, device=cairo_pdf, family="Arial")
    tibble()
})

object_contrasts.df %>% #filter(std_type == "replicate") %>%
    dplyr::left_join(dplyr::select(bait_checks.df, object_id, obj_bait_full_id = bait_full_id, obj_organism = bait_organism)) %>%
    dplyr::mutate(p_value_compressed = 10^(-sapply(-log10(p_value), mlog10_pvalue_compress)),
                  p_value_capped = pmax(1E-20, p_value),
                  p_value_range = if_else(p_value <= 1E-7, "high", "low"),
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
        ggsave(filename = file.path(analysis_path, 'plots', str_c(ms_folder,'_', fit_version),
                                    str_c("volcanos_contrasts_", sel_object_contrast_thresholds.df$std_type[[1]], modelobj_suffix),
                                    paste0(project_id, '_', fit_version, '_volcano_',
                                           sel_object_contrast_thresholds.df$contrast[[1]], '.pdf')),
               device=cairo_pdf, plot = p, width=15, height=18, device=cairo_pdf, family="Arial")
        tibble()
    })

filter(msdata_full$proteins, str_detect(gene_name, "(^|;)EXO")) %>%
    dplyr::inner_join(msdata_full$protein2object) %>% View()

# plot peptide heatmaps
candidate_genes <- c("ATL2", "MAVS", "AIFM1", "SHC1", "NELFB", "ILF3", "C14orf166", "PRAF2", "UNC93B1",
                     "EIF4H", "G3BP1", "PIP4K2C", "HLA-E", "TNFAIP2", "MARC1", "JAK1",
                     "TNFRSF10B", "TNFRSF10A", "TNFRSF21", "IFITM10", "SQSTM1")
sel_objects.df <- dplyr::semi_join(modelobjs_df,
                                   dplyr::filter(msdata_full$proteins, gene_name %in% candidate_genes) %>%
                                   dplyr::inner_join(msdata_full[[str_c("protein2", modelobj)]]) %>%
                                   dplyr::select(!!modelobj_idcol))
sel_objects.df <- dplyr::semi_join(modelobjs_df,
                                   dplyr::filter(msdata_full$proteins, str_detect(gene_name, "(^|;)(CD320)($|;)")) %>%
                                   dplyr::inner_join(msdata_full[[str_c("protein2", modelobj)]]) %>%
                                   dplyr::select(!!modelobj_idcol))
sel_objects.df <- dplyr::semi_join(modelobjs_df, dplyr::select(iactions_4graphml_pre.df, !!modelobj_idcol := object_id))
sel_objects.df <- dplyr::filter(modelobjs_df, is_viral)

sel_pepmodstates.df <- dplyr::inner_join(sel_objects.df, msdata_full[[str_c(modelobj, "2pepmod")]]) %>%
    dplyr::filter(is_specific) %>%
    dplyr::inner_join(msdata_full$pepmodstates) %>%
    dplyr::inner_join(select(msdata_full$pepmods, pepmod_id, peptide_id, mod_seq=seq)) %>%
    dplyr::inner_join(select(msdata_full$peptides, peptide_id, unmod_seq=seq)) %>%
    dplyr::select(object_id, pepmod_id, majority_protein_acs, protac_label, gene_label, gene_names,
                  mod_seq, unmod_seq, charge, pepmodstate_id)
if (exists("fit_stats") && has_name(fit_stats, "subobjects")) {
    sel_pepmodstates.df <- dplyr::left_join(sel_pepmodstates.df,
                                             dplyr::select(filter(fit_stats$subobjects, var=="suo_shift_unscaled"),
                                                           pepmodstate_id, pms_median_log2 = median_log2)) %>%
        dplyr::arrange(object_id, desc(coalesce(pms_median_log2, min(pms_median_log2, na.rm=TRUE))), unmod_seq, mod_seq, charge)
} else {
    sel_pepmodstates.df <- dplyr::arrange(sel_pepmodstates.df, object_id, unmod_seq, mod_seq, charge) %>%
        mutate(pms_median_log2 = NA_real_)
}
sel_pepmodstates.df <- dplyr::mutate(sel_pepmodstates.df,
                                     # FIXME remove pepmod_id once mod_seq is really mod_seq
                                     pepmodstate_ext = paste0(mod_seq, ".", charge, " (", pepmod_id, ")", if_else(is.na(pms_median_log2), "", "*"))
                                                              %>% factor(., levels=.))

sel_pepmod_intens.df <- dplyr::inner_join(dplyr::select(sel_pepmodstates.df, object_id, protac_label, gene_label, gene_names,
                                                        pepmod_id, pepmodstate_id, charge, pepmodstate_ext, pms_median_log2) %>%
                                                 dplyr::distinct(),
                                             msdata_full$pepmodstate_intensities) %>%
    dplyr::inner_join(distinct(select(msdata$msruns, msrun, batch, replicate, bait_type, bait_full_id, bait_id, condition)) %>%
                      dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
                      dplyr::arrange(bait_type, bait_id, bait_full_id, batch, replicate) %>%
                      dplyr::mutate(msrun.2 = factor(msrun, levels=msrun))) %>%
    dplyr::mutate(intensity_norm = exp(-total_msrun_shift)*intensity,
                  msrun = msrun.2,
                  msrun.2 = NULL) %>%
    dplyr::group_by(pepmodstate_id) %>%
    dplyr::mutate(intensity_norm_trunc = pmin(intensity_norm, quantile(intensity_norm, 0.95, na.rm=TRUE))) %>%
    dplyr::ungroup()

sel_pepmods.df <- dplyr::group_by(sel_pepmod_intens.df, pepmod_id) %>%
    dplyr::summarize(n_pepmod_quants = sum(!is.na(intensity)),
                     median_quant = median(intensity, na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(n_pepmod_quants), desc(median_quant), pepmod_id)

group_by(sel_pepmod_intens.df, object_id) %>% do({
    shown_pepmod_intens.df <- mutate(., !!modelobj_idcol := object_id)#sel_pepmod_intens.df
    gene_name <- str_remove(shown_pepmod_intens.df$gene_label[[1]], "\\.\\.\\.$")
    obj_id <- shown_pepmod_intens.df$object_id[[1]]
    is_viral <- inner_join(shown_pepmod_intens.df, modelobjs_df)$is_viral[[1]]
    if (is_viral) {
        bait_df <- filter(bait_checks.df, object_id == obj_id)
        if (nrow(bait_df) > 0) {
            gene_name <- bait_df$bait_full_id[[1]]
        }
    }
    message("Plotting ", gene_name)
p <- ggplot(shown_pepmod_intens.df) +
    geom_tile(aes(x=msrun, y=pepmodstate_ext,
                  fill=intensity_norm, color=ident_type), size=0.5, width=0.85, height=0.85) +
    theme_bw_ast(base_family = "", base_size = 10) +
    theme(axis.text.x = element_text(angle = -90, hjust=0, vjust=0)) +
    facet_grid(object_id + gene_label ~ ., scales = "free_y", space="free_y") +
    guides(color=guide_legend("ident_type", override.aes = list(fill=NA, size=2))) +
    ggtitle(str_c(gene_name,  " (pg_id=", obj_id,
                  ", ac=", shown_pepmod_intens.df$protac_label[[1]], ") peptide map"),
            subtitle = semi_join(modelobjs_df, dplyr::select(shown_pepmod_intens.df, !!modelobj_idcol))$protein_description[[1]]) +
    scale_fill_distiller(na.value="#00000000", type="div", palette = "Spectral") +
    scale_color_manual(na.value="#00000000",
                       values=c("MULTI-MSMS"="black", "MULTI-MATCH-MSMS"="khaki",
                                "MSMS"="cornflowerblue", "MULTI-SECPEP"="firebrick",
                                "MULTI-MATCH"="gray"))
  fname <- file.path(analysis_path, "plots", str_c(ms_folder, "_", data_version), str_c("peptide_heatmaps", if_else(is_viral, "/viral", "")),
                     paste0(project_id, "_", ms_folder, '_', data_version, "_pepmod_heatmap_", gene_name, "_", obj_id, ".pdf"))
  message(fname)
ggsave(filename = fname,
       plot = p, width=30, height=3 + n_distinct(shown_pepmod_intens.df$object_id) + min(20, 0.1*n_distinct(shown_pepmod_intens.df$pepmodstate_id)),
       device=cairo_pdf, family="Arial")
    tibble()
})

# bait expression
bait_intensities.df = dplyr::select(bait_checks.df, bait_full_id, bait_id, bait_orgcode, bait_object_id=object_id) %>%
    dplyr::inner_join(select(conditions.df, condition, bait_type, bait_full_id)) %>%
    dplyr::left_join(select(filter(fit_stats$iactions, var=="iaction_labu_replCI"),
                            condition, object_id, median_log2, ends_with("%"), bait_full_id = object_label)) %>%
    dplyr::arrange(bait_type, bait_id, bait_orgcode) %>%
    dplyr::mutate(bait_full_id = factor(bait_full_id, levels=unique(bait_full_id)),
                  bait_orgcode = if_else(bait_type == "control", "control", as.character(bait_orgcode)),
                  median_log2 = median_log2 + global_labu_shift) %>%
    dplyr::mutate_at(vars(ends_with("%")), ~./log(2)+global_labu_shift) %>%
    dplyr::filter(bait_full_id != "Ctrl_Gaussia_luci")

bait_orgcode_palette <- c("SARS2"="red", "CVHSA"="orange", "CVHNL"="darkgreen", "CVH22"="darkgreen", "HUMAN"="darkblue", "control"="gray")
p <- ggplot(bait_intensities.df,
       aes(x = bait_full_id, group = bait_full_id, color=bait_orgcode, fill=bait_orgcode)) +
    geom_boxplot(aes(middle=median_log2,
                     lower=`25%`, upper=`75%`, ymin = `2.5%`, ymax=`97.5%`), stat="identity", alpha=0.5) +
    scale_color_manual(values=bait_orgcode_palette) +
    scale_fill_manual(values=bait_orgcode_palette) +
    theme_bw_ast() +
    theme(axis.text.x = element_text(angle=-45, hjust=0))
ggsave(filename = file.path(analysis_path, "plots", str_c(ms_folder, "_", data_version),
                            paste0(project_id, "_", ms_folder, '_', fit_version, "_bait_expression.pdf")),
       plot = p, width=15, height=8, device=cairo_pdf, family="Arial")
