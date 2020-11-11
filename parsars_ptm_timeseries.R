source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)

sel_std_type <- "median"
sel_objects.df <- dplyr::filter(modelobjs_df, is_viral)
sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(object_label, str_c(c("GlyGly_VIM_K334_M1", "GlyGly_EDF1_K93_M1",
                                                                  "GlyGly_EGFR_K1182_M1", "GlyGly_EGFR_K708_M1", "GlyGly_EGFR_K754_M1",
                                                                  "GlyGly_EGFR_K737_M1", "GlyGly_EGFR_K867_M1", "GlyGly_EGFR_K739_M1",
                                                                  "GlyGly_GABARAP_K13_M.?", "GlyGly_LAMTOR1_K20_M.?", "GlyGly_MAP1LC3A_K42_M.?",
                                                                  "GlyGly_VPS33A_K278_M.?", "GlyGly_VAMP8_K47_M.?"), collapse="|")))

sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(object_label, str_c(c("Phospho_EGFR_S991_M1",
                                                                  "Phospho_VIM_S5_M1", "Phospho_VIM_S34_M1", "Phospho_VIM_S56_M1", "Phospho_VIM_S72_M1",
                                                                  "Phospho_RAB7A_S72_M1", "Phospho_EGFR_T693_M1", "Phospho_SPAG9_T217_M1", "Phospho_PIK3R4_T927_M1",
                                                                  "Phospho_JUN_S63_M1", "Phospho_JUN_S73_M1", "Phospho_NCK2_Y110_M1",
                                                                  "Phospho_ITGB4_S1356_M1", "Phospho_SOS1_S1134_M1", "Phospho_MAPKAPK2_T334_M1",
                                                                  "Phospho_ACE2_Y781_M.?", "Phospho_ACE2_S787_M.?", "Phospho_VIM_S420_M.?",
                                                                  "Phospho_VIM_S34_M.?", "Phospho_SQSTM1_S272_M.?"
                                                                  ), collapse="|")))
sel_objects.df <- dplyr::semi_join(modelobjs_df, object_contrasts.df) #%>% dplyr::filter(object_id >= 12010)
sel_objects.df <- dplyr::semi_join(modelobjs_df,
                                   dplyr::distinct(
                                       bind_rows(
                                           #dplyr::select(dplyr::filter(object_effects.df, #std_type==sel_std_type & 
                                           #                            is_hit_nomschecks & effect_type == "treatmentXtimepoint"), object_id),
                                           dplyr::select(dplyr::filter(object_contrasts.df, #std_type==sel_std_type & 
                                                                       is_hit_nomschecks & contrast_kind == "treatment_vs_treatment"), object_id)
                                       )))

treatment_palette <- c(mock="gray", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02")
fp_treatment_palette <- c(mock="lightgray", SARS_CoV2 = "#DCDE71", SARS_CoV = "#BA6386")

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
            dplyr::mutate_at(vars("mean", ends_with("%")),
                             list(~exp(. + global_labu_shift)))
        sel_obj_msdata.df <- dplyr::select(sel_obj.df, object_id) %>%
            dplyr::inner_join(fit_stats$observations) %>%
            dplyr::filter(var == "obs_labu") %>%
            dplyr::mutate_at(vars("mean", ends_with("%")),
                             list(~exp(. + global_labu_shift))) %>%
            dplyr::inner_join(dplyr::select(msdata_full$msruns, msrun, condition, treatment, timepoint_num)) %>%
            dplyr::mutate(intensity_norm_scaled = `50%`)
        if (nrow(sel_obj_iactions.df) > 0) {
            h <- 3
            p <-
                ggplot(data=sel_obj_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
                geom_ribbon(aes(x = timepoint_num, ymin = `2.5%`, ymax=`97.5%`),
                            alpha=0.5, fill=NA, stat = "identity", linetype = "dotted", size=0.25) +
                geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                            alpha=0.5, stat = "identity", size=0.5) +
                geom_path(aes(x = timepoint_num, y = `50%`), alpha=0.5, size=1, stat="identity") +
                #geom_point(data=sel_obj_msdata.df,
                #           aes(y = intensity_norm_scaled),
                #           position = position_jitter(width = 0.75, height = 0), size=0.5, show.legend=FALSE) +
                theme_bw_ast(base_family = "", base_size = 12) +
                theme(panel.border = element_rect(color="lightgray", size = 1), panel.background = element_blank(), axis.line = element_blank()) +
                scale_x_continuous("Time, h.p.i.", breaks=unique(msdata$msruns$timepoint_num), minor_breaks = NULL,
                                   labels=scales::label_number(accuracy = 1)) +
                scale_color_manual("Treatment", values=treatment_palette, guide="none") +
                scale_fill_manual("Treatment", values=treatment_palette, guide="none") +
                scale_y_log10("PTM Intensity", n.breaks=8, minor_breaks = NULL, labels=scales::label_scientific())
            if (exists("fp.env")) {
                sel_fp_iactions.df <- semi_join(fp.env$fit_stats$iactions, semi_join(ptm2protregroup.df, sel_obj.df),
                                                by=c("protregroup_id"="fp_protregroup_id")) %>%
                    filter(var == if_else(sel_std_type == "median", "iaction_labu", "iaction_labu_replCI")) %>%
                    dplyr::inner_join(conditions.df)
                if (nrow(sel_fp_iactions.df) > 0L) {
                    sel_fp_iactions.df <- dplyr::mutate_at(sel_fp_iactions.df, vars("mean", ends_with("%")),
                                                list(~exp(. + fp.env$global_labu_shift)))
                    fp_plot <- ggplot(data=sel_fp_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
                        geom_ribbon(aes(x = timepoint_num, ymin = `2.5%`, ymax=`97.5%`),
                                    fill=NA, stat = "identity", linetype = "dotted", size=0.25) +
                        geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                                    fill=NA, stat = "identity", size=0.5) +
                        geom_ribbon(aes(x = timepoint_num, ymin = `25%`, ymax=`75%`),
                                    alpha=0.3, color=NA, stat = "identity", size=0.5) +
                        geom_path(aes(x = timepoint_num, y = `50%`), size=0.5, stat="identity") +
                        theme_bw_ast(base_family = "", base_size = 12) +
                        theme(panel.border = element_rect(color="lightgray", size = 1), panel.background = element_blank(), axis.line = element_blank()) +
                        scale_x_continuous("Time, h.p.i.", breaks=unique(msdata_full$msruns$timepoint_num),
                                           minor_breaks = NULL, labels=scales::label_number(accuracy = 1)) +
                        scale_color_manual(values=treatment_palette, guide="none") +
                        scale_fill_manual(values=treatment_palette, guide="none") +
                        scale_y_log10("Protein Intensity", n.breaks=5, labels=scales::label_scientific(), minor_breaks = NULL,
                                      expand = expansion(mult = 0.25, add = 0.25))
                    p <- ggpubr::ggarrange(p, fp_plot, ncol=1, align="h", common.legend=TRUE, heights=c(0.5, 0.5))
                    h <- h + 3
                }
            }
            p <- p + ggtitle(str_c(sel_obj.df$object_label, " timecourse"),
                             subtitle=str_c(sel_obj.df$protein_description, " (npepmodstates=", sel_obj.df$n_pepmodstates, ")"))
            plot_path <- file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                   str_c("timecourse4paper_", sel_ptm_type, "_", sel_std_type,
                                         modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
            if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
            ggsave(p, file = file.path(plot_path, str_c(project_id, "_", msfolder, '_', fit_version, "_",
                                                        str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], ".pdf")),
                   width=3, height=h, device = cairo_pdf)
        }
        tibble()
    })
