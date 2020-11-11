source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)

sel_std_type <- "median"
sel_objects.df <- dplyr::filter(modelobjs_df, is_viral)
sel_objects.df <- dplyr::filter(modelobjs_df, str_detect(object_label, str_c(c("EDF1", "SQSTM1", "RAB7A", "EGFR", "SPAG9", "PIK3R4", "FN1", "SERPINE1",
                                                                  "ACE2", "GABARAP", "LAMTOR1", "MAP1LC3A", "VPS33A", "JUN",
                                                                  "NCK2", "ITGB4", "SOS1", "MAPKAPK2", "SMAD3", "APOB"), collapse="|")))

treatment_palette <- c(mock="gray", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02")

sel_objects.df %>%
group_by(object_id) %>% do({
    sel_obj.df <- .
    message("Plotting ", sel_obj.df$object_label, " time course")
    sel_var <- if_else(sel_std_type == "median", "iaction_labu", "iaction_labu_replCI")
    sel_obj_iactions.df <- dplyr::inner_join(dplyr::select(sel_obj.df, object_id),
                                             dplyr::filter(fit_stats$iactions, var == sel_var)) %>%
        dplyr::inner_join(conditions.df) %>%
        dplyr::mutate_at(vars(mean, ends_with("%")),
                         list(~exp(. + global_labu_shift)))
    sel_obj_msdata.df <- dplyr::select(sel_obj.df, object_id) %>%
        dplyr::inner_join(fit_stats$observations) %>%
        dplyr::filter(var == "obs_labu") %>%
        dplyr::mutate_at(vars(mean, ends_with("%")),
                         list(~exp(. + global_labu_shift))) %>%
        dplyr::inner_join(dplyr::select(msdata_full$msruns, msrun, condition, treatment, timepoint_num)) %>%
        dplyr::mutate(intensity_norm_scaled = `50%`)
    if (nrow(sel_obj_iactions.df) > 0) {
        p <- ggplot(data=sel_obj_iactions.df, aes(x = timepoint_num, color=treatment, fill=treatment)) +
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
            scale_x_continuous("Time, h.p.i.", breaks=unique(msdata$msruns$timepoint_num), minor_breaks = NULL, labels=scales::label_number(accuracy = 1)) +
            scale_color_manual("Treatment", values=treatment_palette, guide="none") +
            scale_fill_manual("Treatment", values=treatment_palette, guide="none") +
            scale_y_log10("Protein Intensity", n.breaks=5, minor_breaks = NULL, labels=scales::label_scientific())
        plot_path <- file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                               str_c("timecourse4paper_", sel_std_type,
                                     modelobj_suffix, if_else(sel_obj.df$is_viral[[1]], "/viral", "")))
        if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
        ggsave(p, file = file.path(plot_path, str_c(project_id, "_", msfolder, '_', fit_version, "_",
                                                    str_replace(sel_obj.df$object_label[[1]], "/", "-"), "_", sel_obj.df$object_id[[1]], ".pdf")),
               width=3, height=3, device = cairo_pdf)
    }
    tibble()
})

sel_objects.df <- dplyr::filter(modelobjs_df, gene_names %in% c("S", "N")) %>%
    dplyr::mutate(organism_short = case_when(organism == "Severe acute respiratory syndrome coronavirus 2" ~ "SARS_CoV2",
                                             organism == "Human SARS coronavirus" ~ "SARS_CoV",
                                             TRUE ~ organism))
sel_var <- if_else(sel_std_type == "median", "iaction_labu", "iaction_labu_replCI")
sel_obj_iactions.df <- dplyr::inner_join(dplyr::select(sel_objects.df, object_id, organism_short),
                                         dplyr::filter(fit_stats$iactions, var == sel_var)) %>%
    dplyr::inner_join(conditions.df) %>%
    dplyr::mutate_at(vars(mean, ends_with("%")),
                     list(~exp(. + global_labu_shift))) %>%
    dplyr::filter(treatment == organism_short)
sel_obj_msdata.df <- dplyr::select(sel_objects.df, object_id) %>%
    dplyr::inner_join(fit_stats$observations) %>%
    dplyr::filter(var == "obs_labu") %>%
    dplyr::mutate_at(vars(mean, ends_with("%")),
                     list(~exp(. + global_labu_shift))) %>%
    dplyr::inner_join(dplyr::select(msdata_full$msruns, msrun, condition, treatment, timepoint_num)) %>%
    dplyr::mutate(intensity_norm_scaled = `50%`)
p <- ggplot(data=sel_obj_iactions.df, aes(x = timepoint_num, group=object_id, color=treatment, fill=treatment)) +
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
    scale_x_continuous("Time, h.p.i.", breaks=unique(msdata$msruns$timepoint_num), minor_breaks = NULL, labels=scales::label_number(accuracy = 1)) +
    scale_color_manual("Treatment", values=treatment_palette, guide="none") +
    scale_fill_manual("Treatment", values=treatment_palette, guide="none")# +
    scale_y_log10("Protein Intensity", n.breaks=5, minor_breaks = NULL, labels=scales::label_scientific())#
p
plot_path <- file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                       str_c("timecourse4paper_", sel_std_type))
if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
ggsave(p, file = file.path(plot_path, str_c(project_id, "_", msfolder, '_', fit_version, "_S_N_linear_growth.pdf")),
       width=3, height=3, device = cairo_pdf)
