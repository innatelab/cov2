project_id <- 'cov2'
plot_version <- "20201029c"
message('Project ID=', project_id)

datasets_info <- list(
    ubi = list(
        msfolder = "snaut_parsars_ptm_20200907",
        subdataset = "ubi",
        data_version = "20201012",
        fit_version = "20201012",
        ptm_type = "GlyGly",
        ptm_extractor_version = "20201012"
    ),
    phospho = list(
        msfolder = "snaut_parsars_phospho_20201005",
        subdataset = NA,
        data_version = "20201012",
        fit_version = "20201012",
        ptm_type = "Phospho",
        ptm_extractor_version = "20201012"
    )
)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

require(dplyr)
require(readr)
require(stringr)
require(rlang)
require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)

source(file.path(misc_scripts_path, 'ggplot_ext.R'))
source(file.path(project_scripts_path, 'cov2_plots_common.R'))

ptm_envs <- lapply(datasets_info, function(dsinfo) {
    message("Loading ", dsinfo$msfolder, "...")
    ptm_env <- new_environment()
    load(file.path(scratch_path, str_c(project_id, '_msdata_full_', dsinfo$msfolder, '_', dsinfo$data_version, '.RData')), envir=ptm_env)
    load(file.path(scratch_path, str_c(project_id, '_msglm_data_', dsinfo$msfolder, '_', dsinfo$data_version, '.RData')), envir=ptm_env)
    load(file.path(scratch_path, str_c(project_id, '_msglm_fit_', dsinfo$msfolder, '_', dsinfo$fit_version, '.RData')), envir=ptm_env)
    return (ptm_env)
})

for (ds in names(datasets_info)) {
    dsinfo <- datasets_info[[ds]]
    ptm_envs[[ds]]$ptmns_grouped.df <- read_tsv(file.path(data_path, dsinfo$msfolder, str_c("ptm_extractor_", dsinfo$data_version), "ptmns_grouped.txt.gz"),
                             col_types=cols(
                                 .default = col_guess(),
                                 ptmn_id = "i", ptm_id = "i", nselptms = "i",
                                 ptmgroup_id = "i", ptmngroup_id = "i",
                                 ptm_pos = "i", ptm_AA_seq = "c", ptm_is_known = col_logical(),
                                 protein_rank = "i",
                                 ptmid_is_reference = col_logical()
                             )) %>%
    dplyr::group_by(ptmgroup_id) %>%
    dplyr::mutate(nptmngroups = n_distinct(ptmngroup_id),
                  msfolder = dsinfo$msfolder) %>%
    dplyr::group_by(ptm_id) %>%
    dplyr::mutate(nptmns = n_distinct(ptmn_id)) %>%
    dplyr::ungroup()
}
ptmns_grouped.df <- bind_rows(lapply(names(datasets_info), function(ds) ptm_envs[[ds]]$ptmns_grouped.df))

msdata_full <- set_names(lapply(names(ptm_envs$ubi$msdata_full), function(dfname) {
    bind_rows(lapply(names(ptm_envs), function(dsname) {
        env <- ptm_envs[[dsname]]
        df <- env$msdata_full[[dfname]]
        if (has_name(df, "ptmn_id") && !has_name(df, "ptm_type")) {
            message("Adding ptm_type to ", dfname)
            df <- dplyr::left_join(df, dplyr::select(env$msdata_full$ptmns, ptmn_id, ptm_type))
        }
        if (!is.na(datasets_info[[dsname]]$subdataset) && has_name(df, "dataset")) {
            message("Filtering for dataset in ", dfname)
            df <- dplyr::filter(df, dataset == datasets_info[[dsname]]$subdataset)
        }
        if (dfname %in% c("pepmodstates", "msruns", "pepmodstate_intensities", "ptmn2pepmodstate", "ptmns", "ptmn_stats", "ptm2gene")) {
            message("Adding msfolder to ", dfname)
            df$msfolder <- datasets_info[[dsname]]$msfolder
        }
        return(df)
    }))
}), names(ptm_envs$ubi$msdata_full))

object_iactions.df <- bind_rows(lapply(names(ptm_envs), function(dsname) {
    dplyr::mutate(ptm_envs[[dsname]]$fit_stats$iactions,
                  msfolder = ptm_envs[[dsname]]$msfolder)
}))

conditions.df <- bind_rows(lapply(names(ptm_envs), function(dsname) {
    dplyr::mutate(ptm_envs[[dsname]]$conditions.df,
                  msfolder = ptm_envs[[dsname]]$msfolder)
}))

global_labu_shift <- ptm_envs$ubi$global_pepmodstate_labu_shift
msfolder <- datasets_info$ubi$msfolder

proteins.df <- ptm_envs$ubi$msdata_full$proteins %>%
    dplyr::mutate(organism_short = case_when(organism == "Severe acute respiratory syndrome coronavirus 2 OX=2697049" ~ "SARS_CoV2",
                                             organism == "Human SARS coronavirus OX=694009" ~ "SARS_CoV",
                                             TRUE ~ organism))

ptm_extractor_folder <- str_c("ptm_extractor_", datasets_info$ubi$ptm_extractor_version)
viral_agn.df <- read_tsv(file.path(data_path, datasets_info$ubi$msfolder, ptm_extractor_folder,
                         str_c("viral_agn_", datasets_info$ubi$ptm_extractor_version, ".txt"))) %>%
    left_join(dplyr::select(proteins.df, protein_ac, is_viral, is_contaminant, genename, organism_short))
viral_ptm2ptm.df <- bind_rows(lapply(datasets_info, function(dsinfo) {
    read_tsv(file.path(data_path, dsinfo$msfolder, str_c("ptm_extractor_", dsinfo$ptm_extractor_version),
                       str_c("viral_ptm2ptm_", dsinfo$ptm_extractor_version, ".txt"))) %>%
        dplyr::filter(ptm_type == dsinfo$ptm_type) %>%
        dplyr::mutate(msfolder = dsinfo$msfolder)
    })) %>%
    left_join(dplyr::select(proteins.df, protein_ac, is_viral, is_contaminant, genename, organism_short)) %>%
    left_join(dplyr::rename_all(dplyr::select(proteins.df, protein_ac, is_viral, is_contaminant, genename), ~ str_c("hom_", .))) %>%
    mutate(protein_ac_pair = if_else(protein_ac < hom_protein_ac,
                                     str_c(protein_ac, "-", hom_protein_ac),
                                     str_c(hom_protein_ac, "-", protein_ac)),
           genename_short = str_remove(genename, "SARS_CoV\\d?_"),
           hom_genename_short = str_remove(hom_genename, "SARS_CoV\\d?_"),
           genename_pair = if_else(is.na(hom_genename_short) | (genename_short == hom_genename_short),
                                   genename_short,
                                   if_else(genename_short < hom_genename_short,
                                           str_c(genename_short, "_", hom_genename_short),
                                           str_c(hom_genename_short, "_", genename_short))),
           ptm_correct = (ptm_type == "GlyGly" & ptm_AA_seq %in% "K") |
               (ptm_type == "Phospho" & ptm_AA_seq %in% c("S","T","Y")),
           hom_ptm_correct = (ptm_type == "GlyGly" & hom_ptm_AA %in% "K") |
               (ptm_type == "Phospho" & hom_ptm_AA %in% c("S","T","Y")),
           hom_organism_short = case_when(hom_organism == "Severe acute respiratory syndrome coronavirus 2 OX=2697049" ~ "SARS_CoV2",
                                          hom_organism == "Human SARS coronavirus OX=694009" ~ "SARS_CoV",
                                          TRUE ~ hom_organism))

perseus_report.df <- read_tsv(file.path(data_path, datasets_info$phospho$msfolder,
                                       "cov2_snaut_parsars_phospho_20201005_contrasts_report_20201012_long_psp.txt"),
                             comment="#", guess_max=5000)

viral_perseus_report.df <- read_tsv(file.path(data_path, datasets_info$phospho$msfolder, ptm_extractor_folder,
                                        "viral_aa_all_20201012_ext_PSP.txt"), col_types = cols("ptm_pos" = "i"))

viral_perseus_report.df2 <- left_join(dplyr::mutate(viral_perseus_report.df, row_ix=row_number()),
          dplyr::transmute(ptmns_grouped.df, ptmgroup_id, protein_ac, ptm_site=str_c(ptm_AA_seq, ptm_pos))) %>%
    dplyr::arrange(row_ix, desc(!is.na(ptmgroup_id))) %>%
    dplyr::group_by(row_ix) %>%
    dplyr::filter(row_number() == 1L) %>%
    dplyr::ungroup() %>%
    dplyr::select(-row_ix)

ptm_qvalue_ident_max <- 1E-3
ptm_qvalue_quant_max <- 1E-2
ptm_locprob_ident_min <- 0.75
ptm_locprob_quant_min <- 0.5

ptmnXtreatment_stats.df <- dplyr::select(msdata_full$ptmns, msfolder, is_viral, ptm_type, ptmn_id, ptm_id) %>% dplyr::filter(is_viral) %>%
    dplyr::left_join(msdata_full$msruns) %>%
    dplyr::left_join(msdata_full$ptmn2pepmodstate) %>%
    dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
    dplyr::left_join(dplyr::select(msdata_full$ptmn_locprobs, evidence_id, ptmn_id, dataset, ptm_locprob)) %>%
    dplyr::group_by(ptm_type, msfolder, ptm_id, ptmn_id, treatment) %>%
    dplyr::summarise(ptm_qvalue_min = min(psm_qvalue, na.rm=TRUE),
                     ptm_locprob_max = max(ptm_locprob, na.rm=TRUE),
                     n_idented_and_localized  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_qvalue, 1) <= ptm_qvalue_ident_max) &
                                                                                             (coalesce(ptm_locprob, 0) >= ptm_locprob_ident_min)]),
                     n_quanted  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_qvalue, 1) <= ptm_qvalue_quant_max) &
                                                                               (coalesce(ptm_locprob, 0) >= ptm_locprob_quant_min)]),
                     .groups="drop")

ptmnXcondition_stats.df <- dplyr::select(msdata_full$ptmns, msfolder, is_viral, ptm_type, ptmn_id, ptm_id) %>% dplyr::filter(is_viral) %>%
    dplyr::left_join(msdata_full$msruns) %>%
    dplyr::left_join(msdata_full$ptmn2pepmodstate) %>%
    dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
    dplyr::left_join(dplyr::select(msdata_full$ptmn_locprobs, evidence_id, ptmn_id, dataset, ptm_locprob)) %>%
    dplyr::group_by(ptm_type, msfolder, ptm_id, ptmn_id, condition) %>%
    dplyr::summarise(ptm_qvalue_min = min(psm_qvalue, na.rm=TRUE),
                     ptm_locprob_max = max(ptm_locprob, na.rm=TRUE),
                     n_idented_and_localized  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_qvalue, 1) <= ptm_qvalue_ident_max) &
                                                                                                 (coalesce(ptm_locprob, 0) >= ptm_locprob_ident_min)]),
                     n_quanted  = n_distinct(str_c(msrun, " ", pepmodstate_id)[(coalesce(psm_qvalue, 1) <= ptm_qvalue_quant_max) &
                                                                                   (coalesce(ptm_locprob, 0) >= ptm_locprob_quant_min)]),
                     .groups="drop")

ptm_status_levels <- c("N/A", "potential", "low conf.", "observed")

viral_ptmn_aligned.df <- select(viral_ptm2ptm.df, -agn_match_ratio) %>% distinct() %>%
    dplyr::mutate(treatment = organism_short, hom_treatment = hom_organism_short) %>%
    dplyr::left_join(distinct(select(msdata_full$ptm2gene, msfolder, ptm_type, ptm_id, protein_ac, is_viral, is_contaminant, ptm_pos, ptm_AA_seq))) %>%
    dplyr::left_join(dplyr::select(ptmns_grouped.df, msfolder, ptmn_id, ptm_id, ptmgroup_id, ptmngroup_id)) %>%
    dplyr::left_join(distinct(select(msdata_full$ptm2gene, msfolder, ptm_type, hom_ptm_id=ptm_id, hom_protein_ac=protein_ac, hom_ptm_pos=ptm_pos, hom_ptm_AA=ptm_AA_seq))) %>%
    dplyr::left_join(select(msdata_full$ptmns, msfolder, ptm_type, ptmn_id, ptm_id)) %>%
    dplyr::left_join(dplyr::select(ptmns_grouped.df, msfolder, hom_ptm_id=ptm_id, hom_ptmgroup_id=ptmgroup_id) %>% distinct()) %>%
    #dplyr::left_join(select(ptmns.df, ptm_type, hom_ptmn_id=ptmn_id, hom_ptm_id=ptm_id)) %>%
    dplyr::left_join(dplyr::select(ptmnXtreatment_stats.df, msfolder, ptm_type, treatment, ptmn_id, n_idented_and_localized, n_quanted,
                                   pms_qvalue_min=ptm_qvalue_min, ptm_locprob_max)) %>%
    dplyr::left_join(dplyr::select(ptmnXtreatment_stats.df, msfolder, ptm_type, hom_treatment=treatment, hom_ptm_id=ptm_id,
                                   hom_n_idented_and_localized = n_idented_and_localized, hom_n_quanted = n_quanted,
                                   hom_pms_qvalue_min=ptm_qvalue_min, hom_ptm_locprob_max=ptm_locprob_max) %>%
                         dplyr::group_by(msfolder, ptm_type, hom_treatment, hom_ptm_id) %>% # best stats for hom_ptm_id (across multiplicity)
                         dplyr::summarise(hom_pms_qvalue_min=min(hom_pms_qvalue_min, na.rm=TRUE),
                                          hom_ptm_locprob_max=max(hom_ptm_locprob_max, na.rm=TRUE),
                                          hom_n_idented_and_localized = max(hom_n_idented_and_localized, na.rm=TRUE),
                                          hom_n_quanted = max(hom_n_quanted, na.rm=TRUE),
                                          .groups="drop")) %>%
    dplyr::mutate(is_observed = coalesce(n_idented_and_localized, 0) > 0 &
                                coalesce(ptm_locprob_max, 0) >= ptm_locprob_min &
                                coalesce(pms_qvalue_min, 1) <= ptm_qvalue_max,
                  is_observed_lowconf = coalesce(n_quanted) > 0,
                  ptm_status = case_when(is_observed ~ "observed",
                                         is_observed_lowconf ~ "low conf.",
                                         ptm_correct ~ "potential",
                                         TRUE ~ "N/A") %>% factor(ordered=TRUE, levels=ptm_status_levels)) %>%
    dplyr::filter(pmin(coalesce(pms_qvalue_min, 1), coalesce(hom_pms_qvalue_min, 1)) <= ptm_qvalue_max &
                  pmax(coalesce(ptm_locprob_max, 0), coalesce(hom_ptm_locprob_max, 0)) >= ptm_locprob_min &
                  pmax(coalesce(n_idented_and_localized, 0), coalesce(hom_n_idented_and_localized, 0)) > 0)
viral_ptmn_aligned_ext.df <- dplyr::left_join(viral_ptmn_aligned.df,
                                              dplyr::filter(dplyr::select(viral_perseus_report.df, protein_ac, ptm_pos, flanking_15AAs, Motifs), !is.na(flanking_15AAs)) %>% dplyr::distinct())
    
write_tsv(viral_ptmn_aligned_ext.df, file.path(analysis_path, "reports", str_c("viral_ptms_aligned_", plot_version, ".txt")), na = "")

viral_ptmn_aligned_new2old.df <- bind_rows(
    left_join(select(viral_ptmn_aligned.df, ptm_type, ptm_id, ptmn_id, genename, protein_ac, ptm_pos, ptm_status),
              select(viral_ptmn_aligned_old.df, ptm_type, ptm_id, ptmn_id_old=ptmn_id, ptm_status_old=ptm_status,
                     genename_old=genename, protein_ac_old=protein_ac, ptm_pos_old=ptm_pos)),
    right_join(select(viral_ptmn_aligned.df, ptm_type, ptm_id, ptmn_id, genename, protein_ac, ptm_pos, ptm_status),
               select(viral_ptmn_aligned_old.df, ptm_type, ptm_id, ptmn_id_old=ptmn_id, ptm_status_old=ptm_status,
                      genename_old=genename, protein_ac_old=protein_ac, ptm_pos_old=ptm_pos))) %>%
    dplyr::distinct() %>% dplyr::filter(!is.na(ptm_id))
View(filter(viral_ptmn_aligned_new2old.df, coalesce(ptm_status, "") != coalesce(ptm_status_old, "")))
write_tsv(filter(viral_ptmn_aligned_new2old.df, coalesce(ptm_status, "") != coalesce(ptm_status_old, "")),
          file.path(analysis_path, "reports", str_c("viral_ptms_aligned_new2old_diff_", plot_version, ".txt")), na = "")

viral_ptmgroup2ptmgroup.df <- bind_rows(
dplyr::left_join(dplyr::select(dplyr::filter(viral_ptmn_aligned.df, organism_short == "SARS_CoV2"),
                               ptm_type, msfolder, genename_pair, protein_ac_pair, CoV2_ptm_id = ptm_id, CoV2_ptm_status = ptm_status, SARS_ptm_id = hom_ptm_id),
                 dplyr::select(dplyr::filter(viral_ptmn_aligned.df, organism_short == "SARS_CoV"),
                               ptm_type, msfolder, genename_pair, protein_ac_pair, SARS_ptm_id = ptm_id, SARS_ptm_status = ptm_status, CoV2_ptm_id = hom_ptm_id)),
dplyr::right_join(dplyr::select(dplyr::filter(viral_ptmn_aligned.df, organism_short == "SARS_CoV2"),
                               ptm_type, msfolder, genename_pair, protein_ac_pair, CoV2_ptm_id = ptm_id, CoV2_ptm_status = ptm_status, SARS_ptm_id = hom_ptm_id),
                 dplyr::select(dplyr::filter(viral_ptmn_aligned.df, organism_short == "SARS_CoV"),
                               ptm_type, msfolder, genename_pair, protein_ac_pair, SARS_ptm_id = ptm_id, SARS_ptm_status = ptm_status, CoV2_ptm_id = hom_ptm_id))
) %>%
    dplyr::left_join(dplyr::select(ptmns_grouped.df, ptm_type, msfolder, CoV2_ptmgroup_id = ptmgroup_id, CoV2_ptm_id = ptm_id)) %>%
    dplyr::left_join(dplyr::select(ptmns_grouped.df, ptm_type, msfolder, SARS_ptmgroup_id = ptmgroup_id, SARS_ptm_id = ptm_id)) %>%
    dplyr::select(-CoV2_ptm_id, -SARS_ptm_id) %>%
    dplyr::mutate_at(vars(CoV2_ptm_status, SARS_ptm_status), ~factor(ifelse(.x == "low conf.", "observed", as.character(.x)),
                                                                     levels=ptm_status_levels, ordered=TRUE)) %>%
    distinct() %>%
    dplyr::group_by(ptm_type, msfolder, genename_pair, protein_ac_pair, CoV2_ptmgroup_id, SARS_ptmgroup_id) %>%
    dplyr::summarise(CoV2_ptm_status = max(CoV2_ptm_status, na.rm=TRUE),
                     SARS_ptm_status = max(SARS_ptm_status, na.rm=TRUE),
                     .groups="drop") %>%
    dplyr::ungroup()

viral_ptmgroup2ptmgroup_total_stats.df <- group_by(viral_ptmgroup2ptmgroup.df, ptm_type, msfolder, CoV2_ptm_status, SARS_ptm_status) %>%
    dplyr::summarise(n = n())
write_tsv(viral_ptmgroup2ptmgroup_total_stats.df,
          file.path(analysis_path, "reports", str_c("viral_ptmgroup2ptmgroup_total_stats_", plot_version, ".txt")))
viral_ptmgroup2ptmgroup_gene_stats.df <- group_by(viral_ptmgroup2ptmgroup.df, ptm_type, msfolder, genename_pair, protein_ac_pair, CoV2_ptm_status, SARS_ptm_status) %>%
    dplyr::summarise(n = n())
write_tsv(viral_ptmgroup2ptmgroup_gene_stats.df,
          file.path(analysis_path, "reports", str_c("viral_ptmgroup2ptmgroup_gene_stats_", plot_version, ".txt")))

viral_ptmn_intensity.df <- dplyr::full_join(viral_ptmn_aligned.df, conditions.df) %>%
    #dplyr::left_join(dplyr::select(ptmnXcondition_stats.df, ptm_type, msfolder, condition,
    #                               ptm_cond_qvalue_min = ptm_qvalue_min,
    #                               ptm_cond_locprob_max = ptm_locprob_max,
    #                               n_cond_idented_and_localized = n_idented_and_localized)) %>%
    dplyr::left_join(dplyr::select(filter(object_iactions.df, var == "iaction_labu"), ptm_type, ptmn_id, nselptms, condition, median_log2)) %>%
    dplyr::mutate(ptm_id = if_else(is.na(ptm_id), -hom_ptm_id, ptm_id)) %>% # assign some ptm_id for unobserved ptm
    dplyr::filter(treatment == organism_short) %>%
    # show only the highest ptmn_id per ptm_id
    dplyr::group_by(protein_ac, ptm_pos, condition) %>% dplyr::filter(row_number(-coalesce(100 + median_log2, 0)) == 1L) %>% dplyr::ungroup() %>%
    dplyr::mutate(timepoint_label = factor(str_c(timepoint, "h p.i."), levels=str_c(sort(unique(conditions.df$timepoint)), "h p.i.")),
                  shown_ptm_pos = if_else(is.na(agn_ptm_pos), ptm_pos, agn_ptm_pos),
                  shown_median_log2 = pmax(if_else(!is.na(median_log2) & ptm_status %in% c("observed", "low conf."),
                                                   global_labu_shift/log(2) + median_log2, 3.5), 3.5),
                  shown_median_log2 = if_else(organism_short == "SARS_CoV2", shown_median_log2, -shown_median_log2),
                  ptm_label = str_c(ptm_AA_seq, ptm_pos),# if_else(coalesce(nselptms, 1) > 1, str_c(" (M", nselptms,")"), "")),
                  #ptm_status = if_else((ptm_status == "observed" &
                  #                           !((coalesce(ptm_cond_qvalue_min, 1) <= ptm_qvalue_max) &
                  #                             (coalesce(ptm_cond_locprob_max, 0) >= ptm_locprob_min) &
                  #                             (coalesce(n_cond_idented_and_localized, 0) > 0))), "low conf.", ptm_status),
                  ptm_status = as.character(ptm_status),
                  ptm_typeXstatus = case_when(ptm_status == "observed" ~ ptm_type,
                                              ptm_status == "N/A" ~ ptm_status,
                                              TRUE ~ str_c(ptm_type, " (", ptm_status, ")")))

View(filter(viral_ptmn_intensity.df, genename_short == "S" & ptm_pos == 478))

protein_regions.df <- read_tsv(file=file.path(data_path, "cov_cov2-viralprotein_domains_26102020.txt")) %>%
    rename(display_label = display_labbel) %>%
    dplyr::mutate(display_label = display_label == "+",
                  display_pos = display_pos == "+")

protein_agn_regions.df <- dplyr::inner_join(protein_regions.df, dplyr::select(viral_agn.df, protein_ac, start_pos=ptm_pos, agn_start_pos=agn_ptm_pos)) %>%
    dplyr::inner_join(dplyr::select(viral_agn.df, protein_ac, end_pos=ptm_pos, agn_end_pos=agn_ptm_pos)) %>%
    dplyr::left_join(dplyr::select(proteins.df, protein_ac, is_viral, is_contaminant, genename, organism_short)) %>%
    dplyr::filter(organism_short == "SARS_CoV2") %>% # domain are shared between the viruses, so leave just one
    dplyr::arrange(protein_ac, domain, start_pos, end_pos, agn_start_pos, rev(agn_end_pos)) %>%
    dplyr::group_by(protein_ac, domain, start_pos, end_pos) %>%
    dplyr::filter(row_number() == 1L) %>% dplyr::ungroup() %>%
    dplyr::mutate(agn_len = agn_end_pos - agn_start_pos) %>%
    dplyr::arrange(protein_ac, agn_len)

ptm_palette <- c("Phospho" = "#5e268f", "GlyGly" = "black", "N/A" = "gray") # to make GlyGly different from SARS
ptm_fill_palette <- c("Phospho" = "#5e268f", "GlyGly" = "#bf1c2c", "N/A" = "gray")
ptm_typeXstatus_shape_palette <- c("Phospho"=22L, "Phospho (low conf.)"=22L, "Phospho (potential)" = 22L,
                                   "GlyGly"=21L, "GlyGly (low conf.)"=21L, "GlyGly (potential)" = 21L, "N/A" = 4L)
ptm_typeXstatus_fill_palette <- c("Phospho"="#5e268f", "Phospho (low conf.)"="#c29ae5", "Phospho (potential)" = "white",
                                  "GlyGly"="#bf1c2c", "GlyGly (low conf.)"="#ee9099", "GlyGly (potential)" = "white", "N/A" = NA)
ptm_status_width_palette <- c("observed"=0.5, "low conf."=0.5, "potential" = 0.5, "N/A" = 0.5)
domain_type_palette <- c("localization" = "lemonchiffon", "function" = "palegreen3")
domain_type_color_palette <- c("localisation" = "lemonchiffon4", "function" = "palegreen4")
treatment_palette <- c(mock="gray", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02")

viral_ptmn_intensity.df %>% dplyr::filter(timepoint == 24) %>%
group_by(genename_pair, agn_len) %>%
    group_walk(~{
        sel_iactions.df <- mutate(.x, y_start=if_else(shown_median_log2 > 0, 2, -2))
        viral_gene <- .y$genename_pair[[1]]
        genenames_df <- semi_join(viral_ptm2ptm.df, .y) %>% dplyr::select(genename, hom_genename) %>% dplyr::distinct()
        genenames <- c(genenames_df$genename, genenames_df$hom_genename)
        genenames <- sort(unique(genenames[!is.na(genenames)]), decreasing = TRUE)
        sel_regions.df <- dplyr::semi_join(protein_agn_regions.df, dplyr::select(.x, protein_ac)) %>% dplyr::filter(display_pos) %>%
            dplyr::mutate(x_center = 0.5*(agn_start_pos+agn_end_pos),
                          y_min = -2,#if_else(organism_short == "SARS_CoV2", 0, -3),
                          y_max = 2,#if_else(organism_short == "SARS_CoV2", 3, 0),
                          y_label = 1.75)#if_else(organism_short == "SARS_CoV2", 1.5, -1.5))
        message("Plotting PTMs of ", viral_gene)
        p <- ggplot(sel_iactions.df, aes(x = agn_ptm_pos, y = shown_median_log2)) + #, color = organism_short)) +
            #annotate(geom="segment", color="darkgray", size = 1.5, x=0, y=0, xend=.y$agn_len, yend=0) +
            #geom_rect(data = sel_regions.df, aes(xmin=agn_start_pos, xmax=agn_end_pos, x=x_center,
            #                                     y=y_label, ymin=y_min, ymax=y_max,
            #                                     fill=domain_type, color=domain_type), alpha = 0.4) +
            annotate(geom="rect", color="darkgray", fill="gray", size = 0.25, xmin=0, xmax=.y$agn_len, ymin=-2, ymax=2) +
            geom_segment(data = sel_regions.df,
                         aes(x=agn_start_pos, xend=agn_start_pos, y=y_min, yend=y_max, color=domain_type), size=0.25) +
            geom_segment(data = sel_regions.df,
                         aes(x=agn_end_pos, xend=agn_end_pos, y=y_min, yend=y_max, color=domain_type), size=0.25) +
            geom_text_repel(data=sel_regions.df,
                            aes(x=agn_start_pos, y=y_min, label=domain_short, color=domain_type),
                            nudge_y=-1, min.segment.length=0.2, segment.alpha=0.5, segment.size=0.25, vjust=0, size=4) +
            geom_text_repel(data=sel_regions.df,
                            aes(x=agn_end_pos, y=y_max, label=domain_short, color=domain_type),
                            nudge_y=1, min.segment.length=0.2, segment.alpha=0.5, segment.size=0.25, vjust=0, size=4) +
            scale_fill_manual("Domain", values = domain_type_palette) +
            scale_color_manual("Domain", values = domain_type_color_palette) +
            new_scale_color() +
            geom_segment(aes(x = agn_ptm_pos, xend = agn_ptm_pos, y = y_start, yend = shown_median_log2,
                             color=treatment)) + #, size=ptm_status)) +
            geom_point(data=sel_iactions.df, aes(color=treatment, y=y_start), shape=15, size=1) +
            scale_color_manual("Virus", values = treatment_palette) +
            new_scale_color() + new_scale_fill() +
            geom_point(aes(shape=ptm_typeXstatus, fill=ptm_typeXstatus), color="gray", size=2.5) +
            #scale_color_manual("Virus", values = ptm_typeXstatus_color_palette) + new_scale_color() +
            geom_text_repel(data=filter(sel_iactions.df, treatment == "SARS_CoV2" & ptm_status != "N/A"),
                            aes(label = ptm_label, color=ptm_type),
                            min.segment.length=0.2, segment.alpha=0.5, segment.size=0.25, nudge_y = 1, size=4) +
            geom_text_repel(data=filter(sel_iactions.df, treatment == "SARS_CoV" & ptm_status != "N/A"),
                            aes(label = ptm_label, color=ptm_type), 
                            min.segment.length=0.2, segment.alpha=0.5, segment.size=0.25, nudge_y = -1, size=4) +
            scale_color_manual("PTM", values = ptm_palette) +
            scale_fill_manual("PTM", values = ptm_typeXstatus_fill_palette) +
            scale_shape_manual("PTM", values = ptm_typeXstatus_shape_palette) +
            #scale_size_manual("PTM", values = ptm_status_width_palette) +
            scale_y_continuous(expression(log[2](Intensity), parse=TRUE)) +
            scale_x_continuous("Aligned Position", limits=c(0L, .y$agn_len[[1]])) +
            facet_grid(timepoint_label ~ .) +
            theme_bw_ast(base_family = "") +
            theme(panel.border = element_blank(), panel.grid = element_blank(), axis.line = element_blank()) +
            ggtitle(str_c("PTMs on ", str_c(genenames, collapse="/"), " proteins"))
        plot_path <- file.path(analysis_path, "plots", str_c(datasets_info$ubi$msfolder, '_', datasets_info$ubi$fit_version),
                               str_c("viral_ptm_", plot_version))
        if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
        ggsave(p, file = file.path(plot_path, str_c(project_id, "_", datasets_info$ubi$msfolder, '_', datasets_info$ubi$fit_version, "_viral_",
                                                    viral_gene, "_", plot_version, #"_qvalue_1E-3", 
                                                    ".pdf")),
               width=4+nrow(sel_iactions.df)/15, height=4, device = cairo_pdf, family="Arial")
    })

