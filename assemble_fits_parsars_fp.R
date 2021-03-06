# assembly of CoV AP-MS data
#
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200830"
fit_version <- "20200830"
msfolder <- 'snaut_parsars_fp_20200829'
message("Assembling fit results for project ", project_id,
        " (dataset v", data_version, ", fit v", fit_version, ")")

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))
require(msglm)

mop.max_nprocesses <- 32
mop.nprocesses <- 32
source(file.path(pipeline_scripts_path, 'init_cluster.R'))

require(dplyr)
require(stringr)
require(maxquantUtils)

message('Loading data...')
load(file.path(scratch_path, str_c(project_id, '_msdata_full_', msfolder, '_', data_version, '.RData')))
load(file.path(scratch_path, str_c(project_id, '_msglm_data_', msfolder, '_', data_version, '.RData')))

modelobj <- "protregroup"
quantobj <- "pepmodstate"
#modelobj <- "protgroup"
#quantobj <- "protgroup"
source(file.path(project_scripts_path, 'setup_modelobj.R'))

message('Loading MSGLM model fit results...')
strip_samples <- TRUE

fit_path <- file.path(scratch_path, str_c(project_id, '_', msfolder, '_msglm', modelobj_suffix))
fit_files <- list.files(fit_path, str_c(project_id, '_', msfolder, '_msglm', modelobj_suffix, 
                                        '_', fit_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- tibble(filename = as.character(fit_files)) %>%
  tidyr::extract(filename, "chunk_id", ".+_(\\d+).RData$", convert=TRUE, remove=FALSE) %>%
  dplyr::arrange(chunk_id) %>%
  dplyr::mutate(object_id = modelobjs_df[[modelobj_idcol]][chunk_id])
id_range_breaks <- which(c(fit_files.df$chunk_id[-1] - 1L, nrow(modelobjs_df)) != 
                         c(fit_files.df$chunk_id[-length(fit_files.df$chunk_id)], nrow(modelobjs_df)))
fit_files.df[sort(c(id_range_breaks, id_range_breaks+1L)), ]
#fit_files.df <- dplyr::filter(fit_files.df, protgroup_id <= 520)

write_lines(setdiff(1:nrow(modelobjs_df), 
                    fit_files.df$chunk_id),
            file=file.path(scratch_path, str_c(project_id, "_", msfolder, '_', fit_version, "_pending_chunk_ids")))

#load(file.path(fit_path, fit_files[1]))

if (!is.null(mop.cluster)) {
  clusterEvalQ(mop.cluster, library(dplyr))
  clusterEvalQ(mop.cluster, library(msglm))
  clusterExport(mop.cluster, varlist=c('fit_path', 'fit_files.df', 'process_msglm_chunk'))
  fit_reports <- clusterApplyLB(mop.cluster, seq_along(fit_files.df$chunk_id), process_msglm_chunk)
  clusterEvalQ(mop.cluster, gc())
} else {
  fit_reports <- lapply(seq_along(fit_files.df$chunk_id), process_msglm_chunk)
}
names(fit_reports) <- sapply(fit_reports, function(report) paste0(report$results_info$fit_version, '_', report$model_data$objects$object_id[1]))

fit_stats <- lapply(names(fit_reports[[1]]$msglm_results), join_msglm_reports, fit_reports, 'stats')
names(fit_stats) <- names(fit_reports[[1]]$msglm_results)

fit_contrasts <- lapply(names(fit_reports[[1]]$msglm_results), join_msglm_reports, fit_reports, 'contrast_stats')
names(fit_contrasts) <- names(fit_reports[[1]]$msglm_results)

rm(fit_reports)

if (modelobj == "protgroup") {
# FIXME stats should be quantobj-dependent
iactions.df <- tidyr::expand(msdata_full$protgroup_intensities, protgroup_id, msrun) %>%
  dplyr::left_join(msdata_full$protgroup_intensities) %>%
  dplyr::inner_join(msdata$msruns) %>%
  dplyr::mutate(is_quanted = !is.na(intensity),
                is_idented = replace_na(ident_type == "By MS/MS", FALSE)) %>%
  dplyr::group_by(condition, protgroup_id) %>%
  dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_quanted]),
                   nmsruns_idented = n_distinct(msrun[is_idented])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(object_id = protgroup_id)

modelobjs_df <- dplyr::mutate(modelobjs_df,
                              is_msvalid_object = npepmods_unique >= 2L)#(nprotgroups_sharing_proteins == 1 || nproteins_have_razor > 0))
} else if (modelobj == "protregroup") {
iactions.df <- expand(msdata_full$pepmodstate_intensities, msrun, pepmod_id) %>%
  dplyr::inner_join(msdata$msruns) %>%
  dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::mutate(is_quanted = !is.na(intensity),
                is_idented = is_quanted) %>%
  dplyr::left_join(filter(msdata$protregroup2pepmod, is_specific)) %>%
  dplyr::group_by(condition, protregroup_id) %>%
  dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_idented]), # count msruns
                   nmsruns_idented = n_distinct(msrun[is_quanted])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(object_id = protregroup_id)

modelobjs_df <- dplyr::mutate(modelobjs_df,
                              is_msvalid_object = npepmods_unique >= 2L)
}

pre_object_effects.df <- dplyr::inner_join(iactions.df, conditionXeffect.df) %>%
  dplyr::group_by(object_id, effect) %>%
  dplyr::summarise(has_quanted = any(!is.na(nmsruns_quanted)),
                   nmsruns_quanted_min = min(nmsruns_quanted, na.rm=TRUE),
                   nmsruns_quanted_max = max(nmsruns_quanted, na.rm=TRUE),
                   has_idented = any(!is.na(nmsruns_idented)),
                   nmsruns_idented_min = min(nmsruns_idented, na.rm=TRUE),
                   nmsruns_idented_max = max(nmsruns_idented, na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(nmsruns_quanted_min = if_else(has_quanted, nmsruns_quanted_min, 0L),
                nmsruns_quanted_max = if_else(has_quanted, nmsruns_quanted_max, 0L),
                nmsruns_idented_min = if_else(has_idented, nmsruns_idented_min, 0L),
                nmsruns_idented_max = if_else(has_idented, nmsruns_idented_max, 0L),
                has_quanted = NULL, has_idented = NULL)

pre_object_contrasts.df <- dplyr::inner_join(iactions.df, conditionXmetacondition.df) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::mutate(is_lhs = weight > 0) %>%
  dplyr::group_by(object_id, contrast, is_lhs) %>%
  dplyr::summarise(has_quanted = any(!is.na(nmsruns_quanted)),
                   nmsruns_quanted_min = min(nmsruns_quanted, na.rm=TRUE),
                   nmsruns_quanted_max = max(nmsruns_quanted, na.rm=TRUE),
                   has_idented = any(!is.na(nmsruns_idented)),
                   nmsruns_idented_min = min(nmsruns_idented, na.rm=TRUE),
                   nmsruns_idented_max = max(nmsruns_idented, na.rm=TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(nmsruns_quanted_min = if_else(has_quanted, nmsruns_quanted_min, 0L),
                nmsruns_quanted_max = if_else(has_quanted, nmsruns_quanted_max, 0L),
                nmsruns_idented_min = if_else(has_idented, nmsruns_idented_min, 0L),
                nmsruns_idented_max = if_else(has_idented, nmsruns_idented_max, 0L)) %>%
  dplyr::group_by(object_id, contrast) %>%
  dplyr::summarise(nmsruns_quanted_lhs_min = nmsruns_quanted_min[is_lhs],
                   nmsruns_quanted_lhs_max = nmsruns_quanted_max[is_lhs],
                   nmsruns_idented_lhs_min = nmsruns_idented_min[is_lhs],
                   nmsruns_idented_lhs_max = nmsruns_idented_max[is_lhs],
                   nmsruns_quanted_rhs_min = nmsruns_quanted_min[!is_lhs],
                   nmsruns_quanted_rhs_max = nmsruns_quanted_max[!is_lhs],
                   nmsruns_idented_rhs_min = nmsruns_idented_min[!is_lhs],
                   nmsruns_idented_rhs_max = nmsruns_idented_max[!is_lhs]) %>%
  dplyr::ungroup()

object_effects.df <- pre_object_effects.df %>% dplyr::inner_join(fit_stats$object_effects) %>%
  dplyr::left_join(select(modelobjs_df, object_id, is_msvalid_object)) %>%
  dplyr::filter(var %in% c('obj_effect', 'obj_effect_replCI')) %>%
  dplyr::mutate(std_type = if_else(str_detect(var, "_replCI$"), "replicate", "median")) %>%
  dplyr::group_by(std_type, var, effect) %>%
  dplyr::mutate(p_value_adj = p.adjust(p_value, method = "BY")) %>%
  dplyr::ungroup()

object_effects_thresholds.df <- select(object_effects.df, effect, std_type) %>%
  distinct() %>%
  mutate(p_value_threshold = case_when(TRUE ~ 0.001),
         median_log2_threshold = case_when(TRUE ~ 0.25))

object_effects.df <- object_effects.df %>%
  select(-any_of(c("p_value_threshold", "median_log2_threshold"))) %>%
  left_join(object_effects_thresholds.df) %>%
  dplyr::mutate(is_signif = p_value <= p_value_threshold & abs(median_log2) >= median_log2_threshold,
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse,
                is_hit = is_hit_nomschecks & (nmsruns_quanted_min>=1) & (nmsruns_quanted_max>=2) & is_msvalid_object,
                change = if_else(is_signif, if_else(median_log2 < 0, "-", "+"), "."))

object_effect_stats.df <- dplyr::group_by(object_effects.df, effect, std_type) %>%
  dplyr::summarise(p_value_001 = quantile(p_value, 0.001),
                   p_value_01 = quantile(p_value, 0.01),
                   p_value_05 = quantile(p_value, 0.05),
                   median_log2_abs_50 = quantile(abs(median_log2[p_value <= 0.1]), 0.5),
                   median_log2_abs_95 = quantile(abs(median_log2[p_value <= 0.1]), 0.95),
                   median_log2_abs_99 = quantile(abs(median_log2[p_value <= 0.1]), 0.99),
                   n_hits = sum(is_hit_nomschecks, na.rm = TRUE),
                   n_plus = sum(change == "+"),
                   n_minus = sum(change == "-")) %>%
  dplyr::ungroup()
View(filter(object_effect_stats.df, std_type == "median") %>% dplyr::arrange(desc(n_hits)))

object_contrasts.df <- dplyr::inner_join(pre_object_contrasts.df, fit_contrasts$iactions) %>%
  dplyr::rename(contrast_offset_log2 = offset) %>%
  dplyr::left_join(select(modelobjs_df, object_id, is_msvalid_object)) %>%
  dplyr::filter(var %in% c('iaction_labu', 'iaction_labu_replCI')) %>%
  dplyr::mutate(std_type = if_else(str_detect(var, "_replCI$"), "replicate", "median")) %>%
  dplyr::inner_join(dplyr::select(contrastXmetacondition.df, contrast, metacondition, contrast_type, condition_role)) %>%
  dplyr::mutate(p_value = pmin(prob_nonpos, prob_nonneg) * if_else(contrast_type == "comparison", 2, 1)) %>%
  dplyr::group_by(std_type, var, contrast) %>%
  dplyr::mutate(p_value_adj = pmin(p.adjust(c(prob_nonpos, prob_nonneg), method = "BY")[1:n()],
                                   p.adjust(c(prob_nonneg, prob_nonpos), method = "BY")[1:n()])) %>%
  dplyr::ungroup()

object_contrasts_thresholds.df <- select(object_contrasts.df, contrast, contrast_type, std_type) %>%
  distinct() %>%
  dplyr::left_join(dplyr::select(contrasts.df, contrast, contrast_kind, contrast_offset_log2 = offset)) %>%
  mutate(p_value_threshold = case_when(TRUE ~ 0.001),
         p_value_threshold_lesser = case_when(TRUE ~ 0.01),
         median_log2_threshold = case_when(TRUE ~ 0.25),
         median_log2_threshold_lesser = case_when(TRUE ~ 0.125))

object_contrasts.df <- object_contrasts.df %>%
  select(-any_of(c("p_value_threshold", "median_log2_threshold", "median_log2_max"))) %>%
  left_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = p_value <= p_value_threshold & abs(median_log2) >= median_log2_threshold,
                is_signif_lesser = p_value <= p_value_threshold_lesser & abs(median_log2) >= median_log2_threshold_lesser,
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse,
                is_hit = is_hit_nomschecks & (pmax(nmsruns_quanted_lhs_max, nmsruns_quanted_rhs_max)>=2) &
                  is_msvalid_object,
                change = if_else(is_signif, if_else(median_log2 < 0, "-", "+"), "."),
                is_specific_virus_lhs = str_detect(treatment_lhs, "SARS"),
                is_specific_virus_rhs = str_detect(treatment_rhs, "SARS"),
                is_mock_rhs = treatment_rhs == "mock") %>%
  dplyr::group_by(std_type, object_id, timepoint_lhs, timepoint_rhs, contrast_kind) %>%
  # declare _vs_mock a hit if it's a hit or it's significant with less stringent threshold
  # and a strong hit in another virus treatment, while the difference between the viruses is not significant
  dplyr::mutate(is_hit_composed = is_hit | if_else(contrast_kind == "treatment_vs_treatment",
                                                   is_signif_lesser & is_mock_rhs &
                                                   any(is_hit[is_specific_virus_lhs & is_mock_rhs]) & 
                                                   any(!is_signif[is_specific_virus_rhs]),
                                                   FALSE),
                composed_hit_treatment = if_else(contrast_kind == "treatment_vs_treatment" & is_hit_composed & is_specific_virus_lhs & is_mock_rhs,
                                                 treatment_lhs, NA_integer_),
                composed_hit_treatments = str_c(composed_hit_treatment[!is.na(composed_hit_treatment)], collapse="+")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(composed_hit_type = if_else(contrast_kind != "treatment_vs_treatment", NA_character_,
                                            if_else(str_detect(composed_hit_treatments, fixed("+")), "shared",
                                                    if_else(composed_hit_treatments != "", composed_hit_treatments, "none"))),
                change = if_else(is_hit_composed, if_else(median_log2 < 0, "-", "+"), "."))
object_contrast_stats.df <- dplyr::group_by(object_contrasts.df, contrast, contrast_type, contrast_kind, std_type) %>%
  dplyr::summarise(p_value_001 = quantile(p_value, 0.001),
                   p_value_01 = quantile(p_value, 0.01),
                   p_value_05 = quantile(p_value, 0.05),
                   median_log2_abs_50 = quantile(abs(median_log2[p_value <= 0.1]), 0.5),
                   median_log2_abs_95 = quantile(abs(median_log2[p_value <= 0.1]), 0.95),
                   median_log2_abs_99 = quantile(abs(median_log2[p_value <= 0.1]), 0.99),
                   n_hits = sum(is_hit_nomschecks, na.rm = TRUE),
                   n_plus = sum(change == "+"),
                   n_minus = sum(change == "-")) %>%
  dplyr::ungroup()
View(filter(object_contrast_stats.df, std_type == "median" & contrast_kind == "treatment_vs_treatment") %>% dplyr::arrange(desc(n_hits)))

object_effects_wide.df <- pivot_wider(object_effects.df,
                                      id_cols = c("std_type", "object_id", "object_label", "majority_protein_acs", "gene_names"),
                                      names_from = "effect", names_sep = ".",
                                      values_from = c("median_log2", "mean_log2", "sd_log2", "p_value", "is_hit", "change"))

object_contrasts_wide.df <- pivot_wider(object_contrasts.df,
                                        id_cols = c("std_type", "object_id", "object_label", "majority_protein_acs", "gene_names"),
                                        names_from = "contrast", names_sep = ".",
                                        values_from = c("median_log2", "mean_log2", "sd_log2", "p_value", "is_hit", "change"))

rfit_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_fit_', msfolder, '_', fit_version, '.RData'))
results_info <- list(project_id = project_id, msfolder=msfolder,
                     data_version = data_version, fit_version = fit_version,
                     modelobj = modelobj, quantobj = quantobj)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts,
     object_effects.df, object_contrasts.df,
     object_effects_wide.df, object_contrasts_wide.df,
     object_contrasts_thresholds.df,
     file = rfit_filepath)
message('Done.')

report_cols <- c("object_label", "object_id", "gene_names",
                 "majority_protein_acs", "protein_descriptions",
                 "is_contaminant", "is_viral")

pre_object_contrasts_report.df <- filter(object_contrasts.df, contrast_kind=="treatment_vs_treatment" & treatment_lhs != "infected") %>%
  dplyr::select(object_id, std_type, contrast, timepoint=timepoint_lhs,
                median_log2, mean_log2, sd_log2, any_of(c("prob_nonpos", "prob_nonneg", "p_value")),
                is_signif, is_hit_nomschecks, is_hit, change, is_hit_composed, composed_hit_type)

objects4report.df <- dplyr::select(modelobjs_df, any_of(report_cols)) %>%
  dplyr::semi_join(dplyr::select(object_contrasts.df, object_id))

object_contrasts_report.df <- objects4report.df %>%
  dplyr::left_join(pivot_wider(pre_object_contrasts_report.df, c(std_type, object_id),
                               names_from = "contrast", values_from = c("is_hit", "is_hit_composed", "composed_hit_type", "change", "median_log2", "p_value", "sd_log2"),
                               names_sep=".")) %>%
  dplyr::select(any_of(report_cols), std_type,
                ends_with("SARS_CoV2@6h_vs_mock@6h"), ends_with("SARS_CoV@6h_vs_mock@6h"), ends_with("SARS_CoV2@6h_vs_SARS_CoV@6h"),
                ends_with("SARS_CoV2@12h_vs_mock@12h"), ends_with("SARS_CoV@12h_vs_mock@12h"), ends_with("SARS_CoV2@12h_vs_SARS_CoV@12h"),
                ends_with("SARS_CoV2@24h_vs_mock@24h"), ends_with("SARS_CoV@24h_vs_mock@24h"), ends_with("SARS_CoV2@24h_vs_SARS_CoV@24h")) %>%
  # composed_hit_type is the same for all comparison of the same timepoint
  dplyr::select(-matches("composed_hit_type.+SARS_CoV2@\\d+h_vs_(SARS_CoV|mock)@")) %>%#, #-matches("is_hit_composed\\.SARS_CoV2.+_vs_SARS_CoV")) %>%
  dplyr::rename_at(vars(matches("composed_hit_type")), ~str_remove(.x, "SARS_CoV@\\d+h_vs_mock@")) %>%
  dplyr::arrange(gene_names, majority_protein_acs, std_type)

write_tsv(filter(object_contrasts_report.df, std_type == "replicate") %>% dplyr::select(-std_type, -object_id),
          file.path(analysis_path, "reports", paste0(project_id, '_', msfolder, '_contrasts_report_', fit_version, '_replicate_wide.txt')))
write_tsv(filter(object_contrasts_report.df, std_type == "median") %>% dplyr::select(-std_type, -object_id),
          file.path(analysis_path, "reports", paste0(project_id, '_', msfolder, '_contrasts_report_', fit_version, '_median_wide.txt')))
write.xlsx(filter(object_contrasts_report.df, std_type == "replicate") %>% dplyr::select(-std_type, -object_id),
           file.path(analysis_path, "reports", paste0(project_id, '_', msfolder, '_contrasts_report_', fit_version, '_replicate_wide.xlsx')))
write.xlsx(filter(object_contrasts_report.df, std_type == "median") %>% dplyr::select(-std_type, -object_id),
           file.path(analysis_path, "reports", paste0(project_id, '_', msfolder, '_contrasts_report_', fit_version, '_median_wide.xlsx')))

object_contrasts_long_report.df <- objects4report.df %>%
  dplyr::left_join(pre_object_contrasts_report.df) %>%
  dplyr::left_join(dplyr::select(contrasts.df, contrast, treatment_lhs, treatment_rhs)) %>%
  dplyr::arrange(gene_names, majority_protein_acs, timepoint, treatment_rhs, treatment_lhs, contrast, std_type)

write_tsv(filter(object_contrasts_long_report.df) %>% dplyr::select(-object_id),
          file.path(analysis_path, "reports", paste0(project_id, '_', msfolder, '_contrasts_report_', fit_version, '_long.txt')))

