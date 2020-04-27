# assembly of CoV AP-MS data
#
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200423"
fit_version <- "20200425"
ms_folder <- 'cov2timecourse_phospho_dia_20200423'
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
load(file.path(scratch_path, str_c(project_id, '_msglm_data_', ms_folder, '_', data_version, '.RData')))
load(file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData')))

modelobj <- "ptmgroup"
quantobj <- "ptmgroup"
source(file.path(project_scripts_path, 'setup_modelobj.R'))

message('Loading MSGLM model fit results...')
strip_samples <- TRUE

fit_path <- file.path(scratch_path, str_c(project_id, '_', ms_folder, '_msglm', modelobj_suffix))
fit_files <- list.files(fit_path, str_c(project_id, '_', ms_folder, '_msglm', modelobj_suffix, 
                                        '_', fit_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- tibble(filename = as.character(fit_files)) %>%
    mutate(chunk_id = as.integer(str_split_fixed(str_remove(filename, ".RData$"), fixed('_'), 9)[, 9])) %>%
  dplyr::arrange(chunk_id) %>%
  dplyr::mutate(object_id = modelobjs_df[[modelobj_idcol]][chunk_id])
id_range_breaks <- which(c(fit_files.df$chunk_id[-1] - 1L, nrow(modelobjs_df)) != 
                         c(fit_files.df$chunk_id[-length(fit_files.df$chunk_id)], nrow(modelobjs_df)))
fit_files.df[sort(c(id_range_breaks, id_range_breaks+1L)), ]
#fit_files.df <- dplyr::filter(fit_files.df, protgroup_id <= 520)

write_lines(setdiff(1:nrow(modelobjs_df), 
                    fit_files.df$chunk_id),
            path=file.path(scratch_path, str_c(project_id, "_", ms_folder, '_', fit_version, "_pending_chunk_ids")))

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

if (modelobj == "ptmgroup") {
  iactions.df <- msdata_full$ptmgroup_intensities_all %>%
    dplyr::inner_join(msdata$msruns) %>%
    dplyr::mutate(is_quanted = !is.na(intensity),
                  is_idented = replace_na(ident_type == "By MS/MS", FALSE)) %>%
    dplyr::group_by(condition, ptmgroup_id) %>%
    dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_quanted]),
                     nmsruns_idented = n_distinct(msrun[is_idented])) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(object_id = ptmgroup_id)
  
  modelobjs_df <- dplyr::mutate(modelobjs_df,
                                is_msvalid_object = TRUE)#(nprotgroups_sharing_proteins == 1 || nproteins_have_razor > 0))
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
  mutate(p_value_threshold = 0.005,
         median_log2_threshold = 0.5,
         median_log2_max = 15)

object_effects.df <- object_effects.df %>%
  select(-any_of(c("p_value_threshold", "median_log2_threshold", "median_log2_max"))) %>%
  left_join(object_effects_thresholds.df) %>%
  dplyr::mutate(mean_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, mean_log2)),
                median_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, median_log2)),
                is_signif = p_value <= p_value_threshold & abs(median_log2) >= median_log2_threshold,
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

object_contrasts.df <- dplyr::inner_join(pre_object_contrasts.df, fit_contrasts$iactions) %>%
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
  mutate(p_value_threshold = 0.005,
         median_log2_threshold = 0.5,
         median_log2_max = 15)

object_contrasts.df <- object_contrasts.df %>%
  select(-any_of(c("p_value_threshold", "median_log2_threshold", "median_log2_max"))) %>%
  left_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(mean_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, mean_log2)),
                median_log2_trunc = pmax(-median_log2_max, pmin(median_log2_max, median_log2)),
                is_signif = p_value <= p_value_threshold & abs(median_log2) >= median_log2_threshold,
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse &
                  ((contrast_type == "comparison") | (median_log2 >= 0.0)),
                is_hit = is_hit_nomschecks & ((nmsruns_quanted_lhs_max>=2) & (nmsruns_quanted_rhs_max>=2)) &
                  is_msvalid_object,
                change = if_else(is_signif, if_else(median_log2 < 0, "-", "+"), "."))

object_contrast_stats.df <- dplyr::group_by(object_contrasts.df, contrast, contrast_type, std_type) %>%
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

object_effects_wide.df <- pivot_wider(object_effects.df,
                                      id_cols = c("std_type", "object_id", "object_label", "majority_protein_acs", "gene_names"),
                                      names_from = "effect", names_sep = ".",
                                      values_from = c("median_log2", "mean_log2", "sd_log2", "p_value", "is_hit", "change"))

object_contrasts_wide.df <- pivot_wider(object_contrasts.df,
                                        id_cols = c("std_type", "object_id", "object_label", "majority_protein_acs", "gene_names"),
                                        names_from = "contrast", names_sep = ".",
                                        values_from = c("median_log2", "mean_log2", "sd_log2", "p_value", "is_hit", "change"))

rfit_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_fit_', ms_folder, '_', fit_version, '.RData'))
results_info <- list(project_id = project_id, ms_folder=ms_folder,
                     data_version = data_version, fit_version = fit_version,
                     modelobj = modelobj, quantobj = quantobj)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts,
     object_effects.df, object_contrasts.df,
     object_effects_wide.df, object_contrasts_wide.df,
     object_contrasts_thresholds.df,
     file = rfit_filepath)
message('Done.')
