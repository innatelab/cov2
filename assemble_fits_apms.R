# assembly of CoV AP-MS data
#
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200515"
fit_version <- "20200515"
ms_folder <- 'mq_apms_20200510'
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

source(file.path(project_scripts_path, 'setup_modelobj.R'))

message('Loading MSGLM model fit results...')
strip_samples <- TRUE

fit_path <- file.path(scratch_path, str_c(project_id, '_', ms_folder, '_msglm', modelobj_suffix))
fit_files <- list.files(fit_path, str_c(project_id, '_', ms_folder, '_msglm', modelobj_suffix, 
                                        '_', fit_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- tibble(filename = as.character(fit_files)) %>%
    mutate(chunk_id = as.integer(str_split_fixed(str_remove(filename, ".RData$"), fixed('_'), 8)[, 8])) %>%
  dplyr::arrange(chunk_id) %>%
  dplyr::mutate(object_id = modelobjs_df[[modelobj_idcol]][chunk_id])
id_range_breaks <- which(c(fit_files.df$chunk_id[-1] - 1L, nrow(modelobjs_df)) != 
                         c(fit_files.df$chunk_id[-length(fit_files.df$chunk_id)], nrow(modelobjs_df)))
fit_files.df[sort(c(id_range_breaks, id_range_breaks+1L)), ]
#fit_files.df <- dplyr::filter(fit_files.df, protgroup_id <= 520)

#write_lines(setdiff(1:nrow(modelobjs_df), fit_files.df$chunk_id),
#            path=file.path(scratch_path, str_c(project_id, "_", ms_folder, '_', fit_version, "_pending_chunk_ids")))

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
iactions.df <- msdata_full$protgroup_intensities_all %>%
  dplyr::left_join(msdata$protgroup_idents) %>%
  dplyr::inner_join(msdata$msruns) %>%
  dplyr::mutate(is_quanted = !is.na(intensity),
                is_idented = replace_na(ident_type == "By MS/MS", FALSE)) %>%
  dplyr::group_by(condition, protgroup_id) %>%
  dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_quanted]),
                   nmsruns_idented = n_distinct(msrun[is_idented])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(object_id = protgroup_id)

modelobjs_df <- dplyr::mutate(modelobjs_df,
                              is_msvalid_object = (nprotgroups_sharing_proteins == 1 || nproteins_have_razor > 0))
} else if (modelobj == "protregroup") {
iactions.df <- tidyr::expand(msdata_full$pepmodstate_intensities, nesting(pepmod_id, pepmodstate_id), msrun) %>%
  dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::mutate(is_quanted = !is.na(intensity)) %>%
  dplyr::inner_join(msdata$msruns) %>%
  dplyr::left_join(msdata$pepmodstates) %>%
  dplyr::left_join(msdata$protregroup2pepmod) %>%
  dplyr::filter(is_specific) %>%
  dplyr::group_by(condition, protregroup_id) %>%
  dplyr::summarize(nmsruns_quanted = n_distinct(msrun[is_idented]), # count msruns
                   nmsruns_idented = n_distinct(msrun[is_quanted])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(object_id = protregroup_id)

modelobjs_df <- dplyr::mutate(modelobjs_df,
                              is_msvalid_object = TRUE)
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
  dplyr::group_by(object_id, contrast_type, contrast, metacondition, condition_role, is_lhs) %>%
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
  dplyr::group_by(object_id, contrast_type, contrast) %>%
  dplyr::summarise(nmsruns_quanted_lhs_min = nmsruns_quanted_min[is_lhs],
                   nmsruns_quanted_lhs_max = nmsruns_quanted_max[is_lhs],
                   nmsruns_idented_lhs_min = nmsruns_idented_min[is_lhs],
                   nmsruns_idented_lhs_max = nmsruns_idented_max[is_lhs],
                   nmsruns_quanted_rhs_min = nmsruns_quanted_min[!is_lhs],
                   nmsruns_quanted_rhs_max = nmsruns_quanted_max[!is_lhs],
                   nmsruns_idented_rhs_min = nmsruns_idented_min[!is_lhs],
                   nmsruns_idented_rhs_max = nmsruns_idented_max[!is_lhs]) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(filter(contrastXcondition.df, weight > 0),
                                 bait_full_id, contrast)) %>%
  dplyr::left_join(dplyr::select(filter(contrastXcondition.df, contrast_type == "comparison" & weight < 0),
                                 bait_full_id_rhs=bait_full_id, contrast))

object_effects.df <- pre_object_effects.df %>% dplyr::inner_join(fit_stats$object_effects) %>%
  dplyr::left_join(select(modelobjs_df, object_id, is_msvalid_object)) %>%
  dplyr::filter(var %in% c('obj_effect', 'obj_effect_replCI')) %>%
  dplyr::mutate(std_type = if_else(str_detect(var, "_replCI$"), "replicate", "median")) %>%
  dplyr::mutate(trunc_mean_log2 = pmax(-5, pmin(5, mean_log2 - prior_mean_log2)) + prior_mean_log2,
                trunc_median_log2 = pmax(-5, pmin(5, median_log2 - prior_mean_log2)) + prior_mean_log2) %>%
  dplyr::group_by(std_type, var, effect) %>%
  dplyr::mutate(p_value_adj = p.adjust(p_value, method = "BY")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(is_signif = (p_value <= 1E-2) & (abs(median_log2 - prior_mean_log2) >= 1.0),
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse, 
                is_hit = is_hit_nomschecks &
                         (nmsruns_idented_max>2) & (nmsruns_quanted_max>2) & (is_msvalid_object),
                change = if_else(is_signif, if_else(mean_log2 < prior_mean_log2, "-", "+"), "."))

object_effects_stats.df <- dplyr::group_by(object_effects.df, effect, std_type) %>%
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

View(filter(fit_stats$object_batch_effects, var=="obj_batch_effect") %>%
  group_by(object_id) %>%
  mutate(has_batcheffects = any(abs(median_log2) >= 1 & p_value <= 0.05)) %>%
  ungroup())

object_contrasts.df <- dplyr::inner_join(pre_object_contrasts.df, fit_contrasts$iactions) %>%
  dplyr::left_join(select(modelobjs_df, object_id, is_msvalid_object)) %>%
  dplyr::left_join(filter(fit_stats$object_batch_effects, var=="obj_batch_effect") %>%
                   group_by(object_id) %>%
                   summarise(has_batcheffects = any(abs(median_log2) >= 1 & p_value <= 0.05)) %>%
                   ungroup()) %>%
  dplyr::filter(var %in% c('iaction_labu', 'iaction_labu_replCI')) %>%
  dplyr::mutate(std_type = if_else(str_detect(var, "_replCI$"), "replicate", "median")) %>%
  dplyr::mutate(p_value = pmin(prob_nonpos, prob_nonneg) * if_else(contrast_type == "comparison", 2, 1)) %>%
  dplyr::group_by(std_type, var, contrast) %>%
  dplyr::mutate(p_value_adj = pmin(p.adjust(c(prob_nonpos, prob_nonneg), method = "BY")[1:n()],
                                   p.adjust(c(prob_nonneg, prob_nonpos), method = "BY")[1:n()])) %>%
  dplyr::ungroup()

weak_bait_ids <- c("ORF3b", "ORF4", "ORF6",
                   "ORF7a", "ORF8", "ORF8a", "ORF9b",
                   "E", "N", "S",
                   "NSP1", "NSP2", "NSP3_macroD",
                   "NSP4", "NSP7", "NSP8", "NSP9", "NSP10",
                   "NSP12", "NSP14", "NSP15", "NSP16"
                   )

object_contrasts_thresholds.df <- select(object_contrasts.df, contrast, contrast_type, std_type, bait_full_id) %>%
  distinct() %>%
  dplyr::left_join(select(conditions.df, bait_full_id, bait_id)) %>%
  mutate(p_value_threshold = case_when(contrast_type == "filter" & bait_id %in% weak_bait_ids ~ 0.001,
                                       TRUE ~ 0.001),
         median_log2_threshold = case_when(bait_id %in% weak_bait_ids ~ 2,
                                           TRUE ~ 4))

object_contrasts.df <- object_contrasts.df %>%
  select(-any_of(c("p_value_threshold", "median_log2_threshold", "median_log2_max"))) %>%
  left_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = p_value <= p_value_threshold & abs(median_log2) >= median_log2_threshold,
                is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse & #!has_batcheffects &
                  ((contrast_type == "comparison") | (median_log2 >= 0.0)),
                is_hit = is_hit_nomschecks &
                  if_else(median_log2 > 0,
                          (nmsruns_idented_lhs_max>2) & (nmsruns_quanted_lhs_max>2),
                          (nmsruns_idented_rhs_max>2) & (nmsruns_quanted_rhs_max>2)) &
                  is_msvalid_object,
                change = if_else(is_signif, if_else(median_log2 < 0, "-", "+"), "."))
# in bait_vs_bait comparison, consider only the hits of bait_vs_controls comparison
object_contrasts_vs_controls_hits.df <- filter(object_contrasts.df, str_detect(contrast, "_vs_others") & is_hit_nomschecks)
object_contrasts.df <- dplyr::select(object_contrasts.df, -any_of(c("is_hit_nomschecks_lhs", "is_hit_nomschecks_rhs"))) %>%
  dplyr::left_join(dplyr::select(object_contrasts_vs_controls_hits.df, std_type, bait_full_id, object_id, is_hit_nomschecks_lhs=is_hit_nomschecks)) %>%
  dplyr::left_join(dplyr::select(object_contrasts_vs_controls_hits.df, std_type, bait_full_id_rhs = bait_full_id, object_id, is_hit_nomschecks_rhs=is_hit_nomschecks)) %>%
  dplyr::mutate(is_hit_nomschecks = is_hit_nomschecks &
                    ((contrast_type != "comparison") |
                     (median_log2 > 0 & coalesce(is_hit_nomschecks_lhs, FALSE)) | (median_log2 < 0 & coalesce(is_hit_nomschecks_rhs, FALSE))),
                is_hit = is_hit & is_hit_nomschecks)
object_batch_contrasts.df <- filter(object_contrasts.df, str_detect(contrast, "_vs_B\\d+_(others|controls)")) %>%
  transmute(std_type, object_id, contrast_batch = contrast,
            contrast = str_replace(contrast, "_vs_B\\d+_", "_vs_"),
            median_log2_batch = median_log2, p_value_batch = p_value,
            is_hit_nomschecks_batch = is_hit_nomschecks,
            is_signif_batch = is_signif, is_hit_batch = is_hit)
object_contrasts.df <- left_join(object_contrasts.df,
                                 select(filter(object_batch_contrasts.df, std_type=="median"), -std_type)) %>%
  mutate(is_hit_nomschecks = is_hit_nomschecks & coalesce(is_hit_nomschecks_batch, TRUE),
         is_hit = is_hit & is_hit_nomschecks)

require(readxl)

msruns_seq.df <- read_tsv(file.path(data_path, data_info$ms_folder, "SARS_COV2_MS_samples_sequence_20200519.txt")) %>%
  left_join(dplyr::select(msdata$msruns, raw_file, condition, bait_full_id, bait_homid, msrun))

conditions_seq.df <- group_by(msruns_seq.df, condition, bait_full_id, bait_homid) %>%
  summarise(msrun1st_ix = min(msrun_ix),
            batch = sort(unique(batch))[[1]],
            lc_column = sort(unique(lc_column))[[1]],
            n_batches = n_distinct(batch),
            n_lc_columns = n_distinct(lc_column)) %>%
  ungroup() %>%
  arrange(msrun1st_ix) %>%
  mutate(condition_ix = row_number())

condition_neighbours.df <- inner_join(conditions_seq.df,
                                      dplyr::select(conditions_seq.df,
                                                    batch, lc_column,
                                                    msrun1st_ix_carryover = msrun1st_ix,
                                                    condition_ix_carryover = condition_ix,
                                                    condition_carryover = condition,
                                                    bait_full_id_carryover = bait_full_id,
                                                    bait_homid_carryover = bait_homid)) %>%
  mutate(nmsruns_carryover = msrun1st_ix - msrun1st_ix_carryover) %>%
  filter(condition != condition_carryover &
         str_remove(bait_homid_carryover, "\\?+$") != str_remove(bait_homid, "\\?+$") &
           between(nmsruns_carryover, 0, 8))

contrast_neighbours.df <- condition_neighbours.df %>%
  inner_join(select(filter(contrastXcondition.df, weight > 0),
                    condition, contrast) %>% unique()) %>%
  inner_join(select(filter(contrastXcondition.df, weight > 0),
                    condition_carryover = condition, contrast_carryover = contrast) %>% unique()) %>%
  inner_join(transmute(filter(contrastXmetacondition.df, contrast_type == "filter" & weight < 0),
                       contrast,
                       metacondition_type_rhs = str_remove(metacondition, "_[^_]+$"))) %>%
  inner_join(transmute(filter(contrastXmetacondition.df, contrast_type == "filter" & weight < 0),
                       contrast_carryover=contrast,
                       metacondition_type_rhs = str_remove(metacondition, "_[^_]+$"))) %>%
  select(contrast, contrast_carryover, nmsruns_carryover) %>%
  distinct()

object_contrasts_carryover.df <- inner_join(select(object_contrasts.df, object_id, std_type,
                  contrast_carryover=contrast, is_signif_carryover=is_signif, is_hit_nomschecks_carryover = is_hit_nomschecks,
                  median_log2_carryover=median_log2, p_value_carryover = p_value),
           contrast_neighbours.df) %>%
  group_by(object_id, std_type, contrast) %>%
  filter(row_number(p_value_carryover) == 1L) %>% # select most significant carryover
  arrange(object_id, std_type, contrast) %>%
  ungroup()

object_contrasts.df <- left_join(object_contrasts.df, object_contrasts_carryover.df) %>%
  mutate(is_carryover = is_hit_nomschecks &
         (coalesce(p_value_carryover, 1.0) <= 0.01 & coalesce(median_log2_carryover, -Inf) >= 1.0) &
         median_log2 <= coalesce(median_log2_carryover, -Inf),
         is_hit = is_hit & !is_carryover)

object_contrast_stats.df <- dplyr::group_by(object_contrasts.df, contrast, bait_full_id, contrast_type, std_type) %>%
  dplyr::summarise(p_value_001 = quantile(p_value, 0.001),
                   p_value_01 = quantile(p_value, 0.01),
                   p_value_05 = quantile(p_value, 0.05),
                   median_log2_abs_50 = quantile(abs(median_log2[p_value <= 0.1]), 0.5),
                   median_log2_abs_95 = quantile(abs(median_log2[p_value <= 0.1]), 0.95),
                   median_log2_abs_99 = quantile(abs(median_log2[p_value <= 0.1]), 0.99),
                   n_hits = sum(is_hit_nomschecks),
                   n_carryover = sum(is_carryover)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(msdata$msruns, bait_full_id, batch) %>%
                   dplyr::group_by(bait_full_id) %>% dplyr::summarise(batches = str_c(unique(batch), collapse=" ")) %>% dplyr::ungroup())
View(filter(object_contrast_stats.df, std_type == "median" & str_detect(contrast, "_vs_others")) %>% dplyr::arrange(desc(batches), desc(n_hits)))
View(filter(object_contrast_stats.df, std_type == "replicate" & str_detect(contrast, "_vs_others")) %>% dplyr::arrange(desc(batches), desc(n_hits)))
View(filter(object_contrast_stats.df, std_type == "replicate" & str_detect(contrast, "ACE2_vs_B4_others")) %>% dplyr::arrange(desc(batches), desc(n_hits)))

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
