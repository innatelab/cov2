# assembly of CoV AP-MS data
#
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200329"
fit_version <- "20200329"
mq_folder <- 'mq_apms_20200329'
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
load(file.path(scratch_path, str_c(project_id, '_msglm_data_', mq_folder, '_', data_version, '.RData')))
load(file.path(scratch_path, str_c(project_id, '_msdata_full_', mq_folder, '_', data_version, '.RData')))

message('Loading MSGLM model fit results...')
strip_samples <- TRUE

modelobj <- "protgroup"
quantobj <- "protgroup"# "pepmodstate"
modelobjs_df <- msdata[[str_c(modelobj, "s")]]
modelobj_idcol <- str_c(modelobj, "_id")
# FIXME should be done in prepare_data
modelobjs_df$object_id <- modelobjs_df[[modelobj_idcol]]
modelobjs_df$object_label <- modelobjs_df[[str_c(modelobj, "_label")]]

chunk_suffix <- case_when(modelobj == "protgroup" ~ "_pg",
                          modelobj == "protregroup" ~ "_prg",
                          TRUE ~ NA_character_)

fit_path <- file.path(scratch_path, paste0(project_id, '_msglm_', mq_folder, chunk_suffix))
fit_files <- list.files(fit_path, paste0(project_id, '_msglm', chunk_suffix, '_', fit_version, '_\\d+\\.RData'))
message('Found ', length(fit_files), ' model file(s)')
fit_files.df <- tibble(filename = as.character(fit_files)) %>%
    mutate(chunk_id = as.integer(str_split_fixed(str_remove(filename, ".RData$"), fixed('_'), 5)[, 5])) %>%
  dplyr::arrange(chunk_id) %>%
  dplyr::mutate(object_id = modelobjs_df[[modelobj_idcol]][chunk_id])
id_range_breaks <- which(c(fit_files.df$chunk_id[-1] - 1L, nrow(modelobjs_df)) != 
                         c(fit_files.df$chunk_id[-length(fit_files.df$chunk_id)], nrow(modelobjs_df)))
fit_files.df[sort(c(id_range_breaks, id_range_breaks+1L)), ]
#fit_files.df <- dplyr::filter(fit_files.df, protgroup_id <= 520)

#write_lines(setdiff(1:nrow(modelobjs_df), fit_files.df$chunk_id),
#            path=file.path(scratch_path, str_c(project_id, "_", fit_version, "_pending_chunk_ids")))

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
iactions.df <- tidyr::expand(msdata_full$protgroup_intensities, msrun, protgroup_id) %>%
  dplyr::left_join(msdata_full$protgroup_intensities) %>%
  dplyr::left_join(msdata_full$protgroup_idents) %>%
  dplyr::inner_join(msdata$msruns) %>%
  dplyr::group_by(condition, protgroup_id) %>%
  dplyr::summarize(n_quanted = n_distinct(msrun[!is.na(intensity)]),
                   n_idented = n_distinct(msrun[ident_type == "By MS/MS"])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(object_id = protgroup_id)
} else if (modelobj == "protregroup") {
iactions.df <- msdata_full$pepmodstate_intensities %>%
  dplyr::mutate(is_quanted = !is.na(intensity.L) | !is.na(intensity.M) | !is.na(intensity.H)) %>%
  dplyr::inner_join(msdata$mschannels) %>%
  dplyr::left_join(msdata$pepmodstates) %>%
  dplyr::left_join(msdata$protregroup2pepmod) %>%
  dplyr::left_join(msdata$protregroup2protgroup) %>%
  dplyr::left_join(supXcond.df) %>%
  dplyr::group_by(condition, protregroup_id) %>%
  dplyr::summarize(n_quanted = n_distinct(msrun[is_idented]), # count msruns
                   n_idented = n_distinct(msrun[is_quanted])) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(object_id = protregroup_id)
}

pre_object_effects.df <- dplyr::inner_join(iactions.df, conditionXeffect.df) %>%
  dplyr::group_by(object_id, effect) %>%
  dplyr::summarise(n_quanted_min = min(n_quanted),
                   n_quanted_max = max(n_quanted),
                   n_idented_min = min(n_idented),
                   n_idented_max = max(n_idented)) %>%
  dplyr::ungroup()

pre_object_contrasts.df <- dplyr::inner_join(iactions.df, conditionXmetacondition.df) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::group_by(object_id, contrast) %>%
  dplyr::summarise(n_quanted_min = replace_na(min(n_quanted[condition_role == "signal"]), 0L),
                   n_quanted_max = replace_na(max(n_quanted[condition_role == "signal"]), 0L),
                   n_idented_min = replace_na(min(n_idented[condition_role == "signal"]), 0L),
                   n_idented_max = replace_na(max(n_idented[condition_role == "signal"]), 0L)) %>%
  dplyr::ungroup()

object_effects.df <- pre_object_effects.df %>% dplyr::inner_join(fit_stats$object_effects) %>%
  dplyr::filter(var %in% c('obj_effect', 'obj_effect_replCI')) %>%
  dplyr::mutate(std_type = if_else(str_detect(var, "_replCI$"), "replicate", "median")) %>%
  dplyr::mutate(trunc_mean_log2 = pmax(-5, pmin(5, mean_log2 - prior_mean_log2)) + prior_mean_log2,
                trunc_median_log2 = pmax(-5, pmin(5, median_log2 - prior_mean_log2)) + prior_mean_log2) %>%
  dplyr::group_by(std_type, var, effect) %>%
  dplyr::mutate(p_value_adj = p.adjust(p_value, method = "BY")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(is_signif = p_value <= 1E-2,
                is_hit = (n_idented_max>0) & (n_quanted_max>0) & is_signif & abs(median_log2 - prior_mean_log2) >= 1.0,
                change = if_else(is_signif, if_else(mean_log2 < prior_mean_log2, "-", "+"), "."))

object_contrasts.df <- dplyr::inner_join(pre_object_contrasts.df, fit_contrasts$iactions) %>%
  dplyr::ungroup() %>%
  dplyr::filter(var %in% c('iaction_labu', 'iaction_labu_replCI')) %>%
  dplyr::mutate(std_type = if_else(str_detect(var, "_replCI$"), "replicate", "median")) %>%
  dplyr::inner_join(dplyr::select(contrastXmetacondition.df, contrast, metacondition, contrast_type, condition_role)) %>%
  dplyr::mutate(trunc_mean = pmax(-10, pmin(10, mean)),
                trunc_median_log2 = pmax(-10, pmin(10, median_log2)),
                p_value = pmin(prob_nonpos, prob_nonneg)) %>%
  dplyr::group_by(std_type, var, contrast) %>%
  dplyr::mutate(p_value_adj = pmin(p.adjust(c(prob_nonpos, prob_nonneg), method = "BY")[1:n()],
                                   p.adjust(c(prob_nonneg, prob_nonpos), method = "BY")[1:n()])) %>%
  dplyr::ungroup()

weak_bait_ids <- c("SARS_CoV_ORF6",
                   "SARS_CoV2_ORF6",
                   "SARS_CoV2_ORF7a",
                   "SARS_CoV_ORF8",
                   "SARS_CoV_ORF9b",
                   "SARS_CoV2_N",
                   "SARS_CoV2_E",
                   "SARS_CoV2_NSP3_macroD",
                   "SARS_CoV2_NSP4",
                   "SARS_CoV2_NSP7",
                   "SARS_CoV2_NSP15",
                   "SARS_CoV2_NSP16")

object_contrasts_thresholds.df <- select(object_contrasts.df, contrast, contrast_type, std_type) %>%
  distinct() %>%
  mutate(bait_full_id = str_remove(contrast, "_vs_.+"),
         p_value_threshold = case_when(bait_full_id %in% weak_bait_ids ~ 0.01,
                                       TRUE ~ 0.001),
         median_log2_threshold = case_when(bait_full_id %in% weak_bait_ids ~ 1,
                                           TRUE ~ 2))


object_contrasts.df <- object_contrasts.df %>%
  select(-any_of(c("p_value_threshold", "median_log2_threshold"))) %>%
  left_join(object_contrasts_thresholds.df) %>%
  dplyr::mutate(is_signif = p_value <= p_value_threshold & abs(median_log2) >= median_log2_threshold,
                is_hit = is_signif & !is_contaminant & !is_reverse &
                  ((contrast_type == "comparison") | (median_log2 >= 0.0)) &
                  (n_quanted_max>0) & (n_idented_max>0),
                change = if_else(is_hit, if_else(median_log2 < 0, "-", "+"), "."))

object_effects_wide.df <- pivot_wider(object_effects.df,
                                      id_cols = c("std_type", "object_id", "object_label", "majority_protein_acs", "gene_names"),
                                      names_from = "effect", names_sep = ".",
                                      values_from = c("median_log2", "mean_log2", "sd_log2", "p_value", "is_hit", "change"))

object_contrasts_wide.df <- pivot_wider(object_contrasts.df,
                                        id_cols = c("std_type", "object_id", "object_label", "majority_protein_acs", "gene_names"),
                                        names_from = "contrast", names_sep = ".",
                                        values_from = c("median_log2", "mean_log2", "sd_log2", "p_value", "is_hit", "change"))

rfit_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_fit_', mq_folder, '_', fit_version, '.RData'))
results_info <- list(project_id = project_id, mq_folder=mq_folder,
                     data_version = data_version, fit_version = fit_version,
                     modelobj = modelobj, quantobj = quantobj)
message('Saving full analysis results to ', rfit_filepath, '...')
save(results_info, fit_stats, fit_contrasts,
     object_effects.df, object_contrasts.df,
     object_effects_wide.df, object_contrasts_wide.df,
     object_contrasts_thresholds.df,
     file = rfit_filepath)
message('Done.')
