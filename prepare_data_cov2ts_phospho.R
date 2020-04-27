# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200427"
fit_version <- "20200427"
ms_folder <- 'cov2timecourse_phospho_dia_20200423'
message('Dataset version is ', data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, 'matrix_utils.R'))
source(file.path(misc_scripts_path, "ggplot_ext.R"))

require(maxquantUtils)
require(msglm)
require(dplyr)
require(rjson)
require(stringr)
require(readr)
require(pheatmap)
require(tidyr)
#require(broom)

msdata_path <- file.path(data_path, ms_folder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  ms_folder = ms_folder,
                  instr_calib_ptm_filename = "instr_QX7_intensity_ptmgroup_calib_cov2_20200427_borg.json",
                  quant_type = "intensity", quant_col_prefix = "intensity",
                  ptm_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib_ptmgroup <- fromJSON(file = file.path(data_path, data_info$instr_calib_ptm_filename))$instr_calib
instr_calib <- instr_calib_ptmgroup

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/cov_baits_20200415.fasta")),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants.fasta"))
)

bad_msruns <- c()

msdata.wide <- read_tsv(file.path(data_path, ms_folder, "COV2_DIA_phospho_0.75probablity_no normalization.txt"),
                          col_types = c('T: PTM_0_num'='c'))
msdata.wide <- rename_at(msdata.wide, vars(matches("^\\w:\\s")), ~str_remove(., "^\\w:\\s"))

intensity_col_prefix <- "20200417_QX7_MaTa_SA_A549_"
msdata_full <- list(
  ptmgroups = select(msdata.wide, -starts_with(intensity_col_prefix)) %>%
    rename(ptmgroup_id = PTM_collapse_key,
           gene_names = PG.Genes,
           protein_descriptions = PG.ProteinDescriptions,
           protein_names = PG.ProteinNames,
           majority_protein_acs = PG.UniProtIds),
  ptmgroup_intensities = pivot_longer(select(msdata.wide, ptmgroup_id = PTM_collapse_key, starts_with(intensity_col_prefix)),
                                      cols = starts_with(intensity_col_prefix),
                                      names_to = "msrun", names_prefix = intensity_col_prefix, values_to = "log2_intensity") %>%
    filter(is.finite(log2_intensity)) %>%
    mutate(intensity = 2^log2_intensity)
)

msdata_full$ptms <- tidyr::separate_rows(dplyr::select(msdata_full$ptmgroups, ptmgroup_id, protein_names, majority_protein_acs, PEP.PeptidePosition),
                                         majority_protein_acs, sep=';') %>%
  dplyr::rename(protein_ac = majority_protein_acs)

msdata_full$msruns = distinct(select(msdata_full$ptmgroup_intensities, msrun)) %>%
  extract(msrun, c("treatment", "timepoint", "replicate"), "(.+)_(\\d+)h(?:pi)?_(\\d+)$", remove = FALSE) %>%
  mutate(treatment = relevel(factor(treatment), "mock"),
         timepoint_num = parse_integer(timepoint),
         replicate = parse_integer(replicate),
         timepoint = factor(timepoint_num),
         is_used = !(msrun %in% bad_msruns)) %>%
  arrange(treatment, timepoint, replicate) %>%
  mutate(msrun = factor(msrun, levels=msrun),
         condition = str_c(treatment, "_", timepoint, "h"),
         condition = factor(condition, levels=unique(condition)))

msdata_full$ptmgroup_intensities <- msdata_full$ptmgroup_intensities %>%
  dplyr::mutate(ident_type = factor(if_else(TRUE, "By MS/MS", "By matching"), levels = c("By MS/MS", "By matching")),
                msrun = factor(msrun, levels=levels(msdata_full$msruns$msrun)))
# msdata_full <- append_ptmgroups_info(msdata_full, msdata.wide,
#                                       proteins_info = dplyr::bind_rows(
#                                         dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE, is_expected_stable=FALSE),
#                                         dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE,
#                                                       # define stable proteins based on our expectations of their biology and how their timeseries look like
#                                                       is_expected_stable = str_detect(gene_name, "^(M?RP[LS]\\d+|GAPDH|ACT[ABNR]\\d+|TUB[AB]\\d?\\w?|CCT\\d?\\w?)$")),
#                                         dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE, is_expected_stable=FALSE)),
#                                       import_columns = c("is_viral", "is_contaminant", "is_expected_stable"))
# msdata_full$proteins <- mutate(msdata_full$proteins,
#                                protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

msdata_full$ptmgroups <- dplyr::mutate(msdata_full$ptmgroups,
                                       is_reverse = EG.IsDecoy,
                                       is_contaminant = FALSE, # FIXME
                                       is_viral = FALSE, # FIXME
                                       gene_label = strlist_label2(gene_names),
                                       protac_label = strlist_label2(majority_protein_acs), # FIXME wrong order
                                       protein_label = strlist_label2(protein_names),
                                       protein_description = strlist_label2(protein_descriptions),
                                       ptmgroup_label = ptmgroup_id)

# condition = treatment X timepoint
conditions.df <- dplyr::select(msdata_full$msruns, condition, treatment, timepoint, timepoint_num) %>%
  dplyr::distinct()
conditions.df <- dplyr::left_join(conditions.df,
                                  dplyr::select(conditions.df, treatment, timepoint_after=timepoint_num)) %>% # exclude 3h because it's the reference
  dplyr::mutate(is_after = timepoint_num >= timepoint_after) %>%
  tidyr::pivot_wider(all_of(c("condition", "treatment", "timepoint", "timepoint_num")),
                     names_prefix = "after", names_from = timepoint_after, values_from = is_after) %>%
  rename_at(vars(starts_with("after")), ~str_c(., "h"))

msdata_full$msrun_stats <- msrun_statistics(msdata_full, obj="ptmgroup") %>%
  dplyr::mutate(na_ratio = n_missing / n)

msdata <- msdata_full[c('ptmgroup_intensities',
                        'msruns', 'ptmgroups')]#, 'protein2ptmgroup')]
msdata$ptmgroup_intensities <- semi_join(msdata$ptmgroup_intensities, filter(msdata$msruns, is_used))

# setup experimental design matrices
conditionXeffect_orig.mtx <- model.matrix(
  formula(str_c("~ 1 + (", str_c(str_subset(colnames(conditions.df), "^after\\d+h"), collapse =" + "), ") * treatment")),
  conditions.df
)
conditionXeffect.mtx <- conditionXeffect_orig.mtx[, str_subset(colnames(conditionXeffect_orig.mtx),
                                                               str_c("^(\\(Intercept\\)|treatment[^:]+)$"), negate=TRUE)]
dimnames(conditionXeffect.mtx) <- list(condition = as.character(conditions.df$condition),
                                       effect = colnames(conditionXeffect.mtx))

pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  mutate(treatment = effect_factor(effect, 'treatment',
                                   levels(conditions.df$treatment), NA),
         timepoint = effect_factor(effect, 'after', str_c(unique(conditions.df$timepoint_num), "hTRUE"), NA) %>%
           as.character() %>% str_remove_all("^after|hTRUE$") %>% factor(levels = levels(conditions.df$timepoint)))

effects.df <- left_join(effects.df,
                        as_tibble(as.table(conditionXeffect.mtx)) %>%
                          filter(n != 0) %>%
                          mutate(timepoint = effect_factor(effect, 'after', str_c(unique(conditions.df$timepoint_num), "hTRUE"), "initial"),
                                 is_after = timepoint != "initial") %>%
                          filter(!is.na(timepoint)) %>%
                          pivot_wider("effect", names_prefix="after", names_from ="timepoint", values_from=is_after, values_fn = list(is_after=any))) %>%
  mutate_at(vars(starts_with("after")), ~coalesce(., FALSE)) %>%
  mutate(effect_type = case_when(!is.na(timepoint) & !is.na(treatment) ~ "treatmentXtimepoint",
                                 !is.na(timepoint) & is.na(treatment) ~ "timepoint",
                                 !is.na(treatment) ~ "treatment",
                                 TRUE ~ NA_character_),
         effect_label = case_when(effect_type == "treatmentXtimepoint" ~ str_c(treatment, "+", timepoint, "h"),
                                  effect_type == "treatment" ~ as.character(treatment),
                                  effect_type == "timepoint" ~ str_c(timepoint, "h"),
                                  TRUE ~ NA_character_),
         prior_mean = 0,
         prior_tau = case_when(effect_type == "treatmentXtimepoint" ~ 0.25,
                               effect_type == "treatment" ~ 0.5,
                               effect_type == "timepoint" ~ 1.0,
                               TRUE ~ NA_real_),
         is_positive = FALSE)

conditionXeffect.df <- conditionXeffect_frame(conditionXeffect.mtx, effects.df)
inv_conditionXeffect.mtx <- frame2matrix(conditionXeffect.df,
                                         "condition", "effect", "cond_w",
                                         rows = rownames(conditionXeffect.mtx),
                                         cols = colnames(conditionXeffect.mtx))
pheatmap(inv_conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_inv_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

msrunXreplEffect.mtx <- replicate_effects_matrix(
  msdata$msruns,#mutate(msdata$msruns, batch_condition=str_c("B", batch, "_", condition)),
  replicate_col = "replicate", condition_col = "condition")
pheatmap(msrunXreplEffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_msruns_", ms_folder, "_", fit_version, ".pdf")),
         width = 16, height = 20)

msrunXreplEffect.df <- as_tibble(as.table(msrunXreplEffect.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)

#TODO 0h?
#compound_metaconditions <- c("allTreatments_initialTimepoint")
compound_metaconditions <- c()
all_metaconditions <- c(levels(conditions.df$condition), compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in levels(conditions.df$condition)) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
}
#conditionXmetacondition.mtx[filter(conditions.df, timepoint_num == min(timepoint_num))$condition, "allTreatments_initialTimepoint"] <- TRUE
pheatmap(ifelse(conditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_metaconditions_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)

all_contrasts <- character()

contrasts.df <- bind_rows(
  # all treatment pairs at each timepoint
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_rhs = treatment, timepoint_lhs = timepoint))) %>%
    filter(as.integer(treatment_lhs) > as.integer(treatment_rhs)) %>%
    mutate(timepoint_rhs = timepoint_lhs,
           contrast_kind = "treatment_vs_treatment"),
  # all timepoints of the same treatment
  left_join(dplyr::mutate(dplyr::select(conditions.df, metacondition_lhs = condition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(conditions.df, metacondition_rhs = condition, treatment_lhs = treatment, timepoint_rhs = timepoint))) %>%
    filter(as.integer(timepoint_lhs) > as.integer(timepoint_rhs)) %>%
    mutate(treatment_rhs = treatment_lhs,
           contrast_kind = "timepoint_vs_timepoint")
) %>%
  mutate(contrast = str_c(treatment_lhs, "@", timepoint_lhs, "h_vs_", treatment_rhs, "@", timepoint_rhs, "h"),
         contrast_type = "comparison")
all_contrasts = contrasts.df$contrast

contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
for (i in 1:nrow(contrasts.df)) {
  contrastXmetacondition.mtx[contrasts.df$contrast[[i]],
                             c(contrasts.df$metacondition_lhs[[i]],
                               contrasts.df$metacondition_rhs[[i]])] <- c(1,-1)
}

contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::rename(weight = n) %>%
  dplyr::mutate(contrast_type = 'comparison',
                condition_role = "signal")

contrastXcondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::select(-n) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition)

pheatmap(contrastXmetacondition.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_contrasts_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 12)

##########
## normalization
msdata4norm.df <- msdata$ptmgroup_intensities %>% ungroup() %>%
  dplyr::filter(ident_type == "By MS/MS" & !is.na(intensity)) %>% 
  dplyr::semi_join(dplyr::filter(msdata$ptmgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::select(ptmgroup_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$ptmgroup_intensities) %>%
  dplyr::inner_join(select(msdata$msruns, msrun, condition, treatment, timepoint)) #%>% #WHY???
#mutate(timepoint=factor(timepoint))

options(mc.cores=8)

## level 1: msruns level (replicates within condition)
msruns_norm <- normalize_experiments(instr_calib, msdata4norm.df,
                                         quant_col = "intensity", obj_col = "ptmgroup_id", mschan_col = "msrun",
                                         cond_col="msrun", condgroup_col='condition',
                                         mcmc.iter = 2000L, stan_method="mcmc",
                                         mcmc.adapt_delta=0.95,
                                         verbose=TRUE, max_objs=1000L)
## level 2: condition level (normalize all treatments at each timepoint)
conditions_norm <- normalize_experiments(instr_calib, msdata4norm.df,
                                   quant_col = "intensity", obj_col = "ptmgroup_id", mschan_col = "msrun",
                                   cond_col="condition", condgroup_col="timepoint", 
                                   mcmc.iter = 2000L, stan_method="mcmc",
                                   mcmc.adapt_delta=0.95,
                                   verbose=TRUE, max_objs=1000L,
                                   mschan_preshifts = msruns_norm,
                                   preshift_col = 'shift')

## level 3: timepoint level (normalize all timepoints using (supposedly) unregulated sites)
### find unreg. proteins

# impute intensities
msdata4norm_imputed.df <- impute_intensities(tidyr::expand(msdata4norm.df, ptmgroup_id, nesting(condition, msrun)) %>%
                                             left_join(msdata4norm.df),
                                             msdata_full$msrun_stats) %>%
  dplyr::mutate(log2_intensity_imputed = log2(intensity_imputed))

msdata4norm_imputed.mtx <- matrix(msdata4norm_imputed.df$log2_intensity_imputed,
                                  nrow = n_distinct(msdata4norm_imputed.df$ptmgroup_id),
                                  dimnames = list(ptmgroup = unique(msdata4norm_imputed.df$ptmgroup_id),
                                                  msrun = unique(msdata4norm_imputed.df$msrun)))
msdata4norm_imputed_pg.mtx <- msdata4norm_imputed.mtx - rowMedians(msdata4norm_imputed.mtx)
# generate ptms umap
ptmgroup4norm_umap = umap(msdata4norm_imputed_pg.mtx, n_components=3,
                          n_neighbors=20, init="laplacian",
                          min_dist=0.2, metric = "euclidean")
ptmgroup4norm_umap_center <- colMedians(ptmgroup4norm_umap)

ptmgroupXcond_stats.df <- msdata4norm_imputed.df %>% 
  group_by(ptmgroup_id, condition) %>% 
  summarise(quant_freq = sum(!is.na(intensity))/n()) %>% 
  ungroup()
ptmgroup4norm_stats.df <- as_tibble(ptmgroup4norm_umap, .name_repair="minimal") %>%
  set_names(c("x_3d", "y_3d", "z_3d")) %>%
  mutate(ptmgroup_id = unique(ptmgroup_intensities4pca.df$ptmgroup_id),
         center_dist = sqrt(abs(x_3d - ptmgroup4norm_umap_center[[1]])^2 +
                            abs(y_3d - ptmgroup4norm_umap_center[[2]])^2 +
                            abs(z_3d - ptmgroup4norm_umap_center[[3]])^2)) %>%
  left_join(group_by(ptmgroupXcond_stats.df, ptmgroup_id) %>%
            summarise(quant_freq_min = min(quant_freq)) %>%
            ungroup())

protgroup_stats.df <- left_join(select(msdata_full$ptmgroups, ptmgroup_id, gene_names),
                                ptmgroup_stats.df) %>%
  group_by(gene_names) %>%
  summarise(n_ptms = n_distinct(ptmgroup_id),
            n_ptms_solid = n_distinct(ptmgroup_id[quant_freq_min >= 0.75])) %>%
  ungroup()

### find stable phosphosites (close to umap center and with enough identification for each condition)
View(filter(ptmgroup4norm_stats.df, center_dist <= 4 & quant_freq_min >= 0.75))

msdata4timepoint_norm.df <- semi_join(msdata4norm.df, filter(ptmgroup4norm_stats.df, center_dist <= 4 & quant_freq_min >= 0.75))

timepoints_norm <- normalize_experiments(instr_calib, msdata4timepoint_norm.df,
                                         quant_col = "intensity", obj_col = "ptmgroup_id", mschan_col = "msrun",
                                         cond_col='timepoint', condgroup_col=NULL, 
                                         mcmc.iter = 2000L, stan_method="mcmc", mcmc.adapt_delta=0.95,
                                         verbose=TRUE, max_objs=500L,
                                         mschan_preshifts = left_join(select(msruns_norm, msrun, condition, msrun_shift=shift),
                                                                      select(conditions_norm, condition, cond_shift=shift)) %>%
                                             mutate(shift_sum = msrun_shift + cond_shift),
                                         preshift_col = 'shift_sum')

## assemble all normalization levels
msruns_hnorm <- list(msruns = msruns_norm,
                     conditions = conditions_norm,
                     timepoints = timepoints_norm)

#total_msrun_shifts.df <- select(msruns_lev1norm, msrun, condition, msrun_shift = shift) %>%
#  mutate(total_msrun_shift = msrun_shift)

total_msrun_shifts.df <- select(msruns_hnorm$msruns, msrun, condition, msrun_shift = shift) %>%
  left_join(select(msruns_hnorm$conditions, condition, timepoint, condition_shift = shift)) %>%
  left_join(select(msruns_hnorm$timepoints, timepoint, timepoint_shift = shift)) %>%
  mutate(total_msrun_shift = msrun_shift + condition_shift + timepoint_shift,
         msrun = factor(msrun, levels=levels(msdata$msrun$msrun)))

## apply normalization
msdata$ptmgroup_intensities <- dplyr::left_join(dplyr::select(msdata$ptmgroup_intensities, -any_of(c("msrun_shift", "total_msrun_shift"))),
                                                dplyr::select(total_msrun_shifts.df, msrun, msrun_shift, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) 
#dplyr::select(-total_msrun_shift)

global_ptmgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns, TRUE) %>%
                                                 dplyr::select(msrun) %>% dplyr::distinct() %>%
                                                 dplyr::inner_join(msdata$ptmgroup_intensities) %>% .$intensity), na.rm=TRUE)

############################
# batch effects
# no batch effects so far
msrunXbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXreplEffect.mtx),
                                     batch_effect = c())

batch_effects.df <- tibble(batch_effect = character(0),
                           is_positive = logical(0),
                           prior_mean = double(0))

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', ms_folder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df,
     conditionXeffect.mtx, inv_conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     instr_calib_ptmgroup,
     global_ptmgroup_labu_shift,
     msruns_hnorm, total_msrun_shifts.df,
     msrunXreplEffect.mtx,
     batch_effects.df, msrunXbatchEffect.mtx,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full,
     file = rfulldata_filepath)

message('Done.')

set.seed(1232)
msdata_full$ptmgroup_intensities_all <- tidyr::expand(msdata_full$ptmgroup_intensities,
                                                      ptmgroup_id, msrun) %>%
  left_join(dplyr::select(msdata_full$ptmgroup_intensities, -any_of(c("total_msrun_shift", "msrun_shift")))) %>%
  impute_intensities(msdata_full$msrun_stats) %>%
  #dplyr::mutate(total_msrun_shift=0, msrun_shift=0) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift, msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, ptmgroup_id)

ptmgroup_intensities4pca.df <- msdata_full$ptmgroup_intensities_all %>%
  dplyr::semi_join(dplyr::filter(msdata_full$ptmgroups, !is_reverse & !is_contaminant)) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$msruns, is_used)) %>%
  dplyr::arrange(msrun, ptmgroup_id)

ptmgroup_intensities.mtx <- matrix(log2(ptmgroup_intensities4pca.df$intensity_norm),
                                    nrow = n_distinct(ptmgroup_intensities4pca.df$ptmgroup_id),
                                    dimnames = list(ptmgroup = unique(ptmgroup_intensities4pca.df$ptmgroup_id),
                                                    msrun = unique(ptmgroup_intensities4pca.df$msrun)))
ptmgroup_intensities_imp.mtx <- matrix(log2(ptmgroup_intensities4pca.df$intensity_imputed_norm),
                                        nrow = n_distinct(ptmgroup_intensities4pca.df$ptmgroup_id),
                                        dimnames = list(ptmgroup = unique(ptmgroup_intensities4pca.df$ptmgroup_id),
                                                        msrun = unique(ptmgroup_intensities4pca.df$msrun)))
require(FactoMineR)
require(Matrix)
msrun_intensities_pca <- PCA(ptmgroup_intensities_imp.mtx, graph = FALSE)
colnames(msrun_intensities_pca$svd$V) <- str_c("comp_", 1:ncol(msrun_intensities_pca$svd$V))
msrun_intensities_pca.df <- as_tibble(msrun_intensities_pca$svd$V) %>%
  dplyr::mutate(msrun = rownames(msrun_intensities_pca$var$coord)) %>%
  dplyr::inner_join(msdata_full$msruns)

require(ggrepel)
p <- ggplot(msrun_intensities_pca.df %>% arrange(treatment, timepoint),
            aes(x=comp_1, y=comp_2, color=condition)) +
  geom_point() +
  geom_text_repel(aes(label=msrun), size=3, vjust=-1.1) +
  theme_bw_ast(base_family = "", base_size = 10) #+
ggsave(p, filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_msruns_pca_norm_", data_version, ".pdf")),
       width = 12, height = 14)

require(pheatmap)
msruns_ordered <- as.character(msdata_full$msruns$msrun)
ptmgroup_hclu = hclust(dist(ptmgroup_intensities_imp.mtx))
pheatmap(ptmgroup_intensities.mtx[, msruns_ordered], cluster_cols=FALSE, cluster_rows=ptmgroup_hclu,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.pdf")), width=16, height=80)

viral_mask <- unique(ptmgroup_intensities4pca.df$ptmgroup_id) %in% filter(msdata_full$ptmgroups, is_viral)$ptmgroup_id
viral_intensities.mtx <- ptmgroup_intensities.mtx[viral_mask, msruns_ordered]
rownames(bait_intensities.mtx) <- 
  tibble(ptmgroup_id = unique(ptmgroup_intensities4pca.df$ptmgroup_id)[viral_mask]) %>%
  inner_join(msdata_full$ptmgroups) %>% inner_join(msdata_full$protein2ptmgroup) %>%
  inner_join(mutate(msdata_full$proteins, full_code = str_c(str_replace(organism, "CoV-2", "CoV2"), '-', gene_name))) %>%
  group_by(ptmgroup_id) %>%
  summarise(full_code = str_c(full_code, collapse=';')) %>% ungroup() %>% .$full_code
pheatmap(bait_intensities.mtx, cluster_cols=FALSE, cluster_rows=FALSE,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_viral_heatmap_intensity.pdf")), width=16, height=8)

require(uwot)
require(matrixStats)
ptmgroup_intensities_imp_pg.mtx <- ptmgroup_intensities_imp.mtx - rowMedians(ptmgroup_intensities_imp.mtx)
objects_umap2d = umap(ptmgroup_intensities_imp_pg.mtx, n_components=2,
                 n_neighbors=20, init="laplacian",
                 min_dist=0.2, metric = "euclidean")

objects_umap.df <- as_tibble(objects_umap2d, .name_repair="minimal") %>%
  set_names(c("x_2d", "y_2d")) %>%
  mutate(object_id = unique(ptmgroup_intensities4pca.df$ptmgroup_id)) %>%
  #dplyr::left_join(
  #  as_tibble(objects.umap3d, .name_repair="minimal") %>%
  #    set_names(c("x_3d", "y_3d", "z_3d")) %>%
  #    mutate(object_id = unique(object_intensities4pca.df$object_id))    
  #) %>%
  dplyr::inner_join(select(msdata$ptmgroups, object_id=ptmgroup_id, object_label=ptmgroup_label,
                           majority_protein_acs, gene_names, protein_descriptions,
                           is_viral, #is_expected_stable,
                           is_contaminant, is_reverse)) %>%
  mutate(category = case_when(replace_na(is_viral, FALSE) ~ "viral",
                              is_reverse ~ "reverse",
                              is_contaminant ~ "contaminant",
                              is.na(gene_names) ~ "NA",
                              #is_expected_stable ~ "stable",
                              str_detect(gene_names, "(^|;)(?:HIST[123])?H[123][ABCDEF]") ~ "histone",
                              str_detect(gene_names, "(^|;)PPP[123456]") ~ "phosphatase",
                              str_detect(gene_names, "(^|;)RP[SL]\\d+") ~ "ribosome",
                              str_detect(gene_names, "(^|;)MRP[SL]\\d+") ~ "mito ribosome",
                              str_detect(gene_names, "(^|;)PSM[ABCD]\\d+") ~ "proteasome",
                              str_detect(gene_names, "(^|;)(MAVS|OAS\\d+|IFIT\\d+|PARP\\d+|DDX58|IFI\\d+|ISG\\d+|MX\\d+|TGFB\\d+|UNC93B1|ZC3HAV1)(;|$)") ~ "ISG",
                              TRUE ~ "default"
  ))# %>%
#dplyr::left_join(filter(protregroup_effects.df, effect == "treatmentSC35M"))
# FIXME duplicate point due to N-to-N ptmgroup to protregroup matching
#dplyr::left_join(filter(object_effects_wide.df, std_type == "replicate") %>%
#                  select(-gene_names, -majority_protein_acs))
require(plotly)
# plot umap + selected proteins (color = protein biological function)
umap <- ggplotly(
  ggplot(objects_umap.df,
         aes(x=x_2d, y=y_2d, color=category)) +
    geom_point(aes(text=object_label, size=!(category %in% c('default', 'NA')))) +
    #geom_text(data=filter(objects_umap.df, gene_names %in% c("RPL11", "RPL12")),
    #          aes(label=object_label), vjust=-1.1, color = "black", size = 3) +
    #geom_text_repel(data=filter(objects_umap.df, is_viral),
    #                aes(label=object_label)) +
    scale_size_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5)) +
    scale_color_manual(values = c("stable" = "black", "viral" = "red", "ISG" = "blue", "hit" = "orange",
                                  "ribosome" = "darkgreen", "mito ribosome" = "darkred", "proteasome" = "khaki",
                                  "histone" = "red", "phosphatase" = "violet",
                                  "default" = "gray",
                                  "NA" = "pink", "reverse" = "pink", "contaminant" = "yellow")) +
    theme_bw_ast(base_family = "", base_size = 10),
  tooltip = "text"
)
umap

protgroup_stats.df <- left_join(select(msdata_full$ptmgroups, ptmgroup_id, gene_names),
                                ptmgroup_stats.df) %>%
  group_by(gene_names) %>%
  summarise(n_ptms = n_distinct(ptmgroup_id),
            n_ptms_solid = n_distinct(ptmgroup_id[quant_freq_min >= 0.75])) %>%
  ungroup()

tmp.df <- mutate(objects_umap.df,
       x_2d_median = median(x_2d),
       y_2d_median = median(y_2d),
       center_dist = sqrt(abs(x_2d - x_2d_median)^2 + abs(y_2d - y_2d_median)^2)) %>%
  left_join(protgroup_stats.df)

p <- ggplot(filter(msdata_full$ptmgroups, str_detect(ptmgroup_id, "(^|;)(?:ATXN2)")) %>%#"(^|;)(?:EIF[34]B)")
            semi_join(filter(ptmgroup4norm_stats.df, quant_freq_min >= 0.75)) %>%
            inner_join(bind_rows(transmute(msdata_full$ptmgroup_intensities_all, ptmgroup_id, msrun,
                                           intensity=intensity_imputed*exp(-msrun_shift), is_imputed, is_norm=FALSE),
                                 transmute(msdata_full$ptmgroup_intensities_all, ptmgroup_id, msrun,
                                           intensity=intensity_imputed_norm, is_imputed, is_norm=TRUE))) %>%
            left_join(msdata_full$msruns),
            aes(x = timepoint_num, y = intensity, color=treatment, fill=treatment)) +
  geom_smooth(alpha=0.25) +
  geom_point(aes(shape=is_imputed)) +
  scale_x_continuous(breaks = unique(conditions.df$timepoint_num)) +
  scale_shape_manual(values = c("FALSE"=16L, "TRUE"=1L)) +
  scale_color_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  scale_fill_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  facet_wrap(~ ptmgroup_id + is_norm, ncol = 4, scales = "free_y") +
  theme_bw_ast()
p
ggsave(p, filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_ATXN2_", data_version, ".pdf")),
       width = 20, height = 12, device = cairo_pdf)
