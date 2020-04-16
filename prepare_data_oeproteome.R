# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200411"
fit_version <- "20200411"
ms_folder <- 'spectronaut_oeproteome_20200411'
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

msdata_path <- file.path(data_path, ms_folder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  ms_folder = ms_folder,
                  instr_calib_protgroup_filename = "instr_QX8_intensity_protgroup_calib_cov2_20200411_borg.json",
                  quant_type = "intensity", quant_col_prefix = "intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib_protgroup <- fromJSON(file = file.path(data_path, data_info$instr_calib_protgroup_filename))$instr_calib

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/cov_baits_20200331.fasta")),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants.fasta"))
)

bad_msruns <- c("FPMS_SARS_CoV2_ORF3_3", "FPMS_SARS_CoV2_ORF7a_2")
msdata.wide <- read.Spectronaut.ProteinsReport(file.path(msdata_path, "20200410_COV2_B1_DIA_DDA library_proteinreport.csv"),
                                               import_data = "quantity")
msdata_colgroups <- attr(msdata.wide, "column_groups")

msdata_full <- list(
  protgroups = msdata.wide[, msdata_colgroups$protgroup],
  protgroup_intensities = pivot_longer.Spectronaut.ProtgroupIntensities(msdata.wide)
)
msdata_full$msruns <- dplyr::select(msdata_full$protgroup_intensities, msrun_ix, raw_file) %>% dplyr::distinct() %>%
  dplyr::arrange(msrun_ix) %>%
  dplyr::mutate(msrun_sn = str_remove(str_remove(raw_file, "^20200326_QX8_OzKa_SA_"), "(?:_\\d{3,})?.raw$")) %>%
  tidyr::extract(msrun_sn, c("bait_code", "replicate"), "^COV2_FPMS_([^_]+)_(\\d+)$", remove=FALSE) %>%
  dplyr::mutate(replicate = parse_integer(replicate)) %>%
  left_join(dplyr::select(baits_info.df, bait_type, bait_code, bait_full_id, bait_id, organism, orgcode)) %>%
  mutate(condition = str_c("FPMS_", bait_full_id),
         msrun = str_c(condition, "_", replicate),
         batch = if_else(bait_full_id %in% c("SARS_CoV2_E", "SARS_CoV2_M", "SARS_CoV2_N", "SARS_CoV2_NSP15", "SARS_CoV2_NSP16", "SARS_CoV2_ORF3"), 2L, 1L),
         is_used = !msrun %in% bad_msruns)

msdata_full$protgroup_intensities <- dplyr::select(msdata_full$protgroup_intensities, -raw_file) %>%
  dplyr::mutate(ident_type = factor(if_else(nevidences > 0L, "By MS/MS", "By matching"), levels = c("By MS/MS", "By matching"))) %>%
  left_join(dplyr::select(msdata_full$msruns, msrun_ix, msrun))
msdata_full <- append_protgroups_info(msdata_full, msdata.wide,
                                      proteins_info = dplyr::bind_rows(
                                        dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE),
                                        dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE),
                                        dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE)),
                                      import_columns = c("is_viral", "is_contaminant"))
msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
    is_reverse = FALSE,
    gene_label = strlist_label2(gene_names),
    protac_label = strlist_label2(protein_acs),
    protgroup_label = case_when(!is.na(gene_label) ~ gene_label,
                                !is.na(protac_label) ~ protac_label,
                                TRUE ~ str_c('#', protgroup_id)))

# condition = bait
conditions.df <- dplyr::select(msdata_full$msruns, condition, bait_full_id, bait_id, bait_type, orgcode) %>%
  dplyr::distinct() %>%
  dplyr::arrange(bait_type, bait_id, orgcode) %>%
  dplyr::mutate(condition = factor(condition, levels=condition))
msdata_full$msruns <- left_join(dplyr::select(msdata_full$msruns, -any_of("condition")),
                                dplyr::select(conditions.df, bait_full_id, condition))

msdata <- msdata_full[c('protgroup_intensities',
                        'msruns', 'protgroups', 'protein2protgroup')]
msdata$protgroup_intensities <- semi_join(msdata$protgroup_intensities, filter(msdata$msruns, is_used))

# setup experimental design matrices
conditionXeffect_orig.mtx <- model.matrix(
  ~ 1 + bait_id + bait_id:orgcode,
  mutate(conditions.df, orgcode = if_else(bait_type == "sample", orgcode, factor("CVHSA2", levels=levels(orgcode)))))
conditionXeffect.mtx <- conditionXeffect_orig.mtx[, colSums(abs(conditionXeffect_orig.mtx)) != 0 &
                                                    !str_detect(colnames(conditionXeffect_orig.mtx), "\\(Intercept\\)|bait_idCtrl.+:orgcode")]
dimnames(conditionXeffect.mtx) <- list(condition = conditions.df$condition,
                                       effect = colnames(conditionXeffect.mtx))

pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  dplyr::mutate(orgcode = effect_factor(effect, "orgcode", levels(conditions.df$orgcode), NA),
                bait_id = effect_factor(effect, "bait_id", levels(conditions.df$bait_id), NA),
                is_positive = FALSE,#!is.na(bait_id) & is.na(orgcode),
                prior_mean = 0.0)
effects.df$effect_label <- sapply(1:nrow(effects.df), function(i) {
  comps <- c()
  org <- as.character(effects.df$orgcode[[i]])
  ba <- as.character(effects.df$bait_id[[i]])
  if (!is.na(ba)) comps <- append(comps, ba)
  if (!is.na(org)) comps <- append(comps, org)
  str_c(comps, collapse="+")
})
effects.df <- dplyr::mutate(effects.df,
                            effect = factor(effect, levels = effect),
                            effect_type = case_when(!is.na(bait_id) & !is.na(orgcode) ~ "baitXvirus",
                                                    !is.na(bait_id) ~ "bait"),
                            effect_label = factor(effect_label, levels=effect_label)) %>%
  dplyr::mutate(prior_tau = case_when(effect_type == "baitXvirus" ~ 0.25,
                                      effect_type == "bait" ~ 1.0,
                                      TRUE ~ 5.0))
# prior_mean is calculated after the normalization

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

bait_conditions <- as.character(filter(conditions.df, bait_type == "sample")$condition)
allminus_metaconditions <- paste0("FPMS_allminus_", unique(filter(conditions.df, bait_type == "sample")$bait_id))
compound_metaconditions <- c(allminus_metaconditions, "FPMS_controls")
all_metaconditions <- c(bait_conditions, compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in bait_conditions) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
}
for (cname in allminus_metaconditions) {
  bait <- str_remove(cname, "^FPMS_allminus_")
  conditionXmetacondition.mtx[, cname] <- TRUE
  conditionXmetacondition.mtx[str_detect(rownames(conditionXmetacondition.mtx), str_c("_", bait, "$")), cname] <- FALSE
}
conditionXmetacondition.mtx[dplyr::filter(conditions.df, bait_type == "control")$condition, "FPMS_controls"] <- TRUE
pheatmap(ifelse(conditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(data_path, paste0(project_id, "_metaconditions_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n) %>%
  dplyr::mutate(is_preserved_condition = FALSE)#condition %in% c("Ctrl_NT", "Ctrl_Gaussia_luci"))

contrasts.df <- bind_rows(
  transmute(filter(conditions.df, bait_type == "sample"),
            metacondition_lhs = condition,
            metacondition_rhs = "FPMS_controls",
            contrast_type = "comparison"),
  transmute(filter(conditions.df, bait_type == "sample"),
            metacondition_lhs = condition,
            metacondition_rhs = str_c("FPMS_allminus_", bait_id),
            contrast_type = "comparison"),
  inner_join(
    select(filter(conditions.df, bait_type == "sample"), metacondition_lhs = condition, orgcode_lhs = orgcode, bait_id),
    select(filter(conditions.df, bait_type == "sample"), metacondition_rhs = condition, orgcode_rhs = orgcode, bait_id),
  ) %>% filter(as.integer(orgcode_lhs) < as.integer(orgcode_rhs)) %>%
    mutate(contrast_type = "comparison") %>%
    select(-orgcode_lhs, -orgcode_rhs, -bait_id)) %>%
  mutate(contrast = str_c(metacondition_lhs, "_vs_", ifelse(str_starts(metacondition_rhs, "allminus"), "others", metacondition_rhs)))

all_contrasts <- contrasts.df$contrast
contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
for (i in 1:nrow(contrasts.df)) {
  contrastXmetacondition.mtx[as.character(contrasts.df$contrast[[i]]),
                             c(as.character(contrasts.df$metacondition_lhs[[i]]),
                               as.character(contrasts.df$metacondition_rhs[[i]]))] <- c(1, -1)
}
pheatmap(contrastXmetacondition.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(data_path, paste0(project_id, "_exp_design_contrasts_", ms_folder, "_", fit_version, ".pdf")),
         width = 11, height = 12)

contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::rename(weight = n) %>%
  dplyr::left_join(select(contrasts.df, contrast, contrast_type)) %>%
  dplyr::mutate(condition_role = if_else(contrast_type == "filter" & weight < 0, "background", "signal"))

contrastXcondition.df <- conditionXmetacondition.df %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition)

# prepare the data to use for MS runs normalization
msrun_stats.df <- msdata$protgroup_intensities %>%
  dplyr::group_by(msrun) %>%
  dplyr::summarise(n_pg_idents = sum(!is.na(ident_type) & ident_type == "By MS/MS"),
                   n_pg_quants = sum(!is.na(intensity))) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(msdata$msruns)

# normalize using the intensities
msdata4norm.df <- msdata$protgroup_intensities %>%
  dplyr::filter(npeptides_quanted > 1 & !is.na(intensity)) %>%
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_contaminant & !is_viral)) %>%
  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$protgroup_intensities)

options(mc.cores=8)

# normalize experiments:
# 1) MS replicates for a given bait
# 2) same viral protein of different strains
# 3) all baits together
msruns_hnorm <- multilevel_normalize_experiments(instr_calib_protgroup,
                                                 filter(msdata$msruns, is_used) %>%
                                                   mutate(batch=as.character(batch),
                                                          batch_bait_full_id = str_c("B", batch, "_", bait_full_id),
                                                          batch_bait_id = str_c("B", batch, "_", bait_id)),
                                                 msdata4norm.df,
                                                 quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
                                                 mcmc.iter = 2000L,
                                                 #mcmc.chains = 6,
                                                 verbose=TRUE,
                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=700L, missing_exp.ratio=0.1),
                                                                    bait_full_id = list(cond_col="batch_bait_full_id", max_objs=500L, missing_exp.ratio=0.1),
                                                                    bait_id = list(cond_col="batch_bait_id", max_objs=300L, missing_exp.ratio=0.1),
                                                                    batch = list(cond_col="batch", max_objs=200L, missing_exp.ratio=0.1)
                                                 ))

# ignore all higher levels of normalization
msruns_hnorm$msruns_shifts <- msruns_hnorm$mschannel_shifts

total_msrun_shifts.df <- msruns_hnorm$msruns_shifts

# apply normalization
msdata_full$protgroup_intensities <- dplyr::left_join(msdata_full$protgroup_intensities,
                                                      dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift))

# use only native medium to calculate average intensity
global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns) %>%
                                                 dplyr::select(msrun) %>% dplyr::distinct() %>%
                                                 dplyr::inner_join(msdata_full$protgroup_intensities) %>% pull(intensity)), na.rm=TRUE)

msrunXbatchEffect_orig.mtx <- model.matrix(
  ~ 1 + batch,
  mutate(msdata$msruns, batch = factor(batch)))
msrunXbatchEffect.mtx <- msrunXbatchEffect_orig.mtx[, colnames(msrunXbatchEffect_orig.mtx) != "(Intercept)", drop=FALSE]
dimnames(msrunXbatchEffect.mtx) <- list(msrun = msdata$msruns$msrun,
                                        batch_effect = colnames(msrunXbatchEffect_orig.mtx)[colnames(msrunXbatchEffect_orig.mtx) != "(Intercept)"])
pheatmap(msrunXbatchEffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE)

batch_effects.df <- tibble(batch_effect=colnames(msrunXbatchEffect.mtx),
                           is_positive=FALSE,
                           prior_mean = 0.0)

# no subbatch effects so far
msrunXsubbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXreplEffect.mtx),
                                        subbatch_effect = c())

subbatch_effects.df <- tibble(subbatch_effect=character(0),
                              is_positive=logical(0),
                              prior_mean = double(0))

bait_checks_protgroup.df <- dplyr::left_join(dplyr::select(baits_info.df, bait_full_id, bait_id, bait_orgcode=orgcode, bait_organism=organism,
                                                           protein_ac = used_uniprot_ac),
                                             dplyr::select(msdata_full$proteins, protein_ac, protgroup_id, prot_organism=organism)) %>%
  dplyr::left_join(dplyr::select(dplyr::filter(msdata$protgroup_intensities, ident_type=="By MS/MS"), protgroup_id, msrun)) %>%
  dplyr::left_join(dplyr::select(msdata$msruns, msrun, observing_bait_full_id = bait_full_id)) %>%
  dplyr::arrange(bait_full_id, protgroup_id, msrun) %>%
  dplyr::group_by(bait_full_id, protgroup_id) %>%
  dplyr::mutate(idented_in_msruns = str_c(unique(msrun), collapse=";"),
                idented_in_FP_of = str_c(unique(observing_bait_full_id), collapse=";"),
                prot_organisms = str_c(unique(prot_organism), collapse=';')) %>%
  dplyr::filter(row_number()==1L) %>%
  dplyr::select(-msrun, -observing_bait_full_id) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(idented_in_msruns = if_else(idented_in_msruns == "", NA_character_, idented_in_msruns),
                idented_in_FP_of = if_else(idented_in_FP_of == "", NA_character_, idented_in_FP_of))

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', ms_folder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df,
     conditionXeffect.mtx, inv_conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     instr_calib_protgroup,
     global_protgroup_labu_shift,
     msruns_hnorm, total_msrun_shifts.df,
     msrunXreplEffect.mtx,
     batch_effects.df, msrunXbatchEffect.mtx,
     subbatch_effects.df, msrunXsubbatchEffect.mtx,
     bait_checks_protgroup.df,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full,
     file = rfulldata_filepath)

message('Done.')

msdata_full$mschannels <- mutate(msdata_full$msruns,
                                 mstag = "Sum",
                                 mschannel = msrun,
                                 msrun_mq = msrun)
msdata_full$protgroup_intensities <- mutate(msdata_full$protgroup_intensities,
                                            ident_type="By MS/MS")
msdata_full$protgroup_tagintensities <- mutate(msdata_full$protgroup_intensities,
                                               mschannel = msrun,
                                               msrun_mq = msrun,
                                               mstag = "Sum")
msdata_full$mschannel_stats <- mschannel_statistics(msdata_full)
set.seed(1232)
msdata_full$protgroup_intensities_all <- tidyr::expand(msdata_full$protgroup_intensities,
                                                protgroup_id, msrun) %>%
  left_join(dplyr::select(msdata_full$protgroup_intensities, -any_of("total_msrun_shift"))) %>%
  dplyr::mutate(mstag = "Sum") %>%
  impute_intensities(msdata_full$mschannel_stats) %>%
  dplyr::mutate(total_msrun_shift=0) %>%
  #dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
  #filter(mstag == "L") %>%
  dplyr::semi_join(dplyr::filter(msdata_full$protgroups, !is_reverse & !is_contaminant)) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$msruns, is_used)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_norm),
                                    nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                    dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                    msrun = unique(protgroup_intensities4pca.df$msrun)))
protgroup_intensities_imp.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_imputed_norm),
                                        nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                        dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                        msrun = unique(protgroup_intensities4pca.df$msrun)))
require(FactoMineR)
require(Matrix)
msrun_intensities_pca <- PCA(protgroup_intensities_imp.mtx, graph = FALSE)
msrun_intensities_pca.df <- as.data.frame(msrun_intensities_pca$svd$V)
colnames(msrun_intensities_pca.df) <- paste0("comp_", 1:ncol(msrun_intensities_pca.df))
msrun_intensities_pca.df <- dplyr::mutate(msrun_intensities_pca.df,
                                          msrun = rownames(msrun_intensities_pca$var$coord)) %>%
    dplyr::inner_join(msdata_full$msruns)

require(ggrepel)
p <- ggplot(msrun_intensities_pca.df %>%
            dplyr::mutate(msrun = str_remove(msrun, "FPMS_(?:SARS_)?")),
       aes(x=comp_1, y=comp_2, color=bait_id)) +
    geom_point() +
    geom_text_repel(aes(label=str_remove(str_remove(msrun, "FPMS_SARS_"), "APMS_")), size=3, vjust=-1.1) +
    theme_bw_ast(base_family = "", base_size = 10) #+
ggsave(p, filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_msruns_pca_", data_version, ".pdf")),
       width = 12, height = 14)

require(pheatmap)
msruns_ordered <- filter(msdata_full$msruns, msrun %in% colnames(protgroup_intensities_imp.mtx)) %>%
  dplyr::arrange(bait_type, bait_id, bait_full_id, replicate) %>% pull(msrun)
protgroup_hclu = hclust(dist(protgroup_intensities_imp.mtx))
pheatmap(protgroup_intensities.mtx[, msruns_ordered], cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.pdf")), width=20, height=100)

viral_mask <- unique(protgroup_intensities4pca.df$protgroup_id) %in% filter(msdata_full$protgroups, is_viral)$protgroup_id
viral_intensities.mtx <- protgroup_intensities.mtx[viral_mask, msruns_ordered]
rownames(bait_intensities.mtx) <- 
  tibble(protgroup_id = unique(protgroup_intensities4pca.df$protgroup_id)[viral_mask]) %>%
  inner_join(msdata_full$protgroups) %>% inner_join(msdata_full$protein2protgroup) %>%
  inner_join(mutate(msdata_full$proteins, full_code = str_c(str_replace(organism, "CoV-2", "CoV2"), '-', gene_name))) %>%
  group_by(protgroup_id) %>%
  summarise(full_code = str_c(full_code, collapse=';')) %>% ungroup() %>% .$full_code
pheatmap(bait_intensities.mtx, cluster_cols=FALSE, cluster_rows=FALSE,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_viral_heatmap_intensity.pdf")), width=16, height=8)

require(uwot)
require(matrixStats)
protgroup_intensities_imp_norm.mtx <- protgroup_intensities_imp.mtx - rowMedians(protgroup_intensities_imp.mtx)
pg_umap2d = umap(protgroup_intensities_imp_norm.mtx, n_components=2,
                 n_neighbors=20, init="laplacian",
                 min_dist=0.2, metric = "euclidean")

objects_umap.df <- as_tibble(pg_umap2d, .name_repair="minimal") %>%
  set_names(c("x_2d", "y_2d")) %>%
  mutate(object_id = unique(protgroup_intensities4pca.df$protgroup_id)) %>%
  #dplyr::left_join(
  #  as_tibble(objects.umap3d, .name_repair="minimal") %>%
  #    set_names(c("x_3d", "y_3d", "z_3d")) %>%
  #    mutate(object_id = unique(object_tagintensities4pca.df$object_id))    
  #) %>%
  dplyr::inner_join(select(msdata$protgroups, object_id=protgroup_id, object_label=protgroup_label,
                           majority_protein_acs, gene_names, is_viral, is_contaminant, is_reverse)) %>%
  mutate(category = case_when(replace_na(is_viral, FALSE) ~ "viral",
                              is_reverse ~ "reverse",
                              is_contaminant ~ "contaminant",
                              is.na(gene_names) ~ "NA",
                              str_detect(gene_names, "(^|;)RP[SL]\\d+") ~ "ribosome",
                              str_detect(gene_names, "(^|;)MRP[SL]\\d+") ~ "mito ribosome",
                              str_detect(gene_names, "(^|;)PSM[ABCD]\\d+") ~ "proteasome",
                              str_detect(gene_names, "(^|;)(MAVS|OAS\\d+|IFIT\\d+|PARP\\d+|DDX58|IFI\\d+|ISG\\d+|MX\\d+|TGFB\\d+|UNC93B1|ZC3HAV1)(;|$)") ~ "ISG",
                              TRUE ~ "default"
  ))# %>%
  #dplyr::left_join(filter(protregroup_effects.df, effect == "treatmentSC35M"))
  # FIXME duplicate point due to N-to-N protgroup to protregroup matching
  #dplyr::left_join(filter(object_effects_wide.df, std_type == "replicate") %>%
   #                  select(-gene_names, -majority_protein_acs))
require(plotly)
# plot umap + selected proteins (color = protein biological function)
ggplotly(
  ggplot(objects_umap.df,
         aes(x=x_2d, y=y_2d, color=category)) +
    geom_point(aes(text=object_label, size=!(category %in% c('default', 'NA')))) +
    #geom_text(data=filter(objects_umap.df, gene_names %in% c("RPL11", "RPL12")),
    #          aes(label=object_label), vjust=-1.1, color = "black", size = 3) +
    #geom_text_repel(data=filter(objects_umap.df, is_viral),
    #                aes(label=object_label)) +
    scale_size_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5)) +
    scale_color_manual(values = c("viral" = "red", "ISG" = "blue", "hit" = "orange",
                                  "ribosome" = "darkgreen", "mito ribosome" = "darkred", "proteasome" = "khaki",
                                  "default" = "gray",
                                  "NA" = "pink", "reverse" = "pink", "contaminant" = "yellow")) +
    theme_bw_ast(base_family = "", base_size = 10),
  tooltip = "text"
)
