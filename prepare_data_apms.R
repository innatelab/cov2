# SARS-CoV/CoV-2 loading and preparing AP-MS data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200329"
fit_version <- "20200329"
mq_folder <- 'mq_apms_20200329'
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

mqdata_path <- file.path(data_path, mq_folder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  mq_folder = mq_folder,
                  instr_calib_filename = "instr_protgroup_LFQ_calib_scaturro_qep5calib_20161110_borg.json",
                  quant_type = "LFQ", quant_col_prefix = "LFQ_Intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib <- fromJSON(file = file.path(data_path, data_info$instr_calib_filename))$instr_calib

orgcodes = list("SARS-CoV-2" = "CVHSA2",
                "SARS-CoV" = "CVHSA",
                "SARS-CoV-GZ02" = "CVHSA",
                "HCoV-NL63" = "CVHNL",
                "HCoV-229E" = "CVH22",
                "ZIKV" = "ZIKV",
                "Gaussia" = "9MAXI",
                "HCV" = "9HEPC")

baits_info.df <- read_tsv(file.path(analysis_path, "data", "baits_info.txt"))
baits_info.df <- mutate(baits_info.df,
                        used_uniprot_ac = ifelse(!is.na(uniprot_ac), uniprot_ac, bait_id)) %>%
  rename(organism = virus, bait_full_id = bait_id, bait_id = short_name) %>%
  mutate(bait_id = relevel(factor(bait_id), "Ctrl_NT"),
         bait_type = factor(bait_type, c("sample", "control")),
         orgcode = factor(orgcodes[organism], levels=unique(as.character(orgcodes))),
         protein_name = str_c(if_else(is.na(uniprot_ac), as.character(bait_full_id), uniprot_ac), "_", orgcode))

msruns.df <- read_tsv(file.path(mqdata_path, "combined", "experimentalDesign.txt"),
                      col_names=TRUE, col_types = list(Fraction="i")) %>%
  rename(raw_file=Name, msfraction=Fraction, msrun=Experiment, is_ptm=PTM) %>%
  extract(msrun, c("bait_full_id", "replicate"), c("^APMS_(.+)_(\\d+)$"), remove = FALSE) %>%
  left_join(select(baits_info.df, bait_full_id, bait_id, bait_type, orgcode))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(mqdata_path, "fasta/cov_baits_20200326.fasta")),
  human = read_innate_uniprot_fasta(file.path(mqdata_path, "fasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta"))
)

msdata.wide <- read.MaxQuant.ProteinGroups(file.path(mqdata_path, 'combined/txt'),
                                           file_name = "proteinGroups.txt",
                                           import_data = c(data_info$quant_type, "ident_type"))
mqevidence <- read.MaxQuant.Evidence(file.path(mqdata_path, 'combined/txt'), 
                                     evidence.pepobj = "pepmodstate", correct_ratios = FALSE)
mqevidence$peptides <- read.MaxQuant.Peptides(file.path(mqdata_path, 'combined/txt'),
                                              file_name = "peptides.txt",
                                              import_data = c(data_info$pep_quant_type, "ident_type"))
mqevidence$peaks <- NULL # exclude big data frame

mqrdata_filepath <- file.path(data_path, str_c(project_id, '_mqdata_APMS_', mq_folder, '.RData'))
message('Saving imported MaxQuant MS data to ', mqrdata_filepath, '...')
save(data_info, msruns.df, fasta.dfs,
     msdata.wide, mqevidence, instr_calib,
     file = mqrdata_filepath, compress = "xz")
load(mqrdata_filepath)

msdata_colgroups <- attr(msdata.wide, "column_groups")

msdata_full <- list()

msruns.df <- mutate(msruns.df, is_used = msrun %in% mqevidence$mschannels$msrun)
msdata_full$msruns <- filter(msruns.df, is_used)

# prepare protgroup intensities (wider format: all mstags in one row)
intensity_prespec_df <- tibble(.name = msdata_colgroups$LFQ) %>%
  extract(.name, c("mstag", "msrun"), remove=FALSE,
          str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
  mutate(.value = str_c("intensity.", mstag)) %>%
  dplyr::inner_join(select(msruns.df, msrun, raw_file))

# prepare protgroup intensities (longer format: each mstag on ots own row)
protgroup_intensities_all.df <- tidyr::pivot_longer_spec(
  dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$LFQ),
  mutate(intensity_prespec_df, .value = "intensity")) %>%
  select(-mstag)

protgroup_idents_all.df <- tidyr::pivot_longer(
  dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$ident_type),
  cols = msdata_colgroups$ident_type,
  values_to = "ident_type", names_to = c("msrun"),
  names_prefix = "ident_type\\.") %>%
  dplyr::left_join(select(msruns.df, msrun) %>% distinct())

msdata_full$protgroup_intensities <- protgroup_intensities_all.df %>%
  dplyr::semi_join(select(msdata_full$msruns, msrun, raw_file)) %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::select(-raw_file)
msdata_full$protgroup_idents <- protgroup_idents_all.df %>%
  dplyr::semi_join(select(msdata_full$msruns, msrun)) %>%
  dplyr::filter(!is.na(ident_type))

strlist_label <- function(strs) {
  str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
  sapply(str_split(strs, delim), strlist_label)
}

msdata_full <- append_protgroups_info(msdata_full, msdata.wide,
                                 proteins_info = dplyr::bind_rows(
                                   dplyr::mutate(fasta.dfs$CoV, is_viral = TRUE),
                                   dplyr::mutate(fasta.dfs$human, is_viral = FALSE)),
                                 import_columns = "is_viral")
msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
    gene_label = strlist_label2(gene_names),
    protac_label = strlist_label2(protein_acs),
    protgroup_label = case_when(!is.na(gene_label) ~ gene_label,
                                !is.na(protac_label) ~ protac_label,
                                TRUE ~ str_c('#', protgroup_id)))

msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  select(pepmod_id, pepmodstate_id, msrun, intensity = intensity.Sum, starts_with("ident_type")) %>%
  mutate(is_idented = replace_na(ident_type.MSMS, FALSE) | replace_na(`ident_type.MULTI-MATCH-MSMS`, FALSE) |
                      replace_na(`ident_type.MULTI-MSMS`, FALSE))
  
# condition = bait
conditions.df <- dplyr::select(msdata_full$msruns, bait_full_id, bait_id, bait_type, orgcode) %>%
  dplyr::distinct() %>%
  dplyr::arrange(bait_type, bait_id, orgcode) %>%
  dplyr::mutate(condition = factor(bait_full_id, levels=bait_full_id))
msdata_full$msruns <- left_join(msdata_full$msruns, select(conditions.df, bait_full_id, condition))

msdata_full$peptides <- mqevidence$peptides %>%
  dplyr::left_join(select(msdata_full$proteins, lead_razor_protein_ac = protein_ac, is_viral)) %>%
  mutate(is_viral=replace_na(is_viral, FALSE))
msdata_full$pepmods <- mqevidence$pepmods %>%
  dplyr::left_join(select(msdata_full$peptides, peptide_id, is_viral))
msdata_full$pepmodstates <- mqevidence$pepmodstates

msdata_full$pepmodstate_tagintensities <- dplyr::select(
  msdata_full$pepmodstate_intensities, pepmodstate_id, pepmod_id, msrun, is_idented, intensity)

#msdata$protgroup_idents_wide <- select(msdata.wide, protgroup_id, !!msdata_colgroups$ident_type)
#msdata$protgroup_idents <- reshape(msdata_full$protgroup_idents_wide, idvar = 'protgroup_id',
#                                        direction = "long", timevar="msrun_mq", v.names = "ident_type") %>%
#  dplyr::filter(!is.na(ident_type))

#msdata$protgroup_stats <- inner_join(msdata$protgroup_stats,
#    msdata_full$protgroup_intensities %>%
#    dplyr::inner_join(msdata_full$mschannels) %>%
#    dplyr::group_by(protgroup_id) %>%
#    summarize(n = n(),
#              n_quanted = sum(!is.na(intensity.L) | !is.na(intensity.H)),
#              median_intensity = median(intensity, na.rm=TRUE)) %>%
#    dplyr::ungroup())
#
#msdata$protgroups <- dplyr::mutate(msdata$protgroups,
#                                  is_full_quant = protgroup_id %in% (dplyr::filter(msdata$protgroup_stats, n_quanted==n) %>% .$protgroup_id),
#                                  is_top_quant = protgroup_id %in% (dplyr::filter(msdata$protgroup_stats, percent_rank(-median_intensity) <= 0.01) %>% .$protgroup_id))

msdata <- msdata_full[c('protgroup_intensities', 'protgroup_idents',
                        'msruns', 'protgroups', 'protein2protgroup')]#, 'pepmodstate_intensities', 'pepmodstates', 'pepmods')]

# setup experimental design matrices
conditionXeffect_orig.mtx <- model.matrix(
  ~ 1 + bait_id + bait_id:orgcode,
  mutate(conditions.df, orgcode = if_else(bait_type == "sample", orgcode, factor("CVHSA2", levels=levels(orgcode)))))
conditionXeffect.mtx <- conditionXeffect_orig.mtx[, colSums(abs(conditionXeffect_orig.mtx)) != 0 &
                                                    !str_detect(colnames(conditionXeffect_orig.mtx), "bait_idCtrl.+:orgcode")]
dimnames(conditionXeffect.mtx) <- list(condition = conditions.df$condition,
                                       effect = colnames(conditionXeffect.mtx))

cairo_pdf(filename = file.path(data_path, paste0(project_id, "_exp_design_APMS_B1_", fit_version, ".pdf")),
          width = 8, height = 6)
pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  dplyr::mutate(orgcode = effect_factor(effect, "orgcode", levels(conditions.df$orgcode), NA),
                bait_id = effect_factor(effect, "bait_id", levels(conditions.df$bait_id), NA),
                is_positive = !is.na(bait_id) & is.na(orgcode))
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
cairo_pdf(filename = file.path(data_path, paste0(project_id, "_exp_design_inv_APMS_B1_", fit_version, ".pdf")),
          width = 8, height = 6)
pheatmap(inv_conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

msrunXreplEffect.mtx <- replicate_effects_matrix(
  msdata$msruns,
  replicate_col = "replicate", condition_col = "condition")
cairo_pdf(filename = file.path(data_path, paste0(project_id, "_exp_design_msruns_APMS_B1_", fit_version, ".pdf")),
          width = 14, height = 12)
pheatmap(msrunXreplEffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

msrunXreplEffect.df <- as_tibble(as.table(msrunXreplEffect.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)

bait_conditions <- as.character(filter(conditions.df, bait_type == "sample")$condition)
allminus_metaconditions <- paste0("allminus_", unique(filter(conditions.df, bait_type == "sample")$bait_id))
compound_metaconditions <- c(allminus_metaconditions, "controls")
all_metaconditions <- c(bait_conditions, compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in bait_conditions) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
}
for (cname in allminus_metaconditions) {
  bait <- str_remove(cname, "^allminus_")
  conditionXmetacondition.mtx[, cname] <- TRUE
  conditionXmetacondition.mtx[str_detect(rownames(conditionXmetacondition.mtx), str_c("_", bait, "$")), cname] <- FALSE
}
conditionXmetacondition.mtx[filter(conditions.df, bait_type == "control")$bait_full_id, "controls"] <- TRUE
cairo_pdf(filename = file.path(data_path, paste0(project_id, "_metaconditions_APMS_B1_", fit_version, ".pdf")),
          width = 8, height = 6)
pheatmap(ifelse(conditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)

contrasts.df <- bind_rows(
  transmute(filter(conditions.df, bait_type == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = "controls",
            contrast_type = "filter"),
  transmute(filter(conditions.df, bait_type == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = str_c("allminus_", bait_id),
            contrast_type = "filter"),
  inner_join(
    select(filter(conditions.df, bait_type == "sample"), metacondition_lhs = bait_full_id, orgcode_lhs = orgcode, bait_id),
    select(filter(conditions.df, bait_type == "sample"), metacondition_rhs = bait_full_id, orgcode_rhs = orgcode, bait_id),
  ) %>% filter(as.integer(orgcode_lhs) < as.integer(orgcode_rhs)) %>%
  mutate(contrast_type = "comparison") %>%
  select(-orgcode_lhs, -orgcode_rhs, -bait_id)) %>%
  mutate(contrast = str_c(metacondition_lhs, "_vs_", ifelse(str_starts(metacondition_rhs, "allminus"), "others", metacondition_rhs)))

all_contrasts <- contrasts.df$contrast
contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
for (i in 1:nrow(contrasts.df)) {
  contrastXmetacondition.mtx[contrasts.df$contrast[[i]],
                             c(contrasts.df$metacondition_lhs[[i]],
                               contrasts.df$metacondition_rhs[[i]])] <- c(1, -1)
}
cairo_pdf(filename = file.path(data_path, paste0(project_id, "_contrasts_APMS_B1_", fit_version, ".pdf")),
          width = 8, height = 8)
pheatmap(contrastXmetacondition.mtx, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()

contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::rename(weight = n) %>%
  dplyr::left_join(select(contrasts.df, contrast, contrast_type)) %>%
  dplyr::mutate(condition_role = if_else(contrast_type == "filter" & weight < 0, "background", "signal"))

contrastXcondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::select(-n) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition)

# prepare the data to use for MS runs normalization
msrun_stats.df <- left_join(msdata_full$protgroup_intensities,
                            msdata_full$protgroup_idents) %>%
  dplyr::group_by(msrun) %>%
  dplyr::summarise(n_pg_idents = sum(replace_na(ident_type, "") == "By MS/MS"),
                   n_pg_quants = sum(!is.na(intensity))) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(msdata$msruns)

# normalize using the intensities
msdata4norm.df <- msdata$protgroup_intensities %>%
  dplyr::left_join(msdata$protgroup_idents) %>%
  dplyr::filter(replace_na(ident_type, "") == "By MS/MS" & !is.na(intensity)) %>%
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$protgroup_intensities)

options(mc.cores=8)

# normalize experiments:
# 1) MS replicates for a given bait
# 2) same viral protein of different strains
# 3) all baits together
msruns_hnorm <- multilevel_normalize_experiments(instr_calib,
    filter(msdata$msruns, is_used),
    msdata4norm.df,
    quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
    mcmc.iter = 2000L,
    #mcmc.chains = 6,
    verbose=TRUE,
    norm_levels = list(msrun = list(cond_col = "msrun", max_objs=500L, missing_exp.ratio=0.1),
                       bait_full_id = list(cond_col="bait_full_id", max_objs=500L, missing_exp.ratio=0.1),
                       bait_id = list(cond_col="bait_id", max_objs=200L, missing_exp.ratio=0.1)
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

#msdata$msruns <- dplyr::arrange(msdata$msruns, bait_type, bait_id, orgcode) %>%
#  mutate(msrun_ix = row_number())

msdata_full$mschannels <- mutate(msdata_full$msruns,
                                 mstag = "Sum",
                                 mschannel = msrun,
                                 msrun_mq = msrun)
msdata_full$protgroup_tagintensities <- mutate(msdata_full$protgroup_intensities,
  mschannel = msrun,
  msrun_mq = msrun,
  mstag = "Sum")

msdata_full$mschannel_stats <- mschannel_statistics(msdata_full)
set.seed(1232)
msdata_full$protgroup_intensities_all <- expand(msdata_full$protgroup_intensities,
                                                protgroup_id, msrun) %>%
  left_join(msdata_full$protgroup_intensities) %>%
  dplyr::mutate(mstag = "Sum") %>%
  impute_intensities(msdata_full$mschannel_stats) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
  #filter(mstag == "L") %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_imputed_norm),
                                    nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                    dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                    msrun = unique(protgroup_intensities4pca.df$msrun)))

require(FactoMineR)
msrun_intensities_pca <- PCA(protgroup_intensities.mtx, graph = FALSE)
msrun_intensities_pca.df <- as.data.frame(msrun_intensities_pca$svd$V)
colnames(msrun_intensities_pca.df) <- paste0("comp_", 1:ncol(msrun_intensities_pca.df))
msrun_intensities_pca.df <- dplyr::mutate(msrun_intensities_pca.df,
                                          msrun = rownames(msrun_intensities_pca$var$coord)) %>%
    dplyr::inner_join(msdata$msruns)

require(ggrepel)
cairo_pdf(filename = file.path(data_path, paste0(project_id, "_msruns_pca_APMS_B1_", fit_version, ".pdf")),
          width = 14, height = 14)
ggplot(msrun_intensities_pca.df,
       aes(x=comp_1, y=comp_2, color=bait_id)) +
    geom_point() +
    geom_text_repel(aes(label=str_remove(str_remove(msrun, "APMS_SARS_"), "APMS_")), vjust=-1.1) +
    theme_bw_ast(base_family = "", base_size = 10) #+
dev.off()

# no batch effects so far
msrunXbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXreplEffect.mtx),
                                     batch_effect = c())

batch_effects.df <- tibble(batch_effect=character(0),
                           is_positive=logical(0))

rdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', mq_folder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df,
     conditionXeffect.mtx, inv_conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     instr_calib, global_protgroup_labu_shift,
     msruns_hnorm, total_msrun_shifts.df,
     msrunXreplEffect.mtx,
     batch_effects.df, msrunXbatchEffect.mtx,
     file = rdata_filepath)

rdata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', mq_folder, '_', data_version, '.RData'))
message('Saving full MS data to ', rdata_filepath, '...')
save(data_info, msdata_full,
     #protgroup_stats.df,
     file = rdata_filepath)

message('Done.')
