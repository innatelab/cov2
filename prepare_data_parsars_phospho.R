# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20201010"
fit_version <- "20201010"
msfolder <- 'snaut_parsars_phospho_20201005'
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

msdata_path <- file.path(data_path, msfolder, str_c("ptm_extractor_", data_version))

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  msfolder = msfolder,
                  mscalib_protgroup_filename = "instr_QX7_intensity_protgroup_calib_cov2_20200430.json",
                  mscalib_pepmodstate_filename = "mscalib_EXPL2_intensity_pepmodstate_cov2_20200828.json",
                  quant_type = "intensity", quant_col_prefix = "intensity",
                  qvalue_max=1E-2, qvalue_ident_max=1E-3,
                  locprob_min=0.5, locprob_ident_min=0.75,
                  empty_observation_sigmoid_scale=1/3,
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$mscalib_protgroup_filename, '...')
mscalib_protgroup <- fromJSON(file = file.path(data_path, data_info$mscalib_protgroup_filename))$instr_calib
mscalib_pepmodstate_json <- fromJSON(file = file.path(data_path, data_info$mscalib_pepmodstate_filename))
mscalib_pepmodstate <- mscalib_pepmodstate_json$mscalib

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/SARS_CoV_20200928.fasta")),
  CoV2 = read_innate_uniprot_fasta(file.path(data_path, "msfasta/SARS_CoV2_20200928.fasta")),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/Human_July2019_with_isoforms_only_swissprot.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants_20200405.fasta"))
)

bad_msruns <- c()

msdata_full <- lapply(list(
  msruns = "rawfiles_info.txt",
  proteins = "proteins.txt.gz",
  ptm_locprobs = "ptm_locprobs.txt.gz",
  pepmodstate_intensities = "pms_intensities.txt.gz",
  protgroups = "protgroups.txt.gz",
  peptide2protgroup = "peptide_to_protgroup.txt.gz",
  peptides = "peptides.txt.gz",
  pepmodstates = "pepmodstates.txt.gz",
  protein2protgroup = "protein_to_protgroup.txt.gz",
  peptide2protein = "peptide_to_protein.txt.gz",
  ptm2protein = "ptm_to_protein.txt.gz",
  ptm2gene = "ptm_to_gene.txt.gz",
  ptmns = "ptmns.txt.gz",
  ptmn2pepmodstate = "ptmn_to_pepmodstate.txt.gz"),
  function(fname) {
    message("Reading ", fname, "...")
    read_tsv(file.path(msdata_path, fname))
  }
)

msdata_full$msruns <- msdata_full$msruns %>%
    dplyr::mutate(dataset = factor(dataset, c("phospho", "ubi")),
                  treatment = factor(treatment, c("mock", "SARS_CoV", "SARS_CoV2")),
                  msrun_ix = as.integer(msrun_ix),
                  rawfile_ix = as.integer(rawfile_ix),
                  timepoint_num = as.integer(as.character(timepoint)),
                  timepoint = factor(timepoint),
                  replicate = as.integer(replicate),
                  condition = str_c(treatment, "_", timepoint, "h"),
                  msrun = str_c(dataset, "_", condition, "_", replicate),
                  is_used = !msrun %in% bad_msruns) %>%
  dplyr::mutate(msrun = factor(msrun, levels=msrun),
                condition = factor(condition, levels=unique(condition)))
msdata_full$pepmodstate_intensities <- dplyr::inner_join(msdata_full$pepmodstate_intensities,
                                                         dplyr::select(msdata_full$msruns, dataset, msrun, msrun_ix, rawfile_ix)) %>%
  dplyr::select(-rawfile_ix) %>%
  dplyr::mutate(ident_type = factor(ifelse(coalesce(qvalue, 1) <= data_info$qvalue_max, "MULTI-MSMS", "MULTI-MATCH"),
                                    c("MULTI-MATCH", "MULTI-MSMS")))
msdata_full$ptm_locprobs <- dplyr::inner_join(msdata_full$ptm_locprobs,
                                              dplyr::select(msdata_full$msruns, dataset, msrun, msrun_ix, rawfile_ix)) %>%
  dplyr::select(-rawfile_ix)

# add more properties to ptmns as it's the obect being modeled
msdata_full$ptmns <- msdata_full$ptmns %>%
  dplyr::left_join(dplyr::select(dplyr::filter(msdata_full$ptm2gene, ptm_is_reference), ptm_id, ptm_pos, protein_ac, is_viral, is_contaminant)) %>%
  dplyr::left_join(dplyr::select(msdata_full$proteins, gene_name=genename, protein_ac, protein_code)) %>%
  dplyr::mutate(ptmn_label_no_ptm_type = str_remove(ptmn_label, "^[^_]+_"))

msdata_full$ptmn_locprobs <- dplyr::inner_join(msdata_full$ptm2gene, msdata_full$ptm_locprobs) %>%
    dplyr::select(ptm_id, pepmodstate_id, dataset, msrun_ix, msrun, ptm_locprob) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(msdata_full$ptmn2pepmodstate)

msdata_full$ptmn_stats <- inner_join(msdata_full$ptmns, msdata_full$ptmn2pepmodstate) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmodstate_intensities, pepmodstate_id, dataset, msrun, intensity, qvalue)) %>%
  dplyr::inner_join(msdata_full$ptmn_locprobs) %>%
  dplyr::mutate(is_quanted = !is.na(intensity),
                is_idented = coalesce(qvalue, 1) <= data_info$qvalue_ident_max,
                is_localized = coalesce(ptm_locprob, 0) >= data_info$locprob_ident_min,
                is_valid_quant = is_quanted & coalesce(ptm_locprob, 0) >= data_info$locprob_min &
                                 coalesce(qvalue, 1) <= data_info$qvalue_max) %>%
  dplyr::group_by(ptmn_id, ptmn_label, ptm_id, ptm_type, dataset) %>%
  dplyr::summarise(n_quanted = sum(is_quanted),
                   n_idented = sum(is_idented),
                   n_msruns = n_distinct(msrun),
                   n_pepmodstates = n_distinct(pepmodstate_id),
                   n_localized = sum(is_localized),
                   n_idented_and_localized = sum(is_localized & is_idented),
                   n_valid_quants = sum(is_valid_quant),
                   pms_qvalue_min = min(qvalue, na.rm=TRUE),
                   ptm_locprob_max = max(ptm_locprob, na.rm=TRUE),
                   .groups="drop")

msdata_full$msrun_pepmodstate_stats <- msrun_statistics(msdata_full, obj = "pepmodstate") %>%
  dplyr::mutate(na_ratio = n_missing / n)

all_proteins.df <- dplyr::bind_rows(
    dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE, is_expected_stable=FALSE),
    dplyr::mutate(fasta.dfs$CoV2, is_viral=TRUE, is_contaminant=FALSE, is_expected_stable=FALSE),
    dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE,
                  # define stable proteins based on our expectations of their biology and how their timeseries look like
                  is_expected_stable = str_detect(gene_name, "^(M?RP[LS]\\d+|GAPDH|ACT[ABNR]\\d+|TUB[AB]\\d?\\w?|CCT\\d?\\w?)$")),
    dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE, is_expected_stable=FALSE))
#msdata_full <- append_protgroups_info(msdata_full, pgdata.wide,
#                                      proteins_info = all_proteins.df,
#                                      import_columns = c("is_viral", "is_contaminant", "organism"))
msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
    is_reverse = FALSE,
    gene_label = strlist_label2(gene_names),
    protac_label = strlist_label2(majority_protein_acs),
    protein_label = strlist_label2(protein_names),
    protgroup_label = case_when(#is_viral ~ protein_label,
                                !is.na(gene_label) ~ gene_label,
                                !is.na(protac_label) ~ protac_label,
                                TRUE ~ str_c('#', protgroup_id)))

# condition = treatment X timepoint
conditions.df <- dplyr::select(msdata_full$msruns, condition, treatment, timepoint, timepoint_num) %>%
  dplyr::distinct()
conditions.df <- dplyr::left_join(conditions.df,
    dplyr::select(conditions.df, treatment, timepoint_after=timepoint_num)) %>% # exclude 3h because it's the reference
    dplyr::mutate(is_after = timepoint_num >= timepoint_after) %>%
    tidyr::pivot_wider(all_of(c("condition", "treatment", "timepoint", "timepoint_num")),
                       names_prefix = "after", names_from = timepoint_after, values_from = is_after) %>%
    rename_at(vars(starts_with("after")), ~str_c(., "h")) %>%
    mutate(after0h = NULL, # not needed for exp_design
           infected = treatment != "mock")

msdata <- msdata_full[c('pepmodstate_intensities', 'pepmods', 'pepmodstates',
                        'msruns', 'msrun_pepmodstate_stats', 'proteins',
                        'peptide2protein', 'ptmns', 'ptm2gene', 'ptmn2pepmodstate')]
# keep only Phospho & GlyGly
msdata$ptm2gene <- dplyr::filter(msdata_full$ptm2gene, ptm_type %in% c("Phospho", "GlyGly"))
msdata$ptmns <- dplyr::semi_join(dplyr::filter(dplyr::semi_join(msdata_full$ptmns, msdata$ptm2gene), nselptms <= 3),
                                 dplyr::filter(msdata_full$ptmn_stats, n_valid_quants > 0 & n_idented_and_localized > 0))
msdata$ptm2gene <- dplyr::semi_join(msdata_full$ptm2gene, msdata$ptmns)
msdata$ptmn2pepmodstate <- dplyr::semi_join(msdata_full$ptmn2pepmodstate, msdata$ptmns)
msdata$protein2ptmn <- dplyr::inner_join(msdata$ptmns, msdata$ptm2gene)
msdata$ptmn_locprobs <- semi_join(msdata_full$ptmn_locprobs, msdata$ptmn2pepmodstate)
msdata$pepmodstate_intensities <- msdata_full$pepmodstate_intensities %>% filter(coalesce(qvalue, 1.0) <= data_info$qvalue_max) %>%
                                  semi_join(dplyr::select(filter(msdata$msruns, is_used), msrun)) %>%
                                  semi_join(msdata$pepmodstates) %>%
                                  semi_join(filter(msdata$ptmn_locprobs, coalesce(ptm_locprob, 0) >= data_info$locprob_min))
# setup experimental design matrices
conditionXeffect_orig.mtx <- model.matrix(
  formula(str_c("~ 1 + (", str_c(str_subset(colnames(conditions.df), "^after\\d+h"), collapse =" + "), ") * (infected + treatment)")),
  conditions.df
)
colnames(conditionXeffect_orig.mtx) <- str_replace(colnames(conditionXeffect_orig.mtx), "treatment", "strain")
cov2_effects <- str_subset(colnames(conditionXeffect_orig.mtx), "CoV2$")
conditionXeffect.mtx <- conditionXeffect_orig.mtx
conditionXeffect.mtx[, cov2_effects] <-
    conditionXeffect.mtx[, cov2_effects] - conditionXeffect.mtx[, str_remove(cov2_effects, "2$")]
conditionXeffect.mtx <- conditionXeffect.mtx[, str_detect(colnames(conditionXeffect.mtx),
                                                          "(^strainSARS_CoV2?|^infectedTRUE|CoV)$", negate=TRUE) &
                                             (apply(conditionXeffect.mtx, 2, function(x) min(abs(x))) == 0)] # remove intercept
dimnames(conditionXeffect.mtx) <- list(condition = as.character(conditions.df$condition),
                                       effect = colnames(conditionXeffect.mtx))

plots_path <- file.path(analysis_path, "plots", str_c(msfolder, "_", fit_version))
if (!dir.exists(plots_path)) dir.create(plots_path)

pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(analysis_path, "plots", str_c(msfolder, "_", fit_version),
                              paste0(project_id, "_exp_design_", msfolder, "_", fit_version, ".pdf")),
         width = 6, height = 6)

effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  mutate(infected = effect_factor(effect, 'infected',
                                  c(FALSE, TRUE), NA),
         strain = effect_factor(effect, 'strain',
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
  mutate(effect_type = case_when(!is.na(timepoint) & !is.na(strain) ~ "strainXtimepoint",
                                 !is.na(timepoint) & !is.na(infected) ~ "infectedXtimepoint",
                                 !is.na(timepoint) ~ "timepoint",
                                 TRUE ~ NA_character_),
         effect_label = case_when(effect_type == "strainXtimepoint" ~ str_c(strain, "+", timepoint, "h"),
                                  effect_type == "infectedXtimepoint" ~ str_c("infected+", timepoint, "h"),
                                  effect_type == "timepoint" ~ str_c(timepoint, "h"),
                                  TRUE ~ NA_character_),
         prior_mean = 0,
         prior_tau = case_when(effect_type == "strainXtimepoint" ~ 0.25,
                               effect_type == "infectedXtimepoint" ~ 0.25,
                               effect_type == "timepoint" ~ 0.5,
                               TRUE ~ NA_real_),
         prior_df2 = case_when(effect_type == "strainXtimepoint" ~ 4.0,
                               effect_type == "infectedXtimepoint" ~ 4.0,
                               TRUE ~ 2.0),
         is_positive = FALSE)

conditionXeffect.df <- conditionXeffect_frame(conditionXeffect.mtx, effects.df)

#TODO 0h?
compound_metaconditions <- str_c("infected_", unique(filter(conditions.df, infected)$timepoint), "h")
all_metaconditions <- c(levels(conditions.df$condition), compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in levels(conditions.df$condition)) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
  if (str_detect(cname, "CoV")) {
    conditionXmetacondition.mtx[cname, str_replace(cname, "SARS_CoV2?", "infected")] <- TRUE
  }
}
pheatmap(ifelse(conditionXmetacondition.mtx, 1L, 0L), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plots_path, paste0(project_id, "_metaconditions_", msfolder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)
metaconditions.df <- bind_rows(mutate(conditions.df, metacondition=as.character(condition)),
                               tibble(metacondition=compound_metaconditions) %>%
                               tidyr::extract("metacondition", c("treatment", "timepoint"), "(.+)_(\\d+)h", remove=FALSE) %>%
                               mutate(infected=TRUE,
                                      timepoint_num=parse_integer(timepoint))) %>%
  mutate(treatment = factor(treatment, levels=c(levels(conditions.df$treatment), "infected")))

contrasts.df <- bind_rows(
  # all treatment pairs at each timepoint
  left_join(dplyr::mutate(dplyr::select(metaconditions.df, metacondition_lhs = metacondition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(metaconditions.df, metacondition_rhs = metacondition, treatment_rhs = treatment, timepoint_lhs = timepoint))) %>%
  filter((as.integer(treatment_lhs) > as.integer(treatment_rhs)) &
         !((treatment_lhs == "infected") & (treatment_rhs != "mock"))) %>%
  mutate(timepoint_rhs = timepoint_lhs,
         contrast_kind = "treatment_vs_treatment"),
  # all timepoints of the same treatment
  left_join(dplyr::mutate(dplyr::select(metaconditions.df, metacondition_lhs = metacondition, treatment_lhs = treatment, timepoint_lhs = timepoint)),
            dplyr::mutate(dplyr::select(metaconditions.df, metacondition_rhs = metacondition, treatment_lhs = treatment, timepoint_rhs = timepoint))) %>%
  filter(as.integer(timepoint_lhs) > as.integer(timepoint_rhs)) %>%
  mutate(treatment_rhs = treatment_lhs,
         contrast_kind = "timepoint_vs_timepoint")
) %>%
  mutate(contrast = str_c(treatment_lhs, "@", timepoint_lhs, "h_vs_", treatment_rhs, "@", timepoint_rhs, "h"),
         contrast_type = "comparison",
         offset = 0.0)
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
         filename = file.path(plots_path, paste0(project_id, "_exp_design_contrasts_", msfolder, "_", fit_version, ".pdf")),
         width = 8, height = 12)

##########
## normalization
msdata4norm.df <- msdata_full$pepmodstate_intensities %>% ungroup() %>%
  dplyr::filter(!is.na(intensity_norm) & qvalue <= 1E-4) %>% 
  dplyr::select(pepmodstate_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(dplyr::filter(msdata_full$pepmodstates, !coalesce(is_decoy, FALSE) & str_detect(pepmodstate_seq, "Phospho|GlyGly")) %>% # select phospho or ubi pepmods
                    dplyr::select(pepmodstate_id, pepmodstate_seq, peptide_id)) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$proteins, !is_viral & !is_contaminant) %>%
                   dplyr::inner_join(msdata_full$peptide2protein)) %>%
  # select some representative pepmodstates from each intensity range to avoid bias by the most abundant species
  dplyr::inner_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::group_by(pepmodstate_id) %>%
  dplyr::summarise(intensity_med = median(intensity),
                   normshift_abs_max = max(abs(log2(normfactor)), na.rm=TRUE),
                   .groups="drop") %>%
  dplyr::filter(normshift_abs_max <= 2) %>%
  dplyr::mutate(intensity_med_log2 = log2(intensity_med),
                intensity_bin = as.integer(cut(intensity_med_log2, 50))) %>%
  dplyr::group_by(intensity_bin) %>%
  dplyr::slice_sample(n = 200) %>%
  dplyr::inner_join(msdata_full$pepmodstate_intensities)

options(mc.cores=8)

msruns_hnorm <- multilevel_normalize_experiments(mscalib_pepmodstate,
                                                 dplyr::mutate(msdata$msruns, dataset_condition = str_c(dataset, "_", condition),
                                                               msdata$msruns, dataset_timepoint = str_c(dataset, "_", timepoint)),
                                                 msdata4norm.df,
                                                 quant_col = "intensity_norm",
                                                 #quant_col = "intensity",
                                                 obj_col = "pepmodstate_id", mschan_col = "msrun",
                                                 mcmc.iter = 2000L,
                                                 #mcmc.chains = 6,
                                                 verbose=TRUE,
                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=1000L, missing_exp.ratio=0.1),
                                                                    condition = list(cond_col="dataset_condition", max_objs=500L, missing_exp.ratio=0.2),
                                                                    timepoint = list(cond_col="dataset_timepoint", condgroup_col="dataset", max_objs=500L, missing_exp.ratio=0.3)
                                                 ))
#saveRDS(msruns_hnorm, file.path(scratch_path, str_c(project_id, "_", msfolder, "_", fit_version, ".rds")))
#msruns_hnorm <- readRDS(file.path(scratch_path, str_c(project_id, "_", msfolder, "_20200907.rds")))
total_msrun_shifts.df <- msruns_hnorm$mschannel_shifts

# compare normalizations
norm_plot <- ggplot(inner_join(msdata_full$pepmodstate_intensities, msdata_full$msruns), aes(x=msrun)) +
    geom_violin(aes(y = -log2(normfactor)), draw_quantiles=c(0.25, 0.75), scale="width", fill="lightgray") +
    geom_point(data=inner_join(total_msrun_shifts.df, msdata_full$msruns),
               aes(x = msrun, y=total_msrun_shift/log(2)), size=1.5, color="red") +
    facet_wrap(~ treatment + dataset, ncol=2, scales = "free") +
    theme_bw_ast(base_family = "") +
    theme(axis.text.x = element_text(hjust = 0, angle=-45))
ggsave(norm_plot, filename = file.path(plots_path,
                             paste0(project_id, "_", msfolder, "_normalization_", fit_version, ".pdf")),
       width = 8, height = 12)

## apply normalization
#msdata$protgroup_intensities <- dplyr::left_join(msdata$protgroup_intensities,
#                                                 dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
#  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) 
#dplyr::select(-total_msrun_shift)

global_pepmodstate_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns, TRUE) %>%
                                             dplyr::select(msrun) %>% dplyr::distinct() %>%
                                             dplyr::inner_join(msdata$pepmodstate_intensities) %>% .$intensity), na.rm=TRUE)

pepmodstate_labu_min <- inner_join(msdata$pepmodstate_intensities, total_msrun_shifts.df) %>%
  mutate(intensity_norm = intensity * exp(-total_msrun_shift)) %>%
  .$intensity_norm %>% log() %>%
  quantile(0.001, na.rm=TRUE) - global_pepmodstate_labu_shift - 2

############################
# batch effects
# no batch effects so far
msrunXbatchEffect.mtx <- zero_matrix(msrun = msdata$msruns$msrun,
                                     batch_effect = c())

batch_effects.df <- tibble(batch_effect = character(0),
                           is_positive = logical(0),
                           prior_mean = double(0))

# no subbatch effects so far
msrunXsubbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXbatchEffect.mtx),
                                        subbatch_effect = c())

subbatch_effects.df <- tibble(subbatch_effect=character(0),
                              is_positive=logical(0))

# estimate the effect scales per each msrun
pepmodstate_intensities.df <- msdata_full$pepmodstate_intensities %>%
  dplyr::inner_join(total_msrun_shifts.df) %>%
  dplyr::mutate(intensity_norm = intensity_norm*exp(-total_msrun_shift),
                intensity_norm_log2 = log2(intensity_norm)) %>%
  inner_join(msdata_full$msruns) %>%
  dplyr::group_by(pepmodstate_id) %>%
  dplyr::mutate(intensity_median = median(intensity_norm, na.rm=TRUE),
                intensity_mock0h = median(intensity_norm[condition == "mock_0h"], na.rm=TRUE),
                .groups="drop") %>%
  dplyr::group_by(pepmodstate_id, timepoint) %>%
  dplyr::mutate(intensity_mock = median(intensity_norm[treatment == "mock"], na.rm=TRUE),
                intensity_infected = median(intensity_norm[treatment != "mock"], na.rm=TRUE),
                .groups="drop") %>%
  dplyr::group_by(pepmodstate_id, condition) %>%
  dplyr::mutate(intensity_cond = median(intensity_norm, na.rm=TRUE),
                .groups="drop") %>%
  dplyr::mutate(intensity_cond_log2 = log2(intensity_cond),
                intensity_infected_log2 = log2(intensity_infected),
                intensity_mock_log2 = log2(intensity_mock),
                intensity_mock0h_log2 = log2(intensity_mock0h),
                intensity_median_log2 = log2(intensity_median),
                intensity_cond_vs_mock_log2 = intensity_cond_log2 - intensity_mock_log2,
                intensity_cond_vs_mock0h_log2 = intensity_cond_log2 - intensity_mock0h_log2,
                intensity_infected_vs_mock_log2 = intensity_infected_log2 - intensity_mock_log2,
                intensity_infected_vs_mock0h_log2 = intensity_infected_log2 - intensity_mock0h_log2,
                intensity_vs_mock_log2 = intensity_norm_log2 - intensity_mock_log2,
                intensity_vs_mock0h_log2 = intensity_norm_log2 - intensity_mock0h_log2)

require(broom)

pepmodstate_intensities4glm.df <- pepmodstate_intensities.df %>%
  dplyr::semi_join(dplyr::filter(msdata_full$pepmodstates, nselptms == 1)) %>%
  dplyr::filter(between(intensity_vs_mock_log2, -2, 2) &
                  between(intensity_infected_vs_mock_log2, -2, 2)) %>%
  dplyr::filter(treatment != "mock") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(intensity_infected_vs_mock_bin = as.integer(cut(intensity_infected_vs_mock_log2, 50))) %>%
  dplyr::group_by(condition, intensity_infected_vs_mock_bin) %>%
  dplyr::slice_sample(n = 500) %>%
  dplyr::ungroup()

infection_scale.glm <- glm(intensity_vs_mock_log2 ~ 0 + intensity_infected_vs_mock_log2:msrun,
                           data=pepmodstate_intensities4glm.df)

infection_scales.df <- broom::tidy(infection_scale.glm) %>%
  tidyr::extract(term, c("msrun"), ":msrun(.+)$", remove=FALSE) %>%
  dplyr::inner_join(msdata_full$msruns) %>%
  dplyr::mutate(msrun_scale = estimate,
                msrun_scale_log2 = log2(estimate)) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(msrun_scale_log2_median = median(msrun_scale_log2)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(msrun_scale_norm = 2^(msrun_scale_log2 - msrun_scale_log2_median))

infection_scales_plot <- ggplot(infection_scales.df %>% dplyr::mutate(msrun_short = str_remove(msrun, "phospho_SARS_"))) +
  geom_bar(aes(x = msrun_short, y = msrun_scale_norm, fill=treatment), stat="identity") +
  scale_fill_manual(values=treatment_palette) +
  scale_y_log10() +
  facet_wrap(~ timepoint, scales="free_x", ncol=2) +
  theme_bw_ast() +
  theme(axis.text.x = element_text(angle = -45, hjust=0))
ggsave(infection_scales_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                               str_c(project_id, "_", msfolder, '_', fit_version, "_infection_scales.pdf")),
       width=8, height=8, device = cairo_pdf)

pepmodstate_intensities4mock_glm.df <- pepmodstate_intensities.df %>%
  dplyr::semi_join(dplyr::filter(msdata_full$pepmodstates, nselptms == 1)) %>%
  dplyr::filter(between(intensity_vs_mock0h_log2, -2, 2) &
                  between(intensity_cond_vs_mock0h_log2, -2, 2)) %>%
  dplyr::filter(treatment == "mock") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(intensity_infected_vs_mock_bin = as.integer(cut(intensity_cond_vs_mock0h_log2, 50))) %>%
  dplyr::group_by(condition, intensity_infected_vs_mock_bin) %>%
  dplyr::slice_sample(n = 500) %>%
  dplyr::ungroup()

mock_scale.glm <- glm(intensity_vs_mock0h_log2 ~ 0 + intensity_cond_vs_mock0h_log2:msrun,
                      data=pepmodstate_intensities4mock_glm.df)

mock_scales.df <- broom::tidy(mock_scale.glm) %>%
  tidyr::extract(term, c("msrun"), ":msrun(.+)$", remove=FALSE) %>%
  dplyr::inner_join(msdata_full$msruns) %>%
  dplyr::mutate(msrun_scale = estimate,
                msrun_scale_log2 = log2(estimate)) %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(msrun_scale_log2_median = median(msrun_scale_log2)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(msrun_scale_norm = 2^(msrun_scale_log2 - msrun_scale_log2_median))

mock_scales_plot <- ggplot(mock_scales.df %>% dplyr::mutate(msrun_short = str_remove(msrun, "phospho_")) %>% dplyr::filter(condition != "mock_0h")) +
  geom_bar(aes(x = msrun_short, y = msrun_scale_norm, fill=treatment), stat="identity") +
  scale_fill_manual(values=treatment_palette) +
  scale_y_log10() +
  facet_wrap(~ timepoint, scales="free_x", ncol=2) +
  theme_bw_ast() +
  theme(axis.text.x = element_text(angle = -45, hjust=0))
mock_scales_plot
ggsave(mock_scales_plot, file = file.path(analysis_path, "plots", str_c(msfolder, '_', fit_version),
                                          str_c(project_id, "_", msfolder, '_', fit_version, "_mock_scales.pdf")),
       width=8, height=8, device = cairo_pdf)

msrunXeffect.df <- inner_join(conditionXeffect.df, msdata$msruns) %>%
  dplyr::left_join(bind_rows(dplyr::select(mock_scales.df, msrun, msrun_scale=msrun_scale_norm),
                             dplyr::select(infection_scales.df, msrun, msrun_scale=msrun_scale_norm))) %>%
  dplyr::mutate(scaled_mult = mult * msrun_scale)
msrunXeffect_wide.df <- pivot_wider(bind_rows(msrunXeffect.df, tibble(msrun = setdiff(msdata_full$msruns$msrun, msrunXeffect.df$msrun),
                                                                      effect="tmp_intercept", scaled_mult = 0)),
                                    msrun, values_from = scaled_mult, names_from = effect, values_fill = 0) %>%
  dplyr::select(-tmp_intercept)
msrunXeffect.mtx <- as.matrix(dplyr::select(msrunXeffect_wide.df, -msrun))
rownames(msrunXeffect.mtx) <- msrunXeffect_wide.df$msrun
msrunXeffect.mtx <- msrunXeffect.mtx[msdata_full$msruns$msrun, ]

pheatmap(msrunXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plots_path, paste0(project_id, "_msrun_exp_design_", msfolder, "_", fit_version, ".pdf")),
         width = 8, height = 12)

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', msfolder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df, contrasts.df,
     conditionXeffect.mtx, conditionXeffect.df,
     msrunXeffect.mtx, msrunXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     mscalib_pepmodstate,
     global_pepmodstate_labu_shift,
     pepmodstate_labu_min,
     total_msrun_shifts.df, #msruns_hnorm, 
     batch_effects.df, msrunXbatchEffect.mtx,
     subbatch_effects.df, msrunXsubbatchEffect.mtx,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', msfolder, '_', data_version, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full,
     file = rfulldata_filepath)

message('Done.')

set.seed(1232)
msdata_full$protgroup_intensities_all <- tidyr::expand(msdata_full$protgroup_intensities,
                                                protgroup_id, msrun) %>%
  left_join(dplyr::select(msdata_full$protgroup_intensities, -any_of("total_msrun_shift"))) %>%
  impute_intensities(msdata_full$msrun_stats) %>%
  #dplyr::mutate(total_msrun_shift=0) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift), # FIXME -shift
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
  dplyr::semi_join(dplyr::select(dplyr::filter(msdata_full$protgroups, !is_reverse & !is_contaminant), protgroup_id)) %>%
  dplyr::semi_join(dplyr::select(dplyr::filter(msdata_full$msruns, is_used), msrun)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_norm),
                                    nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                    dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                    msrun = unique(protgroup_intensities4pca.df$msrun)))
protgroup_intensities_imp.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_imputed_norm),
                                        nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                        dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                        msrun = unique(protgroup_intensities4pca.df$msrun)))

# PCA of msruns
msrun_intensities_pca <- stats::prcomp(protgroup_intensities_imp.mtx, scale. = TRUE)
msrun_intensities_pca.df <- as_tibble(msrun_intensities_pca$rotation,
                                      rownames="msrun") %>%
  dplyr::inner_join(msdata$msruns)

treatment_palette <- c(mock="gray", SARS_CoV2 = "goldenrod", SARS_CoV = "brown")
require(ggrepel)
p <- ggplot(msrun_intensities_pca.df,
            aes(x=PC1, y=PC2,
                color=treatment
            )) +
  geom_point() +
  geom_text_repel(aes(label=str_remove(msrun, "SARS_")), vjust=-1.1,
                  box.padding = 0.1, force = 0.5) +
  theme_bw_ast(base_family = "", base_size = 6) +
  scale_color_manual(values=treatment_palette)
ggsave(p, filename = file.path(analysis_path, 'plots', str_c(msfolder, "_", fit_version),
                               paste0(project_id, "_msruns_pca_", msfolder, "_", fit_version, ".pdf")),
       width = 10, height = 10, device = cairo_pdf)

require(pheatmap)
msrun_ordered <- as.character(filter(msdata_full$msruns, !is.na(rawfile_ix))$msrun)
protgroup_hclu = hclust(dist(protgroup_intensities_imp.mtx))
pheatmap(protgroup_intensities.mtx[, msruns_ordered], cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", str_c(msfolder, "_", data_version),
                          paste0(project_id, "_", msfolder, "_", data_version, "_heatmap_intensity.png")), width=16, height=80)

viral_protgroup_stats.df <- protgroup_intensities4pca.df %>%
  dplyr::mutate(row_ix = match(protgroup_id, unique(protgroup_id))) %>%
  dplyr::inner_join(dplyr::select(msdata_full$protgroups, protgroup_label, protgroup_id, organism, is_viral)) %>%
  dplyr::filter(is_viral) %>%
  group_by(organism, is_viral, protgroup_id, protgroup_label, row_ix) %>%
  dplyr::summarise(n_quants = sum(!is.na(intensity_norm)),
                   log_intensity_median = median(log(intensity_norm), na.rm=TRUE),
                   .groups="drop") %>%
  dplyr::arrange(organism, desc(log_intensity_median), desc(n_quants))

viral_intensities.mtx <- protgroup_intensities.mtx[viral_protgroup_stats.df$row_ix, msruns_ordered]
rownames(viral_intensities.mtx) <- viral_protgroup_stats.df$protgroup_label
pheatmap(viral_intensities.mtx, cluster_cols=FALSE, cluster_rows=FALSE,
         file = file.path(analysis_path, "plots", str_c(msfolder, "_", data_version),
                          paste0(project_id, "_", msfolder, "_", data_version, "_viral_heatmap_intensity.pdf")), width=12, height=8)

msdata_full$pepmodstate_intensities_all <- tidyr::expand(msdata_full$pepmodstate_intensities,
                                                         pepmodstate_id, msrun) %>%
  left_join(dplyr::select(msdata_full$pepmodstate_intensities, pepmodstate_id, msrun, intensity, intensity_norm, qvalue)) %>%
  impute_intensities(msdata_full$msrun_stats) %>%
  #dplyr::mutate(total_msrun_shift=0) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift), # FIXME -shift
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::left_join(dplyr::filter(dplyr::select(msdata_full$pepmodstates, pepmod_id, pepmodstate_id, is_pg_specific), is_pg_specific)) %>%
  dplyr::left_join(dplyr::select(msdata_full$pepmods, pepmod_id, peptide_id)) %>%
  dplyr::left_join(msdata_full$peptide2protgroup) %>%
  dplyr::arrange(msrun, protgroup_id, pepmodstate_id)

pepmodstate_intensities4pca.df <- msdata_full$pepmodstate_intensities_all %>%
  dplyr::semi_join(dplyr::select(dplyr::filter(msdata_full$protgroups, !is_reverse & !is_contaminant), protgroup_id)) %>%
  dplyr::semi_join(dplyr::select(dplyr::filter(msdata_full$msruns, is_used), msrun)) %>%
  dplyr::arrange(msrun, protgroup_id, pepmodstate_id)

pepmodstate_intensities.mtx <- matrix(log2(if_else(pepmodstate_intensities4pca.df$qvalue <= 1E+4,
                                                   pepmodstate_intensities4pca.df$intensity_norm,
                                                   NA_real_)),
                                      nrow = n_distinct(pepmodstate_intensities4pca.df$pepmodstate_id),
                                      dimnames = list(pepmodstate = unique(pepmodstate_intensities4pca.df$pepmodstate_id),
                                                      msrun = unique(pepmodstate_intensities4pca.df$msrun)))
pepmodstate_intensities_imp.mtx <- matrix(log2(pepmodstate_intensities4pca.df$intensity_imputed_norm),
                                        nrow = n_distinct(pepmodstate_intensities4pca.df$pepmodstate_id),
                                        dimnames = list(pepmodsate = unique(pepmodstate_intensities4pca.df$pepmodstate_id),
                                                        msrun = unique(pepmodstate_intensities4pca.df$msrun)))
msruns_hclu = hclust(dist(t(pepmodstate_intensities_imp.mtx)))

viral_pepmodstate_stats.df <- pepmodstate_intensities4pca.df %>%
  dplyr::mutate(row_ix = match(pepmodstate_id, unique(pepmodstate_id))) %>%
  dplyr::semi_join(dplyr::select(msdata_full$protgroups, protgroup_id, is_viral) %>% dplyr::filter(is_viral)) %>%
  group_by(protgroup_id, pepmod_id, pepmodstate_id, row_ix) %>%
  dplyr::summarise(n_quants = sum(!is.na(intensity_norm)),
                   log_intensity_median = median(log(intensity_norm), na.rm=TRUE),
                   .groups="drop") %>%
  dplyr::inner_join(dplyr::select(msdata_full$protgroups, protgroup_id, protgroup_label, organism, is_viral)) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, charge)) %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmods, pepmod_id, peptide_seq)) %>%
  mutate(row_label = str_c(protgroup_label, ": ", peptide_seq, ".", charge)) %>%
  dplyr::arrange(organism, protgroup_id, desc(log_intensity_median), desc(n_quants))

msrunix_ordered <- as.character(filter(msdata_full$msruns, !is.na(rawfile_ix)) %>% dplyr::arrange(msrun_ix) %>% .$msrun)
viral_intensities.mtx <- pepmodstate_intensities.mtx[viral_pepmodstate_stats.df$row_ix, msrun_ordered]
rownames(viral_intensities.mtx) <- viral_pepmodstate_stats.df$row_label
pheatmap(viral_intensities.mtx, cluster_cols=FALSE, cluster_rows=FALSE,
         file = file.path(analysis_path, "plots", str_c(msfolder, "_", data_version),
                          str_c(project_id, "_", msfolder, "_", data_version, "_viral_pepmod_heatmap_intensity.pdf")), width=14, height=60)

require(uwot)
require(matrixStats)
protgroup_intensities_imp_pg.mtx <- protgroup_intensities_imp.mtx - rowMedians(protgroup_intensities_imp.mtx[, str_detect(colnames(protgroup_intensities_imp.mtx), "mock_0h_")], na.rm = TRUE)
pg_umap2d = umap(protgroup_intensities_imp_pg.mtx, n_components=2,
                 n_neighbors=20, init="laplacian",
                 min_dist=0.2, metric = "euclidean")

objects_umap.df <- as_tibble(pg_umap2d, .name_repair="minimal") %>%
  set_names(c("x_2d", "y_2d")) %>%
  mutate(object_id = unique(protgroup_intensities4pca.df$protgroup_id)) %>%
  #dplyr::left_join(
  #  as_tibble(objects.umap3d, .name_repair="minimal") %>%
  #    set_names(c("x_3d", "y_3d", "z_3d")) %>%
  #    mutate(object_id = unique(object_intensities4pca.df$object_id))    
  #) %>%
  dplyr::inner_join(select(msdata$protgroups, object_id=protgroup_id, object_label=protgroup_label,
                           majority_protein_acs, gene_names, protein_descriptions,
                           is_viral, is_contaminant, is_reverse)) %>%
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
  # FIXME duplicate point due to N-to-N protgroup to protregroup matching
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

treatmentXprotein_palette <- c("SARS_CoV2" = "chocolate", "SARS_COV2_NSP8" = "darkorange", "SARS_COV2_N" = "orange",
                               "mock_S" = "gray", "mock_NSP8" = "gray", "mock_N" = "gray")

shown_msdata.df <- filter(msdata_full$protgroups, is_viral) %>% #str_detect(gene_names, "(^|;)(SERPIN)")) %>%
    inner_join(msdata_full$protgroup_intensities_all) %>%
    left_join(msdata_full$msruns) %>%
    mutate(treatmentXprotein = str_c(treatment, "_", gene_names))

p <- ggplot(shown_msdata.df,
       aes(x = timepoint_num, y = intensity_imputed,
           color=treatment, fill=treatment)) +
  geom_smooth(alpha=0.25) +
  geom_point(aes(shape=is_imputed)) +
  ylab("Normalized MS Intensity") +
  scale_x_continuous("timepoint", breaks = unique(msdata_full$msruns$timepoint_num)) +
  scale_shape_manual(values = c("FALSE"=16L, "TRUE"=1L)) +
  scale_color_manual(values=treatment_palette) +
  scale_fill_manual(values=treatment_palette) +
  facet_wrap(~ protgroup_label, ncol = 3, scales = "free_y") +
  theme_bw_ast()
p
ggsave(p, filename = file.path(analysis_path, "plots", str_c(msfolder, "_", data_version),
                               paste0(project_id, "_", msfolder, "_viral_proteins_", data_version, ".pdf")),
       width = 15, height = 30, device = cairo_pdf)

p <- ggplot(shown_msdata.df,
            aes(x = timepoint_num, y = intensity_imputed,
                color=treatmentXprotein, fill=treatmentXprotein)) +
    geom_smooth(alpha=0.4) +
    ylab("Intensity") +
    scale_x_continuous("time", breaks = unique(msdata_full$msruns$timepoint_num)) +
    scale_color_manual(values=treatmentXprotein_palette, guide="none") +
    scale_fill_manual(values=treatmentXprotein_palette, guide="none") +
    theme_bw_ast() +
    theme(panel.grid = element_blank(), axis.line = element_blank(),
          panel.border = element_blank(),
          axis.line.x.bottom = element_line(color = "black"),
          axis.line.y.left = element_line(color = "black"))
p
ggsave(p, filename = file.path(analysis_path, "plots", msfolder, paste0(project_id, "_", msfolder, "_SARS_CoV2_proteins_simple_", data_version, ".pdf")),
       width = 3, height = 3, device = cairo_pdf)

shared_cov_peptides.df <- filter(msdata_full$proteins, is_viral) %>%
  dplyr::select(protein_ac, protein_name, gene_name, taxon_id) %>%
  dplyr::inner_join(msdata_full$peptide2protein) %>%
  dplyr::group_by(peptide_id) %>%
  dplyr::filter(any(str_detect(protein_name, "_CVHSA")) & any(str_detect(protein_name, "_SARS2"))) %>%
  dplyr::ungroup()

shared_cov_msdata.df <- dplyr::select(shared_cov_peptides.df, gene_name, peptide_id) %>%
  dplyr::distinct() %>%
  dplyr::inner_join(dplyr::select(msdata_full$pepmods, peptide_id, pepmod_id)) %>%
  dplyr::inner_join(msdata_full$pepmodstate_intensities) %>%
  dplyr::inner_join(dplyr::select(msdata_full$msruns, msrun, treatment, timepoint, timepoint_num, replicate))

shared_cov_msdata_wide.df <- pivot_wider(shared_cov_msdata.df, c(gene_name, pepmodstate_id, peptide_id, timepoint, timepoint_num, replicate),
                                         names_from = treatment, values_from = c(intensity_norm, qvalue)) %>%
  dplyr::mutate(log2_ratio.CoV2_vs_SARS = log2(intensity_norm_SARS_CoV2) - log2(intensity_norm_SARS_CoV))

p <- ggplot(dplyr::filter(shared_cov_msdata_wide.df, abs(log2_ratio.CoV2_vs_SARS) <= 5),
            aes(x = timepoint_num, y = log2_ratio.CoV2_vs_SARS)) +
  geom_smooth(alpha=0.5, fill="gray", color="firebrick") +
  geom_point(position = position_jitter(width=1, height=0), size=0.75, alpha=0.5) +
  scale_x_continuous(breaks = unique(msdata_full$msruns$timepoint_num)) +
  facet_wrap(~ gene_name, ncol = 3, scales = "free_y") +
  theme_bw_ast()
p
ggsave(p, filename = file.path(analysis_path, "plots", str_c(msfolder, "_", data_version),
                               paste0(project_id, "_", msfolder, "_shared_viral_peptides_ratio_", data_version, ".pdf")),
       width = 8, height = 10, device = cairo_pdf)
