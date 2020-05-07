# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200506"
fit_version <- "20200506"
ms_folder <- 'cov2timecourse_dia_20200423'
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
                  instr_calib_protgroup_filename = "instr_QX7_intensity_protgroup_calib_cov2_20200430.json",
                  instr_calib_pepmod_filename = "instr_QX7_intensity_pepmod_calib_cov2_20200506.json",
                  quant_type = "intensity", quant_col_prefix = "intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib_protgroup <- fromJSON(file = file.path(data_path, data_info$instr_calib_protgroup_filename))$instr_calib
instr_calib_pepmod <- fromJSON(file = file.path(data_path, data_info$instr_calib_pepmod_filename))$instr_calib
instr_calib <- instr_calib_protgroup

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/cov_baits_20200415.fasta")),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants.fasta"))
)

bad_msruns <- c()
#filter out msruns with >35% missing values in msdata4norm.df
#msdata4lev1norm.df <- semi_join(msdata4norm.df, filter(msdata_full$msrun_stats, na_ratio <= 0.35))

pgdata.wide <- read.Spectronaut.ProteinsReport(file.path(msdata_path, "20200422_193114_20200422 COVID proteome NO normalization_Report.txt"),
                                               import_data = "quantity", delim='\t')
pgdata_colgroups <- attr(pgdata.wide, "column_groups")
pmsdata.wide <- read.Spectronaut.PepmodstatesReport(file.path(msdata_path, "COVID proteome normalization_peptide_Report.txt"),
                                                    import_data = "quantity", delim='\t')
pmsdata_colgroups <- attr(pmsdata.wide, "column_groups")

msdata_full <- list(
    protgroups = pgdata.wide[pgdata_colgroups$protgroup] %>%
        dplyr::mutate(protgroup_id = row_number()),
    protgroup_intensities = filter(pivot_longer.Spectronaut.ProtgroupIntensities(msdata.wide), is.finite(intensity)),
    pepmodstate_intensities = filter(pivot_longer.Spectronaut.PepmodstateIntensities(pmsdata.wide), is.finite(intensity)),
    pepmodstates = dplyr::select(pmsdata.wide[pmsdata_colgroups$pepmodstate], protgroup_sn_id,
                                 pepmod_seq, charge, peptide_seq, q_value) %>%
        dplyr::mutate(pepmodstate_id = row_number() - 1L)
)
msdata_full$pepmods <- dplyr::select(msdata_full$pepmodstates, pepmod_seq, peptide_seq) %>%
    distinct() %>%
    mutate(pepmod_id = row_number() - 1L)
msdata_full$peptides <- dplyr::select(msdata_full$pepmods, peptide_seq) %>%
    distinct() %>%
    mutate(peptide_id = row_number() - 1L)
msdata_full$pepmods <- left_join(msdata_full$pepmods, dplyr::select(msdata_full$peptides, peptide_id, peptide_seq)) %>%
    dplyr::arrange(pepmod_id)
msdata_full$pepmodstates <- left_join(msdata_full$pepmodstates,
                                      dplyr::select(msdata_full$pepmods, pepmod_id, pepmod_seq)) %>%
    dplyr::select(-pepmod_seq, -peptide_seq) %>%
    dplyr::left_join(dplyr::select(msdata_full$protgroups, protgroup_sn_id, protgroup_id)) %>%
    #dplyr::mutate(is_used = !is.na(protgroup_id)) %>%
    dplyr::arrange(pepmodstate_id)
msdata_full$peptide2protein <- dplyr::distinct(dplyr::select(pmsdata.wide, majority_protein_acs, protgroup_sn_id, peptide_seq)) %>%
    expand_collapsed(collapsed_col = "majority_protein_acs", separated_col = "protein_ac", extra_cols = "peptide_seq") %>%
    dplyr::select(-majority_protein_acs) %>%
    dplyr::left_join(msdata_full$peptides)

msdata_full$msruns <- dplyr::select(msdata_full$protgroup_intensities, msrun_ix, raw_file) %>% dplyr::distinct() %>%
  dplyr::mutate(msrun = str_remove(str_remove(raw_file, "^20200418_QX7_MaTa_SA_proteome_A549_"), "(?:_\\d{3,})?.raw$")) %>% 
  tidyr::extract(msrun, c("treatment", "timepoint", "replicate"), "(.*)_(\\d+)hpi_(\\d+)$", remove=FALSE) %>%
  dplyr::mutate(replicate = parse_integer(replicate), timepoint_num = parse_double(timepoint)) %>%
  #left_join(dplyr::select(baits_info.df, bait_type, bait_code, bait_full_id, bait_id, organism, orgcode)) %>%
  mutate(condition = str_remove(msrun, '_\\d+$'),
         treatment = relevel(factor(treatment), "mock"),
         timepoint = factor(timepoint_num),
         #batch = batch,
         is_used = !msrun %in% bad_msruns) %>%
  dplyr::arrange(treatment, timepoint, replicate) %>%
  dplyr::mutate(msrun = factor(msrun, levels=msrun),
                condition = factor(condition, levels=unique(condition)))

msdata_full$protgroup_intensities <- dplyr::select(msdata_full$protgroup_intensities, -raw_file) %>%
  dplyr::mutate(ident_type = factor(if_else(nevidences > 0L, "By MS/MS", "By matching"), levels = c("By MS/MS", "By matching"))) %>%
  left_join(dplyr::select(msdata_full$msruns, msrun_ix, msrun)) %>%
  dplyr::select(-msrun_ix)
msdata_full$pepmodstate_intensities <- dplyr::select(msdata_full$pepmodstate_intensities, -raw_file) %>%
    left_join(dplyr::select(msdata_full$msruns, msrun_ix, msrun)) %>%
    left_join(dplyr::select(msdata_full$pepmodstates, protgroup_id, pepmodstate_id, pepmod_id)) %>%
    dplyr::select(-msrun_ix)
msdata_full$pepmod_intensities <- dplyr::select(msdata_full$pepmodstate_intensities, -pepmodstate_id) %>%
    dplyr::distinct()

all_proteins.df <- dplyr::bind_rows(
    dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE, is_expected_stable=FALSE),
    dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE,
                  # define stable proteins based on our expectations of their biology and how their timeseries look like
                  is_expected_stable = str_detect(gene_name, "^(M?RP[LS]\\d+|GAPDH|ACT[ABNR]\\d+|TUB[AB]\\d?\\w?|CCT\\d?\\w?)$")),
    dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE, is_expected_stable=FALSE))
msdata_full <- append_protgroups_info(msdata_full, pgdata.wide,
                                      proteins_info = all_proteins.df,
                                      import_columns = c("is_viral", "is_contaminant", "is_expected_stable"))
msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
    is_reverse = FALSE,
    gene_label = strlist_label2(gene_names),
    protac_label = strlist_label2(protein_acs),
    protein_label = strlist_label2(protein_names),
    protgroup_label = case_when(is_viral ~ protein_label,
                                !is.na(gene_label) ~ gene_label,
                                !is.na(protac_label) ~ protac_label,
                                TRUE ~ str_c('#', protgroup_id)))

# redefine protein groups (protregroups) considering only peptides that are clicked
pepmods.df <- dplyr::left_join(pmsdata.wide, dplyr::select(msdata_full$pepmods, pepmod_id, pepmod_seq)) %>%
    dplyr::select(pepmod_id, majority_protein_acs, protein_acs, pepmod_seq, peptide_seq, q_value) %>%
    dplyr::mutate(is_used = coalesce(q_value, 1.0) <= 0.01)
save(file = file.path(msdata_path, str_c(project_id, "_", ms_folder, '_', data_version, "_pepmods.RData")),
     pepmods.df, all_proteins.df)
# .. run protregroup_cov2ts_proteome.jl
msdata_full$protregroups <- read_tsv(file.path(data_path, ms_folder,
                                               str_c(project_id, "_", ms_folder, '_', data_version, "_protregroups_acs.txt")),
                                     col_types = list(protregroup_id = "i")) %>%
  dplyr::mutate(is_contaminant = str_detect(majority_protein_acs, "(;|^)CON__"),
                is_reverse = str_detect(majority_protein_acs, "(;|^)REV__"))

msdata_full$protein2protregroup <- dplyr::select(msdata_full$protregroups, protregroup_id, protein_ac=majority_protein_acs) %>%
  separate_rows(protein_ac, sep=fixed(";"), convert=TRUE) %>%
  dplyr::mutate(is_majority = TRUE) %>%
  dplyr::group_by(protregroup_id) %>%
  dplyr::mutate(protein_ac_rank = row_number()) %>%
  dplyr::ungroup()

msdata_full$protregroup2pepmod <- bind_rows(
  select(msdata_full$protregroups, protregroup_id, pepmod_id=spec_pepmod_ids) %>%
    separate_rows(pepmod_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = TRUE),
  select(msdata_full$protregroups, protregroup_id, pepmod_id=pepmod_ids) %>%
    separate_rows(pepmod_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = FALSE)) %>%
  dplyr::group_by(protregroup_id, pepmod_id) %>%
  dplyr::summarise(is_specific = any(is_specific)) %>%
  dplyr::ungroup()

msdata_full$protregroups <- dplyr::inner_join(msdata_full$protregroups,
  dplyr::inner_join(msdata_full$protregroups, msdata_full$protein2protregroup) %>%
  dplyr::inner_join(msdata_full$proteins) %>%
  dplyr::arrange(protregroup_id, protein_ac_rank) %>%
  dplyr::group_by(protregroup_id) %>%
  dplyr::summarise(gene_names = str_c(gene_name, collapse=';'),
                   gene_label = strlist_label(gene_name),
                   protein_names = str_c(protein_name, collapse=';'),
                   protein_label = strlist_label(protein_name),
                   protein_descriptions = str_c(protein_description, collapse=';'),
                   protein_description = strlist_label(protein_description),
                   organism = str_c(unique(organism), collapse=';'),
                   is_viral = any(coalesce(is_viral, FALSE))) %>%
  dplyr::ungroup()) %>%
  dplyr::mutate(protac_label = sapply(str_split(majority_protein_acs, fixed(";")), strlist_label),
                protregroup_label = case_when(is_viral ~ protein_label,
                                              !is.na(gene_label) ~ gene_label,
                                              !is.na(protac_label) ~ protac_label,
                                              TRUE ~ str_c('#', protregroup_id))) %>%
  dplyr::left_join(dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmods) %>%
                   dplyr::group_by(protregroup_id) %>%
                   dplyr::summarise(npeptides = n_distinct(peptide_id),
                                    npepmods = n_distinct(pepmod_id),
                                    npeptides_unique = n_distinct(peptide_id[is_specific]),
                                    npepmods_unique = n_distinct(pepmod_id[is_specific])) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(npeptides_unique_razor = npeptides_unique,
                                 npeptides_razor = 0L,
                                 npepmods_unique_razor = npepmods_unique,
                                 npepmods_razor = 0L))
# prepare protgroup intensities (wider format: all mstags in one row)
intensity_prespec_df <- tibble(.name = msdata_colgroups$LFQ) %>%
  extract(.name, c("mstag", "msrun"), remove=FALSE,
          str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
  mutate(.value = str_c("intensity.", mstag)) %>%
  dplyr::inner_join(select(msruns.df, msrun, raw_file))

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmodstate_intensities) %>%
  dplyr::mutate(is_idented = !is.na(intensity)) %>%
  dplyr::group_by(msrun, protregroup_id) %>%
  dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                   nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                   npepmods_idented = sum(is_idented),
                   nspecpepmods_idented = sum(is_idented & is_specific)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))

# condition = treatment X timepoint
conditions.df <- dplyr::select(msdata_full$msruns, condition, treatment, timepoint, timepoint_num) %>%
  dplyr::distinct()
conditions.df <- dplyr::left_join(conditions.df,
    dplyr::select(conditions.df, treatment, timepoint_after=timepoint_num)) %>% # exclude 3h because it's the reference
    dplyr::mutate(is_after = timepoint_num >= timepoint_after) %>%
    tidyr::pivot_wider(all_of(c("condition", "treatment", "timepoint", "timepoint_num")),
                       names_prefix = "after", names_from = timepoint_after, values_from = is_after) %>%
    rename_at(vars(starts_with("after")), ~str_c(., "h"))

msdata_full$msrun_stats <- msrun_statistics(msdata_full) %>%
  dplyr::mutate(na_ratio = n_missing / n)
  
msdata <- msdata_full[c('protgroup_intensities', 'pepmod_intensities', 'pepmods',
                        'msruns', 'protgroups', 'protregroups', 'protein2protregroup',
                        'protein2pepmod', 'protregroup2pepmod')]
msdata$protgroup_intensities <- semi_join(msdata$protgroup_intensities, filter(msdata$msruns, is_used))
msdata$pepmod_intensities <- semi_join(msdata$pepmod_intensities, filter(msdata$msruns, is_used))

# setup experimental design matrices
conditionXeffect_orig.mtx <- model.matrix(
  formula(str_c("~ 1 + (", str_c(str_subset(colnames(conditions.df), "^after\\d+h"), collapse =" + "), ") * treatment")),
  conditions.df
)
conditionXeffect.mtx <- conditionXeffect_orig.mtx[, str_detect(colnames(conditionXeffect_orig.mtx),
                                                               "^treatment[^:]+$", negate=TRUE) &
                                                    apply(conditionXeffect_orig.mtx, 2, min) == 0] # remove intercept and afterT0
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
         prior_tau = case_when(effect_type == "treatmentXtimepoint" ~ 0.5,
                               effect_type == "treatment" ~ 0.5,
                               effect_type == "timepoint" ~ 1.0,
                               TRUE ~ NA_real_),
         prior_df2 = case_when(effect_type == "treatmentXtimepoint" ~ 2.0,
                               TRUE ~ 1.0),
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
msdata4norm.df <- msdata$protgroup_intensities %>% ungroup() %>%
  dplyr::filter(ident_type == "By MS/MS" & !is.na(intensity)) %>% 
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$protgroup_intensities) %>%
  dplyr::inner_join(select(msdata$msruns, msrun, condition, treatment, timepoint)) #%>% #WHY???
  #mutate(timepoint=factor(timepoint))

options(mc.cores=8)

msruns_hnorm <- multilevel_normalize_experiments(instr_calib_protgroup, msdata$msruns,
                                                 msdata4norm.df,
                                                 quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
                                                 mcmc.iter = 2000L,
                                                 #mcmc.chains = 6,
                                                 verbose=TRUE,
                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=2000L, missing_exp.ratio=0.1),
                                                                    condition = list(cond_col="condition", max_objs=2000L, missing_exp.ratio=0.1),
                                                                    timepoint = list(cond_col="timepoint", max_objs=1000L, missing_exp.ratio=0.3)
                                                 ))

#total_msrun_shifts.df <- msruns_hnorm$mschannel_shifts
total_msrun_shifts.df <- transmute(msdata$msruns, msrun, total_msrun_shift=0.0) # no normalization

## apply normalization
msdata$protgroup_intensities <- dplyr::left_join(msdata$protgroup_intensities,
                                                 dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) 
#dplyr::select(-total_msrun_shift)

global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns, TRUE) %>%
                                               dplyr::select(msrun) %>% dplyr::distinct() %>%
                                               dplyr::inner_join(msdata$protgroup_intensities) %>% .$intensity), na.rm=TRUE)
global_pepmod_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns, TRUE) %>%
                                        dplyr::select(msrun) %>% dplyr::distinct() %>%
                                        dplyr::inner_join(msdata$pepmod_intensities) %>% .$intensity), na.rm=TRUE)

############################
# batch effects
# no batch effects so far
msrunXbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXreplEffect.mtx),
                                     batch_effect = c())

batch_effects.df <- tibble(batch_effect = character(0),
                           is_positive = logical(0),
                           prior_mean = double(0))

# no subbatch effects so far
msrunXsubbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXreplEffect.mtx),
                                        subbatch_effect = c())

subbatch_effects.df <- tibble(subbatch_effect=character(0),
                              is_positive=logical(0))

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', ms_folder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df, contrasts.df,
     conditionXeffect.mtx, inv_conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     instr_calib_protgroup, instr_calib_pepmod,
     global_protgroup_labu_shift, global_pepmod_labu_shift,
     total_msrun_shifts.df, #msruns_hnorm, 
     msrunXreplEffect.mtx,
     batch_effects.df, msrunXbatchEffect.mtx,
     subbatch_effects.df, msrunXsubbatchEffect.mtx,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData'))
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
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
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
protgroup_hclu = hclust(dist(protgroup_intensities_imp.mtx))
pheatmap(protgroup_intensities.mtx[, msruns_ordered], cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.pdf")), width=16, height=80)

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
protgroup_intensities_imp_pg.mtx <- protgroup_intensities_imp.mtx - rowMedians(protgroup_intensities_imp.mtx)
pg_umap2d = umap(protgroup_intensities_imp_norm.mtx, n_components=2,
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
                           is_viral, is_expected_stable, is_contaminant, is_reverse)) %>%
  mutate(category = case_when(replace_na(is_viral, FALSE) ~ "viral",
                              is_reverse ~ "reverse",
                              is_contaminant ~ "contaminant",
                              is.na(gene_names) ~ "NA",
                              is_expected_stable ~ "stable",
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

p <- ggplot(filter(msdata_full$protgroups, str_detect(gene_names, "(^|;)(?:TUB[AB]\\d?)")) %>%
       inner_join(msdata_full$protgroup_intensities_all) %>%
       left_join(msdata_full$msruns) %>%
       aes(x = timepoint_num, y = intensity_imputed_norm, color=treatment, fill=treatment)) +
  geom_smooth(alpha=0.25) +
  geom_point(aes(shape=is_imputed)) +
  scale_x_continuous(breaks = unique(conditions.df$timepoint_num)) +
  scale_shape_manual(values = c("FALSE"=16L, "TRUE"=1L)) +
  scale_color_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  scale_fill_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  facet_wrap(~ gene_names, ncol = 3, scales = "free_y") +
  theme_bw_ast()
p
ggsave(p, filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_LDH_", data_version, ".pdf")),
       width = 16, height = 12, device = cairo_pdf)
