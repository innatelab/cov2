# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200829"
fit_version <- "20200829"
msfolder <- 'snaut_parsars_fp_20200829'
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

msdata_path <- file.path(data_path, msfolder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  msfolder = msfolder,
                  mscalib_protgroup_filename = "instr_QX7_intensity_protgroup_calib_cov2_20200430.json",
                  mscalib_pepmod_filename = "mscalib_EXPL2_intensity_pepmodstate_cov2_20200828.json",
                  quant_type = "intensity", quant_col_prefix = "intensity", qvalue_max=1E-4,
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$mscalib_protgroup_filename, '...')
mscalib_protgroup <- fromJSON(file = file.path(data_path, data_info$mscalib_protgroup_filename))$instr_calib
mscalib_pepmod <- fromJSON(file = file.path(data_path, data_info$mscalib_pepmod_filename))$instr_calib
mscalib <- mscalib_protgroup

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/SARS_CoV_20200817.fasta")),
  CoV2 = read_innate_uniprot_fasta(file.path(data_path, "msfasta/SARS_CoV2_20200817.fasta")),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants_20200405.fasta"))
)

bad_msruns <- c()
#filter out msruns with >35% missing values in msdata4norm.df
#msdata4lev1norm.df <- semi_join(msdata4norm.df, filter(msdata_full$msrun_stats, na_ratio <= 0.35))

pgdata.wide <- read.Spectronaut.ProteinsReport(file.path(msdata_path, "Protein_PG_normalized_precursor cutoff 0.001_Report.txt"),
                                               import_data = "quantity", delim='\t')
pgdata_colgroups <- attr(pgdata.wide, "column_groups")

pmsdata.wide <- read.Spectronaut.PepmodstatesReport(file.path(msdata_path, "Protein_EG_normalized_precursor cutoff 0.001_Report.txt"),
                                                    import_data = c("eg_quantity", "eg_qvalue", "eg_snratio", "eg_normfactor"), delim='\t')
pmsdata_colgroups <- attr(pmsdata.wide, "column_groups")

msdata_full <- list(
    protgroups = pgdata.wide[pgdata_colgroups$protgroup],
    pepmodstates = dplyr::select(pmsdata.wide[pmsdata_colgroups$pepmodstate],
                                 pepmodstate_id, pepmod_seq, charge, peptide_seq, is_pg_specific),
    protgroup_intensities = pivot_longer.Spectronaut.Metrics(pgdata.wide, c("protgroup_sn_id"), c("quantity", "runstats")) %>%
          dplyr::distinct() %>% dplyr::filter(!is.na(intensity)),
    pepmodstate_intensities = pivot_longer.Spectronaut.Metrics(pmsdata.wide, c("pepmodstate_id"), c("eg_quantity", "eg_qvalue", "eg_snratio", "eg_normfactor")) %>%
          dplyr::distinct() %>% dplyr::filter(!is.na(intensity))
)
msdata_full$msruns <- read_tsv(file.path(msdata_path, "msruns_info.txt"), col_types = c(msrun_ix="i", sample_ix="i")) %>%
    dplyr::left_join(dplyr::select(msdata_full$pepmodstate_intensities, rawfile_ix, rawfile, msrun_normfactor) %>%
                     dplyr::filter(!is.na(msrun_normfactor)) %>% dplyr::distinct() %>%
                     tidyr::extract(rawfile, "msrun_ix", "_(\\d+)(?:_\\d+)?\\.raw", remove=FALSE) %>%
                     dplyr::mutate(msrun_ix = parse_integer(msrun_ix))) %>%
    dplyr::arrange(msrun_ix) %>%
    dplyr::mutate(treatment = factor(treatment, c("mock", "SARS_CoV", "SARS_CoV2")),
                  timepoint_num = timepoint,
                  timepoint = factor(timepoint),
                  replicate = as_integer(replicate),
                  condition = str_c(treatment, "_", timepoint, "h"),
                  msrun = str_c(condition, "_", replicate),
                  is_used = !msrun %in% bad_msruns,
                  msrun_mult = 1/msrun_normfactor,
                  msrun_shift = log(msrun_mult)) %>%
  dplyr::arrange(treatment, timepoint_num, replicate) %>%
  dplyr::mutate(#msrun = factor(msrun, levels=msrun),
                condition = factor(condition, levels=unique(condition)))

msdata_full$protgroup_intensities <- dplyr::left_join(msdata_full$protgroup_intensities,
                                               dplyr::select(msdata_full$msruns, rawfile_ix, msrun, msrun_mult)) %>%
  dplyr::mutate(intensity_norm = intensity,
                intensity = intensity_norm * msrun_mult,
                ident_type = factor("By MS/MS", levels = c("By MS/MS", "By matching"))) %>%
  dplyr::left_join(pivot_longer.Spectronaut.Metrics(pgdata.wide, c("protgroup_id", "protgroup_sn_id"), c("quantity", "pg_stats")) %>%
                   dplyr::rename(intensity.2 = intensity)) %>%
  dplyr::select(-rawfile_ix, -rawfile)
msdata_full$pepmodstate_intensities <- dplyr::mutate(msdata_full$pepmodstate_intensities,
                                             intensity_norm = intensity,
                                             intensity = intensity_norm/msrun_normfactor,
                                             msrun_mult = NULL)
msdata_full$pepmods <- msdata_full$pepmodstate_intensities %>%
    dplyr::left_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, pepmod_seq, peptide_seq, is_pg_specific)) %>%
    dplyr::group_by(pepmod_seq, peptide_seq, is_pg_specific) %>%
    dplyr::summarise(qvalue = min(qvalue, na.rm=TRUE), .groups="drop") %>%
    dplyr::distinct() %>%
    dplyr::mutate(pepmod_id = row_number() - 1L)
msdata_full$peptides <- dplyr::select(msdata_full$pepmods, peptide_seq, is_pg_specific) %>%
    distinct() %>%
    mutate(peptide_id = row_number() - 1L)
msdata_full$pepmods <- left_join(msdata_full$pepmods, dplyr::select(msdata_full$peptides, peptide_id, peptide_seq)) %>%
    dplyr::arrange(pepmod_id)
msdata_full$pepmodstates <- left_join(msdata_full$pepmodstates,
                                      dplyr::select(msdata_full$pepmods, pepmod_id, pepmod_seq)) %>%
    dplyr::select(-pepmod_seq, -peptide_seq) %>%
    dplyr::arrange(pepmodstate_id)
msdata_full$peptide2protein <- dplyr::distinct(dplyr::select(pmsdata.wide, protein_acs, peptide_seq)) %>%
    tidyr::separate_rows(c("protein_acs"), sep=";") %>%
    dplyr::rename(protein_ac = protein_acs) %>%
    dplyr::left_join(dplyr::select(msdata_full$peptides, peptide_seq, peptide_id))
msdata_full$protein2protgroup = dplyr::select(msdata_full$protgroups, protgroup_id, protein_acs) %>%
  tidyr::separate_rows(c("protein_acs"), sep=";") %>%
  dplyr::rename(protein_ac = protein_acs)
msdata_full$peptide2protgroup <- dplyr::inner_join(msdata_full$peptide2protein,
                                                   msdata_full$protein2protgroup) %>%
  dplyr::select(-protein_ac) %>% dplyr::distinct()
msdata_full$pepmodstate_intensities <- dplyr::select(msdata_full$pepmodstate_intensities, -rawfile) %>%
    left_join(dplyr::select(msdata_full$msruns, rawfile_ix, msrun)) %>%
    left_join(dplyr::select(msdata_full$pepmodstates, pepmodstate_id, pepmod_id)) %>%
    dplyr::select(-rawfile_ix)

all_proteins.df <- dplyr::bind_rows(
    dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE, is_expected_stable=FALSE),
    dplyr::mutate(fasta.dfs$CoV2, is_viral=TRUE, is_contaminant=FALSE, is_expected_stable=FALSE),
    dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE,
                  # define stable proteins based on our expectations of their biology and how their timeseries look like
                  is_expected_stable = str_detect(gene_name, "^(M?RP[LS]\\d+|GAPDH|ACT[ABNR]\\d+|TUB[AB]\\d?\\w?|CCT\\d?\\w?)$")),
    dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE, is_expected_stable=FALSE))
msdata_full <- append_protgroups_info(msdata_full, pgdata.wide,
                                      proteins_info = all_proteins.df,
                                      import_columns = c("is_viral", "is_contaminant", "organism"))
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
pepmods.df <- dplyr::left_join(pmsdata.wide, dplyr::select(msdata_full$pepmods, pepmod_id, pepmod_seq, qvalue)) %>%
    dplyr::select(pepmod_id, protein_acs, pepmod_seq, peptide_seq, qvalue) %>%
    dplyr::mutate(pepmod_rank = if_else(coalesce(qvalue, 1.0) <= 0.0001, 1L, -1L))
save(file = file.path(msdata_path, str_c(project_id, "_", msfolder, '_', data_version, "_pepmods.RData")),
     pepmods.df, all_proteins.df)
# .. run protregroup_cov2ts_proteome.jl
msdata_full$protregroups <- read_tsv(file.path(data_path, msfolder,
                                               str_c(project_id, "_", msfolder, '_', data_version, "_protregroups_acs.txt")),
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
intensity_prespec_df <- tibble(.name = pgdata_colgroups$quantity) %>%
  extract(.name, c("mstag", "msrun"), remove=FALSE,
          str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
  mutate(.value = str_c("intensity.", mstag)) %>%
  dplyr::inner_join(select(msdata_full$msruns, msrun, rawfile))

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
    rename_at(vars(starts_with("after")), ~str_c(., "h")) %>%
    mutate(after0h = NULL, # not needed for exp_design
           infected = treatment != "mock")

msdata_full$msrun_stats <- msrun_statistics(msdata_full) %>%
  dplyr::mutate(na_ratio = n_missing / n)
  
msdata <- msdata_full[c('protgroup_intensities', 'pepmodstate_intensities', 'pepmods', 'pepmodstates',
                        'msruns', 'protgroups', 'protregroups', 'protein2protregroup',
                        'peptide2protein', 'protregroup2pepmod')]
msdata$msruns <- filter(msdata$msruns, !is.na(rawfile_ix))
msdata$pepmodstate_intensities$ident_type <- "MULTI-MSMS" # HACK
msdata$msrun_pepmodstate_stats <- msrun_statistics(msdata, obj="pepmodstate")
msdata$protgroup_intensities <- semi_join(msdata$protgroup_intensities, filter(msdata$msruns, is_used))
msdata$pepmodstate_intensities <- filter(semi_join(msdata$pepmodstate_intensities, filter(msdata$msruns, is_used)),
                                         coalesce(qvalue, 1.0) <= data_info$qvalue_max)

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
                               effect_type == "timepoint" ~ 1.0,
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
         filename = file.path(analysis_path, "plots", str_c(msfolder, "_", fit_version),
                              paste0(project_id, "_metaconditions_", msfolder, "_", fit_version, ".pdf")),
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
         filename = file.path(analysis_path, "plots", str_c(msfolder, "_", fit_version),
                              paste0(project_id, "_exp_design_contrasts_", msfolder, "_", fit_version, ".pdf")),
         width = 8, height = 12)

##########
## normalization
#msdata4norm.df <- msdata$protgroup_intensities %>% ungroup() %>%
#  dplyr::filter(ident_type == "By MS/MS" & !is.na(intensity)) %>% 
#  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
#  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
#  dplyr::inner_join(msdata$protgroup_intensities) %>%
#  dplyr::inner_join(select(msdata$msruns, msrun, condition, treatment, timepoint)) #%>% #WHY???
#  #mutate(timepoint=factor(timepoint))

options(mc.cores=8)

#msruns_hnorm <- multilevel_normalize_experiments(instr_calib_protgroup, msdata$msruns,
#                                                 msdata4norm.df,
#                                                 quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
#                                                 mcmc.iter = 2000L,
#                                                 #mcmc.chains = 6,
#                                                 verbose=TRUE,
#                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=2000L, missing_exp.ratio=0.1),
#                                                                    condition = list(cond_col="condition", max_objs=2000L, missing_exp.ratio=0.1),
#                                                                    timepoint = list(cond_col="timepoint", max_objs=1000L, missing_exp.ratio=0.3)
#                                                 ))

#total_msrun_shifts.df <- msruns_hnorm$mschannel_shifts
total_msrun_shifts.df <- transmute(msdata$msruns, msrun, total_msrun_shift=msrun_shift)

## apply normalization
#msdata$protgroup_intensities <- dplyr::left_join(msdata$protgroup_intensities,
#                                                 dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
#  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift)) 
#dplyr::select(-total_msrun_shift)

global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns, TRUE) %>%
                                               dplyr::select(msrun) %>% dplyr::distinct() %>%
                                               dplyr::inner_join(msdata$protgroup_intensities) %>% .$intensity), na.rm=TRUE)
global_pepmodstate_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns, TRUE) %>%
                                             dplyr::select(msrun) %>% dplyr::distinct() %>%
                                             dplyr::inner_join(msdata$pepmodstate_intensities) %>% .$intensity), na.rm=TRUE)

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

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', msfolder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df, contrasts.df,
     conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     instr_calib_protgroup, instr_calib_pepmod,
     global_protgroup_labu_shift, global_pepmodstate_labu_shift,
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
