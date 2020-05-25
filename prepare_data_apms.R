# SARS-CoV/CoV-2 loading and preparing AP-MS data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200515"
fit_version <- "20200515"
ms_folder <- 'mq_apms_20200510'
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

mqdata_path <- file.path(data_path, ms_folder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  ms_folder = ms_folder,
                  instr_calib_protgroup_filename = "instr_protgroup_LFQ_calib_scaturro_qep5calib_20161110_borg.json",
                  #instr_calib_pepmodstate_filename = "instr_QEP5_intensity_pepmodstate_calib_cov2_20200428_borg.json",
                  instr_calib_pepmodstate_filename = "instr_QEP5_intensity_pepmodstate_calib_cov2_20200430.json",
                  #instr_calib_pepmodstate_filename = "instr_pepmod_intensity_raw_calib_laudenbach_pcp_20170128_borg.json",
                  quant_type = "LFQ", quant_col_prefix = "LFQ_Intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib_protgroup <- fromJSON(file = file.path(data_path, data_info$instr_calib_protgroup_filename))$instr_calib
instr_calib_pepmodstate <- fromJSON(file = file.path(data_path, data_info$instr_calib_pepmodstate_filename))$instr_calib

source(file.path(project_scripts_path, 'prepare_data_common.R'))

msruns.df <- read_tsv(file.path(mqdata_path, "combined", "experimentalDesign.txt"),
                      col_names=TRUE, col_types = list(Fraction="i")) %>%
  rename(raw_file=Name, msfraction=Fraction, msrun=Experiment, is_ptm=PTM) %>%
  extract(msrun, c("sample_type", "batch", "bait_full_id", "replicate", "tech_replicate"),
          c("^([^_]+)_B(\\d+)_(.+)_(\\d+)(?:#(\\d+))?$"), remove = FALSE) %>%
  mutate(bait_full_id_tmp = str_remove(bait_full_id, "\\?+$")) %>%
  left_join(select(baits_info.df, bait_full_id_tmp=bait_full_id, bait_id, bait_homid, bait_kind, orgcode)) %>%
  select(-bait_full_id_tmp) %>%
  mutate(sample_type = relevel(factor(sample_type), "FPMS"),
         # FIXME remove for the next iteration
         batch = if_else(bait_full_id == "HUMAN_ACE2", "4", batch),
         # FIXME remove for the next iteration
         bait_full_id = if_else(bait_full_id == "HCoV_ORF4" & batch==2L, "HCoV_ORF4?", as.character(bait_full_id)),
         # fix bait id for the candidates to throw away, homid is not touched
         bait_id = str_c(bait_id, str_remove(bait_full_id, "^[^?]+")),
         bait_full_id = if_else(batch == "2" & str_detect(bait_full_id, "SARS_CoV_NSP(1|3_macroD)\\?"),
                                str_c("Ctrl_", bait_full_id), as.character(bait_full_id))) %>%
  # declare AP runs with bait sequence problems as controls
  dplyr::mutate(bait_kind = factor(if_else(str_detect(bait_full_id, "Ctrl_.+\\?"),
                                           "control", as.character(bait_kind)),
                                   levels=levels(bait_kind)),
                bait_id = relevel(factor(bait_id), "Ctrl_NT"),
                bait_full_id = relevel(factor(bait_full_id), "Ctrl_NT"))
baits_info.df$bait_id = factor(baits_info.df$bait_id, levels=levels(msruns.df$bait_id))
baits_info.df$bait_full_id = factor(baits_info.df$bait_full_id, levels=levels(msruns.df$bait_full_id))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "cov_baits_20200415.fasta")) %>%
    dplyr::select(-protein_name) %>%
    left_join(select(baits_info.df, protein_ac=uniprot_ac, protein_name=bait_full_id)),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta"))
)

msdata.wide <- read.MaxQuant.ProteinGroups(file.path(mqdata_path, 'combined/txt'),
                                           file_name = "proteinGroups.txt",
                                           import_data = c(data_info$quant_type, "ident_type", "ms2_count"))
mqevidence <- read.MaxQuant.Evidence(file.path(mqdata_path, 'combined/txt'), 
                                     evidence.pepobj = "pepmodstate", correct_ratios = FALSE)
mqevidence$peptides <- read.MaxQuant.Peptides(file.path(mqdata_path, 'combined/txt'),
                                              file_name = "peptides.txt",
                                              import_data = c("ident_type"))
mqevidence$peaks <- NULL # exclude big data frame

mqrdata_filepath <- file.path(scratch_path, str_c(project_id, '_mqdata_APMS_', ms_folder, '.RData'))
message('Saving imported MaxQuant MS data to ', mqrdata_filepath, '...')
save(data_info, msruns.df, fasta.dfs,
     msdata.wide, mqevidence, instr_calib_pepmodstate, instr_calib_protgroup,
     file = mqrdata_filepath, compress = "xz")
load(mqrdata_filepath)

msdata_colgroups <- attr(msdata.wide, "column_groups")

msdata_full <- list()

msruns.df <- mutate(msruns.df, is_used = msrun %in% mqevidence$mschannels$msrun)
msdata_full$msruns <- filter(msruns.df, is_used)

msdata_full <- append_protgroups_info(msdata_full, msdata.wide,
                                      proteins_info = dplyr::bind_rows(
                                        dplyr::mutate(fasta.dfs$CoV, is_viral = TRUE),
                                        dplyr::mutate(fasta.dfs$human, is_viral = FALSE)),
                                      import_columns = c("is_viral", "organism"))

msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

pacnoiso_stats.df <- left_join(msdata_full$protein2protgroup,
                               select(msdata_full$proteins, protein_ac, protein_ac_noiso)) %>%
  #filter(is_majority) %>%
  group_by(protein_ac_noiso) %>%
  summarise(n_noiso_pgs_have_razor = sum(npeptides_razor > 0),
            n_noiso_pgs = n_distinct(protgroup_id)) %>%
  ungroup()

msdata_full$proteins <- left_join(msdata_full$proteins, pacnoiso_stats.df)

pg_razor_stats.df <- filter(msdata_full$protein2protgroup, is_majority) %>%
  left_join(select(msdata_full$proteins, protein_ac, n_noiso_pgs, n_noiso_pgs_have_razor)) %>%
  group_by(protgroup_id) %>%
  summarise(nproteins_have_razor = sum(npeptides_razor > 0),
            nprotgroups_sharing_proteins = max(n_noiso_pgs)) %>%
  dplyr::ungroup()

msdata_full$protgroups <- left_join(msdata_full$protgroups, pg_razor_stats.df)
msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
                                        gene_label = strlist_label2(gene_names),
                                        protac_label = strlist_label2(protein_acs),
                                        protein_label = strlist_label2(protein_names),
                                        protgroup_label = case_when(is_viral ~ protein_label,
                                                                    !is.na(gene_label) ~ gene_label,
                                                                    !is.na(protac_label) ~ protac_label,
                                                                    TRUE ~ str_c('#', protgroup_id)))

msdata_full$peptides <- mqevidence$peptides %>%
  dplyr::left_join(select(msdata_full$proteins, lead_razor_protein_ac = protein_ac, is_viral)) %>%
  mutate(is_viral=replace_na(is_viral, FALSE))
apms_pepmod_ids <- filter(mqevidence$pepmodstate_intensities, ident_type %in% c("MULTI-MSMS", "MSMS")) %>%
  semi_join(filter(msruns.df, sample_type == "APMS")) %>%
  pull(pepmod_id) %>% unique()
msdata_full$pepmods <- mqevidence$pepmods %>%
  dplyr::left_join(select(msdata_full$peptides, peptide_id, is_viral)) %>%
  dplyr::mutate(is_in_apms = pepmod_id %in% apms_pepmod_ids,
                pepmod_rank = if_else(is_in_apms, 1L, -1L)) # TODO update if proteome would be co-quanted 
msdata_full$pepmodstates <- mqevidence$pepmodstates

# redefine protein groups (protregroups)
pepmods.df <- dplyr::select(msdata_full$pepmods, pepmod_id, protgroup_ids, protein_acs, lead_protein_acs, seq, modifs, charges, is_reverse, pepmod_rank)
proteins.df <- msdata_full$proteins
save(file = file.path(mqdata_path, str_c(project_id, "_", ms_folder, '_', data_version, "_pepmods.RData")),
     pepmods.df, proteins.df)
# .. run protregroup_apms.jl
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

# prepare protgroup intensities (longer format: each mstag on ots own row)
protgroup_intensities_all.df <- tidyr::pivot_longer_spec(
  dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$LFQ),
  mutate(intensity_prespec_df, .value = "intensity")) %>%
  select(-mstag)

protgroup_idents_all.df <- tidyr::pivot_longer(
  dplyr::select(msdata.wide, protgroup_id, !!msdata_colgroups$ident_type, !!msdata_colgroups$ms2_count),
  cols = c(msdata_colgroups$ident_type, msdata_colgroups$ms2_count),
  names_to = c(".value", "msrun"),
  names_pattern = "(ident_type|ms2_count)\\.(.*)") %>%
  dplyr::left_join(select(msruns.df, msrun) %>% distinct())
if (!has_name(protgroup_idents_all.df, "ident_type")) {
  # when matching disabled no ident_type is written
  protgroup_idents_all.df <- mutate(protgroup_idents_all.df,
                                    ident_type = factor(if_else(ms2_count > 0, "By MS/MS", NA_character_),
                                                        levels=c("By matching", "By MS/MS")))
}
msdata_full$protgroup_intensities <- protgroup_intensities_all.df %>%
  dplyr::semi_join(select(msdata_full$msruns, msrun, raw_file)) %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::select(-raw_file)
msdata_full$protgroup_idents <- protgroup_idents_all.df %>%
  dplyr::semi_join(select(msdata_full$msruns, msrun)) %>%
  dplyr::filter(!is.na(ident_type))

msdata_full$pepmodstate_intensities <- mqevidence$pepmodstate_intensities %>%
  select(pepmod_id, pepmodstate_id, msrun, intensity = intensity.Sum, ident_type) %>%
  dplyr::filter(!is.na(intensity)) %>%
  dplyr::mutate(is_idented = str_detect(ident_type, "MSMS"))

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmodstate_intensities) %>%
  dplyr::group_by(msrun, protregroup_id) %>%
  dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                   nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                   npepmods_idented = sum(is_idented),
                   nspecpepmods_idented = sum(is_idented & is_specific)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))

# condition = bait
conditions.df <- dplyr::select(msdata_full$msruns, sample_type, bait_full_id, bait_id, bait_homid, bait_kind, orgcode) %>%
  dplyr::distinct() %>%
  dplyr::arrange(sample_type, bait_kind, bait_homid, bait_id, orgcode) %>%
  dplyr::mutate(condition = str_c(sample_type, "_", bait_full_id),
                condition = factor(condition, levels=condition))
msdata_full$msruns <- left_join(msdata_full$msruns, select(conditions.df, sample_type, bait_full_id, condition))

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

msdata <- msdata_full[c('protgroup_intensities', 'protgroup_idents', 'protregroup_idents',
                        'msruns', 'protgroups', 'protein2protgroup',
                        'pepmodstate_intensities',
                        'pepmodstates', 'pepmods',
                        'protregroups', 'protein2protregroup', 'protregroup2pepmod')]
msdata$msrun_pepmodstate_stats <- msrun_statistics(msdata, obj="pepmodstate")

# setup experimental design matrices
require(Matrix)
require(purrr)
# not grouping baits by bait_homid, since we are only interested in SARS_vs_SARS2 comparisons
bait_exp_designs <- filter(conditions.df, sample_type=="APMS") %>%
# temporary append "?" to avoid grouping together alternative samples
group_by(bait_id) %>%
mutate(n_homologs = n()) %>%
ungroup() %>%
mutate(has_homologs = n_homologs > 1) %>%
group_by(has_homologs) %>% nest() %>%
mutate(mtx = map(data, function(conds.df){
   if (any(conds.df$n_homologs > 1)) {
      res <- mutate(conds.df, bait_id_copy = bait_id) %>%
      group_by(bait_id_copy) %>% nest() %>%
      mutate(mtx = map(data, function(bait_conds.df){
          bait_conds.df <- mutate(bait_conds.df, orgcode = relevel(factor(as.character(orgcode)), "SARS2"),
                                  effect = str_c("bait_id", bait_id,
                                                 if_else(orgcode == "SARS2", "", str_c(":orgcode", orgcode)))) %>%
            arrange(desc(orgcode))
          res <- cbind(matrix(rep.int(1, nrow(bait_conds.df)), ncol = 1),
                       -contr.sum(nrow(bait_conds.df)))
          dimnames(res) <- list(condition = bait_conds.df$condition,
                                effect = bait_conds.df$effect)
          return(res)
      }))
      mtx <- as.matrix(bdiag(res$mtx))
      dimnames(mtx) <- list(condition = flatten_chr(map(res$mtx, rownames)),
                            effect = flatten_chr(map(res$mtx, colnames)))
      return(mtx)
   } else {
      mtx <- model.matrix(~ 1 + bait_id, conds.df)
      dimnames(mtx) <- list(condition = conds.df$condition,
                            effect = colnames(mtx))
      return(mtx)
   }
}))

bait_exp_design.mtx <- as.matrix(bdiag(bait_exp_designs$mtx))
dimnames(bait_exp_design.mtx) <-
  list(condition = flatten_chr(map(bait_exp_designs$mtx, rownames)),
       effect = flatten_chr(map(bait_exp_designs$mtx, colnames)))

conditionXsample_type.mtx <- model.matrix(~ 1 + sample_type, data = conditions.df)
conditionXeffect_orig.mtx <- cbind(conditionXsample_type.mtx,
                                   rbind(zero_matrix(condition = as.character(filter(conditions.df, sample_type != "APMS")$condition),
                                                     effect = flatten_chr(map(bait_exp_designs$mtx, colnames))),
                                         bait_exp_design.mtx[as.character(filter(conditions.df, sample_type == "APMS")$condition), ]))
dimnames(conditionXeffect_orig.mtx) <-
    list(condition = conditions.df$condition,
         effect = c(colnames(conditionXsample_type.mtx), flatten_chr(map(bait_exp_designs$mtx, colnames))))
conditionXeffect.mtx <- conditionXeffect_orig.mtx[as.character(conditions.df$condition),
                                                  colSums(abs(conditionXeffect_orig.mtx)) != 0 &
                                                    !str_detect(colnames(conditionXeffect_orig.mtx), "\\(Intercept\\)|bait_idCtrl.+:orgcode")]

pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(analysis_path, 'plots', str_c(ms_folder, "_", fit_version),
                              paste0(project_id, "_exp_design_", ms_folder, "_", fit_version, ".pdf")),
         width = 14, height = 14)

effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  dplyr::mutate(sample_type = effect_factor(effect, "sample_type", levels(conditions.df$sample_type), NA),
                orgcode = effect_factor(effect, "orgcode", levels(conditions.df$orgcode), NA),
                bait_id = effect_factor(effect, "bait_id", levels(conditions.df$bait_id), NA),
                is_positive = FALSE,#!is.na(bait_id) & is.na(orgcode),
                prior_mean = 0.0)
effects.df$effect_label <- sapply(1:nrow(effects.df), function(i) {
  comps <- c()
  st <- as.character(effects.df$sample_type[[i]])
  org <- as.character(effects.df$orgcode[[i]])
  ba <- as.character(effects.df$bait_id[[i]])
  if (!is.na(st)) comps <- append(comps, st)
  if (!is.na(ba)) comps <- append(comps, ba)
  if (!is.na(org)) comps <- append(comps, org)
  str_c(comps, collapse="+")
})
effects.df <- dplyr::mutate(effects.df,
                            effect = factor(effect, levels = effect),
                            effect_type = case_when(!is.na(sample_type) ~ "sample_type",
                                                    !is.na(bait_id) & !is.na(orgcode) ~ "baitXvirus",
                                                    !is.na(bait_id) ~ "bait"),
                            effect_label = factor(effect_label, levels=effect_label)) %>%
  dplyr::mutate(prior_tau = case_when(effect_type == "baitXvirus" ~ 1.0,
                                      effect_type == "bait" ~ 1.0,
                                      effect_type == "sample_type" ~ 5.0,
                                      TRUE ~ 5.0))
# prior_mean for baitXvirus is calculated as the shift between the homologous baits at the normalization step

conditionXeffect.df <- conditionXeffect_frame(conditionXeffect.mtx, effects.df)

batch_complements.df <- distinct(select(filter(msdata$msruns, sample_type=="APMS" & bait_kind == "sample"), batch, bait_homid)) %>%
  dplyr::mutate(allminus_metacondition = paste0("B", batch, "_allminus_", bait_homid))

bait_conditions <- as.character(filter(conditions.df, sample_type == "APMS" & bait_kind == "sample")$bait_full_id)
allminus_metaconditions <- paste0("allminus_", unique(filter(conditions.df, sample_type == "APMS" & bait_kind == "sample")$bait_homid))
compound_metaconditions <- c(allminus_metaconditions, "controls",
                             batch_complements.df$allminus_metacondition)
all_metaconditions <- c("FPMS", bait_conditions, compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in bait_conditions) {
  conditionXmetacondition.mtx[str_c("APMS_", cname), cname] <- TRUE
}
for (cname in allminus_metaconditions) {
  hombait <- str_remove(cname, "^allminus_")
  minus_conditions.df <- filter(conditions.df, sample_type == "APMS" & ((!str_detect(bait_id, "\\?$") & bait_homid != hombait) | bait_kind == "control"))
  conditionXmetacondition.mtx[as.character(minus_conditions.df$condition), cname] <- TRUE
}
# don't include controls because they are not batch-specific
for (i in 1:nrow(batch_complements.df)) {
  minus_conditions.df <- filter(msdata$msruns, batch==batch_complements.df$batch[[i]] & 
                                sample_type=="APMS" & bait_kind=="sample" & bait_homid != batch_complements.df$bait_homid[[i]]) %>% select(condition) %>% distinct()
  conditionXmetacondition.mtx[as.character(minus_conditions.df$condition), batch_complements.df$allminus_metacondition[[i]]] <- TRUE
}

# exclude S from the background of ACE2 and vice versa
conditionXmetacondition.mtx[as.character(filter(conditions.df, bait_id == "S" & bait_kind == "sample" & sample_type == "APMS")$condition),
                            str_detect(colnames(conditionXmetacondition.mtx), "allminus_ACE2")] <- FALSE
conditionXmetacondition.mtx["APMS_HUMAN_ACE2", str_detect(colnames(conditionXmetacondition.mtx), "allminus_S")] <- FALSE

conditionXmetacondition.mtx[as.character(filter(conditions.df, sample_type=="APMS" & bait_kind == "control")$condition), "controls"] <- TRUE
conditionXmetacondition.mtx[as.character(filter(conditions.df, sample_type=="FPMS")$condition), "FPMS"] <- TRUE

pheatmap(ifelse(conditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, 'plots', str_c(ms_folder, "_", fit_version),
                              paste0(project_id, "_metaconditions_", ms_folder, "_", fit_version, ".pdf")),
         width = 20, height = 12)

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n) %>%
  dplyr::left_join(dplyr::select(conditions.df, condition, bait_full_id, bait_id, bait_homid))

contrasts.df <- bind_rows(
  transmute(filter(conditions.df, sample_type=="APMS" & bait_kind == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = "controls",
            contrast_type = "filter"),
  transmute(filter(conditions.df, sample_type=="APMS" & bait_kind == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = str_c("allminus_", bait_homid),
            contrast_type = "filter"),
  inner_join(
    transmute(filter(conditions.df, sample_type == "APMS" & bait_kind == "sample"),
              metacondition_lhs = bait_full_id, orgcode_lhs = orgcode, bait_homid,
              bait_version = if_else(str_detect(bait_full_id, "\\?$"), "0", "1")),
    transmute(filter(conditions.df, sample_type == "APMS" & bait_kind == "sample"),
              metacondition_rhs = bait_full_id, orgcode_rhs = orgcode, bait_homid,
              bait_version = if_else(str_detect(bait_full_id, "\\?$"), "0", "1"))
  ) %>% filter(as.integer(orgcode_lhs) < as.integer(orgcode_rhs)) %>%
  mutate(contrast_type = "comparison") %>%
  select(-orgcode_lhs, -orgcode_rhs, -bait_homid),
  transmute(filter(msdata$msruns, sample_type=="APMS" & bait_kind == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = str_c("B", batch, "_allminus_", bait_homid),
            contrast_type = "filter") %>% distinct()
) %>%
mutate(contrast = str_c(metacondition_lhs, "_vs_", str_replace(metacondition_rhs, "allminus.+", "others")),
       metacondition_lhs = factor(metacondition_lhs, levels=all_metaconditions),
       metacondition_rhs = factor(metacondition_rhs, levels=all_metaconditions))

all_contrasts <- contrasts.df$contrast
contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
for (i in 1:nrow(contrasts.df)) {
  contrastXmetacondition.mtx[contrasts.df$contrast[[i]],
                             c(contrasts.df$metacondition_lhs[[i]],
                               contrasts.df$metacondition_rhs[[i]])] <- c(1, -1)
}
pheatmap(contrastXmetacondition.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, 'plots', str_c(ms_folder, "_", fit_version),
                              paste0(project_id, "_exp_design_contrasts_", ms_folder, "_", fit_version, ".pdf")),
         width = 18, height = 30)

contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::rename(weight = n) %>%
  dplyr::left_join(select(contrasts.df, contrast, contrast_type)) %>%
  dplyr::mutate(condition_role = if_else(contrast_type == "filter" & weight < 0, "background", "signal"),
                # define the rule to dynamically select control baits for _vs_others comparison based on the protein abundance
                # select the baits with abundance within 50%-90% percentiles
                quantile_min = if_else(condition_role == "background" & str_detect(contrast, "_vs_(B\\d+_)?others"), 0.5, NA_real_),
                quantile_max = if_else(condition_role == "background" & str_detect(contrast, "_vs_(B\\d+_)?others"), 0.9, NA_real_))

contrastXcondition.df <- dplyr::select(contrastXmetacondition.df, -starts_with("quantile")) %>%
  dplyr::inner_join(conditionXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition) %>%
  # don't remove from the background if top-enriched:
  #   a) the negative controls;
  #   b) the proteins 
  dplyr::mutate(is_preserved_condition = (weight < 0) & (contrast_type == "filter") &
                ((condition %in% c("APMS_Ctrl_NT", "APMS_Ctrl_Gaussia_luci", "APMS_Ctrl_SARS_CoV_NSP1?", "APMS_Ctrl_SARS_CoV_NSP3_macroD?")) |
                   # baits that have too many things in common, but when we dynamically exclude them,
                   # the remaining baits don't provide the background enough to remove the shared proteins
                 (str_detect(contrast, "SARS_CoV2?_M_vs_(B\\d+_)?others")) &
                    bait_homid %in% c("ORF3")) | #, "ORF7b")) |
                 (str_detect(contrast, str_c("(", str_c(str_replace(filter(conditions.df, bait_homid=="ORF3")$bait_full_id, fixed("?"), fixed("\\?")), collapse="|"),
                                             ")_vs_(B\\d+_)?others")) &
                    bait_homid %in% c("M")) | #, "ORF7b")) |
                 (str_detect(contrast, "SARS_CoV2?_ORF7b_vs_(B\\d+_)?others") &
                    bait_homid %in% c()))#"M", "ORF3"))))

# prepare the data to use for MS runs normalization
msrun_stats.df <- left_join(msdata_full$protgroup_intensities,
                            msdata_full$protgroup_idents) %>%
  dplyr::group_by(msrun) %>%
  dplyr::summarise(n_pg_idents = sum(!is.na(ident_type) & ident_type == "By MS/MS"),
                   n_pg_quants = sum(!is.na(intensity))) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(msdata$msruns)

# normalize using the intensities
msdata4norm.df <- msdata$pepmodstate_intensities %>%
  dplyr::filter(is_idented & !is.na(intensity)) %>%
  dplyr::semi_join(dplyr::filter(msdata$pepmods, !is_reverse & !is_contaminant & !is_viral) %>%
                   dplyr::inner_join(dplyr::select(msdata$pepmodstates, pepmodstate_id, pepmod_id))) %>%
  dplyr::select(pepmodstate_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$pepmodstate_intensities)

options(mc.cores=8)

# normalize experiments:
# 1) MS replicates for a given bait
# 2) all baits of the same batch
# 3) all batches together (keeping FPMS/APMS separate)
msruns_hnorm <- multilevel_normalize_experiments(instr_calib_pepmodstate,
    filter(msdata$msruns, is_used) %>%
    mutate(sample_type_batch = str_c(sample_type, "_B", batch),
           sample_type_batch_bait_full_id = str_c(sample_type_batch, "_", bait_full_id)),
    msdata4norm.df,
    quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
    mcmc.iter = 1000L,
    #mcmc.chains = 6,
    verbose=TRUE,
    norm_levels = list(msrun = list(cond_col = "msrun", max_objs=1000L, missing_exp.ratio=0.1),
                       bait_full_id = list(cond_col="sample_type_batch_bait_full_id", max_objs=500L, missing_exp.ratio=0.2),
                       batch = list(cond_col="sample_type_batch", condgroup_col="sample_type", max_objs=300L, missing_exp.ratio=0.2)
                       ))
msruns_hnorm$msruns_shifts <- msruns_hnorm$mschannel_shifts
total_msrun_shifts.df <- msruns_hnorm$msruns_shifts
setdiff(msruns.df$msrun, total_msrun_shifts.df$msrun)

# normalize (calculate prior_mean of similarity effect) homologous baits using shared specific interactions
prev_apms.env <- new.env(parent=baseenv())
load(file.path(scratch_path, paste0(project_id, '_msglm_fit_', 'mq_apms_20200427', '_', '20200503', '.RData')), envir=prev_apms.env)
load(file.path(scratch_path, paste0(project_id, '_msglm_data_', 'mq_apms_20200427', '_', '20200503', '.RData')), envir=prev_apms.env)
iactions4norm.df <- filter(prev_apms.env$object_contrasts.df, std_type == "median" & str_detect(contrast, "_vs_others")) %>%
    group_by(bait_id, contrast) %>%
    dplyr::filter(is_hit_nomschecks | row_number(prob_nonpos) <= 30) %>%
    ungroup()
msdata4hombait_norm.df <- inner_join(iactions4norm.df,
                                     semi_join(filter(msdata$msruns, is_used & sample_type=="APMS"),
                                               filter(prev_apms.env$msdata$msruns), by=c("raw_file", "bait_full_id"))) %>%
    inner_join(filter(msdata$protregroup2pepmod, is_specific)) %>%
    select(msrun, condition, bait_id, bait_full_id, pepmod_id) %>%
    inner_join(msdata4norm.df) %>%
    dplyr::select(condition, bait_id, bait_full_id, msrun, pepmodstate_id, intensity) %>%
    distinct()

hombait_norm <- normalize_experiments(instr_calib_pepmodstate, msdata4hombait_norm.df,
                                      quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
                                      cond_col='condition', condgroup_col='bait_id', 
                                      mcmc.iter = 2000L, stan_method="mcmc", mcmc.adapt_delta=0.9,
                                      verbose=TRUE, max_objs=1000L,
                                      mschan_preshifts = msruns_hnorm$msruns_shifts, mschan_shift.min = -2,
                                      nmschan_ratio.min = 0.5,
                                      preshift_col = 'total_msrun_shift')

baitXorgcode_effect_priors.df <- inner_join(hombait_norm, conditions.df) %>%
filter(abs(shift) >= shift_sd & n_objects >= 10) %>%
inner_join(dplyr::select(baits_info.df, bait_full_id, orgcode)) %>%
group_by(bait_id) %>%
mutate(has_cov2 = any(orgcode == "SARS2"),
       n_hombaits = n(),
       effect = str_c("bait_id", bait_id, if_else(orgcode == "SARS2", "", str_c(":orgcode", orgcode))),
       effect_prior = 0.5*(shift[orgcode == "SARS2"] - shift)) %>%
ungroup()

# correct priors of baitXvirus effects where possible
effects.df <- left_join(effects.df,
                        dplyr::filter(baitXorgcode_effect_priors.df, n_hombaits > 1 & has_cov2 & orgcode != "SARS2") %>%
                        dplyr::transmute(effect,# = factor(effect, levels=levels(effects.df$effect)),
                                         new_prior_mean = effect_prior)) %>%
  mutate(prior_mean = if_else(!is.na(new_prior_mean), new_prior_mean, prior_mean)) %>%
  select(-new_prior_mean) %>%
  mutate(effect = factor(effect, levels=effect))

# apply normalization
msdata_full$protgroup_intensities <- dplyr::left_join(msdata_full$protgroup_intensities,
                                                    dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift))

# use only native medium to calculate average intensity
global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns) %>%
                                               dplyr::select(msrun) %>% dplyr::distinct() %>%
                                               dplyr::inner_join(msdata_full$protgroup_intensities) %>% pull(intensity)), na.rm=TRUE)
global_pepmodstate_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns) %>%
                                                 dplyr::select(msrun) %>% dplyr::distinct() %>%
                                                 dplyr::inner_join(msdata_full$pepmodstate_intensities) %>% pull(intensity)), na.rm=TRUE)

#msdata$msruns <- dplyr::arrange(msdata$msruns, bait_kind, bait_id, orgcode) %>%
#  mutate(msrun_ix = row_number())

msdata_full$msrun_stats <- msrun_statistics(msdata_full)
set.seed(1232)
msdata_full$protgroup_intensities_all <- tidyr::expand(msdata_full$protgroup_intensities,
                                                protgroup_id, msrun) %>%
  left_join(dplyr::select(msdata_full$protgroup_intensities, -total_msrun_shift)) %>%
  dplyr::mutate(mstag = "Sum") %>%
  impute_intensities(msdata_full$msrun_stats) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
  #filter(mstag == "L") %>%
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca_agg.df <- inner_join(protgroup_intensities4pca.df,
                                               dplyr::select(msdata$msruns, msrun, batch, condition)) %>%
  dplyr::mutate(batch_condition = str_c("B", batch, "_", condition)) %>%
  dplyr::group_by(batch_condition, protgroup_id) %>%
  dplyr::summarise(intensity_imputed_norm = median(intensity_imputed_norm)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(log2_intensity_imputed_norm = log2(intensity_imputed_norm))

protgroup_intensities.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_norm),
                                    nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                    dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                    msrun = unique(protgroup_intensities4pca.df$msrun)))
protgroup_intensities_imp.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_imputed_norm),
                                    nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                    dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                    msrun = unique(protgroup_intensities4pca.df$msrun)))
protgroup_intensities_imp_agg.mtx <- matrix(log2(protgroup_intensities4pca_agg.df$intensity_imputed_norm),
                                        nrow = n_distinct(protgroup_intensities4pca_agg.df$protgroup_id),
                                        dimnames = list(protgroup = unique(protgroup_intensities4pca_agg.df$protgroup_id),
                                                        batch_condition = unique(protgroup_intensities4pca_agg.df$batch_condition)))

# PCA of msruns
msrun_intensities_pca <- stats::prcomp(protgroup_intensities_imp.mtx, scale. = TRUE)
msrun_intensities_pca.df <- as_tibble(msrun_intensities_pca$rotation,
                                      rownames="msrun") %>%
  dplyr::inner_join(msdata$msruns)
# UMAP of msruns
require(uwot)
msrun_intensities_umap2d = umap(t(protgroup_intensities_imp.mtx), n_components=2,
                 n_neighbors=10, init="laplacian",
                 min_dist=0.2, metric = "euclidean")
msrun_intensities_umap2d.df <- as_tibble(msrun_intensities_umap2d, .name_repair="minimal") %>%
  set_names(c("PC1", "PC2")) %>%
  mutate(msrun = colnames(protgroup_intensities_imp.mtx)) %>%
  dplyr::inner_join(msdata$msruns)

require(ggrepel)
p <- ggplot(msrun_intensities_umap2d.df,
       aes(x=PC1, y=PC2, color=bait_id)) +
    geom_point() +
    geom_text_repel(aes(label=str_remove(str_remove(msrun, "SARS_"), "APMS_")), vjust=-1.1,
                    box.padding = 0.1, force = 0.5) +
    theme_bw_ast(base_family = "", base_size = 6) +
    scale_color_discrete(guide="none")
ggsave(p, filename = file.path(analysis_path, 'plots', str_c(ms_folder, "_", fit_version),
                            paste0(project_id, "_msruns_umap_", ms_folder, "_", fit_version, ".pdf")),
       width = 30, height = 30, device = cairo_pdf)

require(pheatmap)
condition_hclu = hclust(dist(t(protgroup_intensities_imp_agg.mtx)))
msruns_ordered <- filter(msdata$msruns, msrun %in% colnames(protgroup_intensities_imp.mtx)) %>%
  mutate(batch_condition = factor(str_c("B", batch, "_", condition), levels=condition_hclu$labels[condition_hclu$order])) %>%
  dplyr::arrange(sample_type, batch_condition, replicate) %>%
  #dplyr::arrange(bait_kind, bait_id, batch, bait_full_id, replicate) %>%
  pull(msrun)
protgroup_hclu = hclust(dist(protgroup_intensities_imp.mtx))
pheatmap(log2(protgroup_intensities.mtx[, msruns_ordered]), cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version),
                          paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.pdf")), width=40, height=100)
pheatmap(log2(protgroup_intensities.mtx[, msruns_ordered]), cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version),
                          paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.png")), width=50, height=60)

msrunXbatchEffect_orig.mtx <- model.matrix(
  ~ 1 + batch + batch:sample_type,
  mutate(msdata$msruns, batch = factor(batch)))
batch_effects_mask <- str_detect(colnames(msrunXbatchEffect_orig.mtx), "\\(Intercept\\)|batch1", negate = TRUE)
msrunXbatchEffect.mtx <- msrunXbatchEffect_orig.mtx[, batch_effects_mask, drop=FALSE]
dimnames(msrunXbatchEffect.mtx) <- list(msrun = msdata$msruns$msrun,
                                        batch_effect = colnames(msrunXbatchEffect_orig.mtx)[batch_effects_mask])
pheatmap(msrunXbatchEffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
        file = file.path(analysis_path, "plots", str_c(ms_folder, "_", fit_version),
                         paste0(project_id, "_", ms_folder, "_", fit_version, "_batch_effects.pdf")), width=12, height=40)

batch_effects.df <- tibble(batch_effect=colnames(msrunXbatchEffect.mtx),
                           is_positive=FALSE,
                           prior_mean = 0.0)

# no subbatch effects so far
msrunXsubbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXbatchEffect.mtx),
                                        subbatch_effect = c())

subbatch_effects.df <- tibble(subbatch_effect=character(0),
                           is_positive=logical(0))

bait_checks_protgroup.df <- dplyr::left_join(dplyr::select(baits_info.df, bait_full_id, bait_id, bait_homid, bait_orgcode=orgcode, bait_organism=organism,
                                                           protein_ac = used_uniprot_ac),
                                             dplyr::select(msdata_full$proteins, protein_ac, protgroup_id, prot_organism=organism)) %>%
  dplyr::left_join(dplyr::select(dplyr::filter(msdata$protgroup_idents, ident_type=="By MS/MS"), protgroup_id, msrun)) %>%
  dplyr::left_join(dplyr::select(msdata$msruns, msrun, sample_type, observing_bait_full_id = bait_full_id)) %>%
  dplyr::arrange(bait_full_id, protgroup_id, desc(sample_type), msrun) %>%
  dplyr::group_by(bait_full_id, protgroup_id) %>%
  dplyr::mutate(idented_in_msruns = str_c(unique(msrun), collapse=";"),
                idented_in_AP_of = str_c(unique(observing_bait_full_id[sample_type=="APMS"]), collapse=";"),
                prot_organisms = str_c(unique(prot_organism), collapse=';')) %>%
  dplyr::filter(row_number()==1L) %>%
  dplyr::select(-msrun, -sample_type, -observing_bait_full_id) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(idented_in_msruns = if_else(idented_in_msruns == "", NA_character_, idented_in_msruns),
                idented_in_AP_of = if_else(idented_in_AP_of == "", NA_character_, idented_in_AP_of))

bait_checks_protregroup.df <- dplyr::select(baits_info.df, bait_full_id, bait_id, bait_homid,
                                            bait_orgcode=orgcode, bait_organism=organism) %>%
  # FIXME using protein_ac from fasta, because there's discrepancy between fasta and baits_info
  dplyr::left_join(dplyr::select(fasta.dfs$CoV, bait_id=gene_name, bait_organism=organism, protein_ac)) %>%
  dplyr::left_join(dplyr::select(filter(msdata_full$protein2protregroup, is_majority), protein_ac, protregroup_id)) %>%
  dplyr::left_join(dplyr::select(msdata_full$proteins, protein_ac, prot_organism=organism)) %>%
  dplyr::left_join(dplyr::select(dplyr::filter(msdata$protregroup_idents, ident_type=="By MS/MS"), protregroup_id, msrun)) %>%
  dplyr::left_join(dplyr::select(msdata$msruns, msrun, sample_type, observing_bait_full_id = bait_full_id)) %>%
  dplyr::arrange(bait_full_id, protregroup_id, desc(sample_type), msrun) %>%
  dplyr::group_by(bait_full_id, protregroup_id) %>%
  dplyr::mutate(idented_in_msruns = str_c(unique(msrun), collapse=";"),
                idented_in_AP_of = str_c(unique(observing_bait_full_id[sample_type=="APMS"]), collapse=";"),
                prot_organisms = str_c(unique(prot_organism), collapse=';')) %>%
  dplyr::filter(row_number()==1L) %>%
  dplyr::select(-msrun, -sample_type, -observing_bait_full_id) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(idented_in_msruns = if_else(idented_in_msruns == "", NA_character_, idented_in_msruns),
                idented_in_AP_of = if_else(idented_in_AP_of == "", NA_character_, idented_in_AP_of))

pepmodstate_labu_min <- inner_join(msdata$pepmodstate_intensities, total_msrun_shifts.df) %>%
  mutate(intensity_norm = intensity * exp(-total_msrun_shift)) %>%
  .$intensity_norm %>% log() %>%
  quantile(0.001, na.rm=TRUE) - global_pepmodstate_labu_shift - 5

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', ms_folder, '_', fit_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df,
     conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     instr_calib_protgroup, instr_calib_pepmodstate,
     global_protgroup_labu_shift, global_pepmodstate_labu_shift, pepmodstate_labu_min,
     msruns_hnorm, total_msrun_shifts.df,
     batch_effects.df, msrunXbatchEffect.mtx,
     subbatch_effects.df, msrunXsubbatchEffect.mtx,
     bait_checks_protgroup.df, bait_checks_protregroup.df,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full, fasta.dfs,
     #protgroup_stats.df,
     file = rfulldata_filepath)

message('Done.')
