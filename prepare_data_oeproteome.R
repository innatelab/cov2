# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200923"
fit_version <- "20200923"
msfolder <- 'snaut_oefp_20200923'
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
                  mscalib_protgroup_filename = "instr_QX8_intensity_protgroup_calib_cov2_20200519.json",
                  mscalib_pepmodstate_filename = "mscalib_QX8_intensity_pepmodstate_cov2_20200923.json",
                  quant_type = "intensity", quant_col_prefix = "intensity",
                  pep_quant_type = "intensity",
                  qvalue_max = 0.001,
                  empty_observation_sigmoid_scale=1/3)

message('Loading MS instrument calibration data from ', data_info$mscalib_pepmodstate_filename, '...')
mscalib_protgroup <- fromJSON(file = file.path(data_path, data_info$mscalib_protgroup_filename))$instr_calib
mscalib_pepmodstate <- fromJSON(file = file.path(data_path, data_info$mscalib_pepmodstate_filename))$mscalib

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/cov_baits_20200415.fasta")),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/Human_July2019_with_isoforms_only_swissprot.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants_20200405.fasta"))
)

require(readxl)
require(tidyr)

samples_to_use_wide.df <- read_tsv(file.path(data_path, "samples_to_use.txt"))
samples_to_use.df <- pivot_longer(samples_to_use_wide.df, cols = starts_with("batch."),
                                  names_prefix = "batch.", names_to="batch", values_to = "is_measured") %>%
  filter(coalesce(is_measured, "") == "x") %>%
  mutate(sample = str_c("B", batch, "_", str_remove(sample, "B\\d+_")),
         is_used = str_detect(organism, "^Ctrl") | batch == coalesce(batch_used, 0))

msruns_info.df <- read_xlsx(file.path(msdata_path, "SARS_COV2_MS_samples_library.xlsx"),
                            sheet = "samples", range="A1:E146") %>%
  rename(msrun_sn = `ms run name`) %>%
  extract(sample, c('bait_code', 'replicate'), "([^_]+)_(\\d+)", remove = FALSE, convert=TRUE) %>%
  full_join(filter(tidyr::expand(., nesting(bait_code, batch), replicate), replicate <= 4)) %>%
  mutate(is_control = bait_code %in% LETTERS,
         batch = parse_integer(str_remove(batch, "^B"))) %>%
  # fix expansion issues
  filter(!(is_control & ((batch %in% c(2, 4)) & replicate <= 4))) %>%
  mutate(batch = parse_integer(str_remove(batch, "^B")),
         is_wrong_replicate = ((batch == 3) & is_control) | ((batch == 2) & (bait_code == "D") & (replicate == 5)),
         replicate = if_else(is_wrong_replicate, (batch-1L)*4L + replicate, replicate),
         sample = if_else(is.na(sample) | is_wrong_replicate, str_c(bait_code, "_", replicate), sample),
         msrun_sn = str_remove(msrun_sn, "^_SA_"),
         msrun_sn = if_else(is.na(msrun_sn) | is_wrong_replicate, str_c("COV2_FPMS_", sample), msrun_sn),
         is_batch5 = batch == 5) %>%
  arrange(batch, bait_code, replicate) %>%
  dplyr::group_by(bait_code) %>%
  dplyr::mutate(is_older_batch = !is_control & batch < max(batch)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(transmute(samples_to_use.df,
                             bait_code=case_when(bait_code %in% c("31","32") & batch=="4" ~ str_c(bait_code, ".2"),
                                                 TRUE ~ bait_code),
                             batch=parse_integer(batch), is_used))

# fix duplication in the proteins report
#msdata_orig.wide <- read_tsv(file.path(msdata_path, "oe_proteome_unnormalized_20200519.txt"))
#msdata_orig_noeg.wide <- distinct(select(msdata_orig.wide, -ends_with("EG.Qvalue")))
#write_tsv(path = file.path(msdata_path, "oe_proteome_unnormalized_20200519_dedup.txt"),
#          distinct(msdata_orig_noeg.wide))

setdiff(msruns_info.df$msrun_sn, msdata_full$msruns$msrun_sn)
setdiff(msdata_full$msruns$msrun_sn, msruns_info.df$msrun_sn)

bad_msruns <- c("FPMS_B5_Ctrl_ZIKV_Capsid_4", "FPMS_B2_SARS_CoV_ORF8b_2",
                "FPMS_B1_SARS_CoV2_ORF3_3", "FPMS_B1_SARS_CoV2_ORF7a_2",
                "FPMS_B2_SARS_CoV2_NSP1_3", "FPMS_B3_SARS_CoV2_NSP2_2",
                "FPMS_B1_SARS_CoV_ORF6_1", "FPMS_B2_Ctrl_ZIKV_Capsid_5",
                "FPMS_B5_SARS_CoV2_NPS14_1")

pgdata.wide <- read.Spectronaut.ProteinsReport(file.path(msdata_path, "PG_Report.xls"),
                                               import_data = "quantity", delim='\t')
pgdata_colgroups <- attr(pgdata.wide, "column_groups")

pmsdata.wide <- read.Spectronaut.PepmodstatesReport(file.path(msdata_path, "EG_Report.xls"),
                                                    import_data = c("eg_quantity", "eg_qvalue", "eg_snratio", "eg_normfactor"), delim='\t')
pmsdata_colgroups <- attr(pmsdata.wide, "column_groups")

msdata_full <- list(
    protgroups = pgdata.wide[, pgdata_colgroups$protgroup],
    pepmodstates = dplyr::select(pmsdata.wide[, pmsdata_colgroups$pepmodstate],
                                 pepmodstate_id, pepmod_seq, charge, peptide_seq, is_pg_specific),
    protgroup_intensities = pivot_longer.Spectronaut.Metrics(pgdata.wide, c("protgroup_sn_id"), c("quantity", "runstats", "pg_normfactor")) %>%
          dplyr::distinct() %>% dplyr::filter(!is.na(intensity)),
    pepmodstate_intensities = pivot_longer.Spectronaut.Metrics(pmsdata.wide, c("pepmodstate_id"), c("eg_quantity", "eg_qvalue", "eg_snratio", "eg_normfactor")) %>%
          dplyr::distinct() %>% dplyr::filter(!is.na(intensity))
)
msdata_full$msruns <- dplyr::select(msdata_full$pepmodstate_intensities, rawfile_ix, rawfile, msrun_normfactor) %>% dplyr::distinct() %>%
  dplyr::arrange(rawfile_ix) %>%
  tidyr::extract(rawfile, c("msdate", "instrument", "msrun_sn"), "^(\\d+)_([^_]+)_OzKa_SA_(.+)\\.raw", remove=FALSE) %>%
  dplyr::mutate(msrun_sn = str_remove(msrun_sn, "(?:_\\d{3,})$") %>%
                str_replace("31_2_", "31.2_") %>% str_replace("32_2_", "32.2_") %>%
                str_replace("([A-G])(\\d)$", "\\1_\\2")) %>%
  tidyr::extract(msrun_sn, c("batch5", "bait_code_orig", "replicate"), "^COV2_FPMS_(B5_)?([^_]+)_(\\d+)$", remove=FALSE) %>%
  dplyr::mutate(bait_code = str_remove(bait_code_orig, "\\.\\d+$")) %>%
  dplyr::mutate(is_batch5 = (coalesce(batch5, "") != "") | (msdate == '20200426'),
                replicate = parse_integer(replicate)) %>%
  left_join(dplyr::select(baits_info.df, bait_kind, bait_code, bait_full_id, bait_id, bait_homid, organism, orgcode)) %>%
  left_join(dplyr::select(msruns_info.df, msrun_sn, batch, is_batch5, is_used)) %>%
  dplyr::mutate(is_used = coalesce(is_used, FALSE), # removes LMix
         bait_full_id = str_c(bait_full_id, ifelse(!is_used, "?", "")),
         batch = case_when(is.na(batch) & bait_code_orig == "31.2" ~ 4L,
                           TRUE ~ batch),
         # fix bait id for the candidates to throw away, homid is not touched
         bait_id = str_c(bait_id, ifelse(!is_used, "?", "")),
         condition = str_c("FPMS_", bait_full_id),
         msrun = str_c("FPMS_B", batch, "_", bait_full_id, "_", replicate),
         is_used = is_used & !(msrun %in% bad_msruns))
write_tsv(select(msdata_full$msruns, rawfile, msrun, is_used),
          path = file.path(data_path, msfolder, str_c(msfolder, "_", data_version, ".txt")))
msdata_full$msruns <- dplyr::filter(msdata_full$msruns, is_used)

msdata_full$protgroups <- mutate(msdata_full$protgroups,
                                 npeptides = NA_integer_,
                                 npeptides_unique = NA_integer_,
                                 npeptides_unique_razor = NA_integer_)
msdata_full$protgroup_intensities <- dplyr::left_join(msdata_full$protgroup_intensities,
                                               dplyr::select(msdata_full$msruns, rawfile_ix, msrun, msrun_normfactor)) %>%
  dplyr::inner_join(dplyr::select(msdata_full$protgroups, protgroup_sn_id, protgroup_id)) %>%
  dplyr::mutate(intensity_norm = intensity,
                intensity = intensity_norm / msrun_normfactor,
                ident_type = factor("By MS/MS", levels = c("By MS/MS", "By matching"))) %>% #factor(if_else(nevidences > 0L, "By MS/MS", "By matching"), levels = c("By MS/MS", "By matching")))
  dplyr::select(-rawfile_ix, -rawfile, -protgroup_sn_id)

msdata_full$pepmodstate_intensities <- dplyr::mutate(msdata_full$pepmodstate_intensities,
                                             intensity_norm = intensity,
                                             intensity = intensity_norm/msrun_normfactor,
                                             msrun_normfactor = NULL)
msdata_full$pepmods <- dplyr::select(msdata_full$pepmodstates, pepmodstate_id, pepmod_seq, peptide_seq, is_pg_specific) %>%
    dplyr::left_join(msdata_full$pepmodstate_intensities) %>%
    dplyr::group_by(pepmod_seq, peptide_seq, is_pg_specific) %>%
    dplyr::summarise(qvalue = min(qvalue, na.rm=TRUE), .groups="drop") %>%
    dplyr::distinct() %>%
    dplyr::mutate(pepmod_id = row_number() - 1L)
msdata_full$peptides <- dplyr::group_by(msdata_full$pepmods, peptide_seq, is_pg_specific) %>%
    dplyr::summarise(qvalue = min(qvalue, na.rm=TRUE),
                     n_pepmods = n_distinct(pepmod_id),
                     .groups="drop") %>%
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
    dplyr::select(-rawfile_ix) %>%
    dplyr::mutate(ident_type = "MULTI-MSMS") # HACK

all_proteins.df <- dplyr::bind_rows(
    dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE),
    dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE),
    dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE))

msdata_full <- append_protgroups_info(msdata_full, pgdata.wide,
                                      proteins_info = all_proteins.df,
                                      import_columns = c("is_viral", "is_contaminant", "organism"))
msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
    is_reverse = FALSE,
    is_contaminant = coalesce(is_contaminant, TRUE),
    is_viral = coalesce(is_viral, FALSE),
    is_used = q_value <= 0.001,
    gene_label = strlist_label2(gene_names),
    protac_label = strlist_label2(protein_acs),
    protein_label = strlist_label2(protein_names),
    protgroup_label = case_when(is_viral ~ protein_label,
                                !is.na(gene_label) ~ gene_label,
                                !is.na(protac_label) ~ protac_label,
                                TRUE ~ str_c('#', protgroup_id)))

msdata_full$msrun_protgroup_stats <- msrun_statistics(msdata_full) %>%
  dplyr::mutate(na_ratio = n_missing / n)
msdata_full$msrun_pepmodstate_stats <- msrun_statistics(msdata_full, obj="pepmodstate")

# prepare complete intensities matrix for early PCA
set.seed(1232)
msdata_full$protgroup_intensities_all <- tidyr::expand(semi_join(msdata_full$protgroup_intensities, filter(msdata_full$msruns, is_used)),
                                                       protgroup_id, msrun) %>%
  left_join(msdata_full$protgroup_intensities) %>%
  dplyr::mutate(mstag = "Sum") %>%
  impute_intensities(msdata_full$msrun_protgroup_stats) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::arrange(msrun, protgroup_id)

pg_intensities_imp.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_imputed),
                                 nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                 dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                 msrun = unique(protgroup_intensities4pca.df$msrun)))
# PCA of msruns
require(uwot)
pg_intensities_umap = umap(pg_intensities_imp.mtx / rowMedians(pg_intensities_imp.mtx), n_components=3,
                           n_neighbors=20, init="laplacian",
                           min_dist=0.2, metric = "euclidean")

pg_umap.df <- as_tibble(pg_intensities_umap, .name_repair="minimal") %>%
  set_names(c("x_3d", "y_3d", "z_3d")) %>%
  mutate(protgroup_id = parse_integer(rownames(pg_intensities_imp.mtx))) %>%
  dplyr::inner_join(select(msdata$protgroups, protgroup_id, protgroup_label,
                           majority_protein_acs, gene_names, protein_descriptions,
                           is_viral, is_contaminant, is_reverse))
# find the batch-specific point on umap
pg_umap_batch_center.df <- dplyr::filter(pg_umap.df, str_detect(gene_names, "^(GALK1|PNP|PGLS)$")) %>%
  dplyr::summarise(x_3d_center = median(x_3d),
                   y_3d_center = median(y_3d),
                   z_3d_center = median(z_3d),
                   .groups="drop")
pg_umap.df <- mutate(pg_umap.df,
                     batch_distance = sqrt((x_3d - pg_umap_batch_center.df$x_3d_center)^2 +
                                           (y_3d - pg_umap_batch_center.df$y_3d_center)^2 +
                                           (z_3d - pg_umap_batch_center.df$z_3d_center)^2),
                     is_batch_affected = batch_distance <= 2)

# redefine protein groups (protregroups)
peptides.df <- dplyr::select(msdata_full$peptides, peptide_id, peptide_seq, qvalue) %>%
    dplyr::left_join(dplyr::transmute(dplyr::inner_join(dplyr::filter(msdata_full$pepmodstate_intensities,
                                                                      coalesce(qvalue, 1.0) <= data_info$qvalue_max),
                                                        msdata_full$msruns) %>%
                                      dplyr::inner_join(msdata_full$pepmods),
                                      peptide_id = peptide_id,
                                      has_data = TRUE) %>% distinct()) %>%
    dplyr::left_join(dplyr::select(pmsdata.wide, peptide_seq, protein_acs) %>% dplyr::distinct()) %>%
    dplyr::mutate(peptide_rank = if_else(coalesce(has_data, FALSE), 1L, -1L))
save(file = file.path(msdata_path, str_c(project_id, "_", msfolder, '_', data_version, "_peptides.RData")),
     peptides.df, all_proteins.df)
# .. run protregroup_oefp_proteome.jl
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

msdata_full$protregroup2peptide <- bind_rows(
  select(msdata_full$protregroups, protregroup_id, peptide_id=spec_peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = TRUE),
  select(msdata_full$protregroups, protregroup_id, peptide_id=peptide_ids) %>%
    separate_rows(peptide_id, sep=fixed(";"), convert=TRUE) %>%
    mutate(is_specific = FALSE)) %>%
  dplyr::group_by(protregroup_id, peptide_id) %>%
  dplyr::summarise(is_specific = any(is_specific)) %>%
  dplyr::ungroup()

msdata_full$protregroup2pepmod <- inner_join(msdata_full$protregroup2peptide,
                                             msdata_full$pepmods) %>%
  dplyr::select(-peptide_id) %>% dplyr::distinct()
  
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
  tidyr::extract(.name, c("mstag", "msrun"), remove=FALSE,
          str_c("^", data_info$quant_col_prefix, "\\.(\\S+)\\s(\\S+)")) %>%
  mutate(.value = str_c("intensity.", mstag)) %>%
  dplyr::inner_join(select(msdata_full$msruns, msrun, rawfile))

msdata_full$protregroup_idents <- dplyr::inner_join(msdata_full$protregroup2pepmod, msdata_full$pepmodstate_intensities) %>%
  dplyr::mutate(is_idented = !is.na(intensity) & coalesce(qvalue, 1.0) <= data_info$qvalue_max) %>%
  dplyr::group_by(msrun, protregroup_id) %>%
  dplyr::summarise(npepmods_quanted = sum(!is.na(intensity)),
                   nspecpepmods_quanted = sum(!is.na(intensity) & is_specific),
                   npepmods_idented = sum(is_idented),
                   nspecpepmods_idented = sum(is_idented & is_specific)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ident_type = factor(if_else(nspecpepmods_idented > 0, "By MS/MS", "By matching")))

msdata <- msdata_full[c('protgroup_intensities', 'pepmodstate_intensities', 'pepmods', 'pepmodstates',
                        'msruns', 'protgroups', 'protregroups', 'protein2protregroup',
                        'peptide2protein', 'protregroup2pepmod',
                        'msrun_pepmodstate_stats', 'msrun_protgroup_stats')]
msdata$msruns <- filter(msdata$msruns, is_used)
msdata$protgroup_intensities <- filter(semi_join(msdata$protgroup_intensities, msdata$msruns), !is.na(intensity))

# condition = bait
conditions.df <- dplyr::select(msdata$msruns, condition, bait_full_id, bait_id, bait_homid, bait_kind, orgcode) %>%
  dplyr::distinct() %>%
  dplyr::arrange(bait_kind, bait_homid, bait_id, orgcode) %>%
  dplyr::mutate(condition = relevel(factor(condition, levels=condition), "FPMS_Ctrl_NT"),
                bait_id = relevel(factor(bait_id, levels=unique(bait_id)), "Ctrl_NT"))
msdata$msruns <- mutate(msdata$msruns,
                        condition = factor(condition, levels=levels(conditions.df$condition)))

# setup experimental design matrices
require(Matrix)
require(purrr)
# not grouping baits by bait_homid, since we are only interested in SARS_vs_SARS2 comparisons
bait_exp_designs <- conditions.df %>%
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
            arrange(orgcode)
          res <- cbind(matrix(rep.int(1, nrow(bait_conds.df)), ncol = 1),
                       -contr.sum(nrow(bait_conds.df)))
          dimnames(res) <- list(condition = as.character(bait_conds.df$condition),
                                effect = bait_conds.df$effect)
          return(res)
      }))
      mtx <- as.matrix(bdiag(res$mtx))
      dimnames(mtx) <- list(condition = flatten_chr(map(res$mtx, rownames)),
                            effect = flatten_chr(map(res$mtx, colnames)))
      return(mtx)
   } else {
      conds.df <- mutate(conds.df, bait_id=relevel(factor(bait_id), "Ctrl_NT"))
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

conditionXeffect_orig.mtx <- bait_exp_design.mtx
conditionXeffect.mtx <- conditionXeffect_orig.mtx[as.character(conditions.df$condition),
                                                  colSums(abs(conditionXeffect_orig.mtx)) != 0 &
                                                    !str_detect(colnames(conditionXeffect_orig.mtx), "\\(Intercept\\)|bait_idCtrl.+:orgcode")]
dimnames(conditionXeffect.mtx) <- list(condition = as.character(conditions.df$condition),
                                       effect = colnames(conditionXeffect.mtx))
plots_path <- file.path(analysis_path, 'plots', str_c(msfolder, "_", fit_version))
if (!dir.exists(plots_path)) { dir.create(plots_path) }
pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(plots_path, paste0(project_id, "_exp_design_", msfolder, "_", fit_version, ".pdf")),
         width = 14, height = 14)

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
  dplyr::mutate(prior_tau = case_when(effect_type == "baitXvirus" ~ 0.1,
                                      effect_type == "bait" ~ 0.1,
                                      TRUE ~ 5.0))
# prior_mean for baitXvirus is calculated as the shift between the homologous baits at the normalization step?

conditionXeffect.df <- conditionXeffect_frame(conditionXeffect.mtx, effects.df)

batch_complements.df <- distinct(select(filter(msdata$msruns, bait_kind == "sample"), batch, bait_homid)) %>%
  dplyr::mutate(allminus_metacondition = paste0("B", batch, "_allminus_", bait_homid))

bait_conditions <- as.character(filter(conditions.df, bait_kind == "sample")$bait_full_id)
allminus_metaconditions <- paste0("allminus_", unique(filter(conditions.df, bait_kind == "sample")$bait_homid))
compound_metaconditions <- c(allminus_metaconditions, "controls",
                             batch_complements.df$allminus_metacondition)
all_metaconditions <- c(bait_conditions, compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in bait_conditions) {
  conditionXmetacondition.mtx[str_c("FPMS_", cname), cname] <- TRUE
}
for (cname in allminus_metaconditions) {
  hombait <- str_remove(cname, "^allminus_")
  minus_conditions.df <- filter(conditions.df, ((!str_detect(bait_id, "\\?$") & bait_homid != hombait) | bait_kind == "control"))
  conditionXmetacondition.mtx[as.character(minus_conditions.df$condition), cname] <- TRUE
}
# don't include controls because they are not batch-specific
for (i in 1:nrow(batch_complements.df)) {
  minus_conditions.df <- filter(msdata$msruns, batch==batch_complements.df$batch[[i]] & 
                                bait_kind=="sample" & bait_homid != batch_complements.df$bait_homid[[i]]) %>% select(condition) %>% distinct()
  conditionXmetacondition.mtx[as.character(minus_conditions.df$condition), batch_complements.df$allminus_metacondition[[i]]] <- TRUE
}

# don't include viral controls into "controls" metacondition
conditionXmetacondition.mtx[as.character(filter(conditions.df, bait_full_id %in% c("Ctrl_NT", "Ctrl_Gaussia_luci"))$condition), "controls"] <- TRUE
pheatmap(ifelse(conditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(plots_path, paste0(project_id, "_metaconditions_", msfolder, "_", fit_version, ".pdf")),
         width = 22, height = 16)

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n) %>%
  dplyr::left_join(dplyr::select(conditions.df, condition, bait_full_id, bait_id, bait_homid))

hombait_contrasts.df <-   inner_join(
  transmute(filter(conditions.df, bait_kind == "sample"),
            metacondition_lhs = bait_full_id, orgcode_lhs = orgcode, bait_homid,
            bait_version = if_else(str_detect(bait_full_id, "\\?$"), "0", "1")),
  transmute(filter(conditions.df, bait_kind == "sample"),
            metacondition_rhs = bait_full_id, orgcode_rhs = orgcode, bait_homid,
            bait_version = if_else(str_detect(bait_full_id, "\\?$"), "0", "1"))
) %>% filter(as.integer(orgcode_lhs) < as.integer(orgcode_rhs)) %>%
  mutate(contrast_type = "comparison") %>%
  select(-orgcode_lhs, -orgcode_rhs, -bait_homid)

contrasts.df <- bind_rows(
  transmute(filter(conditions.df, bait_kind == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = "controls",
            contrast_type = "comparison",
            has_offset=FALSE),
  transmute(filter(conditions.df, bait_kind == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = str_c("allminus_", bait_homid),
            contrast_type = "comparison",
            has_offset=FALSE),
  transmute(filter(msdata$msruns, bait_kind == "sample"),
            metacondition_lhs = bait_full_id,
            metacondition_rhs = str_c("B", batch, "_allminus_", bait_homid),
            contrast_type = "comparison",
            has_offset=FALSE) %>% distinct(),
  mutate(hombait_contrasts.df, has_offset=FALSE),
  mutate(hombait_contrasts.df, has_offset=TRUE),
) %>%
mutate(contrast = str_c(metacondition_lhs, "_vs_", str_replace(metacondition_rhs, "allminus.+", "others"),
                        ifelse(has_offset, "_corrected", "")),
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
         filename = file.path(plots_path, paste0(project_id, "_exp_design_contrasts_", msfolder, "_", fit_version, ".pdf")),
         width = 18, height = 30)

contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
  dplyr::rename(weight = n) %>%
  dplyr::left_join(dplyr::select(contrasts.df, contrast, contrast_type)) %>%
  dplyr::mutate(condition_role = "signal",
                is_lhs = weight > 0,
                # define the rule to dynamically select control baits for _vs_others comparison based on the protein abundance
                # select the baits with abundance within 50%-90% percentiles
                quantile_min = dplyr::if_else(!is_lhs & str_detect(contrast, "_vs_(B\\d+_)?others"), 0.2, NA_real_),
                quantile_max = dplyr::if_else(!is_lhs & str_detect(contrast, "_vs_(B\\d+_)?others"), 0.8, NA_real_))

contrastXcondition.df <- dplyr::select(contrastXmetacondition.df, -starts_with("quantile")) %>%
  dplyr::inner_join(conditionXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition) %>%
  dplyr::mutate(is_preserved_condition = FALSE) # nothing preserved

# prepare the data to use for MS runs normalization
msrun_stats.df <- msdata$protgroup_intensities %>%
  dplyr::group_by(msrun) %>%
  dplyr::summarise(n_pg_idents = sum(!is.na(ident_type) & ident_type == "By MS/MS"),
                   n_pg_quants = sum(!is.na(intensity))) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(msdata$msruns)

# normalize using the intensities
msdata4norm.df <- msdata$pepmodstate_intensities %>%
  dplyr::filter(#npeptides_quanted > 1 & 
    !is.na(intensity)) %>%
  dplyr::semi_join(dplyr::filter(msdata$protregroups, !is_contaminant & !is_viral & !is_reverse) %>%
                   dplyr::inner_join(dplyr::filter(msdata$protregroup2pepmod, is_specific)) %>%
                   # remove proteins affected by batch effect
                   dplyr::inner_join(dplyr::select(msdata_full$protein2protregroup, protregroup_id, protein_ac)) %>%
                   dplyr::inner_join(msdata_full$protein2protgroup) %>%
                   dplyr::anti_join(dplyr::select(dplyr::filter(pg_umap.df, is_batch_affected), protgroup_id))) %>%
  dplyr::select(pepmodstate_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$pepmodstate_intensities)

options(mc.cores=8)

# normalize experiments:
# 1) MS replicates for a given bait
# 2) all baits of the same batch
# 3) all batches together
msruns_hnorm <- multilevel_normalize_experiments(mscalib_pepmodstate,
                                                 filter(msdata$msruns, is_used) %>%
                                                   mutate(batch=as.character(batch),
                                                          batch_bait_full_id = str_c("B", batch, "_", bait_full_id)),
                                                 msdata4norm.df,
                                                 quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
                                                 mcmc.iter = 1000L,
                                                 #mcmc.chains = 6,
                                                 verbose=TRUE,
                                                 #mschan_preshifts_df = msruns_hnorm$levels$msrun$mschannel_shifts,
                                                 #mschan_preshift_col = "msrun_shift",
                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=1000L, missing_exp.ratio=0.1),
                                                                    bait_full_id = list(cond_col="batch_bait_full_id", max_objs=500L, missing_exp.ratio=0.2),
                                                                    batch = list(cond_col="batch", max_objs=200L, missing_exp.ratio=0.3)
                                                 ))

#save(msruns_hnorm, file=str_c("oeproteome_msruns_hnorm_", fit_version, ".RData"))
#load(file=str_c("oeproteome_msruns_hnorm_", fit_version, ".RData"))

# ignore all higher levels of normalization
msruns_hnorm_old <- msruns_hnorm
msruns_hnorm <- msruns_hnorm_old
msruns_hnorm2 <- msruns_hnorm
msruns_hnorm2$levels$msrun <- msruns_hnorm$levels$msrun
msruns_hnorm2$levels$msrun$level_shifts <- semi_join(msruns_hnorm$levels$msrun$level_shifts, dplyr::select(msdata$msruns, msrun)) %>%
  group_by(batch_bait_full_id) %>%
  mutate(old_shift = shift,
         avg_shift = mean(shift),
         shift = shift - avg_shift) %>%
  ungroup()
msruns_hnorm2$levels$bait_full_id$level_shifts <- left_join(msruns_hnorm$levels$bait_full_id$level_shifts,
                                                            dplyr::select(msruns_hnorm2$levels$msrun$level_shifts, batch_bait_full_id, avg_msrun_shift = avg_shift) %>% distinct()) %>%
  mutate(old_shift = shift,
         shift = shift + avg_msrun_shift) %>%
  group_by(batch) %>%
  mutate(avg_shift = mean(shift),
         shift = shift - avg_shift) %>%
  ungroup()
msruns_hnorm2$levels$batch$level_shifts <- left_join(msruns_hnorm$levels$batch$level_shifts,
                                                     dplyr::select(msruns_hnorm2$levels$bait_full_id$level_shifts, batch, avg_bait_shift = avg_shift) %>% distinct()) %>%
  mutate(old_shift = shift,
         shift = shift + avg_bait_shift) %>%
  mutate(avg_shift = mean(shift),
         shift = shift - avg_shift) %>%
  ungroup()
nrow(msruns_hnorm2$levels$msrun$level_shifts)
View(msruns_hnorm2$levels$msrun$level_shifts)
View(msruns_hnorm2$levels$bait_full_id$level_shifts)
View(msruns_hnorm2$levels$batch$level_shifts)
msruns_hnorm2$msruns_shifts <- left_join(
  dplyr::select(msruns_hnorm2$levels$msrun$level_shifts, msrun, batch_bait_full_id, msrun_shift = shift),
  dplyr::select(msruns_hnorm2$levels$bait_full_id$level_shifts, batch_bait_full_id, batch, bait_full_id_shift = shift)) %>%
left_join(dplyr::select(msruns_hnorm2$levels$batch$level_shifts, batch, batch_shift = shift)) %>%
  dplyr::mutate(total_msrun_shift = msrun_shift + bait_full_id_shift + batch_shift)
msruns_hnorm <- msruns_hnorm2

msruns_hnorm$msruns_shifts <- msruns_hnorm$mschannel_shifts
total_msrun_shifts.df <- mutate(msruns_hnorm$msruns_shifts,
                                total_msrun_shift = msrun_shift + bait_full_id_shift)
setdiff(msruns.df$msrun, total_msrun_shifts.df$msrun)

# normalize (calculate offset for contrasts) homologous baits using shared specific interactions
prev_oeprot.env <- new.env(parent=baseenv())
load(file.path(scratch_path, paste0(project_id, '_msglm_fit_', 'spectronaut_oeproteome_20200527', '_', '20200527', '.RData')), envir=prev_oeprot.env)
load(file.path(scratch_path, paste0(project_id, '_msglm_data_', 'spectronaut_oeproteome_20200527', '_', '20200527', '.RData')), envir=prev_oeprot.env)
iactions4norm.df <- filter(prev_oeprot.env$object_contrasts.df, std_type == "median" & str_detect(contrast, "_vs_controls")) %>%
    dplyr::left_join(dplyr::select(baits_info.df, bait_full_id, bait_id)) %>%
    group_by(bait_full_id, contrast) %>%
    dplyr::filter((abs(median_log2) >= 0.5 & p_value <= 1E-3) | row_number(pmin(prob_nonpos, prob_nonneg)) <= 30) %>%
    ungroup()
msdata4hombait_norm.df <- inner_join(iactions4norm.df,
                                     select(prev_oeprot.env$msdata$msruns, msrun, rawfile=raw_file, bait_full_id)) %>%
    inner_join(select(prev_oeprot.env$msdata$protein2protgroup, object_id=protgroup_id, protein_ac)) %>%
    select(protein_ac, rawfile) %>% distinct() %>%
    inner_join(select(msdata$msruns, msrun, rawfile, condition, bait_id, bait_full_id)) %>%
    inner_join(select(msdata$protein2protregroup, protein_ac, protregroup_id)) %>%
    inner_join(select(dplyr::filter(msdata$protregroup2pepmod, is_specific), protregroup_id, pepmod_id)) %>%
    inner_join(select(msdata$pepmodstates, pepmod_id, pepmodstate_id)) %>%
    inner_join(msdata4norm.df) %>%
    dplyr::select(condition, bait_id, bait_full_id, msrun, protregroup_id, pepmodstate_id, intensity) %>%
    distinct()

hombait_norm <- normalize_experiments(mscalib_pepmodstate, msdata4hombait_norm.df,
                                      quant_col = "intensity", obj_col = "pepmodstate_id", mschan_col = "msrun",
                                      cond_col='condition', condgroup_col='bait_id', 
                                      mcmc.iter = 2000L, stan_method="mcmc", mcmc.adapt_delta=0.9,
                                      verbose=TRUE, max_objs=1000L,
                                      mschan_preshifts = msruns_hnorm$msruns_shifts, mschan_shift.min = -2,
                                      nmschan_ratio.min = 0.5,
                                      preshift_col = 'total_msrun_shift')

hombait_contrast_offsets.df <- inner_join(hombait_norm, conditions.df) %>%
filter(abs(shift) >= shift_sd & n_objects >= 5) %>%
inner_join(dplyr::select(baits_info.df, bait_full_id, orgcode)) %>%
mutate(is_cov2 = orgcode == "SARS2") %>%
group_by(bait_id) %>%
summarise(has_cov2 = any(is_cov2),
          n_hombaits = n(),
          n_objects = min(n_objects),
          contrast = str_c(bait_full_id[is_cov2], "_vs_", bait_full_id[!is_cov2], "_corrected"),
          offset = shift[!is_cov2] - shift[is_cov2], # compensate the average difference in binding affinities,
       ) %>%
ungroup()

# correct contrasts for homologous baits comparison where possible
contrasts.df <- left_join(contrasts.df,
                          dplyr::filter(hombait_contrast_offsets.df, n_hombaits == 2 & has_cov2) %>%
                          dplyr::transmute(contrast, new_offset = offset)) %>%
  mutate(offset = if_else(is.na(new_offset), 0.0, new_offset)) %>%
  select(-new_offset)

# apply normalization
msdata_full$protgroup_intensities <- dplyr::left_join(dplyr::select(msdata_full$protgroup_intensities, -any_of("total_msrun_shift")),
                                                      dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift))

# use only native medium to calculate average intensity
global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns) %>%
                                                 dplyr::select(msrun) %>% dplyr::distinct() %>%
                                                 dplyr::inner_join(msdata_full$protgroup_intensities) %>% pull(intensity)), na.rm=TRUE)
global_pepmodstate_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns) %>%
                                                 dplyr::select(msrun) %>% dplyr::distinct() %>%
                                                 dplyr::inner_join(msdata_full$pepmodstate_intensities) %>% pull(intensity)), na.rm=TRUE)

set.seed(1232)
msdata_full$protgroup_intensities_all <- tidyr::expand(semi_join(msdata_full$protgroup_intensities, filter(msdata_full$msruns, is_used)),
                                                protgroup_id, msrun) %>%
  left_join(dplyr::select(msdata_full$protgroup_intensities, -any_of("total_msrun_shift"))) %>%
  dplyr::mutate(mstag = "Sum") %>%
  impute_intensities(msdata_full$msrun_protgroup_stats) %>%
  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, protgroup_id)

pg_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
  #filter(mstag == "L") %>%
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant)) %>%
  dplyr::arrange(msrun, protgroup_id)

pg_intensities4pca_agg.df <- inner_join(pg_intensities4pca.df,
                                        dplyr::select(msdata$msruns, msrun, batch, condition)) %>%
  dplyr::mutate(batch_condition = str_c("B", batch, "_", condition)) %>%
  dplyr::group_by(batch_condition, protgroup_id) %>%
  dplyr::summarise(intensity_imputed_norm = median(intensity_imputed_norm)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(log2_intensity_imputed_norm = log2(intensity_imputed_norm)) %>%
  dplyr::arrange(batch_condition, protgroup_id)

pg_intensities.mtx <- matrix(log2(pg_intensities4pca.df$intensity_norm),
                             nrow = n_distinct(pg_intensities4pca.df$protgroup_id),
                             dimnames = list(protgroup = unique(pg_intensities4pca.df$protgroup_id),
                                             msrun = unique(pg_intensities4pca.df$msrun)))
pg_intensities_imp.mtx <- matrix(log2(pg_intensities4pca.df$intensity_imputed_norm),
                                 nrow = n_distinct(pg_intensities4pca.df$protgroup_id),
                                 dimnames = list(pg = unique(pg_intensities4pca.df$protgroup_id),
                                                 msrun = unique(pg_intensities4pca.df$msrun)))
pg_intensities_imp_agg.mtx <- matrix(log2(pg_intensities4pca_agg.df$intensity_imputed_norm),
                                     nrow = n_distinct(pg_intensities4pca_agg.df$protgroup_id),
                                     dimnames = list(protgroup = unique(pg_intensities4pca_agg.df$protgroup_id),
                                                     batch_condition = unique(pg_intensities4pca_agg.df$batch_condition)))

# PCA of msruns
msrun_intensities_pca <- stats::prcomp(pg_intensities_imp.mtx, scale. = TRUE)
msrun_intensities_pca.df <- as_tibble(msrun_intensities_pca$rotation[, 1:20],
                                      rownames="msrun") %>%
  dplyr::left_join(dplyr::select(msdata$msruns, -starts_with("PC")))

batch_palette <- c("control" = "gray", "B1" = "deepskyblue", "B2" = "darkorange", "B3" = "darkseagreen", "B4" = "gold", "B5" = "firebrick")
bad_bait_ids <- c("SARS_CoV_ORF9b", "SARS_CoV2_NSP12", "SARS_CoV_NSP12", "SARS_CoV_ORF7a", "HCoV_ORF3", "SARS_CoV2_ORF8")

require(ggrepel)
p <- ggplot(msrun_intensities_pca.df,
            aes(x=PC1, y=PC2,
                color=bait_full_id %in% bad_bait_ids
                #color=str_c("B", batch)
            )) +
  geom_point() +
  geom_text_repel(aes(label=str_remove(str_remove(msrun, "SARS_"), "FPMS_")), vjust=-1.1,
                  box.padding = 0.1, force = 0.5) +
  #scale_color_manual(values=batch_palette) +
  theme_bw_ast(base_family = "", base_size = 6)
ggsave(p, filename = file.path(plots_path, str_c(project_id, "_msruns_pca_", msfolder, "_", fit_version, ".pdf")),
       width = 20, height = 20, device = cairo_pdf)

# UMAP of msruns
require(uwot)
msrun_intensities_umap2d = umap(t(protgroup_intensities_imp.mtx), n_components=2,
                 n_neighbors=10, init="laplacian",
                 min_dist=0.5, metric = "euclidean")
msrun_intensities_umap2d.df <- as_tibble(msrun_intensities_umap2d, .name_repair="minimal") %>%
  set_names(c("x_2d", "y_2d")) %>%
  mutate(msrun = colnames(protgroup_intensities_imp.mtx)) %>%
  dplyr::inner_join(msdata$msruns)

require(ggrepel)
p <- ggplot(msrun_intensities_umap2d.df,
       aes(x=x_2d, y=y_2d,
           color=str_c("B", batch)
           #color=bait_full_id %in% bad_bait_ids
       )) +
    geom_point() +
    geom_text_repel(aes(label=str_remove(str_remove(msrun, "SARS_"), "FPMS_")),
                    box.padding = 0.1, force = 0.5) +
    scale_color_manual(values=batch_palette) +
    theme_bw_ast(base_family = "", base_size = 6)
ggsave(p, filename = file.path(plots_path, str_c(project_id, "_msruns_umap_", msfolder, "_", fit_version, ".pdf")),
       width = 25, height = 25, device = cairo_pdf)

require(matrixStats)
protgroup_intensities_imp_norm.mtx <- protgroup_intensities_imp.mtx / rowMedians(protgroup_intensities_imp.mtx)

require(uwot)
pg_intensities_umap2d = umap(protgroup_intensities_imp_norm.mtx, n_components=2,
                             n_neighbors=20, init="laplacian",
                             min_dist=0.2, metric = "euclidean")

objects_umap.df <- as_tibble(pg_intensities_umap2d, .name_repair="minimal") %>%
  set_names(c("x_2d", "y_2d")) %>%
  mutate(protgroup_id = parse_integer(rownames(protgroup_intensities_imp_norm.mtx)),
         object_id = protgroup_id) %>%
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
                              str_detect(gene_names, "(^|;)(ALDOA|DERA|ENO1|ENO2|ENO3|G6PD|GALK1|GAPDH|GPI|IDH1|KYNU|LDHA|MDH1|ME1|NAMPT|NAXD|NNMT|PDXK|PFKL|PFKM|PFKP|PGAM1|PGD|PGK1|PGLS|PGM1|PGM2|PGM2L1|PKM|PNP|PNPO|PSAT1|RBKS|TALDO1|TKT|TPI1)(;|$)") ~ "pyridine",
                              #str_detect(gene_names, "(^|;)(?:HIST[123])?H[123][ABCDEF]") ~ "histone",
                              #str_detect(gene_names, "(^|;)PPP[123456]") ~ "phosphatase",
                              str_detect(gene_names, "(^|;)RP[SL]\\d+") ~ "ribosome",
                              str_detect(gene_names, "(^|;)MRP[SL]\\d+") ~ "mito ribosome",
                              str_detect(gene_names, "(^|;)PSM[ABCD]\\d+") ~ "proteasome",
                              str_detect(gene_names, "TGF") ~ "TGF",
                              str_detect(gene_names, "(^|;)(MAVS|OAS\\d+|IFIT\\d+|PARP\\d+|DDX58|IFI\\d+|ISG\\d+|MX\\d+|UNC93B1|ZC3HAV1)(;|$)") ~ "ISG",
                              TRUE ~ "default"
  ))# %>%
#dplyr::left_join(filter(protregroup_effects.df, effect == "treatmentSC35M"))
# FIXME duplicate point due to N-to-N protgroup to protregroup matching
#dplyr::left_join(filter(object_effects_wide.df, std_type == "replicate") %>%
#                  select(-gene_names, -majority_protein_acs))
require(plotly)
# plot umap + selected proteins (color = protein biological function)
object_categories_palette <-  c("stable" = "black", "pyridine" = "red", "TGF" = "violet", "viral" = "orange", "ISG" = "blue", "hit" = "orange",
                                "ribosome" = "darkgreen", "mito ribosome" = "darkred", "proteasome" = "khaki",
                                "histone" = "orange", "phosphatase" = "violet",
                                "default" = "gray",
                                "NA" = "pink", "reverse" = "pink", "contaminant" = "yellow")
umap <- ggplotly(
  ggplot(objects_umap.df,
         aes(x=x_2d, y=y_2d, color=category)) +
    geom_point(aes(text=object_label, size=!(category %in% c('default', 'NA')))) +
    #geom_text(data=filter(objects_umap.df, gene_names %in% c("RPL11", "RPL12")),
    #          aes(label=object_label), vjust=-1.1, color = "black", size = 3) +
    #geom_text_repel(data=filter(objects_umap.df, is_viral),
    #                aes(label=object_label)) +
    scale_size_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5)) +
    scale_color_manual(values = object_categories_palette) +
    theme_bw_ast(base_family = "", base_size = 10),
  tooltip = "text"
)
umap

p <- ggplot(filter(msdata_full$protgroups, is_viral) %>% #str_detect(gene_names, "(^|;)(SERPIN)")) %>%
              inner_join(msdata_full$protgroup_intensities_all) %>%
              left_join(msdata_full$msruns),
            aes(x = timepoint_num, y = intensity_imputed_norm,
                color=treatment, fill=treatment)) +
  geom_smooth(alpha=0.25) +
  geom_point(aes(shape=is_imputed)) +
  scale_x_continuous(breaks = unique(msdata_full$msruns$timepoint_num)) +
  scale_shape_manual(values = c("FALSE"=16L, "TRUE"=1L)) +
  scale_color_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  scale_fill_manual(values=c("mock"="gray", "SARS_COV2"="red")) +
  facet_wrap(~ gene_names, ncol = 3, scales = "free_y") +
  theme_bw_ast()
p
ggsave(p, filename = file.path(plots_path, paste0(project_id, "_", msfolder, "_SARS_CoV2_proteins_", data_version, ".pdf")),
       width = 12, height = 7, device = cairo_pdf)


p <- ggplot(objects_umap.df,
            aes(x=x_2d, y=y_2d)) +
  geom_point(aes(color=category), size=0.5) +
  geom_text_repel(data=filter(objects_umap.df, category=="pyridine"),
                  aes(label=object_label), vjust=-1.1,
                  box.padding = 0.1, force = 0.5) +
  theme_bw_ast(base_family = "", base_size = 6) +
  scale_color_manual(values = object_categories_palette)
ggsave(p, filename = file.path(plots_path, paste0(project_id, "_protgroups_umap_", msfolder, "_", fit_version, ".pdf")),
       width = 20, height = 20, device = cairo_pdf)

require(pheatmap)
condition_hclu = hclust(dist(t(protgroup_intensities_imp_agg.mtx)))
msruns_ordered <- filter(msdata$msruns, msrun %in% colnames(protgroup_intensities_imp.mtx)) %>%
  mutate(batch_condition = factor(str_c("B", batch, "_", condition), levels=condition_hclu$labels[condition_hclu$order])) %>%
  dplyr::arrange(batch_condition, replicate) %>%
  #dplyr::arrange(bait_kind, bait_id, batch, bait_full_id, replicate) %>%
  pull(msrun)
protgroup_hclu = hclust(dist(protgroup_intensities_imp.mtx / rowMedians(protgroup_intensities_imp.mtx)))
pheatmap(log2(pmax(pmin(protgroup_intensities.mtx[, msruns_ordered]/ rowMedians(protgroup_intensities_imp.mtx),
                   quantile(protgroup_intensities.mtx/ rowMedians(protgroup_intensities_imp.mtx), 0.999, na.rm=TRUE)),
                   quantile(protgroup_intensities.mtx/ rowMedians(protgroup_intensities_imp.mtx), 0.001, na.rm=TRUE))),
         cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(plots_path, paste0(project_id, "_", msfolder, "_", data_version, "_heatmap_intensity.pdf")), width=40, height=100)
pheatmap(log2(pmax(pmin(protgroup_intensities.mtx[, msruns_ordered]/rowMedians(protgroup_intensities_imp.mtx),
                        quantile(protgroup_intensities.mtx/rowMedians(protgroup_intensities_imp.mtx), 0.99, na.rm=TRUE)),
                   quantile(protgroup_intensities.mtx/ rowMedians(protgroup_intensities_imp.mtx), 0.01, na.rm=TRUE))),
         cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(plots_path, paste0(project_id, "_", msfolder, "_", fit_version, "_heatmap_intensity.png")), width=50, height=60)

protgroup_mask <- rownames(protgroup_intensities_imp.mtx) %in% filter(msdata$protgroups, 
  str_detect(gene_names, "(^|;)(ALDOA|DERA|ENO1|ENO2|ENO3|G6PD|GALK1|GAPDH|GPI|IDH1|KYNU|LDHA|MDH1|ME1|NAMPT|NAXD|NNMT|PDXK|PFKL|PFKM|PFKP|PGAM1|PGD|PGK1|PGLS|PGM1|PGM2|PGM2L1|PKM|PNP|PNPO|PSAT1|RBKS|TALDO1|TKT|TPI1)(;|$)") |
  str_detect(gene_names, "(^|;)(M?RP[LS]\\d+|PSM[ABCD]\\d+)(;|$)")
)$protgroup_id
sum(protgroup_mask)

protgroup_intensities_imp.submtx <- protgroup_intensities_imp.mtx[protgroup_mask, ] / rowMedians(protgroup_intensities_imp.mtx[protgroup_mask, ])
submtx_rownames <- dplyr::select(protgroup_intensities4pca.df, protgroup_id) %>%
  dplyr::distinct() %>% dplyr::left_join(dplyr::select(msdata_full$protgroups, protgroup_id, protgroup_label)) %>%
  dplyr::mutate(is_used = protgroup_mask) %>% dplyr::filter(is_used) %>% .$protgroup_label
rownames(protgroup_intensities_imp.submtx) <- submtx_rownames

protgroup_subhclu = hclust(dist(protgroup_intensities_imp.submtx))
protgroup_intensities.submtx <- protgroup_intensities.mtx[protgroup_mask, msruns_ordered] / rowMedians(protgroup_intensities.mtx[protgroup_mask, ])
rownames(protgroup_intensities.submtx) <- submtx_rownames

pheatmap(log2(pmax(pmin(protgroup_intensities.submtx,
                        quantile(protgroup_intensities.submtx, 0.99, na.rm=TRUE)),
                   quantile(protgroup_intensities.submtx, 0.01, na.rm=TRUE))),
         cluster_cols=FALSE, cluster_rows=protgroup_subhclu,
         file = file.path(plots_path, paste0(project_id, "_", msfolder, "_", data_version, "_subheatmap_intensity.png")), width=50, height=30)

filter(msdata_full$protgroups, str_detect(protgroup_label, "JUN"))

ggplot(filter(msdata_full$protgroups, str_detect(protgroup_label, "^CLUH")) %>%
       inner_join(msdata_full$protgroup_intensities_all) %>%
         inner_join(select(msdata_full$msruns, -organism)) %>%
         mutate(batch_condition = str_c("B", batch, "_", condition),
                batch=str_c("B", batch))) +
  geom_point(aes(x = batch_condition, y = intensity_imputed_norm, shape=is_imputed, color=batch),
             position=position_jitter(width=0.5)) +
  scale_y_log10() +
  scale_shape_manual(values=c("TRUE"=1, "FALSE"=16)) +
  theme_bw_ast() +
  theme(axis.text.x = element_text(angle=-45, vjust=0, hjust=0)) +
  facet_wrap(~ protgroup_label + protgroup_id, ncol=1, scales = "free")

#msdata$msruns <- left_join(msdata$msruns,
#                           as_tibble(msrun_intensities_pca$rotation,
#                                     rownames="msrun") %>% dplyr::select(msrun, PC2))

msrunXbatchEffect_orig.mtx <- model.matrix(
  ~ 1 + batch,
  mutate(msdata$msruns, batch = factor(batch)))
batch_effects_mask <- str_detect(colnames(msrunXbatchEffect_orig.mtx), "\\(Intercept\\)", negate = TRUE)

msrunXprcompEffect.mtx <- matrix(2*msdata$msruns$PC2 / median(abs(msdata$msruns$PC2)),
                                 dimnames = list(msruns = msdata$msruns$msrun,
                                                 effects = "PC2"))
# with PC2
msrunXbatchEffect.mtx <- cbind(msrunXbatchEffect_orig.mtx[, batch_effects_mask, drop=FALSE],
                               msrunXprcompEffect.mtx)
dimnames(msrunXbatchEffect.mtx) <- list(msrun = msdata$msruns$msrun,
                                        batch_effect = c(colnames(msrunXbatchEffect_orig.mtx)[batch_effects_mask], 
                                                         colnames(msrunXprcompEffect.mtx)))
# without PC2
msrunXbatchEffect.mtx <- msrunXbatchEffect_orig.mtx[, batch_effects_mask, drop=FALSE]
dimnames(msrunXbatchEffect.mtx) <- list(msrun = msdata$msruns$msrun,
                                        batch_effect = colnames(msrunXbatchEffect_orig.mtx)[batch_effects_mask])

pheatmap(msrunXbatchEffect.mtx[arrange(msdata$msruns, batch, condition, replicate)$msrun,], cluster_rows=FALSE, cluster_cols=FALSE,
        file = file.path(plots_path, paste0(project_id, "_", msfolder, "_", fit_version, "_batch_effects.pdf")), width=12, height=40)

batch_effects.df <- tibble(batch_effect=colnames(msrunXbatchEffect.mtx),
                           is_positive=str_detect(batch_effect, "^PC"),
                           prior_mean = 0.0)

# no subbatch effects so far
msrunXsubbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXbatchEffect.mtx),
                                        subbatch_effect = c())

subbatch_effects.df <- tibble(subbatch_effect=character(0),
                              is_positive=logical(0))

bait_checks_protgroup.df <- dplyr::left_join(dplyr::select(baits_info.df, bait_full_id, bait_id, bait_homid, bait_orgcode=orgcode, bait_organism=organism,
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

protgroup_labu_min <- inner_join(msdata$protgroup_intensities, total_msrun_shifts.df) %>%
  mutate(intensity_norm = intensity * exp(-total_msrun_shift)) %>%
  .$intensity_norm %>% log() %>%
  quantile(0.001, na.rm=TRUE) - global_protgroup_labu_shift - 5

pepmodstate_labu_min <- inner_join(msdata$pepmodstate_intensities, total_msrun_shifts.df) %>%
  mutate(intensity_norm = intensity * exp(-total_msrun_shift)) %>%
  .$intensity_norm %>% log() %>%
  quantile(0.001, na.rm=TRUE) - global_pepmodstate_labu_shift - 5

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', msfolder, '_', fit_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df, contrasts.df,
     conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     mscalib_protgroup, mscalib_pepmodstate,
     protgroup_labu_min, pepmodstate_labu_min,
     global_protgroup_labu_shift, global_pepmodstate_labu_shift,
     msruns_hnorm, hombait_norm, total_msrun_shifts.df,
     batch_effects.df, msrunXbatchEffect.mtx,
     subbatch_effects.df, msrunXsubbatchEffect.mtx,
     bait_checks_protgroup.df,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', msfolder, '_', data_version, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full, fasta.dfs,
     file = rfulldata_filepath)

message('Done.')
