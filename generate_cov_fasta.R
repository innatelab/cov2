# Loading and preparing PCP data
# 
# Author: Alexey Stukalov
###############################################################################

source('~/R/config.R')

project_id <- 'cov2'

source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(misc_scripts_path, 'fasta_utils.R'))

require(Biostrings)
require(readxl)
require(readr)
require(stringr)
require(dplyr)
require(tidyr)

baits_info.df <- read_xlsx(file.path(analysis_path, "data", "baits_info.xlsx"))
baits_info.df <- mutate(baits_info.df,
                        used_uniprot_ac = ifelse(!is.na(uniprot_ac), uniprot_ac, bait_id),
                        organism_code = case_when(virus == "SARS-CoV-2" ~ "SARS2",
                                                  virus == "SARS-CoV" ~ "CVHSA",
                                                  virus == "SARS-CoV-GZ02" ~ "CVHSA",
                                                  virus == "HCoV-NL63" ~ "CVHNL",
                                                  virus == "HCoV-229E" ~ "CVH22",
                                                  virus == "ZIKV" ~ "ZIKV",
                                                  virus == "Gaussia" ~ "9MAXI",
                                                  virus == "HCV" ~ "9HEPC",
                                                  TRUE ~ "UNKNOWN"),
                        protein_name = str_c(if_else(is.na(uniprot_ac), short_name, uniprot_ac), "_", organism_code),
                        fasta_header = str_c("sp|", used_uniprot_ac, "|", protein_name,
                                             if_else(is.na(Description), " NO DESCRIPTION", str_c(" ", Description)),
                                             " OS=", virus, " GN=", short_name)
)
aaseqs <- baits_info.df$aa_sequence
names(aaseqs) <- baits_info.df$fasta_header

bait_aaseqset <- AAStringSet(aaseqs[!is.na(aaseqs)])
Biostrings::writeXStringSet(bait_aaseqset, file.path(analysis_path, "data", "msfasta", "cov_baits_20200415.fasta"))

exp_design_template.df <- read_tsv(file.path(analysis_path, "data", "mq_apms_20200510", "combined", "experimentalDesignTemplate.txt"))

samples_to_use_wide.df <- read_tsv(file.path(data_path, "samples_to_use.txt"))
samples_to_use.df <- pivot_longer(samples_to_use_wide.df, cols = starts_with("batch."),
                                  names_prefix = "batch.", names_to="batch", values_to = "is_measured") %>%
    filter(coalesce(is_measured, "") == "x") %>%
    mutate(is_used = str_detect(organism, "^Ctrl") | batch == coalesce(batch_used, 0))

exp_design.df <- extract(exp_design_template.df,
                         Name, c("date", "instrument", "user", "msrun_type", "project", "sample_type",
                                 "batch_fname", "calib_sample", "bait_code", "replicate", "tech_replicate", "remes_time"),
                         "^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)(?:_B(\\d+))?(?:_(\\d+uL(?:_\\d+)?))?_([^_]+)_([^_]+)(?:_([^_]{1,2}))?(?:_([^_]{4,}))?$", remove=FALSE) %>%
    left_join(select(baits_info.df, bait_kind=bait_type, bait_code, bait_full_id = bait_id, batch_tbl = ap_batch)) %>%
    mutate(bait_kind = if_else(is.na(bait_kind) & bait_code == "ACE2", "sample", bait_kind),
           bait_full_id = case_when(is.na(bait_full_id) & bait_code == "ACE2" ~ "HUMAN_ACE2",
                                    TRUE ~ bait_full_id),
           replicate = case_when(bait_code %in% c("31", "32") & !is.na(tech_replicate) ~ tech_replicate, # fix strange naming scheme
                                 TRUE ~ replicate),
           replicate = parse_integer(replicate),
           # FIX batch
           batch = case_when(!is.na(batch_fname) ~ batch_fname,
                             date == "20200425" ~ "5",
                             sample_type == "FPMS" & date %in% c("20200423", "20200424") ~ "1", # noise calibration based on B1 proteomes
                             (bait_code %in% "ACE2") & (date == "20200409") ~ "4", # correct wrong batch date in the filename for ACE2_4
                             (bait_code %in% c("31", "32")) & (date == "20200321") ~ "2",
                             (bait_code %in% c("31", "32")) & (date == "20200416") ~ "4",
                             ((bait_kind == "control") | is.na(batch_tbl)) & (date == "20200409") ~ "3",
                             ((bait_kind == "control") | is.na(batch_tbl)) & (date == "20200416") ~ "4",
                             (bait_kind == "control") & (date == "20200321") ~ as.character(ceiling(replicate/4)),
                             TRUE ~ batch_tbl),
           batch_ix = parse_integer(batch),
           replicate_ix = case_when(bait_kind == "control" ~ 1L + (replicate-1L) %% 4L,
                                    TRUE ~ replicate)) %>%
    left_join(select(samples_to_use.df, bait_code, batch, is_used)) %>%
    mutate(bait_full_id_tbl = bait_full_id,
           bait_full_id = if_else(is_used, bait_full_id_tbl,
                                  str_replace(bait_full_id_tbl, "(CoV[^_]?)_(.+)$", "\\1_\\2?"))) %>% #rename earlier baits so that they are separate in the analysis
    group_by(sample_type, bait_full_id, batch_ix, replicate_ix) %>%
    mutate(tech_replicate_ix = if_else(rep.int(n() > 1L, n()), row_number(), rep.int(NA_integer_, n()))) %>%
    ungroup() %>%
    mutate(Experiment = str_c(sample_type, "_B", batch_ix, '_', bait_full_id, "_", replicate_ix,
                              if_else(is.na(tech_replicate_ix), "", str_c("#", tech_replicate_ix))))# %>%

setdiff(samples_to_use.df$sample, exp_design.df$sample)

write_tsv(select(filter(exp_design.df, is_used),
                 Name, Fraction, Experiment, PTM),
          path = file.path(analysis_path, "data", "mq_apms_20200524", "experimentalDesign.txt"), na = "")

bait_mapping.df <- read_tsv(file.path(analysis_path, "data", "BaitsIDMapping.txt"))

zikv_path = "/pool/scratch/astukalov/scaturro_zika"
fasta.dfs <- list(
    human = read_innate_uniprot_fasta(file.path(zikv_path, "fasta/UP000005640_9606.fasta")),
    human_additional = read_innate_uniprot_fasta(file.path(zikv_path, "fasta/UP000005640_9606_additional.fasta")),
    viruses = read_innate_uniprot_fasta(file.path(zikv_path, "fasta/KJ776791.2(ZIKV-HPF2013)_HCV_DV.fasta"))
)

ms_data$proteins <- bind_rows(fasta.dfs)
save(ms_data, file = file.path(zikv_path, "zikv_apms2_msglm_data_20171001_with_proteins.RData"))
