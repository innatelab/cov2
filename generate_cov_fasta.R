# Loading and preparing PCP data
# 
# Author: Alexey Stukalov
###############################################################################

source('~/R/config.R')

project_id <- 'cov2'
version <- "20200326"

source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(misc_scripts_path, 'fasta_utils.R'))

require(Biostrings)
require(readxl)
require(readr)
require(stringr)
require(dplyr)
require(tidyr)

#baits_info.df <- read_xls(file.path(analysis_path, "data", "baits_info.xlsx"))
baits_info.df <- read_tsv(file.path(analysis_path, "data", "baits_info.txt"))
baits_info.df <- mutate(baits_info.df,
                        used_uniprot_ac = ifelse(!is.na(uniprot_ac), uniprot_ac, bait_id),
                        organism_code = case_when(virus == "SARS-CoV-2" ~ "CVHSA2",
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
Biostrings::writeXStringSet(bait_aaseqset, file.path(analysis_path, "data", "cov_baits_20200331.fasta"))

exp_design_template.df <- read_tsv(file.path(analysis_path, "data", "mq_apms_20200409", "experimentalDesignTemplate.txt"))

exp_design.df <- extract(exp_design_template.df,
                         Name, c("date", "instrument", "user", "sample_type", "project", "data_type", "bait_code", "replicate", "tech_replicate"),
                         "^([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)_([^_]+)(?:_([^_]+))?$", remove=FALSE) %>%
    left_join(select(baits_info.df, bait_type, bait_code, bait_id, ap_batch)) %>%
    mutate(ap_batch = case_when(bait_code %in% c("31", "32") ~ "2",
                                bait_type == "control" ~ NA_character_,
                                TRUE ~ ap_batch),
           replicate = parse_integer(replicate),
           batch_ix = case_when(bait_type == "control" ~ as_integer(ceiling(replicate/4)),
                                TRUE ~ parse_integer(ap_batch)),
           replicate_ix = case_when(bait_type == "control" ~ replicate - (batch_ix - 1L)*4L,
                                    TRUE ~ replicate)) %>%
    group_by(bait_id, batch_ix, replicate_ix) %>%
    mutate(tech_replicate_ix = if_else(n() > 1, row_number(), NA_integer_)) %>%
    ungroup() %>%
    mutate(Experiment = str_c("APMS_B", batch_ix, '_', bait_id, "_", replicate_ix,
                              if_else(is.na(tech_replicate_ix), "", str_c("#", tech_replicate_ix)))) %>%
    select(Name, Fraction, Experiment, PTM)

write_tsv(exp_design.df, path = file.path(analysis_path, "data", "mq_apms_20200409", "experimentalDesign.txt"), na = "")

bait_mapping.df <- read_tsv(file.path(analysis_path, "data", "BaitsIDMapping.txt"))

zikv_path = "/pool/scratch/astukalov/scaturro_zika"
fasta.dfs <- list(
    human = read_innate_uniprot_fasta(file.path(zikv_path, "fasta/UP000005640_9606.fasta")),
    human_additional = read_innate_uniprot_fasta(file.path(zikv_path, "fasta/UP000005640_9606_additional.fasta")),
    viruses = read_innate_uniprot_fasta(file.path(zikv_path, "fasta/KJ776791.2(ZIKV-HPF2013)_HCV_DV.fasta"))
)

ms_data$proteins <- bind_rows(fasta.dfs)
save(ms_data, file = file.path(zikv_path, "zikv_apms2_msglm_data_20171001_with_proteins.RData"))
