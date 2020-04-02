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

