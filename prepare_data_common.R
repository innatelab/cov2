require(readxl)

strlist_label <- function(strs) {
    str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
    sapply(str_split(strs, delim), strlist_label)
}

orgcodes = list("SARS-CoV-2" = "SARS2",
                "SARS-CoV" = "CVHSA",
                "SARS-CoV-GZ02" = "CVHSA",
                "HCoV-NL63" = "CVHNL",
                "HCoV-229E" = "CVH22",
                "ZIKV" = "ZIKV",
                "Gaussia" = "9MAXI",
                "HCV" = "9HEPC",
                "HUMAN" = "HUMAN")

baits_info.df <- read_xlsx(file.path(analysis_path, "data", "baits_info.xlsx"))
baits_info.df <- mutate(baits_info.df,
                        used_uniprot_ac = ifelse(!is.na(uniprot_ac), uniprot_ac, bait_id)) %>%
  filter(!is.na(bait_id)) %>%
  rename(organism = virus, bait_full_id = bait_id, bait_id = short_name, bait_kind=bait_type) %>%
  mutate(bait_id = relevel(factor(bait_id), "Ctrl_NT"),
         # bait homology id -- which baits are homologs to each other?
         # sometimes (ORF3) the names of the homologs in different strains is different
         bait_homid = case_when(bait_id %in% c("ORF3", "ORF3a", "ORF3b", "ORF4", "ORF4a") ~ "ORF3",
                                TRUE ~ as.character(bait_id)) %>%
                      factor() %>% relevel("Ctrl_NT"),
         bait_kind = factor(bait_kind, c("sample", "control")),
         orgcode = factor(orgcodes[organism], levels=unique(as.character(orgcodes))),
         protein_name = str_c(if_else(is.na(uniprot_ac), as.character(bait_full_id), uniprot_ac), "_", orgcode))

