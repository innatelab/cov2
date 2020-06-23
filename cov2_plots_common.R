source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)

organism_type_map <- c("Gaussia" = "control", "HCoV-229E" = "HCoV", "HCoV-NL63" = "HCoV", "HCV" = "control",
                       "SARS-CoV" = "SARS-CoV", "SARS-CoV-2" = "SARS-CoV-2", "SARS-CoV-GZ02" = "SARS-CoV",
                       "SARS-CoV-2;SARS-CoV" = "SARS-CoV-2", # FIXME
                       "ZIKV" = "control", "Homo sapiens OX=9606" = "host")
organism_type_palette <- c("SARS-CoV" = "sienna",
                           "SARS-CoV-2" = "orange",
                           "HCoV" = "seagreen",
                           "host" = "black",
                           "control" = "gray")
