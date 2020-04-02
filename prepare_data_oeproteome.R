# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200331"
fit_version <- "20200331"
ms_folder <- 'spectronaut_qc_20200331'
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

msdata_path <- file.path(data_path, ms_folder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  ms_folder = ms_folder,
                  instr_calib_filename = "instr_protgroup_LFQ_calib_scaturro_qep5calib_20161110_borg.json",
                  quant_type = "LFQ", quant_col_prefix = "LFQ_Intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib <- fromJSON(file = file.path(data_path, data_info$instr_calib_filename))$instr_calib

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(msdata_path, "fasta/cov_baits_20200326.fasta")),
  human = read_innate_uniprot_fasta(file.path(msdata_path, "fasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta")),
  contaminants = read_contaminants_fasta(file.path(msdata_path, "fasta/contaminants.fasta"))
)


msdata.wide <- read.Spectronaut.ProteinsReport(file.path(msdata_path, "20200331_115716_COVID19_test_directDIA_protein_Report.csv"),
                                               import_data = "quantity")
msdata_colgroups <- attr(msdata.wide, "column_groups")

msdata_full <- list(
  protgroups = msdata.wide[, msdata_colgroups$protgroup],
  protgroup_intensities = pivot_longer.Spectronaut.ProtgroupIntensities(msdata.wide)
)
msdata_full$msruns <- dplyr::select(msdata_full$protgroup_intensities, msrun_ix, raw_file) %>% dplyr::distinct() %>%
  dplyr::arrange(msrun_ix) %>%
  dplyr::mutate(msrun = str_remove(str_remove(raw_file, "^20200326_QX8_OzKa_"), ".raw$"))
msdata_full$msruns <- dplyr::mutate(msdata_full$msruns,
                                    condition = str_extract(msrun, "SA_COV2_FPMS_\\d+"),
                                    bait_code = str_extract(condition, "\\d+$")) %>%
  dplyr::left_join(select(baits_info.df, bait_code, bait_full_id, bait_id))
msdata_full$protgroup_intensities <- dplyr::select(msdata_full$protgroup_intensities, -raw_file)
msdata_full <- append_protgroups_info(msdata_full, msdata.wide,
                                      proteins_info = dplyr::bind_rows(
                                        dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE),
                                        dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE),
                                        dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE)),
                                      import_columns = c("is_viral", "is_contaminant"))
msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

strlist_label <- function(strs) {
  str_c(strs[[1]], if_else(n_distinct(strs) > 1, '...', ''))
}
strlist_label2 <- function(strs, delim=fixed(';')) {
  sapply(str_split(strs, delim), strlist_label)
}

msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
    is_reverse = FALSE,
    gene_label = strlist_label2(gene_names),
    protac_label = strlist_label2(protein_acs),
    protgroup_label = case_when(!is.na(gene_label) ~ gene_label,
                                !is.na(protac_label) ~ protac_label,
                                TRUE ~ str_c('#', protgroup_id)))

rdata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData'))
message('Saving full MS data to ', rdata_filepath, '...')
save(data_info, msdata_full,
     #protgroup_stats.df,
     file = rdata_filepath)

message('Done.')

require(FactoMineR)
require(Matrix)
dim(protgroup_intensities.mtx)
protgroup_intensities.mtx <- as.matrix(msdata.wide[str_detect(colnames(msdata.wide), "\\.PG\\.Quantity$")])
protgroup_intensities.mtx <- protgroup_intensities.mtx[rowSums(is.na(protgroup_intensities.mtx)) == 0, ]
colnames(protgroup_intensities.mtx) <- msdata_full$msruns$msrun
msrun_intensities_pca <- PCA(log2(protgroup_intensities.mtx), graph = FALSE)
msrun_intensities_pca.df <- as.data.frame(msrun_intensities_pca$svd$V)
colnames(msrun_intensities_pca.df) <- paste0("comp_", 1:ncol(msrun_intensities_pca.df))
msrun_intensities_pca.df <- dplyr::mutate(msrun_intensities_pca.df,
                                          msrun = rownames(msrun_intensities_pca$var$coord)) %>%
    dplyr::inner_join(msdata_full$msruns)

require(ggrepel)
cairo_pdf(filename = file.path(data_path, paste0(project_id, "_msruns_pca_OeProteome_QC_", fit_version, ".pdf")),
          width = 14, height = 14)
ggplot(msrun_intensities_pca.df,
       aes(x=comp_1, y=comp_2, color=bait_id)) +
    geom_point() +
    geom_text_repel(aes(label=str_remove(str_remove(msrun, "APMS_SARS_"), "APMS_")), vjust=-1.1) +
    theme_bw_ast(base_family = "", base_size = 10) #+
dev.off()
