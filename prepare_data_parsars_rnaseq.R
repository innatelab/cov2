# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20201020"
fit_version <- "20201020"
data_folder <- 'parsars_rnaseq_20201020'
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

pubdata_path <- "/pool/pub3rdparty"

# ncbi_genes.df <- read_tsv(file.path(pubdata_path, "ncbi_gene2accession_20201018.gz"), na=c("", "-"),
#                           col_types = cols(
#                               `#tax_id` = "i", GeneID = "i",
#                               status = "c",
#                               RNA_nucleotide_accession.version = "c",
#                               RNA_nucleotide_gi = "c",
#                               protein_accession.version = "c",
#                               protein_gi = "c",
#                               genomic_nucleotide_accession.version = "c",
#                               genomic_nucleotide_gi = "i",
#                               start_position_on_the_genomic_accession = "c",
#                               end_position_on_the_genomic_accession = "c",
#                               orientation = "c",
#                               assembly = "c",
#                               mature_peptide_accession.version = "c",
#                               mature_peptide_gi = "c",
#                               Symbol = "c"
#                           ))
#human_ncbi_genes.df <- dplyr::filter(ncbi_genes.df, `#tax_id` == 9606) %>%
#    tidyr::extract(protein_accession.version, c("protein_ac", "protein_ac_ver"), "^(.+)\\.(\\d+)$", remove=FALSE, convert=TRUE)
#write_tsv(human_ncbi_genes.df, file=file.path(pubdata_path, "ncbi_human_gene2accession_20201018.gz"), na = "-")
human_ncbi_genes.df <- read_tsv(file=file.path(pubdata_path, "ncbi_human_gene2accession_20201018.gz"), na = "-")

fasta.dfs <- list(
    CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/SARS_CoV_20200928.fasta")),
    CoV2 = read_innate_uniprot_fasta(file.path(data_path, "msfasta/SARS_CoV2_20200928.fasta")),
    human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta"))#,
    #contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants_20200405.fasta"))
)

geneid2proteinac.df <- transmute(human_ncbi_genes.df, GeneID, protein_acXver=protein_accession.version) %>%
  dplyr::filter(!is.na(protein_acXver)) %>%
  tidyr::extract("protein_acXver", c("protein_ac", "protein_ac_ver"), "^([^.]+)(?:\\.(\\d+))?", convert=FALSE, remove=FALSE) %>%
  dplyr::select(GeneID, protein_ac) %>% dplyr::distinct() %>%
  dplyr::semi_join(dplyr::select(fasta.dfs$human, protein_ac))

genename2entrez.df <- bind_rows(
    transmute(human_ncbi_genes.df, gene_name=genomic_nucleotide_accession.version, entrez_id=GeneID, genename_source = "NCBI_AC"),
    transmute(human_ncbi_genes.df, gene_name=Symbol, entrez_gene=GeneID, genename_source = "NCBI_Symbol")
) %>% dplyr::filter(!is.na(gene_name)) %>% dplyr::distinct()

hgnc_genes.df <- read_tsv(file.path(pubdata_path, "hgnc_20201018.txt")) %>%
    dplyr::rename(gene_name = `Approved symbol`,
                  protein_ac = `UniProt ID(supplied by UniProt)`,
                  hgnc_id = `HGNC ID`)

ensembl2entrez.df <- read_tsv(file.path(pubdata_path, "ensembl_to_entrez_20201018.txt")) %>%
    dplyr::rename(gene_name = `Gene name`,
                  entrez_id = `NCBI gene (formerly Entrezgene) ID`)
hgnc2entrez.df <- read_tsv(file.path(pubdata_path, "hgnc_to_entrez_20201018.txt")) %>%
    dplyr::rename(gene_name = `Approved symbol`,
                  entrez_id = `NCBI Gene ID`)

genename2entrez.df <- bind_rows(
    dplyr::transmute(ensembl2entrez.df, gene_name, entrez_id, source="Ensembl"),
    dplyr::transmute(hgnc2entrez.df, gene_name, entrez_id, source="HGNC"),
) %>% dplyr::group_by(gene_name) %>%
    dplyr::filter(row_number() == 1L) %>%
    dplyr::ungroup()

treatment_recode <- c("MOCK"="mock", "SARSCoV2"="SARS_CoV2", "SARSCoV"="SARS_CoV")

report_data_path <- file.path(data_path, data_folder)
report_files <- list.files(report_data_path, "^Limma_Combined_.*\\.txt$")

rnaseq_report.df <- bind_rows(lapply(report_files, function(report_file) {
    message("Reading ", report_file)
    read_tsv(file.path(report_data_path, report_file), locale = locale(decimal_mark=",", grouping_mark=""))
}))

geneid2genename.df <- bind_rows(
    transmute(human_ncbi_genes.df, GeneID, GeneName=genomic_nucleotide_accession.version, genename_source = "NCBI_AC"),
    transmute(human_ncbi_genes.df, GeneID, GeneName=Symbol, genename_source = "NCBI_Symbol")
) %>% dplyr::filter(!is.na(GeneName)) %>% dplyr::distinct()

genename2proteinac.df <- dplyr::inner_join(geneid2genename.df, geneid2proteinac.df) %>%
    dplyr::select(-GeneID) %>% dplyr::distinct() %>%
    dplyr::group_by(protein_ac) %>%
    dplyr::mutate(has_genesymbol = any(genename_source == "NCBI_Symbol")) %>%
    dplyr::ungroup()
setdiff(unique(rnaseq_report.df$GeneName), unique(genename2proteinac.df$GeneName))

object_contrasts.df <- rnaseq_report.df %>%
    tidyr::extract(Treatment, c("treatment_lhs", "timepoint_num_lhs"),
                   "(SARSCoV2?|MOCK)(6|12|24)h", remove=FALSE, convert=TRUE) %>%
    tidyr::extract(Reference, c("treatment_rhs", "timepoint_num_rhs"),
                   "(SARSCoV2?|MOCK)(6|12|24)h", remove=FALSE, convert=TRUE) %>%
    dplyr::mutate(treatment_lhs = factor(set_names(treatment_recode[treatment_lhs], NULL), levels=treatment_recode),
                  treatment_rhs = factor(set_names(treatment_recode[treatment_rhs], NULL), levels=treatment_recode),
                  timepoint_lhs = factor(timepoint_num_lhs),
                  timepoint_rhs = factor(timepoint_num_rhs),
                  condition_lhs = str_c(treatment_lhs, "@", timepoint_num_lhs, "h"),
                  condition_rhs = str_c(treatment_rhs, "@", timepoint_num_rhs, "h"),
                  contrast = str_c(condition_lhs, "_vs_", condition_rhs)) %>%
    #dplyr::left_join(dplyr::select(genename2proteinac.df, GeneName, protein_ac)) %>%
    dplyr::rename(gene_name = GeneName, mean_log2 = log2FoldChange, abu_log2_rhs = baseMean, p_value=pvalue, p_value_adj=padj) %>%
    dplyr::select(-Treatment, -Reference, -starts_with("Barcodes"),
                  -starts_with("Within"), -starts_with("Selection"),
                  -Method, -Design, -TreatmentFactor, -NormStrategyIdentifier, -ComparisonName)

object_contrasts_thresholds.df <- select(object_contrasts.df, contrast) %>%
    distinct() %>%
    mutate(contrast_offset_log2 =0.0,
           p_value_threshold = case_when(TRUE ~ 0.05),
           p_value_threshold_lesser = case_when(TRUE ~ 0.001),
           mean_log2_threshold = case_when(TRUE ~ 0.0),
           mean_log2_threshold_lesser = case_when(TRUE ~ 0.1))

object_contrasts.df <- dplyr::left_join(object_contrasts.df, object_contrasts_thresholds.df) %>%
    dplyr::mutate(contrast_type = "comparison",
                  object_id = gene_name,
                  object_label = gene_name,
                  is_signif = p_value_adj <= p_value_threshold & abs(mean_log2) >= mean_log2_threshold,
                  is_signif_lesser = p_value <= p_value_threshold_lesser & abs(mean_log2) >= mean_log2_threshold_lesser,
                  is_hit = is_signif,
                  is_specific_virus_lhs = str_detect(treatment_lhs, "SARS"),
                  is_specific_virus_rhs = str_detect(treatment_rhs, "SARS"),
                  is_mock_rhs = treatment_rhs == "mock") %>%
    dplyr::group_by(object_id, timepoint_lhs, timepoint_rhs) %>%
    # declare _vs_mock a hit if it's a hit or it's significant with less stringent threshold
    # and a strong hit in another virus treatment, while the difference between the viruses is not significant
    dplyr::mutate(is_hit_composed = is_hit | (is_signif_lesser & is_mock_rhs &
                                              any(is_hit[is_specific_virus_lhs & is_mock_rhs]) & 
                                              any(!is_signif[is_specific_virus_rhs])),
                  composed_hit_treatment = if_else(is_hit_composed & is_specific_virus_lhs & is_mock_rhs,
                                                   treatment_lhs, NA_integer_),
                  composed_hit_treatments = str_c(composed_hit_treatment[!is.na(composed_hit_treatment)], collapse="+")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(composed_hit_type = if_else(str_detect(composed_hit_treatments, fixed("+")), "shared",
                                              if_else(composed_hit_treatments != "", composed_hit_treatments, "none")),
                  change = if_else(is_hit_composed, if_else(mean_log2 < 0, "-", "+"), "."))

object_contrasts_wide.df <- pivot_wider(object_contrasts.df,
                                        id_cols = c("object_id", "gene_name"),
                                        names_from = "contrast", names_sep = ".",
                                        values_from = c("mean_log2", "p_value", "is_hit_composed", "change"))

rfit_filepath <- file.path(scratch_path, paste0(project_id, '_', data_folder, '_', fit_version, '.RData'))
results_info <- list(project_id = project_id, data_folder=data_folder,
                     data_version = data_version, fit_version = fit_version)
message('Saving full analysis results to ', rfit_filepath, '...')
save(file=rfit_filepath, results_info,
     rnaseq_report.df, genename2proteinac.df, object_contrasts.df, object_contrasts_thresholds.df)

source(file.path(misc_scripts_path, 'ggplot_ext.R'))

require(Cairo)
require(ggrastr)
require(ggrepel)
require(ggnewscale)
require(ggforce)

object_contrasts_trunc.df <- tibble(contrast = object_contrasts_thresholds.df$contrast,
                                    mean_log2_max = 1.5)

object_contrasts_4show.df <- dplyr::group_by(object_contrasts.df, gene_name, contrast) %>%
    filter(row_number() == 1L) %>%
    dplyr::ungroup()

object_iactions_4show.df <- filter(object_contrasts_4show.df,
                                   treatment_lhs %in% c("SARS_CoV2", "SARS_CoV") &
                                   treatment_rhs %in% c("SARS_CoV2", "SARS_CoV")) %>%
    dplyr::select(contrast, timepoint_lhs, timepoint_rhs,
                  contrast_mean_log2 = mean_log2, contrast_p_value = p_value,
                  is_signif, is_hit, is_hit_composed, composed_hit_type,
                  object_id, object_label) %>%
    dplyr::mutate(condition_lhs = str_remove(contrast, "_vs_.+"),
                  condition_rhs = str_remove(contrast, ".+_vs_")) %>%
    dplyr::mutate(contrast_lhs = str_c(condition_lhs, "_vs_mock@", timepoint_lhs, "h"),
                  contrast_rhs = str_c(condition_rhs, "_vs_mock@", timepoint_rhs, "h")) %>%
    dplyr::left_join(dplyr::select(object_contrasts_trunc.df, contrast, contrast_mean_log2_max = mean_log2_max)) %>%
    dplyr::left_join(dplyr::filter(object_contrasts_4show.df, contrast_type=="comparison") %>%
                     dplyr::transmute(contrast_lhs = contrast, object_id, contrast_lhs_p_value = p_value, contrast_lhs_mean_log2 = mean_log2,
                                      lhs_mean_log2 = abu_log2_rhs + mean_log2,
                                      is_signif_lhs = is_signif, is_hit_lhs = is_hit,
                                      is_hit_composed_lhs = is_hit_composed)) %>%
    dplyr::left_join(dplyr::filter(object_contrasts_4show.df, contrast_type=="comparison") %>%
                     dplyr::transmute(contrast_rhs = contrast, object_id, contrast_rhs_p_value = p_value, contrast_rhs_mean_log2 = mean_log2,
                                      rhs_mean_log2 = abu_log2_rhs + mean_log2,
                                      is_signif_rhs = is_signif, is_hit_rhs = is_hit,
                                      is_hit_composed_rhs = is_hit_composed)) %>%
    dplyr::mutate(is_foreground = coalesce(is_hit_composed_lhs, FALSE) | coalesce(is_hit_composed_rhs, FALSE),
                  # coalesce(is_hit_nomschecks_lhs, FALSE) | coalesce(is_hit_nomschecks_rhs, FALSE),
                  is_hilite = is_foreground,# & is_hit_nomschecks,
                  contrast_lhs_mean_log2_trunc = if_else(is.na(contrast_mean_log2_max), contrast_lhs_mean_log2,
                                                         pmax(pmin(contrast_mean_log2_max, contrast_lhs_mean_log2), -contrast_mean_log2_max)),
                  contrast_rhs_mean_log2_trunc = if_else(is.na(contrast_mean_log2_max), contrast_rhs_mean_log2,
                                                         pmax(pmin(contrast_mean_log2_max, contrast_rhs_mean_log2), -contrast_mean_log2_max)),
                  show_label = is_foreground,# & (is_hit_nomschecks | is_viral),
                  truncation = scatter_truncation(contrast_lhs_mean_log2, contrast_lhs_mean_log2_trunc,
                                                  contrast_rhs_mean_log2, contrast_rhs_mean_log2_trunc,
                                                  composed_hit_type != "none"),
                  truncation_type = point_truncation_type(truncation, is_foreground))

scatter_type <- "contrast"
object_iactions_4show_kde.df <- object_iactions_4show.df %>% filter(!is.na(contrast_lhs_mean_log2) & !is.na(contrast_rhs_mean_log2)) %>%
    group_by(contrast) %>% group_modify(~{
        kde2d4plot(.x, "contrast_lhs_mean_log2_trunc", "contrast_rhs_mean_log2_trunc", n = 400)$density_df
    }) %>% ungroup()

modelobj_suffix <- "gene"
orgcode_palette <- c(HUMAN="black", contaminant="gray", SARS2 = "#F4982A", CVHSA = "#811A02")
composed_hit_type_palette <- c(shared="black", SARS_CoV2 = "#F4982A", SARS_CoV = "#811A02", none="gray")

show_labels <- FALSE
object_iactions_4show.df %>%
    group_by(contrast) %>% do({
        sel_object_contrast.df <- dplyr::filter(., is_hilite | between(percent_rank(contrast_lhs_mean_log2), 0.001, 0.999))
        sel_object_contrast_thresholds.df <- semi_join(object_contrasts_thresholds.df, sel_object_contrast.df)
        message("Plotting ", sel_object_contrast_thresholds.df$contrast[[1]],
                " (", sum(sel_object_contrast.df$show_label), " label(s))")
        manylabels <- sum(sel_object_contrast.df$show_label) > 300
        sel_kde.df <- semi_join(object_iactions_4show_kde.df,
                                dplyr::select(sel_object_contrast.df, contrast)[1, ]) %>%
            dplyr::filter(bin2d_density >= 0.05)
        p <- if (scatter_type == "contrast"){
            ggplot(sel_object_contrast.df,
                   aes(x = contrast_lhs_mean_log2_trunc, y = contrast_rhs_mean_log2_trunc))
        } else {
            ggplot(sel_object_contrast.df,
                   aes(x = lhs_mean_log2, y = rhs_mean_log2))
        }
        p <- p +
            #geom_raster(data=sel_kde.df, aes(fill=bin2d_density), color=NA) +
            stat_contour_filled(data=sel_kde.df, aes(z=bin2d_density, fill=after_stat(level_mid)), bins=20) +
            stat_contour(data=sel_kde.df, aes(z=bin2d_density, color=after_stat(level))) +
            scale_color_gradient("density", low="gray75", high="black", trans=power_trans(0.25), guide=FALSE) +
            scale_fill_gradient("density", low="gray95", high="slategray4", trans=power_trans(0.25)) +
            geom_vline(xintercept=0, size=1, color="dodgerblue4") +
            geom_hline(yintercept=0, size=1, color="dodgerblue4") +
            new_scale_color() +
            #geom_point_rast(data=dplyr::filter(sel_object_contrast.df, !is_foreground),
            #                alpha=0.2, size=0.5, color="darkgray", shape=16L) +
            geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]], color="dodgerblue4", linetype="dashed") +
            geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] - sel_object_contrast_thresholds.df$mean_log2_threshold[[1]],
                        color="dodgerblue4", linetype="dotted", size=0.5) +
            geom_abline(slope=1, intercept=sel_object_contrast_thresholds.df$contrast_offset_log2[[1]] + sel_object_contrast_thresholds.df$mean_log2_threshold[[1]],
                        color="dodgerblue4", linetype="dotted", size=0.5)
        if (show_labels) {
            p <- p +
                geom_text_repel(data=dplyr::filter(sel_object_contrast.df, show_label),
                                aes(label = object_label),
                                size=ifelse(manylabels, 2.5, 3.5),
                                show.legend = FALSE, segment.color = "gray")
        }
        p <- p +
            geom_point(data=dplyr::filter(sel_object_contrast.df, is_foreground),
                       aes(color=composed_hit_type, shape=truncation, size=truncation_type)) +
            scale_color_manual(values=composed_hit_type_palette, na.value="black") +
            #scale_color_manual(values=orgcode_palette, na.value="black", guide="none") +
            scale_shape_manual(values=point_truncation_shape_palette, guide="none") +
            scale_size_manual(values=if_else(manylabels, 0.5, 1.0) * point_truncation_size_palette, guide="none", ) +
            xlab(str_c("log2(fold-change) ", sel_object_contrast.df$contrast_lhs)) +
            ylab(str_c("log2(fold-change) ", sel_object_contrast.df$contrast_rhs)) +
            coord_fixed() +
            theme_bw_ast()
        plot_path <- file.path(analysis_path, 'plots', str_c(data_folder,'_', fit_version),
                               str_c("scatter_", scatter_type, "s_", modelobj_suffix))
        if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
        ggsave(filename = file.path(plot_path,
                                    str_c(project_id, '_', fit_version, "_scatter_", scatter_type, "s_",
                                          str_replace_all(sel_object_contrast_thresholds.df$contrast[[1]], '\\?', 'alt'),
                                          if_else(show_labels, "", "_nolabels"), '.pdf')),
               plot = p, width=16, height=16, device=cairo_pdf, family="Arial")
        tibble()
    })

report_cols <- c("gene_name", "object_id",
                 "is_contaminant", "is_viral")

pre_object_contrasts_report.df <- filter(object_contrasts.df, treatment_lhs != "infected") %>%
    dplyr::select(object_id, contrast, treatment_lhs, treatment_rhs, timepoint=timepoint_lhs,
                  mean_log2, any_of(c("prob_nonpos", "prob_nonneg", "p_value", "p_value_adj")),
                  is_signif, is_hit, change, is_hit_composed, composed_hit_type)

objects4report.df <- dplyr::select(object_contrasts.df, any_of(report_cols)) %>% distinct() %>%
    dplyr::semi_join(dplyr::select(object_contrasts.df, object_id))

object_contrasts_report.df <- objects4report.df %>%
    dplyr::left_join(pivot_wider(pre_object_contrasts_report.df, c(object_id),
                                 names_from = "contrast", values_from = c("is_hit", "is_hit_composed", "composed_hit_type", "change", "mean_log2", "p_value", "p_value_adj"),
                                 names_sep=".")) %>%
    dplyr::select(any_of(report_cols),
                  ends_with("SARS_CoV2@6h_vs_mock@6h"), ends_with("SARS_CoV@6h_vs_mock@6h"), ends_with("SARS_CoV2@6h_vs_SARS_CoV@6h"),
                  ends_with("SARS_CoV2@12h_vs_mock@12h"), ends_with("SARS_CoV@12h_vs_mock@12h"), ends_with("SARS_CoV2@12h_vs_SARS_CoV@12h"),
                  ends_with("SARS_CoV2@24h_vs_mock@24h"), ends_with("SARS_CoV@24h_vs_mock@24h"), ends_with("SARS_CoV2@24h_vs_SARS_CoV@24h")) %>%
    # composed_hit_type is the same for all comparison of the same timepoint
    dplyr::select(-matches("composed_hit_type.+SARS_CoV2@\\d+h_vs_(SARS_CoV|mock)@")) %>%#, #-matches("is_hit_composed\\.SARS_CoV2.+_vs_SARS_CoV")) %>%
    dplyr::rename_at(vars(matches("composed_hit_type")), ~str_remove(.x, "SARS_CoV@\\d+h_vs_mock@")) %>%
    dplyr::arrange(gene_name)

write_tsv(object_contrasts_report.df %>% dplyr::select(-object_id),
          file.path(analysis_path, "reports", paste0(project_id, '_', data_folder, '_contrasts_report_', fit_version, '_wide.txt')))

require(writexl)
write_xlsx(object_contrasts_report.df %>% dplyr::select(-object_id),
           file.path(analysis_path, "reports", paste0(project_id, '_', data_folder, '_contrasts_report_', fit_version, '_wide.xlsx')))

object_contrasts_long_report.df <- objects4report.df %>%
    dplyr::left_join(pre_object_contrasts_report.df) %>%
    #dplyr::left_join(dplyr::select(contrasts.df, contrast, treatment_lhs, treatment_rhs)) %>%
    dplyr::arrange(gene_name, timepoint, treatment_lhs, treatment_rhs, contrast)

write_tsv(filter(object_contrasts_long_report.df) %>% dplyr::select(-object_id),
          file.path(analysis_path, "reports", paste0(project_id, '_', data_folder, '_contrasts_report_', fit_version, '_long.txt')))


