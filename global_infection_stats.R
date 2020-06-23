source(file.path(project_scripts_path, 'cov2_plots_common.R'))

project_id <- 'cov2'
message('Project ID=', project_id)
datasets <- list(
    cov2ts_phospho = list(ms_folder = 'cov2timecourse_phospho_dia_20200423',
                          data_version = "20200428",
                          fit_version = "20200428",
                          label = "Phosphoproteome (DIA)",
                          type = "Phosphoproteome", protocol = "DIA", goal = "timecourse",
                          format = "msglm"),
    cov2ts_proteome = list(ms_folder = 'cov2timecourse_dia_20200423',
                           data_version = "20200429",
                           fit_version = "20200429",
                           label = "Proteome (DIA)", protocol = "DIA", goal = "timecourse",
                           type = "Proteome",
                           format = "msglm"),
    cov2el_proteome = list(ms_folder = 'cov2earlylate_fp_phos_ubi_dda_20200429',
                           data_version = "20200514",
                           fit_version = "20200514",
                           label = "Proteome (DDA)", protocol = "DDA", goal = "PTM crosstalk",
                           type = "Proteome",
                           format = "msglm"),
    cov2ts_rnaseq = list(data_folder = "cov2ts_tx",
                         fit_ver = "20200518",
                         format = "rnaseq",
                         label = "Transcriptome", protocol = NA_character_, goal = "timecourse",
                         type = "Transcriptome"),
    cov2el_phospho = list(ms_folder = "cov2earlylate_fp_phos_ubi_dda_20200601",
                          filename="output_Phospho (STY)Sites",
                          fit_version="20200510",
                          ptm="phospho", prefix="pho",
                          type = "Phosphoproteome",
                          format = "perseus",
                          label = "Phosphoproteome (DDA)", protocol = "DDA", goal = "PTM crosstalk"
                          ),
    cov2el_ubi = list(ms_folder="cov2earlylate_fp_phos_ubi_dda_20200601",
                      filename="output_S0-01_GlyGly (K)Sites",
                      fit_version="20200510",
                      ptm="ubiquitin", prefix="ubi",
                      type = "Ubiquitinome",
                      format = "perseus",
                      label = "Ubiquitinome (DDA)", protocol = "DDA", goal = "PTM crosstalk"
                      )
)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

party3rd_data_path <- file.path(bioinfo_pool_path, "pub3rdparty")

require(rstan)
require(dplyr)
require(tidyselect)
require(stringr)
require(msglm)
require(maxquantUtils)

for (envname in names(datasets)) {
    message("Loading ", envname, " data...")
    if (datasets[[envname]]$format != "msglm") {
        continue
    }
    ds.env <- new.env(parent=baseenv())
    assign(str_c(envname, ".env"), ds.env)
    envinfo <- datasets[[envname]]
    load(file.path(scratch_path, paste0(project_id, '_msglm_data_', envinfo$ms_folder, '_', envinfo$data_version, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msdata_full_', envinfo$ms_folder, '_', envinfo$data_version, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msglm_fit_', envinfo$ms_folder, '_', envinfo$data_version, '.RData')), envir = ds.env)
}
cov2ts_phospho.env$msdata_full$ptmgroup2psitep <- read_tsv(file.path(data_path, datasets$cov2ts_phospho$ms_folder,
                                                                     "COV2_DIA_phospho_0.75probablity_no normalization_psitep_nodata.txt"))
cov2ts_rnaseq.env <- new.env(parent=baseenv())
cov2ts_rnaseq.env$rnaseq_contrasts.df <- cov2ts_rnaseq.env$all_pairwise_comp.df %>%
    # recode ComparisonName into the same format as contrast
    mutate(ComparisonName = str_replace(str_replace(ComparisonName, "SARSCoV2", "SARS_CoV2"), "MOCK", "mock")) %>%
    tidyr::extract(ComparisonName,
                   c("timepoint_lhs", "treatment_lhs",
                     "timepoint_rhs", "treatment_rhs",
                     "timepoint", "treatment"),
                     "(?:(\\d+)h|(SARS_CoV2))_vs_(?:(\\d+)h|(mock))within(?:(\\d+)h|(SARS_CoV2|mock))") %>%
    mutate(timepoint_lhs = if_else(!is.na(timepoint_lhs), timepoint_lhs, timepoint),
           timepoint_rhs = if_else(!is.na(timepoint_rhs), timepoint_rhs, timepoint),
           treatment_lhs = if_else(!is.na(treatment_lhs), treatment_lhs, treatment),
           treatment_rhs = if_else(!is.na(treatment_rhs), treatment_rhs, treatment),
           contrast = str_c(treatment_lhs, "@", timepoint_lhs, "h_vs_",
                            treatment_rhs, "@", timepoint_rhs, "h"))

for (envname in names(datasets)) {
    if (datasets[[envname]]$format != "perseus") {
        continue
    }
    message("Loading perseus ", envname, " data...")
    ds.env <- new.env(parent=baseenv())
    assign(str_c(envname, ".env"), ds.env)
    dsinfo <- perseus_datasets[[envname]]
    
    ds.env$perseus_analysis_orig.df = read_tsv(file.path(data_path, dsinfo$ms_folder,
                                                         "curban_analysis", paste0(dsinfo$filename, ".txt")),
                                               col_names = TRUE, comment="#", quote = "", guess_max = 10000) %>%
        dplyr::rename(ptmgroup_id = id)
    #vars(c("Positions within proteins", "Leading proteins",
    #       "Protein", "Protein names", "Gene names"))
    ds.env$perseus_analysis_long.df <- tidyr::pivot_longer(
        dplyr::select(ds.env$perseus_analysis_orig.df, ptmgroup_id, protgroup_ids = `Protein group IDs`,
                      protein_ac = Protein, ptm_pos = Position, ptm_aa = `Amino acid`, charge = Charge, Multiplicity,
                      gene_names = `Gene names`,
                      contains("Student's T-test")),
        cols=contains("Student's T-test"),
        names_pattern=str_c("^(.+) (", dsinfo$prefix, "_.+)$"),
        names_to=c(".value", "comparison")) %>%
        dplyr::rename(is_signif = `Student's T-test Significant`,
                      delta_log2 = `Student's T-test Difference`,
                      mlog10_p_value = `-Log Student's T-test p-value`,
                      q_value = `Student's T-test q-value`,
                      ttest_stat = `Student's T-test Test statistic`) %>%
        dplyr::filter(!is.na(comparison)) %>%
        tidyr::extract(comparison, c("timepoint"), paste0(dsinfo$prefix, "_SARS_COV2_(\\d+)h_", dsinfo$prefix, "_mock_(?:\\d+)h"), remove=FALSE) %>%
        dplyr::mutate_at(vars(mlog10_p_value, ttest_stat), ~if_else(is.nan(.), NA_real_, .)) %>%
        dplyr::mutate(contrast_lhs = "SARS_COV2",
                      contrast_rhs = "mock",
                      ptm_type = dsinfo$ptm,
                      dataset = envname,
                      timepoint = parse_integer(timepoint),
                      p_value = 10^(-mlog10_p_value),
                      is_signif = coalesce(is_signif == "+", FALSE),
                      change = if_else(is_signif, if_else(delta_log2 > 0, "+", "-"), "."),
                      ptm_label = str_c(gene_names, "_", ptm_aa, ptm_pos)
        )
    # leave the most significant comparison
    ds.env$perseus_analysis_short.df <- dplyr::group_by(ds.env$perseus_analysis_long.df,
                                                        protein_ac, ptm_pos, ptm_aa, timepoint) %>%
        dplyr::filter(row_number(-mlog10_p_value) == 1L) %>%
        dplyr::ungroup()
}

datatype_order = c("Transcriptome", "Proteome", "Phosphoproteome", "Ubiquitinome")
datatype_palette = c(
    Phosphoproteome = "#5e268f",
    Transcriptome = "#383c9b",
    Proteome = "#226430",
    Ubiquitinome = "#bf1c2c"
)

cov2_infection_stats.df <- bind_rows(
    bind_rows(
        mutate(cov2el_phospho.env$perseus_analysis_short.df, dataset = "cov2el_phospho"),
        mutate(cov2el_ubi.env$perseus_analysis_short.df, dataset = "cov2el_ubi")
    ) %>% group_by(dataset, timepoint, change) %>%
    summarise(n_regulated_proteins = n_distinct(unlist(str_split(protgroup_ids, ';'))),
              n_regulated_ptms = n_distinct(ptmgroup_id)),
    bind_rows(
        mutate(cov2ts_phospho.env$object_contrasts.df, dataset = "cov2ts_phospho",
               ptm_label = str_remove(object_label, "_M\\d+")) %>%
        dplyr::filter(std_type == "replicate" & str_detect(contrast, "SARS_.+_vs_mock")) %>%
        extract(contrast, c('timepoint_lhs', 'timepoint_rhs'), 'SARS_COV2@(\\d+)h_vs_mock@(\\d+)h', remove=FALSE) %>%
        mutate(timepoint = parse_integer(timepoint_lhs)) %>%
        separate_rows(majority_protein_acs, sep=';') %>%
        rename(protein_ac = majority_protein_acs) %>%
        left_join(dplyr::filter(cov2ts_proteome.env$msdata_full$protein2protgroup, is_majority)) %>%
        group_by(dataset, timepoint, change) %>%
        summarise(n_regulated_proteins = n_distinct(protgroup_id),
                  n_regulated_ptms = n_distinct(ptm_label)) %>%
        ungroup(),
        mutate(cov2ts_proteome.env$object_contrasts.df, dataset = "cov2ts_proteome") %>%
        dplyr::filter(std_type == "median" & str_detect(contrast, "SARS_.+_vs_mock")) %>%
        tidyr::extract(contrast, c('timepoint_lhs', 'timepoint_rhs'), 'SARS_COV2@(\\d+)h_vs_mock@(\\d+)h', remove=FALSE) %>%
        mutate(timepoint = parse_integer(timepoint_lhs)) %>%
        group_by(dataset, timepoint, change) %>%
        summarise(n_regulated_proteins = n_distinct(protgroup_id),
                  n_regulated_ptms = NA_integer_) %>%
        ungroup(),
        mutate(cov2el_proteome.env$object_contrasts.df, dataset = "cov2el_proteome") %>%
            dplyr::filter(std_type == "replicate" & str_detect(contrast, "SARS_.+_vs_mock")) %>%
            tidyr::extract(contrast, c('timepoint_lhs', 'timepoint_rhs'), 'SARS_COV2@(\\d+)h_vs_mock@(\\d+)h', remove=FALSE) %>%
            mutate(timepoint = parse_integer(timepoint_lhs)) %>%
            group_by(dataset, timepoint, change) %>%
            summarise(n_regulated_proteins = n_distinct(protregroup_id),
                      n_regulated_ptms = NA_integer_) %>%
        ungroup(),
        mutate(cov2ts_rnaseq.env$rnaseq_contrasts.df, dataset = "cov2ts_rnaseq") %>%
        dplyr::filter(treatment_lhs != treatment_rhs & timepoint_lhs == timepoint_rhs) %>%
        mutate(timepoint = parse_integer(timepoint_lhs)) %>%
        group_by(dataset, timepoint, change) %>%
        summarise(n_regulated_proteins = n_distinct(GeneName),
                  n_regulated_ptms = NA_integer_) %>%
        ungroup()
    )
) %>% mutate(dataset_label = sapply(dataset, function(ds) datasets[[ds]]$label),
             dataset_type = factor(sapply(dataset, function(ds) datasets[[ds]]$type),
                                   levels = datatype_order),
             protocol = factor(sapply(dataset, function(ds) datasets[[ds]]$protocol)),
             goal = factor(sapply(dataset, function(ds) datasets[[ds]]$goal),
                           levels=c("timecourse", "PTM crosstalk")),
             n_regulated_proteins_vis = n_regulated_proteins * if_else(change == "+", 1, -1),
             n_regulated_ptms_vis = n_regulated_ptms * if_else(change == "+", 1, -1)) %>%
    dplyr::arrange(dataset_type, timepoint, dataset) %>%
    dplyr::mutate(data_slice = str_c(dataset_label, "@", timepoint, "h"),
                  data_slice = factor(data_slice, levels=unique(data_slice))) %>%
    ungroup()

cov2_infection_stats_long.df = bind_rows(
    dplyr::select(
        dplyr::rename(cov2_infection_stats.df,
                      n_regulated = n_regulated_proteins,
                      n_regulated_vis = n_regulated_proteins_vis) %>%
        dplyr::mutate(entity = "protein"), -n_regulated_ptms, -n_regulated_ptms_vis),
    dplyr::select(
        dplyr::rename(cov2_infection_stats.df,
                      n_regulated = n_regulated_ptms,
                      n_regulated_vis = n_regulated_ptms_vis) %>%
            dplyr::mutate(entity = "ptm"), -n_regulated_proteins, -n_regulated_proteins_vis),
) %>%
mutate(
    n_regulated_vis = n_regulated_vis*if_else(dataset == "cov2ts_phospho", 3, 1)
)

shown_cov2_infection_stats.df <- dplyr::filter(cov2_infection_stats_long.df, change != ".")

p <- ggplot(shown_cov2_infection_stats.df,
       aes(x = factor(timepoint), fill=dataset_type, color=dataset_type, y=n_regulated_vis)) +
    geom_col(aes(alpha=entity), position=position_dodge(preserve="total", width = 0.3), width=1.4, show.legend=FALSE) +
    geom_text(data=filter(shown_cov2_infection_stats.df, n_regulated_vis>0), aes(label=n_regulated),
              size=5, vjust=-0.8, show.legend=FALSE) +
    geom_text(data=filter(shown_cov2_infection_stats.df, n_regulated_vis<0), aes(label=n_regulated),
              size=5, vjust=1.5, show.legend=FALSE) +
    scale_x_discrete('Timepoint') +# breaks=unique(stats.df[['timepoint']])) +
    scale_y_continuous('N regulated', expand = expansion(0.5, 0.1)) +
    scale_color_manual(values = datatype_palette) +
    scale_fill_manual(values = datatype_palette) +
    scale_alpha_manual(values = c("ptm" = 0.25, "protein" = 0.75)) +
    theme_bw_ast(base_size=20) +
    facet_grid(dataset_type ~ goal, scale="free", space = "free_x", switch = "y") +
    theme(panel.grid = element_blank(), axis.line = element_blank(), #plot.margin = element_blank(),
          panel.border = element_rect(color="grey"), panel.spacing.y = unit(5, "pt"))
p
ggsave(file.path(analysis_path, "plots", str_c(project_id, '_infection_nregulated.pdf')), p, width=5, height=6, device=cairo_pdf)
