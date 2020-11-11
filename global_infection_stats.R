project_id <- 'cov2'
plot_version <- "20201027"
message('Project ID=', project_id)
datasets <- list(
    phospho = list(msfolder = 'snaut_parsars_phospho_20201005',
                   data_ver = "20201012",
                   fit_ver = "20201012",
                   ptm_extractor_ver = "20201012",
                   label = "Phosphoproteome",
                   type = "Phosphoproteome", protocol = "DIA", goal = "timecourse",
                   format = "msglm"),
    fp = list(msfolder = 'snaut_parsars_fp_20200829',
              data_ver = "20200830",
              fit_ver = "20200830",
              label = "Proteome", protocol = "DIA", goal = "timecourse",
              type = "Proteome",
              format = "msglm"),
    ubi = list(msfolder="snaut_parsars_ptm_20200907",
               data_ver = "20201012",
               fit_ver="20201012",
               ptm_extractor_ver = "20201012",
               ptm="ubiquitin", prefix="ubi",
               type = "Ubiquitinome",
               format = "msglm",
               label = "Ubiquitinome", protocol = "DIA", goal = "timecourse"
               ),
    rnaseq = list(data_folder = "parsars_rnaseq_20201020",
                  fit_ver = "20201020",
                  format = "limma",
                  label = "Transcriptome", protocol = "RNAseq", goal = "timecourse",
                  type = "Transcriptome")
)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))
source(file.path(project_scripts_path, 'cov2_plots_common.R'))

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
        next
    }
    ds.env <- new.env(parent=baseenv())
    assign(str_c(envname, ".env"), ds.env)
    envinfo <- datasets[[envname]]
    load(file.path(scratch_path, paste0(project_id, '_msglm_data_', envinfo$msfolder, '_', envinfo$data_ver, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msdata_full_', envinfo$msfolder, '_', envinfo$data_ver, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msglm_fit_', envinfo$msfolder, '_', envinfo$fit_ver, '.RData')), envir = ds.env)
    if ("ptm_extractor_ver" %in% names(envinfo)) {
        ds.env$ptmns_grouped.df <- read_tsv(file.path(data_path, envinfo$msfolder, str_c("ptm_extractor_", envinfo$ptm_extractor_ver), "ptmns_grouped.txt.gz"),
                                     col_types=cols(
                                         .default = col_guess(),
                                         ptmn_id = "i", ptm_id = "i", nselptms = "i",
                                         ptmgroup_id = "i", ptmngroup_id = "i",
                                         ptm_pos = "i", ptm_AA_seq = "c", ptm_is_known = col_logical(),
                                         protein_rank = "i",
                                         ptmid_is_reference = col_logical()
                                     )) %>%
            dplyr::group_by(ptmgroup_id) %>%
            dplyr::mutate(nptmngroups = n_distinct(ptmngroup_id)) %>%
            dplyr::group_by(ptm_id) %>%
            dplyr::mutate(nptmns = n_distinct(ptmn_id)) %>%
            dplyr::ungroup()
        ds.env$ptmgroup2protgroup.df <- read_tsv(file.path(data_path, envinfo$msfolder, str_c("ptm2protgroup_", envinfo$ptm_extractor_ver, ".txt.gz")))
    }
}

rnaseq.env <- new.env(parent=baseenv())
load(file.path(scratch_path, paste0(project_id, '_', datasets$rnaseq$data_folder, '_', datasets$rnaseq$fit_ver, '.RData')), envir = rnaseq.env)

datatype_order = c("Transcriptome", "Proteome", "Ubiquitinome", "Phosphoproteome")
datatype_palette = c(
    Phosphoproteome = "#5e268f",
    Transcriptome = "#383c9b",
    Proteome = "#226430",
    Ubiquitinome = "#bf1c2c"
)

sel_std_type <- "median"

cov2_infection_stats.df <- bind_rows(
    bind_rows(
        mutate(phospho.env$object_contrasts.df, dataset = "phospho", entity="PTM", object_id=as_integer(object_id)),
        mutate(ubi.env$object_contrasts.df, dataset = "ubi", entity="PTM", object_id=as_integer(object_id))
    ) %>% dplyr::filter(std_type == sel_std_type & contrast_kind == "treatment_vs_treatment" & treatment_rhs == "mock" & treatment_lhs != "infected") %>%
    group_by(entity, dataset, treatment=treatment_lhs, timepoint=timepoint_rhs, composed_hit_type, change, is_viral, is_contaminant) %>%
    summarise(#n_regulated_proteins = n_distinct(unlist(str_split(protgroup_ids, ';'))),
              n_entities = n_distinct(ptmgroup_id), .groups="drop"),
    mutate(fp.env$object_contrasts.df, dataset = "fp", entity="protein") %>%
    dplyr::filter(std_type == sel_std_type & contrast_kind == "treatment_vs_treatment" & treatment_rhs == "mock" & treatment_lhs != "infected") %>%
    group_by(entity, dataset, treatment=treatment_lhs, timepoint=timepoint_rhs, composed_hit_type, change, is_viral, is_contaminant) %>%
    summarise(n_entities = n_distinct(protregroup_id), .groups="drop"),
    mutate(rnaseq.env$object_contrasts.df, dataset = "rnaseq", entity="gene", is_viral=FALSE, is_contaminant=FALSE) %>%
    dplyr::filter(treatment_rhs == "mock" & treatment_lhs != "infected") %>%
    group_by(entity, dataset, treatment=treatment_lhs, timepoint=timepoint_rhs, composed_hit_type, change, is_viral, is_contaminant) %>%
    summarise(n_entities = n_distinct(object_id), .groups="drop")
) %>% ungroup() %>%
    mutate(n_regulated = if_else(change == ".", 0L, n_entities),
           dataset_label = sapply(dataset, function(ds) datasets[[ds]]$label),
           dataset_type = factor(sapply(dataset, function(ds) datasets[[ds]]$type),
                                 levels = datatype_order),
           protocol = factor(sapply(dataset, function(ds) datasets[[ds]]$protocol)),
           timepoint_num = parse_integer(timepoint),
           timepoint = factor(timepoint_num, ordered = TRUE),
           treatment = factor(treatment, c("mock", "SARS_CoV2", "SARS_CoV")),
           composed_hit_type = factor(composed_hit_type, c("none", "shared", "SARS_CoV2", "SARS_CoV"))) %>%
    dplyr::arrange(dataset_type, timepoint, treatment, dataset, is_contaminant, is_viral, composed_hit_type) %>%
    dplyr::mutate(data_slice = str_c(dataset_label, "@", timepoint, "h"),
                  data_slice = factor(data_slice, levels=unique(data_slice)),
                  timepointXtreatment = str_c(timepoint, "h ", treatment),
                  timepointXtreatment = factor(timepointXtreatment, levels=unique(timepointXtreatment))) %>%
    dplyr::group_by(dataset_type, timepoint, treatment, dataset, is_contaminant, is_viral, change) %>%
    dplyr::mutate(n_regulated_stacked = cumsum(n_regulated)) %>% dplyr::ungroup() %>%
    dplyr::mutate(n_regulated_stacked_prev = n_regulated_stacked - n_regulated,
                  n_regulated_stacked_mid = 0.5*(n_regulated_stacked + n_regulated_stacked_prev)) %>%
    dplyr::mutate_at(vars(starts_with("n_regulated")), list(vis = ~if_else(change == "+", .x, -.x)))
write_tsv(cov2_infection_stats.df, file=file.path(analysis_path, "reports", str_c(project_id, "_nregulated_per_timepoint_", plot_version, ".txt")))

composed_hit_type_fill_pallette <- c("none" = "lightgray", "shared" = "gray80", "SARS_CoV2" = "#F4982A", "SARS_CoV" = "#811A02")
composed_hit_type_pallette <- c("none" = "lightgray", "shared" = "gray", "SARS_CoV2" = "#F4982A", "SARS_CoV" = "#811A02")
shown_cov2_infection_stats.df <- dplyr::filter(cov2_infection_stats.df, change != "." & !is_viral & !is_contaminant) %>%
    dplyr::mutate(composed_hit_type = factor(composed_hit_type, rev(levels(composed_hit_type)))) # to make ggplot dislay in proper order
p <- ggplot(shown_cov2_infection_stats.df,
       aes(x = timepointXtreatment, fill=composed_hit_type, color=composed_hit_type, y=n_regulated_vis)) +
    geom_hline(yintercept = 0, color="black", size=1) +
    geom_col(position=position_stack(), width=0.9, alpha=0.5, show.legend=FALSE) +
    geom_text(data=filter(shown_cov2_infection_stats.df, composed_hit_type != "shared" & n_regulated_vis>0),
              aes(y=n_regulated_stacked_vis, label=n_regulated),
              size=5, vjust=-0.5, show.legend=FALSE) +
    geom_text(data=filter(shown_cov2_infection_stats.df, composed_hit_type != "shared" & n_regulated_vis<0),
              aes(y=n_regulated_stacked_vis, label=n_regulated),
              size=5, vjust=1.5, show.legend=FALSE) +
    geom_text(data=filter(shown_cov2_infection_stats.df, composed_hit_type == "shared"),
              aes(y=n_regulated_stacked_mid_vis, label=n_regulated),
              size=5, nudge_y=-0.5, show.legend=FALSE) +
    scale_x_discrete('Timepoint', labels=c("6h", "", "12h", "", "24h", "", "36h", "")) +# breaks=unique(stats.df[['timepoint']])) +
    scale_y_continuous("", expand = expansion(0.15, 0.1), labels=NULL) +
    scale_color_manual(values = composed_hit_type_pallette) +
    scale_fill_manual(values = composed_hit_type_fill_pallette) +
    #scale_alpha_manual(values = c("ptm" = 0.25, "protein" = 0.75)) +
    theme_bw_ast(base_size=20) +
    facet_grid(dataset_type ~ ., scale="free", space = "free_x", switch = "y") +
    theme(panel.grid = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank(), #plot.margin = element_blank(),
          #axis.text.x = element_text(angle=-45, hjust=0),
          panel.border = element_blank(), panel.spacing.y = unit(5, "pt"))
p
ggsave(file.path(analysis_path, "plots", str_c(project_id, '_infection_nregulated_', plot_version, '.pdf')), p, width=5, height=10, device=cairo_pdf)

cov2_total_stats.df <- bind_rows(
    bind_rows(
        mutate(phospho.env$object_contrasts.df, dataset = "phospho", entity="PTM"),
        mutate(ubi.env$object_contrasts.df, dataset = "ubi", entity="PTM")
    ) %>%
        dplyr::filter(std_type == sel_std_type & contrast_kind == "treatment_vs_treatment" & treatment_rhs == "mock" & treatment_lhs != "infected") %>%
        dplyr::group_by(entity, dataset, is_contaminant) %>%
        summarise(n_regulated = n_distinct(ptmgroup_id[change != "."]),
                  n_total = n_distinct(ptmgroup_id),
                  .groups="drop"),
    bind_rows(
        mutate(phospho.env$object_contrasts.df, dataset = "phospho", entity="protgroup") %>%
            dplyr::inner_join(phospho.env$ptmgroup2protgroup.df),
        mutate(ubi.env$object_contrasts.df, dataset = "ubi", entity="protgroup") %>%
            dplyr::inner_join(ubi.env$ptmgroup2protgroup.df)
    ) %>%
        dplyr::filter(std_type == sel_std_type & contrast_kind == "treatment_vs_treatment" & treatment_rhs == "mock" & treatment_lhs != "infected") %>%
        dplyr::group_by(entity, dataset, is_contaminant) %>%
        summarise(n_regulated = n_distinct(protgroup_id[change != "."]),
                  n_total = n_distinct(protgroup_id),
                  .groups="drop"),
    mutate(fp.env$object_contrasts.df, dataset = "fp", entity="protein") %>%
        dplyr::filter(std_type == sel_std_type & contrast_kind == "treatment_vs_treatment" & treatment_rhs == "mock" & treatment_lhs != "infected") %>%
        group_by(entity, dataset, is_contaminant) %>%
        summarise(n_regulated = n_distinct(protregroup_id[change != "."]),
                  n_total = n_distinct(protregroup_id),
                  .groups="drop"),
    mutate(rnaseq.env$object_contrasts.df, dataset = "rnaseq", entity="gene", is_viral=FALSE, is_contaminant=FALSE) %>%
        dplyr::filter(treatment_rhs == "mock" & treatment_lhs != "infected") %>%
        group_by(entity, dataset, is_contaminant) %>%
        summarise(n_regulated = n_distinct(gene_name[change != "."]),
                  n_total = n_distinct(gene_name),
                  .groups="drop")
)
write_tsv(cov2_total_stats.df, file=file.path(analysis_path, "reports", str_c(project_id, "_ntotal_stats_", plot_version, ".txt")))
