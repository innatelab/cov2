require(readr)
require(dplyr)
require(maxquantUtils)

msdata.wide <- readr::read_tsv("/pool/analysis/astukalov/cov2/data/mq_apms_20200525/combined/txt/proteinGroups.txt",
                               col_types = readr::cols(`Fasta headers` = "c", `id` = "i"),
                               na = maxquantUtils:::MaxQuant_NAs, guess_max = 10000)
msdata_report.df <- dplyr::select(dplyr::filter(msdata.wide, coalesce(Reverse, "") != "+"), `Gene names`, `Majority protein IDs`,
                                  id, `Potential contaminant`, `Number of proteins`,
                                  `Peptide counts (all)`, `Peptide counts (razor+unique)`,
                                  `Peptide counts (unique)`, `Sequence coverage [%]`,
                                  `Unique + razor sequence coverage [%]`, `Unique sequence coverage [%]`,
                                  `Mol. weight [kDa]`, `Sequence lengths`,
                                  `Q-value`, `Score`,
                                  starts_with("Intensity "))
write_tsv(msdata_report.df, "/pool/analysis/astukalov/cov2/data/mq_apms_20200525/proteinGroups_reported.txt")

oemsdata.wide <- read.Spectronaut.ProteinsReport("/pool/analysis/astukalov/cov2/data/spectronaut_oeproteome_20200527/20200527_163029_COVID OE FPMS FINAL not normalied_Protein Group Report_fixed.txt",
                                                 import_data = "quantity", delim='\t')
oemsdata_colgroups <- attr(oemsdata.wide, "column_groups")

oemsdata_full <- list(
    protgroups = oemsdata.wide[, oemsdata_colgroups$protgroup],
    protgroup_intensities = pivot_longer.Spectronaut.ProtgroupIntensities(oemsdata.wide)
)
msruns.df <- read_tsv("/pool/analysis/astukalov/cov2/data/spectronaut_oeproteome_20200527/dia_oeproteome_used_msruns_20200524.txt")

oemsdata_full$msruns <- dplyr::select(oemsdata_full$protgroup_intensities, msrun_ix, raw_file) %>% dplyr::distinct() %>%
    dplyr::left_join(msruns.df) %>%
    dplyr::mutate(msrun = str_replace(str_replace(str_replace(msrun, "CoV2_S", "CoVII_S"), "CoV_S", "CoV2_S"), "CoVII_S", "CoV_S"))

oemsdata_report.df <- left_join(dplyr::select(oemsdata_full$protgroups, protgroup_id, gene_names, majority_protein_acs, protein_names, protein_descriptions, q_value),
                                pivot_wider(dplyr::left_join(oemsdata_full$protgroup_intensities, oemsdata_full$msruns),
                                  protgroup_id, names_from = msrun,
                                  names_prefix = "Intensity.", values_from = intensity)) %>%
    dplyr::select(-protgroup_id)
write_tsv(oemsdata_report.df, "/pool/analysis/astukalov/cov2/data/spectronaut_oeproteome_20200527/spectronaut_oeproteome_20200527_report.txt")

oemsfit_old.env <- new.env(parent=baseenv())
load("/pool/analysis/astukalov/cov2/scratch/cov2_msglm_fit_spectronaut_oeproteome_20200527_20200527.RData", envir=oemsfit_old.env)
oemsfit_new.env <- new.env(parent=baseenv())
load("/pool/analysis/astukalov/cov2/scratch/cov2_msglm_fit_spectronaut_oeproteome_20200527_20200608.RData", envir=oemsfit_new.env)
load("/pool/analysis/astukalov/cov2/scratch/cov2_msglm_data_spectronaut_oeproteome_20200527_20200608.RData", envir=oemsfit_new.env)

apms.env <- new.env(parent=baseenv())
load("/pool/analysis/astukalov/cov2/scratch/cov2_msglm_fit_mq_apms_20200525_20200525.RData", envir=apms.env)
load("/pool/analysis/astukalov/cov2/scratch/cov2_msglm_data_mq_apms_20200525_20200525.RData", envir=apms.env)

oeproteome_batch_contrasts.df <- filter(oemsfit_new.env$object_contrasts.df, std_type == "replicate" & str_detect(contrast, "_vs_B\\d+_others"))

oeproteome_contrasts.df <- filter(oemsfit_new.env$object_contrasts.df, std_type == "replicate" & str_detect(contrast, "_vs_controls")) %>%
    left_join(dplyr::select(oeproteome_batch_contrasts.df, bait_full_id, object_id, std_type, contrast_batch = contrast, is_hit_batch = is_hit_nomschecks,
                            median_log2_batch = median_log2, p_value_batch = p_value)) %>%
    mutate(is_upregulation_new = is_hit_nomschecks & coalesce(median_log2_batch, 0) > 0 & coalesce(median_log2, 0) > 0,
           is_hit_batch = p_value_batch <= 1E-3 & abs(median_log2_batch) >= 0.5,
           is_hit_new = is_hit_nomschecks & is_hit_batch & (coalesce(median_log2_batch, 0) * coalesce(median_log2, 0) > 0))

oeproteome_contrasts.df <- left_join(oeproteome_contrasts.df, oemsfit_old.env$object_contrasts.df,
                                     by=c("contrast", "std_type", "object_id", "bait_full_id"), suffix=c("", ".old")) %>%
    mutate(is_upregulation_old = is_hit_nomschecks.old & coalesce(median_log2.old, 0),
           is_hit_old = is_hit_nomschecks.old)

oeproteome_contrasts_apmsobj.df <- inner_join(oeproteome_contrasts.df,
                                              dplyr::filter(oemsfit_new.env$msdata$protein2protgroup, is_majority)) %>%
    dplyr::select(-protgroup_id, -object_id) %>%
    dplyr::inner_join(dplyr::filter(apms.env$msdata$protein2protregroup, is_majority)) %>%
    dplyr::mutate(object_id = protregroup_id) %>%
    dplyr::group_by(std_type, contrast, object_id) %>%
    dplyr::mutate(majority_protein_acs = str_flatten(unique(unlist(str_split(majority_protein_acs, ';'))), ';'),
                  gene_names = str_flatten(unique(unlist(str_split(gene_names, ';'))), ';')) %>%
    dplyr::filter(row_number(p_value) == 1L) %>%
    dplyr::ungroup()

oeproteome_log2_mindelta <- 4

oeproteome_contrasts_in_apms.df <- filter(apms.env$object_contrasts.df, str_detect(contrast, "_vs_others") & std_type == "replicate" & is_hit) %>%
    dplyr::select(std_type, contrast, bait_full_id, object_id, object_label, median_log2, median_log2_batch, p_value, p_value_batch) %>%
    dplyr::inner_join(dplyr::select(oeproteome_contrasts_apmsobj.df, bait_full_id, object_id, contrast_oeproteome = contrast,
                                    oeproteome_majority_protein_acs = majority_protein_acs,
                                    oeproteome_gene_names = gene_names,
                                    oeproteome_is_hit_new = is_hit_new, oeproteome_is_upregulation_new = is_upregulation_new,
                                    oeproteome_is_hit_old = is_hit_old, oeproteome_is_upregulation_old = is_upregulation_old,
                                    oeproteome_median_log2 = median_log2, oeproteome_median_log2_batch = median_log2_batch,
                                    oeproteome_p_value = p_value, oeproteome_p_value_batch = p_value_batch,
                                    oeproteome_median_log2_old = median_log2.old, oeproteome_p_value_old = p_value.old)) %>%
    mutate(oeproteome_is_upregulation_new = coalesce(oeproteome_is_upregulation_new, FALSE) &
                                            median_log2 > 0 & (median_log2 - coalesce(oeproteome_median_log2, 0)) <= oeproteome_log2_mindelta,
           oeproteome_is_upregulation_old = coalesce(oeproteome_is_upregulation_old, FALSE) &
               median_log2 > 0 & (median_log2 - coalesce(oeproteome_median_log2_old, 0)) <= oeproteome_log2_mindelta,
           oeproteome_oldnew_change = coalesce(oeproteome_is_hit_new, FALSE) != coalesce(oeproteome_is_hit_old, FALSE)) %>%
    group_by(object_id) %>%
    mutate(oeproteome_is_upregulation_new_any = any(coalesce(oeproteome_is_upregulation_new, FALSE)),
           oeproteome_is_upregulation_old_any = any(coalesce(oeproteome_is_upregulation_old, FALSE)),
           oeproteome_is_hit_new_any = any(coalesce(oeproteome_is_hit_new, FALSE)),
           oeproteome_is_hit_old_any = any(coalesce(oeproteome_is_hit_old, FALSE))) %>%
    ungroup() %>%
    mutate(oeproteome_oldnew_change_any = coalesce(oeproteome_is_hit_new_any, FALSE) != coalesce(oeproteome_is_hit_old_any, FALSE))
    
write_tsv(oeproteome_contrasts_in_apms.df, "/pool/analysis/astukalov/cov2/data/spectronaut_oeproteome_20200527/apms_X_oeproteome_20200608_changes.txt")
