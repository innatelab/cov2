project_id <- 'cov2'
message('Project ID=', project_id)
report_version = "20200612"
apms_info <- list(ms_folder = 'mq_apms_20200525',
                  data_version = "20200525",
                  fit_version = "20200525")
datasets <- list(
    apms = apms_info,
    prev_apms = list(ms_folder = 'mq_apms_20200510',
                data_version = "20200515",
                fit_version = "20200515"),
    oeproteome = list(ms_folder = 'spectronaut_oeproteome_20200527',
                      data_version = "20200527",
                      fit_version = "20200608"),
    cov2ts_phospho = list(ms_folder = 'cov2timecourse_phospho_dia_20200423',
                          data_version = "20200428",
                          fit_version = "20200428"),
    cov2ts_proteome = list(ms_folder = 'cov2timecourse_dia_20200423',
                           data_version = "20200429",
                           fit_version = "20200429"),
    cov2el_proteome = list(ms_folder = 'cov2earlylate_fp_phos_ubi_dda_20200429',
                           data_version = "20200514",
                           fit_version = "20200514")
)
perseus_datasets <- list(
    cov2el_phospho = list(ms_folder = "cov2earlylate_fp_phos_ubi_dda_20200601",
                          filename="output_Phospho (STY)Sites",
                          fit_version="20200510",
                          ptm="phospho", prefix="pho"),
    cov2el_ubi = list(ms_folder="cov2earlylate_fp_phos_ubi_dda_20200601",
                      filename="output_S0-01_GlyGly (K)Sites",
                      fit_version="20200510",
                      ptm="ubiquitin", prefix="ubi")
)

message("Assembling fit results for project ", project_id,
        " (APMS dataset v", datasets$apms$data_version, ", fit v", datasets$apms$fit_version, ", ",
          "OE proteome dataset v", datasets$oeproteome$data_version, ", fit v", datasets$oeproteome$fit_version, ")")

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
    if (envname == "apms") {
        ds.env = environment()
    } else {
        ds.env <- new.env(parent=baseenv())
        assign(str_c(envname, ".env"), ds.env)
    }
    envinfo <- datasets[[envname]]
    load(file.path(scratch_path, paste0(project_id, '_msglm_data_', envinfo$ms_folder, '_', envinfo$data_version, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msdata_full_', envinfo$ms_folder, '_', envinfo$data_version, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msglm_fit_', envinfo$ms_folder, '_', envinfo$data_version, '.RData')), envir = ds.env)
}
load(file.path(analysis_path, "networks", paste0(project_id, '_4graph_', datasets$prev_apms$ms_folder, '_', datasets$prev_apms$fit_version, '.RData')), envir = prev_apms.env)
cov2ts_phospho.env$msdata_full$ptmgroup2psitep <- read_tsv(file.path(data_path, datasets$cov2ts_phospho$ms_folder,
                                                                     "COV2_DIA_phospho_0.75probablity_no normalization_psitep_nodata.txt"))
for (envname in names(perseus_datasets)) {
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

dplyr::select(ds.env$perseus_analysis_orig.df, id, matches(str_c("^(.+?) (", dsinfo$prefix, "_.+)\\$")))
View(cov2el_phospho.env$perseus_analysis_long.df)

source(file.path(project_scripts_path, 'setup_modelobj.R'))

bait_checks.df <- get(str_c("bait_checks_", modelobj, ".df"))
bait_checks.df$object_id <- bait_checks.df[[modelobj_idcol]]

#ap_plasmids.df <- read_xlsx(file.path(data_path, "AP_plasmids_20200409.xlsx")) %>%
#    extract(`CONSTRUCT NAME`,
#            c("vector", "terminus", "tag", "gene_name", "organism"),
#            "(EC|p\\w+)[ -]([NC](?:-term)?)(?:[ -](myc|Myc|HA|6Myc|Flag|FLAG|V5|1xHA|2xHA))?[-_ ](RIG-I|[^_ -]+)(?:[-_ (](human|mouse|murine)\\)?)?",
#            remove = FALSE)
#write_tsv(crispr_plasmids.df, file.path(data_path, "AP_plasmids_20200409_with_genename.txt"), na = "")

is_in_prev_apms <- str_c("is_in_", datasets$prev_apms$fit_version, "_apms")
prev_iactions_obj.df <- filter(prev_apms.env$iactions_ex_4graphml.df, str_detect(type, "experiment")) %>%
    dplyr::inner_join(select(filter(prev_apms.env$msdata_full$protein2protregroup, is_majority),
                             protein_ac, dest_object_id = protregroup_id)) %>%
    dplyr::select(bait_full_id, protein_ac, prob_nonpos) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::group_by(bait_full_id, object_id = protregroup_id) %>%
    dplyr::summarize(prev_prob_nonpo = min(prob_nonpos, na.rm=TRUE)) %>%
    dplyr::mutate(!!is_in_prev_apms := TRUE)

oeproteome_batch_contrasts.df <- filter(oeproteome.env$object_contrasts.df, std_type == "median" & str_detect(contrast, "_vs_B\\d+_others"))
oeproteome_effects_obj.df <- filter(oeproteome.env$object_contrasts.df, str_detect(contrast, "_vs_controls") & std_type == "median") %>%
    left_join(dplyr::select(oeproteome_batch_contrasts.df, bait_full_id, object_id, std_type, contrast_batch = contrast,
                            median_log2_batch = median_log2, p_value_batch = p_value)) %>%
    dplyr::mutate(bait_full_id = str_match(contrast, "(.+)_vs_(?:controls|others)")[,2],
                  is_hit_stringent = is_hit,
                  is_hit_batch = coalesce(p_value_batch, 1.0) <= 1E-3 & abs(coalesce(median_log2_batch, 0)) >= 0.25,
                  is_signif = p_value <= 1E-3 & abs(median_log2) >= 0.25,
                  is_hit_nomschecks = is_signif & !is_contaminant & !is_reverse,
                  is_hit = is_hit_nomschecks & is_hit_batch & (coalesce(median_log2_batch, 0) * coalesce(median_log2, 0) > 0)) %>%
    dplyr::inner_join(dplyr::select(filter(oeproteome.env$msdata_full$protein2protgroup, is_majority),
                         protein_ac, object_id = protgroup_id)) %>%
    dplyr::select(bait_full_id, protein_ac,
                  oeproteome_median_log2 = median_log2, oeproteome_p_value = p_value,
                  oeproteome_is_hit = is_hit, oeproteome_is_hit_stringent = is_hit_stringent) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::select(bait_full_id, object_id = protregroup_id, starts_with("oeproteome_")) %>%
    dplyr::arrange(bait_full_id, object_id, desc(oeproteome_is_hit), oeproteome_p_value) %>%
    dplyr::group_by(bait_full_id, object_id) %>%
    dplyr::filter(row_number() == 1L) %>% # pick the most significant if multiple matches
    dplyr::ungroup()

cov2ts_proteome_effects_obj.df <- filter(cov2ts_proteome.env$object_contrasts.df, str_detect(contrast, "SARS.+_vs_mock") & std_type == "replicate") %>%
    dplyr::inner_join(dplyr::select(filter(cov2ts_proteome.env$msdata_full$protein2protgroup, is_majority),
                                    protein_ac, object_id = protgroup_id)) %>%
    dplyr::select(contrast, protein_ac,
                  cov2ts_proteome_median_log2 = median_log2, cov2ts_proteome_p_value = p_value, cov2ts_proteome_is_hit = is_hit_nomschecks) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::mutate(object_id = protregroup_id) %>%
    tidyr::extract(contrast, "timepoint", "mock@(\\d+)h", remove=FALSE) %>%
    dplyr::mutate(timepoint = parse_integer(timepoint)) %>%
    dplyr::group_by(contrast, object_id) %>%
    dplyr::filter(row_number(cov2ts_proteome_p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::group_by(object_id) %>%
    dplyr::summarise(cov2ts_proteome_timepoints = if_else(any(cov2ts_proteome_is_hit),
                                                          str_c(sort(timepoint[cov2ts_proteome_is_hit]), "h",
                                                                if_else(cov2ts_proteome_median_log2[cov2ts_proteome_is_hit] > 0, "+", "-"),
                                                                collapse=" "), NA_character_),
                     cov2ts_proteome_median_log2 = cov2ts_proteome_median_log2[row_number(cov2ts_proteome_p_value) == 1L],
                     cov2ts_proteome_p_value = cov2ts_proteome_p_value[row_number(cov2ts_proteome_p_value) == 1L],
                     cov2ts_proteome_is_hit = any(cov2ts_proteome_is_hit)) %>%
    dplyr::ungroup()

cov2el_proteome_effects_obj.df <- filter(cov2el_proteome.env$object_contrasts.df, str_detect(contrast, "SARS.+_vs_mock") & std_type == "replicate") %>%
    dplyr::inner_join(dplyr::select(filter(cov2el_proteome.env$msdata_full$protein2protregroup, is_majority),
                                    protein_ac, object_id = protregroup_id)) %>%
    dplyr::select(contrast, protein_ac, median_log2, p_value, is_hit = is_hit_nomschecks) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::mutate(object_id = protregroup_id) %>%
    tidyr::extract(contrast, "timepoint", "mock@(\\d+)h", remove=FALSE) %>%
    dplyr::mutate(timepoint = parse_integer(timepoint)) %>%
    dplyr::group_by(contrast, object_id) %>%
    dplyr::filter(row_number(p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::group_by(object_id) %>%
    dplyr::summarise(cov2el_proteome_timepoints = if_else(any(is_hit),
                                                          str_c(sort(timepoint[is_hit]), "h",
                                                                if_else(median_log2[is_hit] > 0, "+", "-"),
                                                                collapse=" "), NA_character_),
                     cov2el_proteome_median_log2 = median_log2[row_number(p_value) == 1L],
                     cov2el_proteome_p_value = p_value[row_number(p_value) == 1L],
                     cov2el_proteome_is_hit = any(is_hit)) %>%
    dplyr::ungroup()

cov2ts_phospho_effects_obj.df <- filter(cov2ts_phospho.env$object_contrasts.df, str_detect(contrast, "SARS.+_vs_mock") & std_type == "replicate") %>%
    dplyr::mutate(ptmgroup_shortid = str_remove(ptmgroup_id, "_M\\d+$")) %>%
    dplyr::inner_join(dplyr::transmute(cov2ts_phospho.env$msdata_full$ptmgroup2psitep,
                                       ptmgroup_shortid, protein_ac, ptm_gene_name=gene_name, ptm_pos=data_ptm_pos, ptm_label=str_c(data_ptm_AA, data_ptm_pos))) %>%
    dplyr::select(contrast, protein_ac, ptm_label, ptm_pos, ptm_gene_name,
                  cov2ts_phospho_median_log2 = median_log2, cov2ts_phospho_p_value = p_value, cov2ts_phospho_is_hit = is_hit_nomschecks) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::inner_join(select(msdata_full$protregroups, protregroup_id, object_label=protregroup_label)) %>%
    tidyr::extract(contrast, "timepoint", "mock@(\\d+)h", remove=FALSE) %>%
    dplyr::mutate(ptm_label = if_else(ptm_gene_name != str_remove(object_label, "\\.\\.\\.$"),
                                      str_c(ptm_gene_name, "_", ptm_label), ptm_label),
                  object_id = protregroup_id,
                  timepoint = parse_integer(timepoint)) %>%
    dplyr::group_by(contrast, object_id, ptm_label) %>%
    dplyr::filter(row_number(cov2ts_phospho_p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::arrange(object_id, ptm_pos, ptm_label, timepoint) %>%
    dplyr::group_by(object_id, ptm_pos, ptm_label) %>%
    dplyr::summarise(timepoints = if_else(any(cov2ts_phospho_is_hit), str_c(timepoint[cov2ts_phospho_is_hit], "h",
                                                                            if_else(cov2ts_phospho_median_log2[cov2ts_phospho_is_hit] > 0, "+", "-"),
                                                                            collapse=" "), NA_character_),
                     cov2ts_phospho_median_log2 = cov2ts_phospho_median_log2[row_number(cov2ts_phospho_p_value) == 1L],
                     cov2ts_phospho_p_value = cov2ts_phospho_p_value[row_number(cov2ts_phospho_p_value) == 1L],
                     cov2ts_phospho_is_hit = any(cov2ts_phospho_is_hit)) %>%
    dplyr::arrange(object_id, ptm_pos, ptm_label) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(cov2ts_phospho_ptms = if_else(any(cov2ts_phospho_is_hit), str_c(ptm_label[cov2ts_phospho_is_hit],
                                                                                     "(", timepoints[cov2ts_phospho_is_hit], ")", collapse=" "), NA_character_),
                     cov2ts_phospho_median_log2 = cov2ts_phospho_median_log2[row_number(cov2ts_phospho_p_value) == 1L],
                     cov2ts_phospho_p_value = cov2ts_phospho_p_value[row_number(cov2ts_phospho_p_value) == 1L],
                     cov2ts_phospho_is_hit = any(cov2ts_phospho_is_hit)) %>%
    dplyr::ungroup()

cov2el_phospho_effects_obj.df <- cov2el_phospho.env$perseus_analysis_short.df %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::inner_join(select(msdata_full$protregroups, protregroup_id, object_label=protregroup_label)) %>%
    dplyr::mutate(object_id = protregroup_id) %>%
    dplyr::group_by(timepoint, object_id, ptm_label) %>%
    dplyr::filter(row_number(p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::arrange(object_id, ptm_pos, ptm_label, timepoint) %>%
    dplyr::group_by(object_id, ptm_pos, ptm_label) %>%
    dplyr::summarise(timepoints = if_else(any(is_signif), str_c(timepoint[is_signif], "h",
                                                                if_else(delta_log2[is_signif] > 0, "+", "-"),
                                                                collapse=" "), NA_character_),
                     delta_log2 = delta_log2[row_number(p_value) == 1L],
                     p_value = p_value[row_number(p_value) == 1L],
                     is_signif = any(is_signif)) %>%
    dplyr::arrange(object_id, ptm_pos, ptm_label) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(cov2el_phospho_ptms = if_else(any(is_signif), str_c(ptm_label[is_signif],
                                                                         "(", timepoints[is_signif], ")", collapse=" "), NA_character_),
                     cov2el_phospho_median_log2 = delta_log2[row_number(p_value) == 1L],
                     cov2el_phospho_p_value = p_value[row_number(p_value) == 1L],
                     cov2el_phospho_is_hit = any(is_signif)) %>%
    dplyr::ungroup()

cov2el_ubi_effects_obj.df <- cov2el_ubi.env$perseus_analysis_short.df %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::inner_join(select(msdata_full$protregroups, protregroup_id, object_label=protregroup_label)) %>%
    dplyr::mutate(object_id = protregroup_id) %>%
    dplyr::group_by(timepoint, object_id, ptm_label) %>%
    dplyr::filter(row_number(p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::arrange(object_id, ptm_pos, ptm_label, timepoint) %>%
    dplyr::group_by(object_id, ptm_pos, ptm_label) %>%
    dplyr::summarise(timepoints = if_else(any(is_signif), str_c(timepoint[is_signif], "h",
                                                                if_else(delta_log2[is_signif] > 0, "+", "-"),
                                                                collapse=" "), NA_character_),
                     delta_log2 = delta_log2[row_number(p_value) == 1L],
                     p_value = p_value[row_number(p_value) == 1L],
                     is_signif = any(is_signif)) %>%
    dplyr::arrange(object_id, ptm_pos, ptm_label) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(cov2el_ubi_ptms = if_else(any(is_signif), str_c(ptm_label[is_signif],
                                                                         "(", timepoints[is_signif], ")", collapse=" "), NA_character_),
                     cov2el_ubi_median_log2 = delta_log2[row_number(p_value) == 1L],
                     cov2el_ubi_p_value = p_value[row_number(p_value) == 1L],
                     cov2el_ubi_is_hit = any(is_signif)) %>%
    dplyr::ungroup()

crispr_plasmids.df <- read_tsv(file.path(data_path, "CRISPR_plasmids_20200410.txt")) %>%
    mutate(organism = if_else(is.na(organism), "human", organism),
           library_plasmid_number = if_else(is.na(library_plasmid_number),
                                            str_c("undef#", cumsum(is.na(library_plasmid_number))), 
                                            library_plasmid_number))

crispr_plasmids_obj.df <- inner_join(filter(crispr_plasmids.df, organism == "human" & !str_detect(Gene_gRNA, "^neg_c")),
                             dplyr::select(msdata_full$proteins, gene_name, protein_ac)) %>%
    dplyr::inner_join(modelobj2protein.df) %>%
    dplyr::select(object_id, library_plasmid_number) %>%
    dplyr::distinct() %>%
    dplyr::arrange(object_id, library_plasmid_number) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(crispr_plasmid_ids = str_c(library_plasmid_number, collapse=" ")) %>%
    dplyr::ungroup()

krogan_apms.df <- read_tsv(file.path(data_path, "published_ppi", "interactorsAll-Krogan-mapped.txt")) %>%
    dplyr::rename(protein_ac = Preys, bait_full_id = Bait, is_hit = significant)

krogan_apms_obj.df <- left_join(krogan_apms.df, modelobj2protein.df) %>%
    dplyr::filter(!is.na(object_id)) %>%
    group_by(bait_full_id, object_id) %>%
    summarise(krogan_is_hit = any(is_hit),
              krogan_MIST = min(MIST),
              krogan_avg_spec = max(AvgSpec),
              krogan_fold_change = max(FoldChange))

virhostnet_confidences = c(low_confidence = "low", literature_mining = "pub", virhost = "virhostnet", high_confidence = "high")

virhostnet_ppi.df <- read_tsv(file.path(data_path, "published_ppi", "SARSCoV2011_VirhostNet_combined_interactions_all.txt")) %>%
    transmute(bait_full_id = bait_id, bait_protein_ac = uniprot_ac, protein_ac = Uniprot_Prey,
              virhostnet_confidence = factor(virhostnet_confidences[confidence],
                                             levels = set_names(virhostnet_confidences, NULL), ordered = TRUE),
              virhostnet_references = Reference)

virhostnet_ppi_obj.df <- left_join(virhostnet_ppi.df, modelobj2protein.df) %>%
    filter(!is.na(object_id)) %>%
    group_by(bait_full_id, object_id) %>%
    summarise(virhostnet_confidence = max(virhostnet_confidence),
              virhostnet_references = str_c(sort(unique(virhostnet_references)), collapse=' '))

# contrasts-based interactions filter
iactions_4graphml_pre.df <- filter(object_contrasts.df, str_detect(contrast, "_vs_others") & std_type == "replicate" &
                                   !(object_id %in% bait_checks.df$object_id))

# how much APMS enrichment should be different from oeproteome fold-change enrichment to exclude that
# interaction is actually upregulation
oeproteome_log2_mindelta <- 4

iactions_4table.df <- filter(object_contrasts.df, str_detect(contrast, "_vs_others") & std_type == "replicate") %>%
    dplyr::mutate(bait_full_id_orig = str_remove(bait_full_id, "\\?$")) %>%
    dplyr::inner_join(dplyr::select(bait_checks.df, bait_full_id_orig=bait_full_id, bait_organism)) %>%
    #left_join(krogan_apms_obj.df) %>%
    #left_join(virhostnet_ppi_obj.df) %>%
    #left_join(crispr_plasmids_obj.df) %>%
    #left_join(prev_iactions_obj.df) %>%
    left_join(oeproteome_effects_obj.df) %>%
    #left_join(hom_object_contrasts.df) %>%
    left_join(cov2ts_proteome_effects_obj.df) %>%
    left_join(cov2ts_phospho_effects_obj.df) %>%
    left_join(cov2el_proteome_effects_obj.df) %>%
    left_join(cov2el_phospho_effects_obj.df) %>%
    left_join(cov2el_ubi_effects_obj.df) %>%
    left_join(dplyr::select(modelobjs_df, object_id, protein_description, any_of(c("npepmods_unique", "npeptides_unique")))) %>%
    dplyr::mutate(oeproteome_is_upregulation = coalesce(oeproteome_is_hit_stringent, FALSE) & coalesce(oeproteome_median_log2, 0) > 0 &
                      median_log2 > 0 & (median_log2 - coalesce(oeproteome_median_log2, 0)) <= oeproteome_log2_mindelta,
                  is_hit = is_hit & !oeproteome_is_upregulation
                  ) %>%
    dplyr::arrange(bait_id, bait_organism, bait_full_id, object_id, p_value) %>%
    dplyr::select(bait_organism, bait_id, bait_full_id, contrast,
                  object_id, object_label, protein_description,
                  is_hit, median_log2, p_value, median_log2_threshold, p_value_threshold,
                  contrast_batch, median_log2_batch, p_value_batch, is_hit_nomschecks_batch,
                  contrast_carryover, median_log2_carryover, p_value_carryover, is_hit_nomschecks_carryover,
                  majority_protein_acs, gene_names, protein_names, is_viral, is_contaminant, is_reverse,
                  any_of(c("npepmods_unique", "npeptides_unique")),
                  nmsruns_quanted = nmsruns_quanted_lhs_max,
                  #matches("is_in_\\d+_apms"), matches("prev_prob"),
                  ends_with("_homology"),
                  starts_with("oeproteome_"),
                  starts_with("cov2ts_proteome"), starts_with("cov2ts_phospho"),
                  starts_with("cov2el_proteome"), starts_with("cov2el_phospho"), starts_with("cov2el_ubi"),
                  #starts_with("krogan_"), 
                  starts_with("virhostnet_")#,
                  #crispr_plasmid_ids
                  ) %>%
    dplyr::select(-oeproteome_is_hit_stringent)
write_tsv(iactions_4table.df,
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_all_', apms_info$ms_folder, '_', report_version, modelobj_suffix, '.txt')),
          na = "")
write_tsv(filter(iactions_4table.df, is_hit),
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_signif_', apms_info$ms_folder, '_', report_version, modelobj_suffix, '.txt')),
          na = "")

hom_object_contrasts.df <- filter(object_contrasts.df, contrast_type == "comparison" & str_detect(contrast, "_corrected") & std_type == "replicate") %>%
    dplyr::left_join(dplyr::select(bait_checks.df, bait_full_id, bait_homid)) %>%
    dplyr::semi_join(dplyr::inner_join(dplyr::select(dplyr::filter(iactions_4table.df, is_hit), bait_full_id, object_id),
                                       dplyr::select(bait_checks.df, bait_full_id, bait_homid)) %>% dplyr::select(bait_homid, object_id)) %>%
    dplyr::select(gene_names, majority_protein_acs, bait_id, contrast, contrast_offset, median_log2, p_value, change, is_hit = is_hit_nomschecks)
write_tsv(hom_object_contrasts.df,
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_signif_homology_comparison_', apms_info$ms_folder, '_', report_version, modelobj_suffix, '.txt')),
          na = "")

bait_labels.df <- distinct(select(iactions_4graphml_pre.df, contrast)) %>%
    dplyr::inner_join(select(filter(contrastXmetacondition.df, condition_role == "signal"), contrast, metacondition)) %>%
    dplyr::inner_join(conditionXmetacondition.df) %>%
    dplyr::inner_join(select(conditions.df, condition, bait_full_id, bait_id, orgcode)) %>%
    dplyr::mutate(bait_full_id_orig = str_remove(bait_full_id, "\\?+$")) %>%
    dplyr::select(-metacondition, -condition)

iactions_4graphml.df <- filter(iactions_4graphml_pre.df, is_hit) %>%
    dplyr::select(-prob_nonneg) %>%
    dplyr::inner_join(bait_labels.df) %>%
    #dplyr::group_by(bait_full_id, bait_id, object_label, object_id) %>%
    #dplyr::summarize(#conditions.vs_background = paste0(condition, collapse=' '),
    #                 prob_nonpos_min.vs_background = prob_nonpos[1],
    #                 median_log2.vs_background = median_log2[1],
    #                 sd_log2.vs_background = sd_log2[1]) %>%
    #dplyr::ungroup() %>%
    left_join(krogan_apms_obj.df) %>%
    left_join(virhostnet_ppi_obj.df) %>%
    left_join(prev_iactions_obj.df) %>%
    left_join(oeproteome_effects_obj.df) %>%
    dplyr::mutate(oeproteome_is_upregulation = coalesce(oeproteome_is_hit, FALSE) & coalesce(oeproteome_median_log2, 0) > 0 &
                      median_log2 > 0 & (median_log2 - coalesce(oeproteome_median_log2, 0)) <= oeproteome_log2_mindelta,
                  is_hit = is_hit & !oeproteome_is_upregulation) %>%
    dplyr::arrange(bait_full_id, bait_id, object_id, prob_nonpos)

special_bait_ids <- c()

# fix bait mapping to object_id
bait2object.df <- left_join(bait_labels.df,
                            rename(bait_checks.df, bait_full_id_orig = bait_full_id, bait_id_orig=bait_id)) %>%
    mutate(is_alternative = str_detect(bait_full_id, "\\?$")) %>%
    arrange(bait_id, bait_full_id_orig, desc(is_alternative), object_id, bait_organism) %>%
    group_by(object_id) %>%
    mutate(new_object_id = if_else(row_number() == 1, object_id, NA_integer_)) %>%
    ungroup() %>%
    mutate(object_id = new_object_id, new_object_id = NULL)

# assign negative object_ids to undetected/special bait ORFs
missed_bait_objects.df <- dplyr::filter(bait2object.df, (is.na(object_id) | bait_full_id %in% special_bait_ids | bait_organism != prot_organism) &
                                        !str_detect(bait_full_id, "Ctrl_") &
                                        bait_full_id %in% msdata$msruns$bait_full_id) %>%
    dplyr::arrange(bait_id, bait_full_id) %>%
    dplyr::mutate(object_id = -row_number()) %>%
    dplyr::transmute(bait_full_id, bait_id, bait_homid,
        object_id, protein_acs=protein_ac, majority_protein_acs=protein_ac,
        gene_names=bait_id, gene_label=bait_id,
        protein_names=str_c(bait_id, '_', bait_orgcode), protac_label=protein_ac,
        fasta_headers = NA_character_, n_proteins=1L,
        score = NA_real_, q_value = NA_real_,
        seqlen=NA_integer_, seqlens = NA_character_, mol_weight_kDA = NA_real_,
        seqcov = 0, unique_razor_seqcov = 0,
        is_contaminant = FALSE, is_reverse = FALSE, is_viral = TRUE, is_full_quant = FALSE, is_top_quant = FALSE,
        is_comp2 = FALSE, object_label = bait_full_id, organism=bait_organism)

bait_object_ids.df <- dplyr::bind_rows(dplyr::inner_join(
        dplyr::select(dplyr::filter(bait2object.df, !(bait_full_id %in% special_bait_ids) &
                                                    !str_detect(bait_full_id, "Ctrl_")),
                      bait_full_id, bait_homid, object_id),
        modelobjs_df),
        missed_bait_objects.df) %>%
    mutate(object_label = bait_full_id) %>%
    dplyr::filter(bait_full_id %in% msdata$msruns$bait_full_id)

iactions_4graphml.df <- iactions_4graphml.df %>%
    #dplyr::left_join(comparisons_4graphml.df) %>%
    #dplyr::left_join(effects_4graphml.df) %>%
    dplyr::inner_join(dplyr::select(modelobjs_df, object_id, object_label)) %>%
    dplyr::mutate(type = 'experiment',
                  weight = sqrt(pmin(100, -log10(prob_nonpos))) *
                           sqrt(pmax(0, median_log2))) %>%
    dplyr::inner_join(dplyr::select(bait_object_ids.df, src_object_id=object_id, bait_full_id)) %>%
    dplyr::mutate(src_object_label = bait_id) %>%
    dplyr::rename(dest_object_id = object_id, dest_object_label = object_label)

objects_4graphml.df <- dplyr::bind_rows(semi_join(modelobjs_df,
                                                  dplyr::select(iactions_4graphml.df, object_id = dest_object_id)) %>%
                                        anti_join(dplyr::select(bait_object_ids.df, object_id)),
                                        dplyr::select(bait_object_ids.df, -bait_full_id, -bait_id)) %>%
    dplyr::left_join(select(bait_object_ids.df, object_id, object_bait_full_id=bait_full_id, object_bait_homid=bait_homid)) %>%
    #dplyr::left_join(dplyr::select(object_batch_effects.df, object_id, batch_effect_p_value = p_value, batch_effect_median_log2 = median_log2)) %>%
    dplyr::mutate(is_detected = object_id %in% iactions_4graphml.df$dest_object_id,
                  is_bait = !is.na(object_bait_full_id),
                  #is_transient_contaminant = !is.na(batch_effect_p_value) & batch_effect_p_value <= 1E-5 & batch_effect_median_log2 > 2.5,
                  is_manual_contaminant = FALSE) %>%
    dplyr::mutate(protein_class = case_when(replace_na(is_viral, FALSE) ~ "viral",
                                            is_contaminant ~ "known_contaminant",
                                            #is_transient_contaminant ~ "transient_contaminant",
                                            replace_na(is_manual_contaminant, FALSE) ~ "manual_contaminant",
                                            is_reverse ~ "reverse",
                                            TRUE ~ NA_character_),
                  exp_role = if_else(is_bait, 'bait', 'prey')) %>%
    left_join(crispr_plasmids_obj.df) %>%
    left_join(filter(semi_join(oeproteome_effects_obj.df, filter(iactions_4graphml_pre.df, is_hit)), oeproteome_is_hit) %>%
              dplyr::arrange(object_id, oeproteome_p_value, bait_full_id) %>%
              dplyr::group_by(object_id) %>%
              dplyr::mutate(oeproteome_bait_full_ids = str_c(bait_full_id, collapse=' '),
                            oeproteome_median_log2_max = max(oeproteome_median_log2)) %>%
              dplyr::filter(row_number(oeproteome_p_value) == 1L) %>% dplyr::ungroup() %>%
              dplyr::select(object_id, starts_with("oeproteome_"))) %>%
    dplyr::mutate(oeproteome_is_hit = coalesce(oeproteome_is_hit, FALSE) & oeproteome_median_log2_max >= 1) %>%
    left_join(cov2ts_proteome_effects_obj.df) %>%
    left_join(cov2ts_phospho_effects_obj.df) %>%
    left_join(cov2el_proteome_effects_obj.df) %>%
    left_join(cov2el_phospho_effects_obj.df) %>%
    left_join(cov2el_ubi_effects_obj.df)

# integrate known PPI
ppi.env <- new.env(parent=baseenv())
load(file.path(party3rd_data_path, paste0('iactions_n_complexes_20191118.RData')), envir = ppi.env)

source(file.path(misc_scripts_path, "network_utils.R"))

all_participants.df <- bind_rows(
    ppi.env$interaction_participants.df %>%
        dplyr::transmute(interaction_id = paste0("i",interaction_ix),
                         participant_role, protein_db, protein_ac,
                         ppi_type = "interaction"),
    ppi.env$complex_participants.df %>%
        dplyr::transmute(interaction_id = paste0("c",complex_ix),
                         participant_role = "member",
                         protein_db = "uniprot", protein_ac = uniprotId,
                         ppi_type = "complex")
) %>% separate_rows(interaction_id, sep="\\|") %>% mutate(participant_ix = row_number())

# map protein groups of the AP-MS networks onto the UniProt AC nodes of the known PPI
ppi_participants.df <- bind_rows(
    expand_protgroup_acs(objects_4graphml.df, "majority_protein_acs", id_col="object_id", ac_col="protein_ac") %>%
        dplyr::mutate(ac_match_rank = 1L),
    expand_protgroup_acs(objects_4graphml.df, "majority_protein_acs", id_col="object_id", ac_col="protein_ac") %>%
        dplyr::mutate(protein_ac = str_remove(protein_ac, "-(?:PRO_)?\\d+$"), ac_match_rank = 2L)
) %>% dplyr::group_by(protein_ac) %>% dplyr::filter(row_number(ac_match_rank) == 1L) %>% dplyr::ungroup() %>%
    dplyr::select(-row_ix, -prot_ix) %>%
    dplyr::inner_join(bind_rows(ppi.env$interaction_participants.df %>%
                                    dplyr::transmute(interaction_id = paste0("i",interaction_ix),
                                                     participant_role, protein_db, protein_ac,
                                                     ppi_type = "interaction"),
                                ppi.env$complex_participants.df %>%
                                    dplyr::transmute(interaction_id = paste0("c",complex_ix),
                                                     participant_role = "member",
                                                     protein_db = "uniprot", protein_ac = uniprotId,
                                                     ppi_type = "complex")))

ppi_confidence_classes.df <- read_tsv(file.path(data_path, "../../vgirault_vzvapms/data/intact_iaction_stats_20191118_weighted_VG.txt")) %>%
    mutate(confidence_class = factor(confidence_class, levels = c("low", "medium", "strong"), ordered = TRUE),
           n_pubmeds_required_strong = 1L,
           n_pubmeds_required_medium = 2L,
           n_pubmeds_required_low = 3L,
           n_pubmeds_required = case_when(confidence_class == "strong" ~ n_pubmeds_required_strong,
                                          confidence_class == "medium" ~ n_pubmeds_required_medium,
                                          confidence_class == "low" ~ n_pubmeds_required_low,
                                          TRUE ~ NA_integer_))

observed_ppi_participants.df <- ppi_participants.df %>%
    dplyr::inner_join(bind_rows(
        dplyr::select(iactions_4graphml.df, object_id = dest_object_id, bait_full_id) %>%
            mutate(bait_full_id=as.character(bait_full_id)),
        dplyr::select(filter(bait_checks.df, !is.na(bait_full_id)), object_id, bait_full_id) %>%
            distinct()) %>% distinct()) %>%
    dplyr::left_join(select(bait_checks.df, bait_object_id = object_id, bait_full_id, bait_id, bait_homid)) %>%
    dplyr::mutate(is_bait = bait_object_id == object_id) %>%
    dplyr::select(-bait_object_id) %>% distinct()

# matrix expand PPIs (pairs should be visible by the same or homologous bait)
# leave only confident interactions
known_ppi_pairs.df <- dplyr::inner_join(observed_ppi_participants.df, observed_ppi_participants.df,
                                        by=c("interaction_id", "bait_homid", "bait_id", 'ppi_type')) %>%
    dplyr::filter(((is_bait.x & !is_bait.y) | ((is_bait.x == is_bait.y) & (object_id.x < object_id.y)))
                  & (protein_ac.x != protein_ac.y)
                  & ((bait_full_id.x == bait_full_id.y) |
                     (str_detect(bait_full_id.x, "SARS_CoV") & str_detect(bait_full_id.y, "SARS_CoV")))
                  ) %>%
    dplyr::left_join(select(ppi.env$interaction_info.df, interaction_ix, pubmed_id, is_negative,
                            interaction_type, interaction_detection_method) %>%
                         mutate(interaction_id = str_c("i", interaction_ix))) %>%  
    dplyr::left_join(select(ppi_confidence_classes.df, -n_iactions)) %>%
    dplyr::rename(role_x = participant_role.x, role_y = participant_role.y) %>%
    dplyr::mutate(is_baitprey = ((role_x == "bait" & role_y == "prey") | (role_x == "prey" & role_y == "bait")),
                  is_neutneut = ((role_x == "neutral component" & role_y == "neutral component") |
                                     (role_x == "unspecified role" & role_y == "unspecified role")),
                  is_fluorescence = ((role_x == "fluorescence donor" & role_y == "fluorescence acceptor") |
                                         (role_x == "fluorescence acceptor" & role_y == "fluorescence donor")),
                  is_luminescence = ((role_x == "luminescence donor" & role_y == "luminescence acceptor") |
                                         (role_x == "luminescence acceptor" & role_y == "luminescence donor")),
                  is_valid_iaction = (!replace_na(is_negative, TRUE) & (is_baitprey | is_neutneut | is_fluorescence | is_luminescence)) |
                      (ppi_type == "complex")) %>%
    dplyr::filter(is_valid_iaction) %>%
    dplyr::group_by(src_object_id = object_id.x, dest_object_id = object_id.y) %>%
    dplyr::summarise(bait_ids = paste0(sort(unique(c(bait_full_id.x, bait_full_id.y))), collapse = ' '),
                     iaction_ids = paste0(unique(interaction_id), collapse = ' '),
                     pubmed_ids = paste0(unique(pubmed_id[!is.na(pubmed_id)]), collapse = ' '),
                     has_complex = any(ppi_type == "complex"),
                     n_baits = n_distinct(c(bait_full_id.x, bait_full_id.y)),
                     n_iactions = n_distinct(interaction_id),
                     n_pubmeds = n_distinct(pubmed_id[!is.na(pubmed_id)]),
                     n_pubmeds_strong = n_distinct(pubmed_id[!is.na(pubmed_id) & confidence_class == "strong"]),
                     n_pubmeds_medium = n_distinct(pubmed_id[!is.na(pubmed_id) & confidence_class == "medium"]),
                     n_pubmeds_weak = n_distinct(pubmed_id[!is.na(pubmed_id) & confidence_class == "weak"]),
                     is_pubmeds_ok_strong = n_pubmeds_strong >= n_pubmeds_required_strong[[1]],
                     is_pubmeds_ok_medium = n_pubmeds_medium >= n_pubmeds_required_medium[[1]],
                     is_pubmeds_ok_low = n_pubmeds_weak >= n_pubmeds_required_low[[1]]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(has_complex | is_pubmeds_ok_strong | is_pubmeds_ok_medium | is_pubmeds_ok_low) %>%
    dplyr::mutate(
        confidence_class_max = factor(case_when(is_pubmeds_ok_strong ~ "strong",
                                                is_pubmeds_ok_medium ~ "medium",
                                                is_pubmeds_ok_low ~ "low",
                                                TRUE ~ NA_character_), levels=levels(ppi_confidence_classes.df$confidence_class)),
        ppi_type = case_when(#has_complex & is_enriched & is_gocc ~ "enriched_complex_gocc",
                                       #has_complex & is_enriched ~ "enriched_complex",
                                       has_complex ~ "complex",
                                       !is.na(confidence_class_max) ~ paste0("ppi_", as.character(confidence_class_max)),
                                       TRUE ~ NA_character_))

# virtual interactions between homologous baits
homology_iactions.df <- dplyr::inner_join(
    dplyr::transmute(filter(objects_4graphml.df, is_bait), object_bait_homid, src_object_id=object_id),
    dplyr::transmute(filter(objects_4graphml.df, is_bait), object_bait_homid, dest_object_id=object_id),
    by="object_bait_homid") %>%
    dplyr::filter(src_object_id < dest_object_id) %>%
    dplyr::select(-object_bait_homid) %>%
    dplyr::mutate(is_homology = TRUE, weight = 50.0)

iactions_ex_4graphml.df <- dplyr::full_join(iactions_4graphml.df,
                                            mutate(known_ppi_pairs.df,
                                                   weight = case_when(str_detect(ppi_type, "complex|virhost") ~ 40.0,
                                                                      str_detect(ppi_type, "ppi") ~ 30.0,
                                                                      TRUE ~ NA_real_))) %>%
    dplyr::full_join(homology_iactions.df) %>%
    dplyr::mutate(type = case_when(!is.na(type) & !is.na(ppi_type) ~ paste0('experiment ', ppi_type),
                                   !is.na(type) & is.na(ppi_type) ~ 'experiment',
                                   !is.na(ppi_type) ~ ppi_type,
                                   replace_na(is_homology, FALSE) ~ "homology",
                                   TRUE ~ NA_character_))

require(RGraphML)

iactions.graphml <- GraphML.generate(objects_4graphml.df,
                                     iactions_ex_4graphml.df,
                                     node_id.column = 'object_id', parent_id.column = NA,
                                     source.column = 'src_object_id', target.column = 'dest_object_id',
                                     node.attrs = c(`Protein Group ID` = 'object_id',
                                                    `Majority ACs` = 'majority_protein_acs',
                                                    `Gene Names` = 'gene_names',
                                                    `Experimental Role` = 'exp_role',
                                                    `Protein Name` = 'protein_label',
                                                    `Protein Description` = 'protein_description',
                                                    `Protein Class` = "protein_class",
                                                    `Seq Length` = "seqlen",
                                                    `Detected` = "is_detected",
                                                    `Organism` = "organism",
                                                    `CRISPR Plasmid IDs` = "crispr_plasmid_ids"
                                     ),
                                     edge.attrs = c(`P-value (vs Background)` = 'prob_nonpos',
                                                    `Enrichment (vs Background)` = 'median_log2',
                                                    `Weight` = 'weight',
                                                    `type` = "type",
                                                     #`Known types` = "known_types",
                                                    `Known Interaction IDs` = "iaction_ids",
                                                    `Krogan hit` = "krogan_is_hit",
                                                    `Krogan MIST score` = "krogan_MIST",
                                                    `Krogan Average SC` = "krogan_avg_spec",
                                                    `Krogan Fold Change` = "krogan_fold_change",
                                                    `VirHostNet Confidence` = "virhostnet_confidence",
                                                    `VirHostNet PubMeds` = "virhostnet_references")
                                     #edge.attrs = c( `P-value (vs Background)` = 'p_value_min.vs_background',
                                     #                 `P-value (WT vs Mock)` = 'p_value.SC35MWT',
                                     #                 `P-value (delNS1 vs Mock)` = 'p_value.SC35MdelNS1',
                                     #                 `Enrichment (vs Background)` = 'median_log2.vs_background',
                                     #                 `Enrichment (WT vs Mock)` = 'median_log2.SC35MWT',
                                     #                 `Enrichment (delNS1 vs Mock)` = 'median_log2.SC35MdelNS1',
                                     #                 `Weight` = 'weight',
                                     #                 `type` = "type" )
)

write(iactions.graphml, file.path(analysis_path, 'networks',
                                  str_c(project_id, '_', datasets$apms$ms_folder, '_', report_version, '.graphml')))
rgraph_filepath <- file.path(analysis_path, 'networks', str_c(project_id, '_4graph_', datasets$apms$ms_folder, '_', report_version, '.RData'))
results_info <- list(project_id = project_id, datasets,
                     modelobj = modelobj, quantobj = quantobj)
object_contrasts_slim.df <- select(object_contrasts.df, contrast, contrast_type, std_type,
                                   bait_full_id, bait_id, object_id,
                                   #conditions_lhs, conditions_rhs,
                                   starts_with("is_hit"), starts_with("is_signif"),
                                   starts_with("prob_"), starts_with("median_log2"), starts_with("p_value"), change)
message('Saving full analysis results to ', rgraph_filepath, '...')
save(results_info, iactions_4table.df, iactions_ex_4graphml.df, objects_4graphml.df,
     object_contrasts_slim.df,
     bait2object.df, missed_bait_objects.df, bait_labels.df,
     file = rgraph_filepath)
