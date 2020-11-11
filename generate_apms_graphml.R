project_id <- 'cov2'
message('Project ID=', project_id)
report_version = "20201104"
apms_info <- list(msfolder = 'mq_apms_20200525',
                  data_ver = "20200525",
                  fit_ver = "20200525")
datasets <- list(
    apms = apms_info,
    prev_apms = list(msfolder = 'mq_apms_20200510',
                data_ver = "20200515",
                fit_ver = "20200515"),
    oeproteome = list(msfolder = 'spectronaut_oeproteome_20200527',
                      data_ver = "20200527",
                      fit_ver = "20200608"),
    phospho = list(msfolder = 'snaut_parsars_phospho_20201005',
                   data_ver = "20201012",
                   fit_ver = "20201012",
                   ptm_extractor_version = "20201012"),
    proteome = list(msfolder = 'snaut_parsars_fp_20200829',
                    data_ver = "20200830",
                    fit_ver = "20200830"),
    ubi = list(msfolder = 'snaut_parsars_ptm_20200907',
               data_ver = "20201012",
               fit_ver = "20201012",
               ptm_extractor_version = "20201012"),
)

cov2_rnaseq = list(data_folder = "parsars_rnaseq_20201020",
                   data_ver = "20201020",
                   fit_ver = "20201020")

message("Assembling fit results for project ", project_id,
        " (APMS dataset v", datasets$apms$data_ver, ", fit v", datasets$apms$fit_ver, ", ",
          "OE proteome dataset v", datasets$oeproteome$data_ver, ", fit v", datasets$oeproteome$fit_ver, ")")

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
    load(file.path(scratch_path, paste0(project_id, '_msglm_data_', envinfo$msfolder, '_', envinfo$data_ver, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msdata_full_', envinfo$msfolder, '_', envinfo$data_ver, '.RData')), envir = ds.env)
    load(file.path(scratch_path, paste0(project_id, '_msglm_fit_', envinfo$msfolder, '_', envinfo$data_ver, '.RData')), envir = ds.env)
}
load(file.path(analysis_path, "networks", paste0(project_id, '_4graph_', datasets$prev_apms$msfolder, '_', datasets$prev_apms$fit_ver, '.RData')), envir = prev_apms.env)
#phospho.env$msdata_full$ptmgroup2psitep <- read_tsv(file.path(data_path, datasets$phospho$msfolder,
#                                                           "COV2_DIA_phospho_0.75probablity_no normalization_psitep_nodata.txt"))
rnaseq.env <- new.env(parent=baseenv())
load(file.path(scratch_path, paste0(project_id, '_', datasets$rnaseq$data_folder, '_', datasets$rnaseq$fit_ver, '.RData')), envir = rnaseq.env)

mscalib_pepmodstate <- instr_calib_pepmodstate
source(file.path(project_scripts_path, 'setup_modelobj.R'))

bait_checks.df <- get(str_c("bait_checks_", modelobj, ".df"))
bait_checks.df$object_id <- bait_checks.df[[modelobj_idcol]]

#ap_plasmids.df <- read_xlsx(file.path(data_path, "AP_plasmids_20200409.xlsx")) %>%
#    extract(`CONSTRUCT NAME`,
#            c("vector", "terminus", "tag", "gene_name", "organism"),
#            "(EC|p\\w+)[ -]([NC](?:-term)?)(?:[ -](myc|Myc|HA|6Myc|Flag|FLAG|V5|1xHA|2xHA))?[-_ ](RIG-I|[^_ -]+)(?:[-_ (](human|mouse|murine)\\)?)?",
#            remove = FALSE)
#write_tsv(crispr_plasmids.df, file.path(data_path, "AP_plasmids_20200409_with_genename.txt"), na = "")

is_in_prev_apms <- str_c("is_in_", datasets$prev_apms$fit_ver, "_apms")
prev_iactions_obj.df <- dplyr::filter(prev_apms.env$iactions_ex_4graphml.df, str_detect(type, "experiment")) %>%
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

proteome_effects_obj.df <- filter(proteome.env$object_contrasts.df, str_detect(contrast, "SARS.+_vs_mock") & std_type == "median") %>%
    dplyr::inner_join(dplyr::select(filter(proteome.env$msdata_full$protein2protgroup, is_majority),
                                    protein_ac, object_id = protgroup_id)) %>%
    dplyr::select(contrast, protein_ac,
                  proteome_median_log2 = median_log2, proteome_p_value = p_value,
                  proteome_is_hit = is_hit_composed, proteome_hit_type = composed_hit_type,
                  proteome_change = change) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::mutate(object_id = protregroup_id) %>%
    tidyr::extract(contrast, "timepoint", "^SARS.+mock@(\\d+)h", remove=FALSE) %>%
    dplyr::mutate(timepoint = parse_integer(timepoint)) %>%
    dplyr::group_by(contrast, object_id, timepoint) %>%
    dplyr::filter(row_number(proteome_p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::ungroup() %>%
    dplyr::arrange(object_id, timepoint, contrast) %>%
    dplyr::mutate(hit_ref = if_else(proteome_is_hit, str_c(timepoint, "h:", proteome_hit_type, "(", proteome_change, ")"), "")) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(proteome_timepoints = if_else(any(proteome_is_hit),
                                                   str_c(unique(hit_ref[proteome_is_hit]), collapse=" "), NA_character_),
                     #proteome_median_log2 = proteome_median_log2[row_number(proteome_p_value) == 1L],
                     #proteome_p_value = proteome_p_value[row_number(proteome_p_value) == 1L],
                     proteome_is_hit = any(proteome_is_hit)) %>%
    dplyr::ungroup()

rnaseq_effects_obj.df <- filter(rnaseq.env$object_contrasts.df, str_detect(contrast, "SARS.+_vs_mock")) %>%
    dplyr::select(contrast, gene_name,
                  rnaseq_mean_log2 = mean_log2, rnaseq_p_value = p_value,
                  rnaseq_is_hit = is_hit_composed, rnaseq_hit_type = composed_hit_type,
                  rnaseq_change = change) %>%
    dplyr::inner_join(dplyr::inner_join(dplyr::filter(msdata_full$protein2protregroup, is_majority),
                                        msdata_full$proteins) %>%
                      dplyr::select(protregroup_id, gene_name) %>% dplyr::distinct()) %>%
    dplyr::mutate(object_id = protregroup_id) %>%
    tidyr::extract(contrast, "timepoint", "^SARS.+mock@(\\d+)h", remove=FALSE) %>%
    dplyr::mutate(timepoint = parse_integer(timepoint)) %>%
    dplyr::group_by(contrast, object_id, timepoint) %>%
    dplyr::filter(row_number(rnaseq_p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::ungroup() %>%
    dplyr::arrange(object_id, timepoint, contrast) %>%
    dplyr::mutate(hit_ref = if_else(rnaseq_is_hit, str_c(timepoint, "h:", rnaseq_hit_type, "(", rnaseq_change, ")"), "")) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(rnaseq_timepoints = if_else(any(rnaseq_is_hit),
                                                   str_c(unique(hit_ref[rnaseq_is_hit]), collapse=" "), NA_character_),
                     #rnaseq_median_log2 = rnaseq_median_log2[row_number(rnaseq_p_value) == 1L],
                     #rnaseq_p_value = rnaseq_p_value[row_number(rnaseq_p_value) == 1L],
                     rnaseq_is_hit = any(rnaseq_is_hit)) %>%
    dplyr::ungroup()

phospho_effects_obj.df <- filter(phospho.env$object_contrasts.df, str_detect(contrast, "^SARS.+_vs_mock") & std_type == "median") %>%
    dplyr::filter(ptmid_is_reference) %>%
    #dplyr::mutate(ptmgroup_shortid = str_remove(ptmgroup_id, "_M\\d+$")) %>%
    #dplyr::inner_join(dplyr::transmute(phospho.env$msdata_full$ptm2gene,
    #                                   ptmgroup_shortid, protein_ac, ptm_gene_name=gene_name, ptm_pos=data_ptm_pos, ptm_label=str_c(data_ptm_AA, data_ptm_pos))) %>%
    dplyr::select(contrast, protein_ac, ptm_AA_seq, ptm_pos, ptm_gene_name=gene_name,
                  phospho_median_log2 = median_log2, phospho_p_value = p_value,
                  phospho_is_hit = is_hit_composed, phospho_hit_type = composed_hit_type,
                  phospho_change = change) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::inner_join(select(msdata_full$protregroups, protregroup_id, object_label=protregroup_label)) %>%
    tidyr::extract(contrast, "timepoint", "^SARS.+mock@(\\d+)h", remove=FALSE) %>%
    dplyr::mutate(ptm_label = str_c(ptm_AA_seq, ptm_pos, "@", ptm_gene_name),
                  object_id = protregroup_id,
                  timepoint = parse_integer(timepoint)) %>%
    dplyr::group_by(contrast, object_id, ptm_label) %>%
    dplyr::filter(row_number(phospho_p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::arrange(object_id, ptm_gene_name, ptm_pos, ptm_label, timepoint, contrast) %>%
    dplyr::mutate(hit_ref = if_else(phospho_is_hit, str_c(timepoint, "h:", phospho_hit_type, "(", phospho_change, ")"), "")) %>%
    dplyr::group_by(object_id, ptm_gene_name, ptm_pos, ptm_label) %>%
    dplyr::summarise(ptm_hits = if_else(any(phospho_is_hit), str_c(unique(hit_ref[phospho_is_hit]), collapse=" "), NA_character_),
                     #phospho_median_log2 = phospho_median_log2[row_number(phospho_p_value) == 1L],
                     #phospho_p_value = phospho_p_value[row_number(phospho_p_value) == 1L],
                     phospho_is_hit = any(phospho_is_hit)) %>%
    dplyr::arrange(object_id, ptm_pos, ptm_label) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(phospho_ptms = if_else(any(phospho_is_hit), str_c(ptm_label[phospho_is_hit],
                                                                       str_c("[", ptm_hits[phospho_is_hit], "]"), collapse=" "), NA_character_),
                     #phospho_median_log2 = phospho_median_log2[row_number(phospho_p_value) == 1L],
                     #phospho_p_value = phospho_p_value[row_number(phospho_p_value) == 1L],
                     phospho_is_hit = any(phospho_is_hit)) %>%
    dplyr::ungroup()

ubi_effects_obj.df <- filter(ubi.env$object_contrasts.df, str_detect(contrast, "^SARS.+_vs_mock") & std_type == "median") %>%
    dplyr::filter(ptmid_is_reference) %>%
    #dplyr::mutate(ptmgroup_shortid = str_remove(ptmgroup_id, "_M\\d+$")) %>%
    #dplyr::inner_join(dplyr::transmute(phospho.env$msdata_full$ptm2gene,
    #                                   ptmgroup_shortid, protein_ac, ptm_gene_name=gene_name, ptm_pos=data_ptm_pos, ptm_label=str_c(data_ptm_AA, data_ptm_pos))) %>%
    dplyr::select(contrast, protein_ac, ptm_AA_seq, ptm_pos, ptm_gene_name=gene_name,
                  ubi_median_log2 = median_log2, ubi_p_value = p_value,
                  ubi_is_hit = is_hit_composed, ubi_hit_type = composed_hit_type,
                  ubi_change = change) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::inner_join(select(msdata_full$protregroups, protregroup_id, object_label=protregroup_label)) %>%
    tidyr::extract(contrast, "timepoint", "^SARS.+mock@(\\d+)h", remove=FALSE) %>%
    dplyr::mutate(ptm_label = str_c(ptm_AA_seq, ptm_pos, "@", ptm_gene_name),
                  object_id = protregroup_id,
                  timepoint = parse_integer(timepoint)) %>%
    dplyr::group_by(contrast, object_id, ptm_label) %>%
    dplyr::filter(row_number(ubi_p_value) == 1L) %>% # pick the most significant if multiple matches
    dplyr::arrange(object_id, ptm_gene_name, ptm_pos, ptm_label, timepoint, contrast) %>%
    dplyr::mutate(hit_ref = if_else(ubi_is_hit, str_c(timepoint, "h:", ubi_hit_type, "(", ubi_change, ")"), "")) %>%
    dplyr::group_by(object_id, ptm_gene_name, ptm_pos, ptm_label) %>%
    dplyr::summarise(ptm_hits = if_else(any(ubi_is_hit), str_c(unique(hit_ref[ubi_is_hit]), collapse=" "), NA_character_),
                     #ubi_median_log2 = ubi_median_log2[row_number(ubi_p_value) == 1L],
                     #ubi_p_value = ubi_p_value[row_number(ubi_p_value) == 1L],
                     ubi_is_hit = any(ubi_is_hit)) %>%
    dplyr::arrange(object_id, ptm_pos, ptm_label) %>%
    dplyr::group_by(object_id) %>%
    dplyr::summarise(ubi_ptms = if_else(any(ubi_is_hit), str_c(ptm_label[ubi_is_hit],
                                                                       str_c("[", ptm_hits[ubi_is_hit], "]"), collapse=" "), NA_character_),
                     #ubi_median_log2 = ubi_median_log2[row_number(ubi_p_value) == 1L],
                     #ubi_p_value = ubi_p_value[row_number(ubi_p_value) == 1L],
                     ubi_is_hit = any(ubi_is_hit)) %>%
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
    left_join(proteome_effects_obj.df) %>%
    left_join(phospho_effects_obj.df) %>%
    left_join(ubi_effects_obj.df) %>%
    left_join(rnaseq_effects_obj.df) %>%
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
                  starts_with("rnaseq"), starts_with("proteome"), starts_with("phospho"), starts_with("ubi"), starts_with("rnaseq"),
                  #starts_with("krogan_"), 
                  starts_with("virhostnet_")#,
                  #crispr_plasmid_ids
                  ) %>%
    dplyr::select(-oeproteome_is_hit_stringent)
write_tsv(iactions_4table.df,
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_all_', apms_info$msfolder, '_', report_version, modelobj_suffix, '.txt')),
          na = "")
write_tsv(filter(iactions_4table.df, is_hit),
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_signif_', apms_info$msfolder, '_', report_version, modelobj_suffix, '.txt')),
          na = "")

require(writexl)

hom_object_contrasts.df <- filter(object_contrasts.df, contrast_type == "comparison" & str_detect(contrast, "_corrected") & std_type == "replicate") %>%
    dplyr::left_join(dplyr::select(bait_checks.df, bait_full_id, bait_homid)) %>%
    dplyr::semi_join(dplyr::inner_join(dplyr::select(dplyr::filter(iactions_4table.df, is_hit), bait_full_id, object_id),
                                       dplyr::select(bait_checks.df, bait_full_id, bait_homid)) %>% dplyr::select(bait_homid, object_id)) %>%
    dplyr::select(gene_names, majority_protein_acs, bait_id, contrast, contrast_offset, median_log2, p_value, change, is_hit = is_hit_nomschecks)
write_tsv(hom_object_contrasts.df,
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_signif_homology_comparison_', apms_info$msfolder, '_', report_version, modelobj_suffix, '.txt')),
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
    left_join(proteome_effects_obj.df) %>%
    left_join(phospho_effects_obj.df) %>%
    left_join(cov2el_proteome_effects_obj.df) %>%
    left_join(cov2el_phospho_effects_obj.df) %>%
    left_join(ubi_effects_obj.df)

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
                                  str_c(project_id, '_', datasets$apms$msfolder, '_', report_version, '.graphml')))
rgraph_filepath <- file.path(analysis_path, 'networks', str_c(project_id, '_4graph_', datasets$apms$msfolder, '_', report_version, '.RData'))
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

