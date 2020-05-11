project_id <- 'cov2'
message('Project ID=', project_id)
datasets <- list(
    apms = list(ms_folder = 'mq_apms_20200427',
                data_version = "20200503",
                fit_version = "20200503"),
    prev_apms = list(ms_folder = 'mq_apms_20200417',
                data_version = "20200420",
                fit_version = "20200420"),
    oeproteome = list(ms_folder = 'spectronaut_oeproteome_20200411',
                      data_version = "20200411",
                      fit_version = "20200411"),
    cov2ts_phospho = list(ms_folder = 'cov2timecourse_phospho_dia_20200423',
                          data_version = "20200428",
                          fit_version = "20200428"),
    cov2ts_proteome = list(ms_folder = 'cov2timecourse_dia_20200423',
                           data_version = "20200429",
                           fit_version = "20200429")
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
    dplyr::select(bait_full_id, protein_ac) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::select(bait_full_id, object_id = protregroup_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(!!is_in_prev_apms := TRUE)

oeproteome_effects_obj.df <- filter(oeproteome.env$object_contrasts.df, str_detect(contrast, "_vs_FPMS_controls") & std_type == "replicate") %>%
    dplyr::mutate(bait_full_id = str_match(contrast, "FPMS_(.+)_vs_FPMS_(controls|others)")[,2]) %>%
    dplyr::inner_join(dplyr::select(filter(oeproteome.env$msdata_full$protein2protgroup, is_majority),
                         protein_ac, object_id = protgroup_id)) %>%
    dplyr::select(bait_full_id, protein_ac,
                  oeproteome_median_log2 = median_log2, oeproteome_p_value = p_value, oeproteome_is_hit = is_hit) %>%
    dplyr::inner_join(filter(msdata_full$protein2protregroup, is_majority)) %>%
    dplyr::select(bait_full_id, object_id = protregroup_id, starts_with("oeproteome_")) %>%
    dplyr::group_by(bait_full_id, object_id) %>%
    dplyr::filter(row_number(oeproteome_p_value) == 1L) %>% # pick the most significant if multiple matches
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
                               !(object_id %in% bait_checks.df$object_id)) %>%
    anti_join(filter(fit_stats$object_batch_effects, var=="obj_batch_effect" & abs(median_log2) >= 1 & p_value <= 0.05) %>%
                  select(object_id) %>% distinct())

iactions_4table.df <- dplyr::inner_join(iactions_4graphml_pre.df,
                                        dplyr::select(bait_checks.df, bait_full_id, bait_id, bait_organism)) %>%
    left_join(krogan_apms_obj.df) %>%
    left_join(virhostnet_ppi_obj.df) %>%
    left_join(crispr_plasmids_obj.df) %>%
    left_join(prev_iactions_obj.df) %>%
    left_join(oeproteome_effects_obj.df) %>%
    left_join(cov2ts_proteome_effects_obj.df) %>%
    left_join(cov2ts_phospho_effects_obj.df) %>%
    left_join(dplyr::select(modelobjs_df, object_id, protein_description, any_of(c("npepmods_unique", "npeptides_unique")))) %>%
    dplyr::arrange(bait_id, bait_organism, bait_full_id, object_id, p_value) %>%
    dplyr::select(bait_organism, bait_id, bait_full_id, contrast,
                  object_label, protein_description,
                  is_hit, median_log2, p_value, median_log2_threshold, p_value_threshold,
                  majority_protein_acs, gene_names, protein_names, is_viral, is_contaminant, is_reverse,
                  any_of(c("npepmods_unique", "npeptides_unique")),
                  nmsruns_quanted = nmsruns_quanted_lhs_max,
                  matches("is_in_\\d+_network"), starts_with("oeproteome_"),
                  starts_with("cov2ts_proteome_"), starts_with("cov2ts_phospho"),
                  starts_with("krogan_"), starts_with("virhostnet_"),
                  crispr_plasmid_ids)

write_tsv(iactions_4table.df,
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_all_', datasets$apms$ms_folder, '_', datasets$apms$fit_version, modelobj_suffix, '.txt')),
          na = "")
write_tsv(filter(iactions_4table.df, is_hit),
          file.path(analysis_path, 'networks',
                    str_c(project_id, '_interactions_signif_', datasets$apms$ms_folder, '_', datasets$apms$fit_version, modelobj_suffix, '.txt')),
          na = "")

bait_labels.df <- distinct(select(iactions_4graphml_pre.df, contrast)) %>%
    dplyr::inner_join(select(filter(contrastXmetacondition.df, condition_role == "signal"), contrast, metacondition)) %>%
    dplyr::inner_join(conditionXmetacondition.df) %>%
    dplyr::inner_join(select(conditions.df, condition, bait_full_id, bait_id, orgcode)) %>%
    dplyr::select(-metacondition, -condition)

iactions_4graphml.df <- filter(iactions_4graphml_pre.df, is_hit) %>%
    dplyr::select(-prob_nonneg) %>%
    dplyr::inner_join(bait_labels.df) %>%
    dplyr::arrange(bait_full_id, bait_id, object_id, prob_nonpos) %>%
    dplyr::group_by(bait_full_id, bait_id, object_label, object_id) %>%
    dplyr::summarize(#conditions.vs_background = paste0(condition, collapse=' '),
                     prob_nonpos_min.vs_background = prob_nonpos[1],
                     median_log2.vs_background = median_log2[1],
                     sd_log2.vs_background = sd_log2[1]) %>%
    dplyr::ungroup() %>%
    left_join(krogan_apms_obj.df) %>%
    left_join(virhostnet_ppi_obj.df) %>%
    left_join(prev_iactions_obj.df) %>%
    left_join(oeproteome_effects_obj.df)

special_bait_ids <- c()

# fix bait mapping to object_id
bait2object.df <- arrange(bait_checks.df, bait_id, bait_full_id, object_id, bait_organism) %>%
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
    dplyr::transmute(bait_full_id, bait_id,
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
                      bait_full_id, object_id),
        modelobjs_df),
        missed_bait_objects.df) %>%
    mutate(object_label = bait_full_id) %>%
    dplyr::filter(bait_full_id %in% msdata$msruns$bait_full_id)

iactions_4graphml.df <- iactions_4graphml.df %>%
    #dplyr::left_join(comparisons_4graphml.df) %>%
    #dplyr::left_join(effects_4graphml.df) %>%
    dplyr::inner_join(dplyr::select(modelobjs_df, object_id, object_label)) %>%
    dplyr::mutate(type = 'experiment',
                  weight = sqrt(pmin(100, -log10(prob_nonpos_min.vs_background))) *
                           sqrt(pmax(0, median_log2.vs_background))) %>%
    dplyr::inner_join(dplyr::select(bait_object_ids.df, src_object_id=object_id, bait_full_id)) %>%
    dplyr::mutate(src_object_label = bait_id) %>%
    dplyr::rename(dest_object_id = object_id, dest_object_label = object_label)

objects_4graphml.df <- dplyr::bind_rows(semi_join(modelobjs_df,
                                                  dplyr::select(iactions_4graphml.df, object_id = dest_object_id)) %>%
                                        anti_join(dplyr::select(bait_object_ids.df, object_id)),
                                        dplyr::select(bait_object_ids.df, -bait_full_id, -bait_id)) %>%
    #dplyr::left_join(dplyr::select(object_batch_effects.df, object_id, batch_effect_p_value = p_value, batch_effect_median_log2 = median_log2)) %>%
    dplyr::mutate(is_bait = object_id %in% bait_object_ids.df$object_id,
                  is_detected = object_id %in% iactions_4graphml.df$dest_object_id,
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
    left_join(cov2ts_phospho_effects_obj.df)

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
    dplyr::left_join(select(bait_checks.df, bait_object_id = object_id, bait_full_id, bait_id)) %>%
    dplyr::mutate(is_bait = bait_object_id == object_id) %>%
    dplyr::select(-bait_object_id) %>% distinct()

# matrix expand PPIs (pairs should be visible by the same bait)
# leave only confident interactions
known_ppi_pairs.df <- dplyr::inner_join(observed_ppi_participants.df, observed_ppi_participants.df,
                                        by=c("interaction_id", "bait_id", 'ppi_type')) %>%
    dplyr::filter(((is_bait.x & !is_bait.y) | ((is_bait.x == is_bait.y) & (object_id.x < object_id.y)))
                  & (protein_ac.x != protein_ac.y)) %>%
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
    dplyr::summarise(bait_ids = paste0(sort(unique(bait_id)), collapse = ' '),
                     iaction_ids = paste0(unique(interaction_id), collapse = ' '),
                     pubmed_ids = paste0(unique(pubmed_id[!is.na(pubmed_id)]), collapse = ' '),
                     has_complex = any(ppi_type == "complex"),
                     n_baits = n_distinct(bait_id),
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
    dplyr::select(filter(objects_4graphml.df, is_bait), gene_label, src_object_id=object_id),
    dplyr::select(filter(objects_4graphml.df, is_bait), gene_label, dest_object_id=object_id),
    by="gene_label") %>%
    dplyr::filter(src_object_id < dest_object_id) %>%
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
                                     edge.attrs = c(`P-value (vs Background)` = 'prob_nonpos_min.vs_background',
                                                    `Enrichment (vs Background)` = 'median_log2.vs_background',
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
                                  str_c(project_id, '_', datasets$apms$ms_folder, '_', datasets$apms$fit_version, '.graphml')))
rgraph_filepath <- file.path(analysis_path, 'networks', str_c(project_id, '_4graph_', datasets$apms$ms_folder, '_', datasets$apms$fit_version, '.RData'))
results_info <- list(project_id = project_id, datasets,
                     modelobj = modelobj, quantobj = quantobj)
message('Saving full analysis results to ', rgraph_filepath, '...')
save(results_info, iactions_4table.df, iactions_ex_4graphml.df, objects_4graphml.df,
     file = rgraph_filepath)
