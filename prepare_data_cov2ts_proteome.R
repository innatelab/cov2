# SARS-CoV/CoV-2 loading and preparing the viral-protein-overexpressed A549 proteome data
# 
# Author: Alexey Stukalov
###############################################################################

project_id <- 'cov2'
message('Project ID=', project_id)
data_version <- "20200423"
fit_version <- "20200423"
ms_folder <- 'cov2timecourse_dia_20200423'
batch <- 'cov2ts'
message('Dataset version is ', data_version)

source("~/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(misc_scripts_path, 'setup_project_paths.R'))

source(file.path(misc_scripts_path, 'fasta_utils.R'))
source(file.path(misc_scripts_path, 'matrix_utils.R'))
source(file.path(misc_scripts_path, "ggplot_ext.R"))

library(devtools)
devtools::load_all('~/projects/R/maxquantutils/')
devtools::load_all('~/projects/R/msglm')
#require(maxquantUtils)
#require(msglm)
require(dplyr)
require(rjson)
require(stringr)
require(readr)
require(pheatmap)
require(tidyr)
require(broom)

msdata_path <- file.path(data_path, ms_folder)

data_info <- list(project_id = project_id,
                  data_ver = data_version, fit_ver = fit_version,
                  ms_folder = ms_folder,
                  instr_calib_protgroup_filename = "instr_QX8_intensity_protgroup_calib_cov2_20200411_borg.json",
                  quant_type = "intensity", quant_col_prefix = "intensity",
                  pep_quant_type = "intensity")

message('Loading MS instrument calibration data from ', data_info$instr_calib_filename, '...')
instr_calib<- fromJSON(file = file.path(data_path, data_info$instr_calib_protgroup_filename))$instr_calib

source(file.path(project_scripts_path, 'prepare_data_common.R'))

fasta.dfs <- list(
  CoV = read_innate_uniprot_fasta(file.path(data_path, "msfasta/cov_baits_20200331.fasta")),
  human = read_innate_uniprot_fasta(file.path(data_path, "msfasta/uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta")),
  contaminants = read_contaminants_fasta(file.path(data_path, "msfasta/contaminants.fasta"))
)

bad_msruns <- c("FPMS_SARS_CoV2_ORF3_3", "FPMS_SARS_CoV2_ORF7a_2")


#TODO:
#maxquantUtils::read.Spectronaut.ProteinsReport() needs to be changed from read_csv to read_delim


msdata.wide <- read.Spectronaut.ProteinsReport(file.path(msdata_path, "20200422_193114_200422 COVID proteome NO normalization_Report.txt"),
                                               import_data = "quantity")
msdata_colgroups <- attr(msdata.wide, "column_groups")

msdata_full <- list(
  protgroups = msdata.wide[, msdata_colgroups$protgroup],
  protgroup_intensities = pivot_longer.Spectronaut.ProtgroupIntensities(msdata.wide)
)
msdata_full$msruns <- dplyr::select(msdata_full$protgroup_intensities, msrun_ix, raw_file) %>% dplyr::distinct() %>%
  dplyr::arrange(msrun_ix) %>%
  dplyr::mutate(msrun_sn = str_remove(str_remove(raw_file, "^20200418_QX7_MaTa_SA_proteome_A549_"), "(?:_\\d{3,})?.raw$")) %>% 
  tidyr::extract(msrun_sn, c("infection","timepoint", "replicate"), "(.*)_(\\d+)hpi_(\\d+)$", remove=FALSE) %>%
  dplyr::mutate(replicate = parse_integer(replicate),timepoint = parse_integer(timepoint)) %>%
  #left_join(dplyr::select(baits_info.df, bait_type, bait_code, bait_full_id, bait_id, organism, orgcode)) %>%
  mutate(condition = str_c(infection,'_' ,timepoint),
         msrun = str_c(condition, "_", replicate),
         infection = factor(infection),
         timepoint= factor(timepoint),
         condition= factor(condition),
         batch = batch,#if_else(bait_full_id %in% c("SARS_CoV2_E", "SARS_CoV2_M", "SARS_CoV2_N", "SARS_CoV2_NSP15", "SARS_CoV2_NSP16", "SARS_CoV2_ORF3"), 2L, 1L),
         is_used = !msrun %in% bad_msruns
         )

msdata_full$protgroup_intensities <- dplyr::select(msdata_full$protgroup_intensities, -raw_file) %>%
  dplyr::mutate(ident_type = factor(if_else(nevidences > 0L, "By MS/MS", "By matching"), levels = c("By MS/MS", "By matching"))) %>%
  left_join(dplyr::select(msdata_full$msruns, msrun_ix, msrun))
msdata_full <- append_protgroups_info(msdata_full, msdata.wide,
                                      proteins_info = dplyr::bind_rows(
                                        dplyr::mutate(fasta.dfs$CoV, is_viral=TRUE, is_contaminant=FALSE),
                                        dplyr::mutate(fasta.dfs$human, is_viral=FALSE, is_contaminant=FALSE),
                                        dplyr::mutate(fasta.dfs$contaminants, is_viral=FALSE, is_contaminant=TRUE)),
                                      import_columns = c("is_viral", "is_contaminant"))
msdata_full$proteins <- mutate(msdata_full$proteins,
                               protein_ac_noiso = str_remove(protein_ac, "-\\d+$"))

msdata_full$protgroups <- dplyr::mutate(msdata_full$protgroups,
    is_reverse = FALSE,
    gene_label = strlist_label2(gene_names),
    protac_label = strlist_label2(protein_acs),
    protgroup_label = case_when(!is.na(gene_label) ~ gene_label,
                                !is.na(protac_label) ~ protac_label,
                                TRUE ~ str_c('#', protgroup_id)))

# condition = bait
conditions.df <- dplyr::select(msdata_full$msruns, condition, infection, timepoint, ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(infection, timepoint) %>%
  dplyr::mutate(infection = factor(infection)) %>% 
  dplyr::mutate(condition = factor(condition, levels=condition))
#msdata_full$msruns <- left_join(dplyr::select(msdata_full$msruns, -any_of("condition")),
#                               dplyr::select(conditions.df, infection, condition))
distinct_cond_timepoints <- as.numeric(as.character(levels(conditions.df$timepoint)))
for (TP in distinct_cond_timepoints) {
  new_col_name <- paste0('after',TP,'h')
  conditions.df <- conditions.df %>% 
    dplyr::mutate(!!sym(new_col_name) := (as.numeric(as.character(timepoint))>=TP ))
}  


msdata <- msdata_full[c('protgroup_intensities',
                        'msruns', 'protgroups', 'protein2protgroup')]
msdata$protgroup_intensities <- semi_join(msdata$protgroup_intensities, filter(msdata$msruns, is_used))

# setup experimental design matrices
##conditionXeffect_orig.mtx <- model.matrix(
##  ~ 1 + bait_id + bait_id:orgcode,
##  mutate(conditions.df, orgcode = if_else(bait_type == "sample", orgcode, factor("CVHSA2", levels=levels(orgcode)))))
##conditionXeffect.mtx <- conditionXeffect_orig.mtx[, colSums(abs(conditionXeffect_orig.mtx)) != 0 &
##                                                    !str_detect(colnames(conditionXeffect_orig.mtx), "\\(Intercept\\)|bait_idCtrl.+:orgcode")]
##dimnames(conditionXeffect.mtx) <- list(condition = conditions.df$condition,
##                                       effect = colnames(conditionXeffect.mtx))
conditionXeffect.mtx <- model.matrix(
  ~ 1 + (after3h + after6h + after12h + after18h + after24h + after30h)*infection,conditions.df
)

dimnames(conditionXeffect.mtx) <- list(condition = conditions.df$condition,
                                       effect = colnames(conditionXeffect.mtx))
conditionXeffect.mtx <- conditionXeffect.mtx[, setdiff(colnames(conditionXeffect.mtx), "(Intercept)")]

pheatmap(conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE, 
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

##effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
##  dplyr::mutate(orgcode = effect_factor(effect, "orgcode", levels(conditions.df$orgcode), NA),
##                bait_id = effect_factor(effect, "bait_id", levels(conditions.df$bait_id), NA),
##                is_positive = FALSE,#!is.na(bait_id) & is.na(orgcode),
##                prior_mean = 0.0)

effects.df <- tibble(effect=colnames(conditionXeffect.mtx)) %>%
  mutate(infection = effect_factor(effect,'infection',
                                   levels(conditions.df$infection),NA))

for (TP in distinct_cond_timepoints) {
  new_col_name <- paste0('after',TP,'h')
  effects.df <- effects.df %>% 
    dplyr::mutate(!!sym(new_col_name) := (str_detect(effect,paste0('after',TP,'h')) ))
}  
effects.df <- effects.df %>% mutate(timepoint = str_extract(effect,"\\d+h"),
                                    timepoint= str_replace(timepoint,"h","") %>% factor(),
                                    is_positive = FALSE,
                                    prior_mean = 0,
                                    prior_tau = case_when(!is.na(timepoint) & is.na(infection) ~ 0.5,
                                                          !is.na(timepoint) & !is.na(infection) ~ 2.0,
                                                          !is.na(infection) ~ 1.0,
                                                          TRUE ~ 1.0
                                    )
                                    )




conditionXeffect.df <- conditionXeffect_frame(conditionXeffect.mtx, effects.df)
inv_conditionXeffect.mtx <- frame2matrix(conditionXeffect.df,
                                         "condition", "effect", "cond_w",
                                         rows = rownames(conditionXeffect.mtx),
                                         cols = colnames(conditionXeffect.mtx))
pheatmap(inv_conditionXeffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_inv_", ms_folder, "_", fit_version, ".pdf")),
         width = 8, height = 6)

msrunXreplEffect.mtx <- replicate_effects_matrix(
  msdata$msruns,#mutate(msdata$msruns, batch_condition=str_c("B", batch, "_", condition)),
  replicate_col = "replicate", condition_col = "condition")
pheatmap(msrunXreplEffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_exp_design_msruns_", ms_folder, "_", fit_version, ".pdf")),
         width = 16, height = 20)

msrunXreplEffect.df <- as_tibble(as.table(msrunXreplEffect.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)

########my code#############
#TODO 0h?
compound_metaconditions <- c("AllVirus0h")
all_metaconditions <- c(levels(conditions.df$condition), compound_metaconditions)
conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
                                            metacondition = all_metaconditions)
for (cname in levels(conditions.df$condition)) {
  conditionXmetacondition.mtx[cname, cname] <- TRUE
}
conditionXmetacondition.mtx[str_detect(rownames(conditionXmetacondition.mtx), "_0"), "AllVirus0h"] <- TRUE

conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
  dplyr::filter(n != 0) %>% dplyr::select(-n)

##
conds <- levels(conditions.df$condition)
mock_conds <- conds[str_detect(conds,'mock')]
virus_conds <- conds[!str_detect(conds,'mock')]
all_timepoints <- c(distinct_cond_timepoints)

all_contrasts <- character()

#phase 1 virus vs mock
for (TP in all_timepoints) {
  shifted_timepoints <- all_timepoints[all_timepoints>=TP]
  contrast_infection <- virus_conds[str_detect(virus_conds, paste(shifted_timepoints,collapse='|'))]
  contrast_mock <- mock_conds[str_detect(mock_conds, paste0('\\D', TP, '\\D?'))]
  all_contrasts <- c(all_contrasts,paste0(contrast_infection,'_VS_',contrast_mock))
  
}

####skip
#phase 2 del vs wt (same virus)
#infection_all <- levels(conditions.df$infection)
#viruses <- infection_all[!infection_all %in% "Mock"]

#for (virus in viruses) {
  #print(virus)
#  virus_del_conds <- virus_conds %>% str_detect(paste0(virus,'_','del')) %>% virus_conds[.]
#  virus_wt_conds <- virus_conds %>% str_detect(paste0(virus,'_','wt')) %>% virus_conds[.]
  #print(virus_del_conds)
  #print("")
#  for (TP in all_timepoints) {
#    shifted_timepoints <- all_timepoints[all_timepoints>=TP]
#    contrast_del <- virus_del_conds[str_detect(virus_del_conds, paste(shifted_timepoints,collapse='|'))]
#    contrast_wt <- virus_wt_conds[str_detect(virus_wt_conds, paste0('\\D', TP, '\\D?'))]
#    all_contrasts <- c(all_contrasts,paste0(contrast_del,'_VS_',contrast_wt))
#  }
#}


contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)

for (contr in all_contrasts) {
  contrastXmetacondition.mtx[contr,str_split(contr,'_VS_')[[1]]] <- c(1,-1)
  
}


contrastXmetacondition.df <- as.data.frame(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(Freq != 0) %>%
  dplyr::rename(weight = Freq) %>%
  dplyr::mutate(contrast_type = 'comparison',
                condition_role = "signal")

contrastXcondition.df <- as.data.frame(as.table(conditionXmetacondition.mtx)) %>% dplyr::filter(Freq != 0) %>%
  dplyr::select(-Freq) %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition)
###########################
##bait_conditions <- as.character(filter(conditions.df, bait_type == "sample")$condition)
##allminus_metaconditions <- paste0("FPMS_allminus_", unique(filter(conditions.df, bait_type == "sample")$bait_id))
##compound_metaconditions <- c(allminus_metaconditions, "FPMS_controls")
##all_metaconditions <- c(bait_conditions, compound_metaconditions)
##conditionXmetacondition.mtx <- false_matrix(condition = levels(conditions.df$condition),
##                                            metacondition = all_metaconditions)
##for (cname in bait_conditions) {
##  conditionXmetacondition.mtx[cname, cname] <- TRUE
##}
##for (cname in allminus_metaconditions) {
##  bait <- str_remove(cname, "^FPMS_allminus_")
##  conditionXmetacondition.mtx[, cname] <- TRUE
## conditionXmetacondition.mtx[str_detect(rownames(conditionXmetacondition.mtx), str_c("_", bait, "$")), cname] <- FALSE
##}
##conditionXmetacondition.mtx[dplyr::filter(conditions.df, bait_type == "control")$condition, "FPMS_controls"] <- TRUE
##pheatmap(ifelse(conditionXmetacondition.mtx, 1.0, 0.0), cluster_rows=FALSE, cluster_cols=FALSE,
##         filename = file.path(data_path, paste0(project_id, "_metaconditions_", ms_folder, "_", fit_version, ".pdf")),
##         width = 8, height = 6)

##conditionXmetacondition.df <- as_tibble(as.table(conditionXmetacondition.mtx)) %>%
##  dplyr::filter(n != 0) %>% dplyr::select(-n) %>%
##  dplyr::mutate(is_preserved_condition = FALSE)#condition %in% c("Ctrl_NT", "Ctrl_Gaussia_luci"))

##contrasts.df <- bind_rows(
##  transmute(filter(conditions.df, bait_type == "sample"),
##            metacondition_lhs = condition,
##            metacondition_rhs = "FPMS_controls",
##            contrast_type = "comparison"),
##  transmute(filter(conditions.df, bait_type == "sample"),
##            metacondition_lhs = condition,
##            metacondition_rhs = str_c("FPMS_allminus_", bait_id),
##            contrast_type = "comparison"),
##  inner_join(
##    select(filter(conditions.df, bait_type == "sample"), metacondition_lhs = condition, orgcode_lhs = orgcode, bait_id),
##    select(filter(conditions.df, bait_type == "sample"), metacondition_rhs = condition, orgcode_rhs = orgcode, bait_id),
##  ) %>% filter(as.integer(orgcode_lhs) < as.integer(orgcode_rhs)) %>%
##    mutate(contrast_type = "comparison") %>%
##    select(-orgcode_lhs, -orgcode_rhs, -bait_id)) %>%
##  mutate(contrast = str_c(metacondition_lhs, "_vs_", ifelse(str_starts(metacondition_rhs, "allminus"), "others", metacondition_rhs)))

##all_contrasts <- contrasts.df$contrast
##contrastXmetacondition.mtx <- zero_matrix(contrast = all_contrasts, metacondition = all_metaconditions)
##for (i in 1:nrow(contrasts.df)) {
##  contrastXmetacondition.mtx[as.character(contrasts.df$contrast[[i]]),
##                             c(as.character(contrasts.df$metacondition_lhs[[i]]),
##                               as.character(contrasts.df$metacondition_rhs[[i]]))] <- c(1, -1)
##}
#########################
pheatmap(contrastXmetacondition.mtx, cluster_rows=FALSE, cluster_cols=FALSE,
         filename = file.path(data_path, paste0(project_id, "_exp_design_contrasts_", ms_folder, "_", fit_version, ".pdf")),
         width = 11, height = 12)

#contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>%
###  dplyr::rename(weight = n) %>%
##  dplyr::left_join(select(contrasts.df, contrast, contrast_type)) %>%
##  dplyr::mutate(condition_role = if_else(contrast_type == "filter" & weight < 0, "background", "signal"))
contrastXmetacondition.df <- as_tibble(as.table(contrastXmetacondition.mtx)) %>% dplyr::filter(n != 0) %>% 
  dplyr::rename(weight = n) %>%
  dplyr::mutate(contrast_type='comparison',
                condition_role= 'signal')

contrastXcondition.df <- conditionXmetacondition.df %>%
  dplyr::inner_join(contrastXmetacondition.df) %>%
  dplyr::arrange(contrast, contrast_type, metacondition, condition)



# prepare the data to use for MS runs normalization
msrun_stats.df <- msdata$protgroup_intensities %>%
  dplyr::group_by(msrun) %>%
  dplyr::summarise(n_pg_idents = sum(!is.na(ident_type) & ident_type == "By MS/MS"),
                   n_pg_quants = sum(!is.na(intensity))) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(msdata$msruns)
##########
##normalization_sabri
msdata4norm.df <- msdata$protgroup_intensities %>% ungroup() %>%
  dplyr::filter(ident_type == "By MS/MS" & !is.na(intensity)) %>% 
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$protgroup_intensities) %>%
  dplyr::inner_join(select(msdata$msruns, msrun,condition,infection,timepoint)) #%>% #WHY???
  #mutate(timepoint=factor(timepoint))
#filter(!str_detect(msrun,'.*_\\d$')) #Temp

#intensities of (non viral, non cont, non rev) protgroups 
msdata4norm_1.df <- protgroup_intensities_all.df %>% ungroup() %>%
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(protgroup_intensities_all.df) %>%
  dplyr::inner_join(select(msdata$msruns, msrun,condition,infection,timepoint))

#check the ratio of NA to MS in all runs
#and filter those with NAs>25% 
runs_na_25 <- msdata4norm_1.df %>% group_by(msrun) %>% summarise(ratio_na =sum(is.na(intensity))/length(intensity)) %>% ungroup() %>% filter(ratio_na<0.25) 
#%>%group_by(cond) %>% summarise(count=length(ratio_na)) %>%  View()

bad_runs <-  msdata4norm_1.df %>% group_by(msrun) %>% summarise(ratio_na =sum(is.na(intensity))/length(intensity)) %>% ungroup() %>% filter(ratio_na>=0.25) 

#filter out msruns with >25% missing values in msdata4norm.df
msdata4norm_good_runs.df <- msdata4norm.df %>% filter(msrun %in% runs_na_25$msrun)




options(mc.cores=8)


# normalization


## level 1: msruns level

msruns_hnorm <- normalize_experiments(instr_calib,
                                      msdata4norm.df,
                                      quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
                                      cond_col="msrun", condgroup_col='condition',
                                      mcmc.iter = 1000L, stan_method="mcmc",
                                      mcmc.adapt_delta=0.95,
                                      verbose=TRUE, max_objs=200L)


## level 2: condition level using the unregulated proteins
### find unreg. proteins



#impute first
msdata$mschannels <- mutate(msdata$msruns, mstag = "L", mschannel = str_c(msrun, "_", mstag))
msdata$protgroup_tagintensities <- mutate(msdata$protgroup_intensities, mstag = "L", mschannel = str_c(msrun, "_", mstag))
protgroup_tagintensities_all.df <- mutate(msdata$protgroup_intensities, mstag = "L", mschannel = str_c(msrun, "_", mstag))
msdata$mschannels <- msdata$mschannels %>% mutate(msrun_mq=msrun)
msdata$mschannel_stats <- mschannel_statistics(msdata)

#impute intensities
protgroup_int_imputed<- impute_intensities(msdata$protgroup_intensities,
                                           msdata$mschannel_stats) %>%
  dplyr::select(-msrun_mq, -msruns_mq) %>%
  dplyr::mutate(log2_intensity_imputed = log2(intensity_imputed))






#msdata4norm_2.df has imputed values and excluding bad runs
msdata4norm_2.df <- msdata$protgroup_intensities %>% ungroup() %>%
  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_reverse & !is_contaminant & !is_viral)) %>%
  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
  dplyr::inner_join(msdata$protgroup_intensities) %>%
  dplyr::inner_join(select(msdata$msruns, msrun,condition,infection,timepoint)) %>% #WHY???
  dplyr::left_join(protgroup_int_imputed %>% select(protgroup_id,msrun,intensity_imputed)) #%>%
  #filter(msrun %in% runs_na_25$msrun)
#filter(!str_detect(msrun,'.*_\\d$')) #Temp



#apply runs normalization to the intensities
runs_norm.df <- dplyr::left_join(msdata4norm_2.df,dplyr::select(msruns_hnorm, msrun, shift)) %>%
  dplyr::mutate(intensity_norm_samples = intensity*exp(-shift),
                intensity_norm_imp_samples = intensity_imputed*exp(-shift))

# normalize the rows (for each protgroup: intensity/median(mock))
normalized_rows <- runs_norm.df %>% group_by(protgroup_id) %>%
  mutate(intensity_norm_row = intensity_norm_samples / median(intensity_norm_samples[condition=="mock_3"], na.rm=TRUE),
         intensity_norm_imp_row = intensity_norm_imp_samples / median(intensity_norm_imp_samples[condition=="mock_3"], na.rm=TRUE)) %>% ungroup()



### find stable prots

#summarize the percentage of non NA intensities per condition per protgroup
perc_na <- normalized_rows %>% 
  group_by(protgroup_id, condition) %>% 
  summarise(ratio_val_int = sum(!is.na(intensity_norm_row))/length(intensity_norm_row)) %>% 
  ungroup()

#join the 2 DF
normed_filtered.df <- normalized_rows %>% left_join(perc_na %>% group_by(protgroup_id) %>% summarise(all_conds_pres=all(ratio_val_int>=0.75)))  

#list protgroups_id with #peptides>4
#protgroups_4peps.list <- msdata.wide %>% filter(n_peptides>=4) %>% pull(protgroup_id)

#lm

#lm_res.df <- normed_filtered.df  %>% filter(protgroup_id %in% protgroups_4peps.list,all_conds_pres) %>%  group_by(protgroup_id)
lm_res.df <- normed_filtered.df  %>% filter(all_conds_pres) %>%  group_by(protgroup_id) %>%
  do({
    df <- .
    #message("Processing ", df$protgroup_id[[1]])
    lmres <- lm(log2(intensity_norm_imp_row) ~ infection * timepoint, df) 
    tidy(lmres)
  })

filtered_lm_res <- lm_res.df %>% filter(term %>% str_detect('.+timepoint*'),term!='(Intercept)') 
# %>% %>% filter(p.value>0.05)

stable_prots <- filtered_lm_res %>% group_by(protgroup_id) %>% 
  summarise(sig=any(p.value <0.05)) %>% 
  filter(sig==FALSE) %>% ungroup() 

housekeep4norm.df <- msdata4norm.df %>% filter(protgroup_id %in% stable_prots$protgroup_id) 
#housekeep4norm.df <- msdata4norm.df %>% filter(protgroup_id %in% c(60:100))



housekeep_hnorm <- normalize_experiments(instr_calib,
                                         housekeep4norm.df,
                                         quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
                                         cond_col='condition', condgroup_col='infection', 
                                         mcmc.iter = 1000L, stan_method="mcmc",
                                         mcmc.adapt_delta=0.95,
                                         verbose=TRUE, max_objs=500L,
                                         mschan_preshifts = msruns_hnorm,
                                         
                                         preshift_col = 'shift')




## level 3: condition @ timepoint=0h

cond_3h_hnorm <- normalize_experiments(instr_calib,
                                       msdata4norm.df %>% filter(timepoint==3),
                                       quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
                                       cond_col="condition", condgroup_col=NULL, 
                                       mcmc.iter = 1000L, stan_method="mcmc",
                                       mcmc.adapt_delta=0.95,
                                       verbose=TRUE, max_objs=500L,
                                       mschan_preshifts = left_join(msruns_hnorm,mutate(housekeep_hnorm,shift2=shift) %>%
                                                                      select(condition,shift2),by='condition') %>% 
                                         mutate(sum_shifts =shift+shift2),
                                       
                                       preshift_col = 'sum_shifts')



## join all levels into a df


msruns_hnorm_all.df <- msruns_hnorm %>% left_join(mutate(housekeep_hnorm,shift2=shift) %>% 
                                                    select(condition,shift2)) %>% 
  left_join(mutate(cond_3h_hnorm, shift3=shift) %>% 
              select(condition,shift3))


#propagate shift3
msruns_hnorm_all.df <- msruns_hnorm_all.df %>% mutate(treatment=str_replace(condition, "_\\d\\d?$", "")) %>% 
  group_by(treatment) %>% mutate(shift3_prop=mean(na.omit(shift3))) %>% ungroup() 



msruns_hnorm_all.df <- msruns_hnorm_all.df %>% mutate(total_msrun_shift=rowSums(msruns_hnorm_all.df %>%   
                                                                                  select(shift,shift2,shift3_prop)   ))



total_msrun_shifts.df <- dplyr::inner_join(msruns_hnorm_all.df,
                                           select(msdata$msruns, msrun)) 


## apply normalization

msdata$protgroup_intensities <- dplyr::left_join(msdata$protgroup_intensities,
                                                 dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  left_join(protgroup_int_imputed ) %>% 
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) 
#dplyr::select(-total_msrun_shift)




global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns, TRUE) %>%
                                                 dplyr::select(msrun) %>% dplyr::distinct() %>%
                                                 dplyr::inner_join(msdata$protgroup_intensities) %>% .$intensity), na.rm=TRUE)



#msdata$protgroup_intensities_all <- impute_intensities(protgroup_tagintensities_all.df,
#                                                       msdata$mschannel_stats) %>%
#  dplyr::select(-mschannel, -mstag) %>%
#  dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
#  dplyr::left_join(dplyr::select(msdata$msruns, msrun, condition)) %>%
#  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
#                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
#                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
#                is_imputed = is.na(intensity)) %>%
#  dplyr::arrange(msrun, protgroup_id)



###############skip###############
# normalize using the intensities
##msdata4norm.df <- msdata$protgroup_intensities %>%
##  dplyr::filter(npeptides_quanted > 1 & !is.na(intensity)) %>%
##  dplyr::semi_join(dplyr::filter(msdata$protgroups, !is_contaminant & !is_viral)) %>%
##  dplyr::select(protgroup_id) %>% dplyr::distinct() %>%
##  dplyr::inner_join(msdata$protgroup_intensities)

##options(mc.cores=8)

# normalize experiments:
# 1) MS replicates for a given bait
# 2) same viral protein of different strains
# 3) all baits together
##msruns_hnorm <- multilevel_normalize_experiments(instr_calib_protgroup,
##                                                 filter(msdata$msruns, is_used) %>%
##                                                   mutate(batch=as.character(batch),
##                                                          batch_bait_full_id = str_c("B", batch, "_", bait_full_id),
##                                                          batch_bait_id = str_c("B", batch, "_", bait_id)),
##                                                 msdata4norm.df,
##                                                 quant_col = "intensity", obj_col = "protgroup_id", mschan_col = "msrun",
##                                                 mcmc.iter = 2000L,
                                                 #mcmc.chains = 6,
##                                                 verbose=TRUE,
##                                                 norm_levels = list(msrun = list(cond_col = "msrun", max_objs=700L, missing_exp.ratio=0.1),
##                                                                    bait_full_id = list(cond_col="batch_bait_full_id", max_objs=500L, missing_exp.ratio=0.1),
##                                                                    bait_id = list(cond_col="batch_bait_id", max_objs=300L, missing_exp.ratio=0.1),
##                                                                    batch = list(cond_col="batch", max_objs=200L, missing_exp.ratio=0.1)
##                                                 ))

# ignore all higher levels of normalization
##msruns_hnorm$msruns_shifts <- msruns_hnorm$mschannel_shifts

##total_msrun_shifts.df <- msruns_hnorm$msruns_shifts

# apply normalization
##msdata_full$protgroup_intensities <- dplyr::left_join(msdata_full$protgroup_intensities,
##                                                      dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
##  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift))

# use only native medium to calculate average intensity
##global_protgroup_labu_shift <- 0.95*median(log(dplyr::filter(msdata$msruns) %>%
##                                                 dplyr::select(msrun) %>% dplyr::distinct() %>%
##                                                 dplyr::inner_join(msdata_full$protgroup_intensities) %>% pull(intensity)), na.rm=TRUE)
############################
#batcheffect
msrunXbatchEffect.mtx <- model.matrix(~ 1,
                                      data=dplyr::mutate(msdata$msruns, batch=factor(batch)))
msrunXbatchEffect.mtx <- msrunXbatchEffect.mtx[, colnames(msrunXbatchEffect.mtx) != "(Intercept)", drop=FALSE]
dimnames(msrunXbatchEffect.mtx) <- list(msrun = msdata$msruns$msrun,
                                        batch_effect = colnames(msrunXbatchEffect.mtx))

batch_effects.df <- tibble(batch_effect = character(0),#colnames(msrunXbatchEffect.mtx),
                           is_positive  = logical(0))

##msrunXbatchEffect_orig.mtx <- model.matrix(
## ~ 1 + batch,
##  mutate(msdata$msruns, batch = factor(batch)))
##msrunXbatchEffect.mtx <- msrunXbatchEffect_orig.mtx[, colnames(msrunXbatchEffect_orig.mtx) != "(Intercept)", drop=FALSE]
##dimnames(msrunXbatchEffect.mtx) <- list(msrun = msdata$msruns$msrun,
##                                       batch_effect = colnames(msrunXbatchEffect_orig.mtx)[colnames(msrunXbatchEffect_orig.mtx) != "(Intercept)"])
##pheatmap(msrunXbatchEffect.mtx, cluster_rows=FALSE, cluster_cols=FALSE)

##batch_effects.df <- tibble(batch_effect=colnames(msrunXbatchEffect.mtx),
##                           is_positive=FALSE,
##                           prior_mean = 0.0)

# no subbatch effects so far
##msrunXsubbatchEffect.mtx <- zero_matrix(msrun = rownames(msrunXreplEffect.mtx),
##                                        subbatch_effect = c())

##subbatch_effects.df <- tibble(subbatch_effect=character(0),
##                              is_positive=logical(0),
##                              prior_mean = double(0))

##bait_checks_protgroup.df <- dplyr::left_join(dplyr::select(baits_info.df, bait_full_id, bait_id, bait_orgcode=orgcode, bait_organism=organism,
##                                                           protein_ac = used_uniprot_ac),
##                                             dplyr::select(msdata_full$proteins, protein_ac, protgroup_id, prot_organism=organism)) %>%
##  dplyr::left_join(dplyr::select(dplyr::filter(msdata$protgroup_intensities, ident_type=="By MS/MS"), protgroup_id, msrun)) %>%
##  dplyr::left_join(dplyr::select(msdata$msruns, msrun, observing_bait_full_id = bait_full_id)) %>%
##  dplyr::arrange(bait_full_id, protgroup_id, msrun) %>%
##  dplyr::group_by(bait_full_id, protgroup_id) %>%
##  dplyr::mutate(idented_in_msruns = str_c(unique(msrun), collapse=";"),
##                idented_in_FP_of = str_c(unique(observing_bait_full_id), collapse=";"),
##                prot_organisms = str_c(unique(prot_organism), collapse=';')) %>%
##  dplyr::filter(row_number()==1L) %>%
##  dplyr::select(-msrun, -observing_bait_full_id) %>%
##  dplyr::ungroup() %>%
##  dplyr::mutate(idented_in_msruns = if_else(idented_in_msruns == "", NA_character_, idented_in_msruns),
##                idented_in_FP_of = if_else(idented_in_FP_of == "", NA_character_, idented_in_FP_of))

rmsglmdata_filepath <- file.path(scratch_path, str_c(project_id, '_msglm_data_', ms_folder, '_', data_version, '.RData'))
message('Saving MS data for MSGLM to ', rmsglmdata_filepath, '...')
save(data_info, msdata,
     conditions.df, effects.df,
     conditionXeffect.mtx, inv_conditionXeffect.mtx, conditionXeffect.df,
     conditionXmetacondition.mtx, conditionXmetacondition.df,
     contrastXmetacondition.mtx, contrastXmetacondition.df, contrastXcondition.df,
     instr_calib,
     global_protgroup_labu_shift,
     msruns_hnorm, total_msrun_shifts.df,
     msrunXreplEffect.mtx,
     batch_effects.df, msrunXbatchEffect.mtx,
     #subbatch_effects.df, msrunXsubbatchEffect.mtx,
     #bait_checks_protgroup.df,
     file = rmsglmdata_filepath)

rfulldata_filepath <- file.path(scratch_path, str_c(project_id, '_msdata_full_', ms_folder, '_', data_version, '.RData'))
message('Saving full MS data to ', rfulldata_filepath, '...')
save(data_info, msdata_full,
     file = rfulldata_filepath)

message('Done.')

msdata_full$mschannels <- mutate(msdata_full$msruns,
                                 mstag = "Sum",
                                 mschannel = msrun,
                                 msrun_mq = msrun)
msdata_full$protgroup_intensities <- mutate(msdata_full$protgroup_intensities,
                                            ident_type="By MS/MS")
msdata_full$protgroup_tagintensities <- mutate(msdata_full$protgroup_intensities,
                                               mschannel = msrun,
                                               msrun_mq = msrun,
                                               mstag = "Sum")
msdata_full$mschannel_stats <- mschannel_statistics(msdata_full)
set.seed(1232)
msdata_full$protgroup_intensities_all <- tidyr::expand(msdata_full$protgroup_intensities,
                                                protgroup_id, msrun) %>%
  left_join(dplyr::select(msdata_full$protgroup_intensities, -any_of("total_msrun_shift"))) %>%
  dplyr::mutate(mstag = "Sum") %>%
  impute_intensities(msdata_full$mschannel_stats) %>%
  dplyr::mutate(total_msrun_shift=0) %>%
  #dplyr::left_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::mutate(intensity_norm = intensity*exp(-total_msrun_shift),
                intensity_imputed_norm = intensity_imputed*exp(-total_msrun_shift),
                log2_intensity_imputed_norm = log2(intensity_imputed_norm),
                is_imputed = is.na(intensity)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities4pca.df <- msdata_full$protgroup_intensities_all %>%
  #filter(mstag == "L") %>%
  dplyr::semi_join(dplyr::filter(msdata_full$protgroups, !is_reverse & !is_contaminant)) %>%
  dplyr::semi_join(dplyr::filter(msdata_full$msruns, is_used)) %>%
  dplyr::arrange(msrun, protgroup_id)

protgroup_intensities.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_norm),
                                    nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                    dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                    msrun = unique(protgroup_intensities4pca.df$msrun)))
protgroup_intensities_imp.mtx <- matrix(log2(protgroup_intensities4pca.df$intensity_imputed_norm),
                                        nrow = n_distinct(protgroup_intensities4pca.df$protgroup_id),
                                        dimnames = list(protgroup = unique(protgroup_intensities4pca.df$protgroup_id),
                                                        msrun = unique(protgroup_intensities4pca.df$msrun)))
require(FactoMineR)
require(Matrix)
msrun_intensities_pca <- PCA(protgroup_intensities_imp.mtx, graph = FALSE)
msrun_intensities_pca.df <- as.data.frame(msrun_intensities_pca$svd$V)
colnames(msrun_intensities_pca.df) <- paste0("comp_", 1:ncol(msrun_intensities_pca.df))
msrun_intensities_pca.df <- dplyr::mutate(msrun_intensities_pca.df,
                                          msrun = rownames(msrun_intensities_pca$var$coord)) %>%
    dplyr::inner_join(msdata_full$msruns)

require(ggrepel)
p <- ggplot(msrun_intensities_pca.df %>% arrange(infection,timepoint) %>% 
            dplyr::mutate(msrun = str_remove(msrun, "FPMS_(?:SARS_)?")),
       aes(x=comp_1, y=comp_2, color=condition)) +
    geom_point() +
    geom_text_repel(aes(label=str_remove(str_remove(msrun, "FPMS_SARS_"), "APMS_")), size=3, vjust=-1.1) +
    theme_bw_ast(base_family = "", base_size = 10) #+
ggsave(p, filename = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_msruns_pca_", data_version, ".pdf")),
       width = 12, height = 14)

require(pheatmap)
msruns_ordered <- filter(msdata_full$msruns, msrun %in% colnames(protgroup_intensities_imp.mtx)) %>%
  dplyr::arrange(infection, replicate) %>% pull(msrun)
protgroup_hclu = hclust(dist(protgroup_intensities_imp.mtx))
pheatmap(protgroup_intensities.mtx[, msruns_ordered], cluster_cols=FALSE, cluster_rows=protgroup_hclu,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_heatmap_intensity.pdf")), width=20, height=100)

viral_mask <- unique(protgroup_intensities4pca.df$protgroup_id) %in% filter(msdata_full$protgroups, is_viral)$protgroup_id
viral_intensities.mtx <- protgroup_intensities.mtx[viral_mask, msruns_ordered]
rownames(bait_intensities.mtx) <- 
  tibble(protgroup_id = unique(protgroup_intensities4pca.df$protgroup_id)[viral_mask]) %>%
  inner_join(msdata_full$protgroups) %>% inner_join(msdata_full$protein2protgroup) %>%
  inner_join(mutate(msdata_full$proteins, full_code = str_c(str_replace(organism, "CoV-2", "CoV2"), '-', gene_name))) %>%
  group_by(protgroup_id) %>%
  summarise(full_code = str_c(full_code, collapse=';')) %>% ungroup() %>% .$full_code
pheatmap(bait_intensities.mtx, cluster_cols=FALSE, cluster_rows=FALSE,
         file = file.path(analysis_path, "plots", ms_folder, paste0(project_id, "_", ms_folder, "_", data_version, "_viral_heatmap_intensity.pdf")), width=16, height=8)

require(uwot)
require(matrixStats)
protgroup_intensities_imp_norm.mtx <- protgroup_intensities_imp.mtx - rowMedians(protgroup_intensities_imp.mtx)
pg_umap2d = umap(protgroup_intensities_imp_norm.mtx, n_components=2,
                 n_neighbors=20, init="laplacian",
                 min_dist=0.2, metric = "euclidean")

objects_umap.df <- as_tibble(pg_umap2d, .name_repair="minimal") %>%
  set_names(c("x_2d", "y_2d")) %>%
  mutate(object_id = unique(protgroup_intensities4pca.df$protgroup_id)) %>%
  #dplyr::left_join(
  #  as_tibble(objects.umap3d, .name_repair="minimal") %>%
  #    set_names(c("x_3d", "y_3d", "z_3d")) %>%
  #    mutate(object_id = unique(object_tagintensities4pca.df$object_id))    
  #) %>%
  dplyr::inner_join(select(msdata$protgroups, object_id=protgroup_id, object_label=protgroup_label,
                           majority_protein_acs, gene_names, is_viral, is_contaminant, is_reverse)) %>%
  mutate(category = case_when(replace_na(is_viral, FALSE) ~ "viral",
                              is_reverse ~ "reverse",
                              is_contaminant ~ "contaminant",
                              is.na(gene_names) ~ "NA",
                              str_detect(gene_names, "(^|;)RP[SL]\\d+") ~ "ribosome",
                              str_detect(gene_names, "(^|;)MRP[SL]\\d+") ~ "mito ribosome",
                              str_detect(gene_names, "(^|;)PSM[ABCD]\\d+") ~ "proteasome",
                              str_detect(gene_names, "(^|;)(MAVS|OAS\\d+|IFIT\\d+|PARP\\d+|DDX58|IFI\\d+|ISG\\d+|MX\\d+|TGFB\\d+|UNC93B1|ZC3HAV1)(;|$)") ~ "ISG",
                              TRUE ~ "default"
  ))# %>%
  #dplyr::left_join(filter(protregroup_effects.df, effect == "treatmentSC35M"))
  # FIXME duplicate point due to N-to-N protgroup to protregroup matching
  #dplyr::left_join(filter(object_effects_wide.df, std_type == "replicate") %>%
   #                  select(-gene_names, -majority_protein_acs))
require(plotly)
# plot umap + selected proteins (color = protein biological function)
umap <- ggplotly(
  ggplot(objects_umap.df,
         aes(x=x_2d, y=y_2d, color=category)) +
    geom_point(aes(text=object_label, size=!(category %in% c('default', 'NA')))) +
    #geom_text(data=filter(objects_umap.df, gene_names %in% c("RPL11", "RPL12")),
    #          aes(label=object_label), vjust=-1.1, color = "black", size = 3) +
    #geom_text_repel(data=filter(objects_umap.df, is_viral),
    #                aes(label=object_label)) +
    scale_size_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5)) +
    scale_color_manual(values = c("viral" = "red", "ISG" = "blue", "hit" = "orange",
                                  "ribosome" = "darkgreen", "mito ribosome" = "darkred", "proteasome" = "khaki",
                                  "default" = "gray",
                                  "NA" = "pink", "reverse" = "pink", "contaminant" = "yellow")) +
    theme_bw_ast(base_family = "", base_size = 10),
  tooltip = "text"
)

