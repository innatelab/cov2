Sys.setenv(TZ='Etc/GMT+1') # issue#612@rstan
#job.args <- c("cov2", "ast_parsars_ptm", "snaut_parsars_ptm_20200907", "20200920", "20200920", "0", "39")
#filter(mutate(modelobjs_df, rn=row_number()), str_detect(ptmn_label, "^GlyGly_SQSTM1_K13_M1"))$rn
if (!exists('job.args')) {
  job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message('Project ID=', project_id)

job_name <- as.character(job.args[[2]])
job_version <- job.args[[5]]
msfolder <- job.args[[3]]
data_version <- job.args[[4]]
fit_version <- job_version
job_id <- as.integer(job.args[[6]])
job_chunk <- as.integer(job.args[[7]])#job_chunk <- 5505L
message('Job ', job_name, '(id=', job_id, '_', job_chunk,
        " data_version=", data_version, " fit_version=", fit_version, " running on ", Sys.info()["nodename"], ")")

#source("~/R/config.R")
source("/projects/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

rdata_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_data_', msfolder, '_', data_version, '.RData'))
message('Loading data from ', rdata_filepath)
load(rdata_filepath)

if (Sys.getenv('SLURM_CPUS_PER_TASK') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) # SLURM way
} else if (Sys.getenv('NSLOTS') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('NSLOTS')) # SGE way
} else {
  mcmc_nchains <- 8
}

require(rlang)
require(dplyr)
require(msglm)
require(rstan)
require(tidyr)
require(stringr)
require(maxquantUtils)

modelobj <- "ptmn"
quantobj <- "pepmodstate"
source(file.path(project_scripts_path, 'setup_modelobj.R'))

sel_object_ids <- modelobjs_df[[modelobj_idcol]][[job_chunk]]
message(sel_object_ids, " ", modelobj, " ID(s): ",
        paste0(sort(unique(dplyr::pull(modelobjs_df[modelobjs_df[[modelobj_idcol]] %in% sel_object_ids, ], object_label))), collapse=' '))
sel_ptm_type <- modelobjs_df$ptm_type[[job_chunk]]
msdata.df <- dplyr::filter(msdata$ptmn2pepmodstate, ptmn_id %in% sel_object_ids) %>%
  dplyr::inner_join(dplyr::select(msdata$pepmodstate_intensities, pepmodstate_id, msrun,
                                  intensity=intensity_norm, qvalue)) %>% # !!! use SN-normalized intensity as the raw intensity
  dplyr::inner_join(dplyr::select(msdata$ptmn_locprobs, ptmn_id, ptm_locprob, pepmodstate_id, msrun)) %>%
  dplyr::inner_join(dplyr::select(msdata$msruns, msrun, condition, treatment, timepoint)) %>%
  dplyr::filter((coalesce(qvalue, 1) <= data_info$qvalue_max) &
                (coalesce(ptm_locprob, 0) >= data_info$locprob_min)) %>% # filter valid PTM localization and identification probabilities
  dplyr::mutate(object_id = ptmn_id)
message('Preparing MS GLM data...')
model_data <- list()
model_data$mschannels <- dplyr::select(msdata$msruns, dataset, condition, msrun, replicate) %>%
  dplyr::filter(dataset == case_when(sel_ptm_type == "GlyGly" ~ "ubi",
                                     sel_ptm_type == "Phospho" ~ "phospho",
                                     TRUE ~ "unknown")) %>%
  #dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::arrange(condition, replicate, msrun) %>% dplyr::distinct() %>%
  dplyr::mutate(mschannel_ix = row_number(),
                msrun_ix = as.integer(factor(msrun, levels=unique(msrun))),
                msproto_ix = 1L,
                zero_msrun_shift = 0L)
experiment_shift_col <- 'zero_msrun_shift' # !!! no normalization shifts since normalized data is used
model_data$mschannels$model_mschannel_shift <- model_data$mschannels[[experiment_shift_col]]
model_data$conditions <- conditions.df %>%
  mutate(condition_ix = row_number())

model_data$interactions <- tidyr::crossing(object_id = sel_object_ids,
                                           condition_ix = model_data$conditions$condition_ix) %>%
    dplyr::inner_join(dplyr::select(model_data$conditions, condition_ix, condition) %>% dplyr::distinct()) %>%
    dplyr::left_join(msdata.df %>% dplyr::select(condition, object_id) %>% dplyr::distinct() %>%
                     dplyr::mutate(is_virtual = FALSE)) %>%
    dplyr::mutate(is_virtual = is.na(is_virtual),
                  iaction_id = paste(condition, object_id))

model_data$interactions <- dplyr::arrange(model_data$interactions, condition_ix, object_id) %>%
  dplyr::mutate(glm_iaction_ix = row_number(),
                glm_object_ix = as.integer(factor(object_id)))

model_data$objects <- dplyr::transmute(model_data$interactions, glm_object_ix, object_id, is_underdefined=FALSE) %>%
  dplyr::distinct() %>% dplyr::arrange(glm_object_ix) %>%
  dplyr::inner_join(dplyr::select(modelobjs_df, ptm_id, ptmn_id, ptmn_label, nselptms, object_id, object_label, ptmn_label_no_ptm_type)) %>%
  dplyr::inner_join(dplyr::select(filter(msdata$ptm2gene, ptm_is_reference), ptm_id, protein_ac, ptm_pos, ptm_AA_seq, ptm_type, contains("is_"))) %>%
  dplyr::inner_join(dplyr::select(msdata$proteins, protein_ac, gene_name=genename, protein_name, contains("is_"))) %>%
  #dplyr::select(-is_fit) %>%
  dplyr::arrange(glm_object_ix)

# arrange pepmodstates by object, by profile cluster and by the number of quantitations
model_data$subobjects <- msdata.df %>%
    dplyr::inner_join(msdata$ptmn2pepmodstate) %>%
    dplyr::group_by(ptmn_id, pepmodstate_id) %>%
    dplyr::summarise(n_quant = sum(!is.na(intensity)),
                     intensity_med = median(intensity, na.rm=TRUE)) %>% ungroup() %>%
    dplyr::inner_join(dplyr::select(model_data$objects, ptmn_id, glm_object_ix)) %>%
    # FIXME cluster per object
    dplyr::inner_join(cluster_msprofiles(msdata.df, msdata$msrun_pepmodstate_stats, obj_col='pepmodstate_id', msrun_col='msrun')) %>%
    dplyr::arrange(glm_object_ix, profile_cluster, desc(n_quant), desc(intensity_med)) %>%
    dplyr::group_by(ptmn_id, glm_object_ix, profile_cluster) %>%
    dplyr::mutate(subobject_group_ix = row_number() %/% 20, # put objects within cluster into groups of 20
                  subobject_local_ix = row_number() %% 20) %>%
    dplyr::ungroup() %>%
    # take the first group of 10 objects from each cluster, then continue with the second group etc
    dplyr::arrange(glm_object_ix, subobject_group_ix, profile_cluster, subobject_local_ix) %>%
    dplyr::mutate(glm_subobject_ix = row_number()) %>%
    dplyr::filter(glm_subobject_ix <= 20) # remove less abundant subobjects of rich objects

# entries for an interaction in all replicate experiments
model_data$observations <- dplyr::inner_join(model_data$interactions, model_data$mschannels) %>%
  dplyr::arrange(glm_object_ix, glm_iaction_ix, mschannel_ix) %>%
  dplyr::mutate(glm_observation_ix = seq_len(n()))

model_data$msdata <- dplyr::inner_join(model_data$observations, model_data$subobjects) %>%
  dplyr::left_join(msdata.df) %>%
  dplyr::arrange(glm_observation_ix, glm_subobject_ix) %>%
  dplyr::mutate(glm_subobservation_ix = seq_len(n()))

model_data$msdata <- mutate(model_data$msdata,
                qdata_ix = if_else(!is.na(intensity), cumsum(!is.na(intensity)), NA_integer_),
                mdata_ix = if_else(is.na(intensity), cumsum(is.na(intensity)), NA_integer_))

model_data <- prepare_effects(model_data, underdefined_iactions=FALSE)

dims_info <- msglm.prepare_dims_info(model_data, object_cols=c('object_id', modelobj_idcol, "object_label", "ptmn_label_no_ptm_type",
                                                               "ptm_type", "ptm_AA_seq", "ptm_pos", "nselptms",
                                                               "protein_ac", "gene_name", "protein_name",
                                                               "is_viral", "is_contaminant"#, "is_decoy"
                                                               ))
# remove unneeded data to free some memory
msdata <- NULL
gc()

msglm.stan_data <- stan.prepare_data(mscalib, model_data,
                                     global_labu_shift = global_labu_shift,
                                     obj_labu_min_scale = 1,
                                     iact_repl_shift_df = 2,
                                     suo_fdr=0.001, reliable_obs_fdr = 0.01, specific_iaction_fdr = 1,
                                     empty_observation_sigmoid_scale = data_info$empty_observation_sigmoid_scale)

message('Running STAN in NUTS mode...')
options(mc.cores=mcmc_nchains)
msglm.stan_fit <- stan.sampling(msglm.stan_data, adapt_delta=0.9, max_treedepth=11L,
                                iter=4000L, chains=mcmc_nchains)

min.iteration <- as.integer(1.25 * msglm.stan_fit@sim$warmup)

msglm_results <- process.stan_fit(msglm.stan_fit, dims_info)

res_prefix <- paste0(project_id, "_", msfolder, "_msglm", modelobj_suffix)
if (!dir.exists(file.path(scratch_path, res_prefix))) {
  dir.create(file.path(scratch_path, res_prefix))
}
rfit_filepath <- file.path(scratch_path, res_prefix, paste0(res_prefix, '_', fit_version, '_', job_chunk, '.RData'))
message('Saving STAN results to ', rfit_filepath, '...')
results_info <- list(project_id = project_id, msfolder = msfolder,
                     data_version = data_version, fit_version = fit_version,
                     job_name = job_name, job_chunk = job_chunk, modelobj = modelobj, quantobj = quantobj)
save(data_info, results_info,
     model_data, msglm.stan_data, msglm_results,
     dims_info, file = rfit_filepath)
message('Done.')
on.exit(unlink(tempdir(), force = TRUE), add=TRUE)
