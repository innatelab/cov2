Sys.setenv(TZ='Etc/GMT+1') # issue#612@rstan
#job.args <- c("cov2", "ast_cov2_msglm", "cov2timecourse_dia_20200423", "20200506", "20200506", "0", "1235")
if (!exists('job.args')) {
  job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message('Project ID=', project_id)

job_name <- as.character(job.args[[2]])
job_version <- job.args[[5]]
ms_folder <- job.args[[3]]
data_version <- job.args[[4]]
fit_version <- job_version
job_id <- as.integer(job.args[[6]])
job_chunk <- as.integer(job.args[[7]])
message('Job ', job_name, '(id=', job_id, '_', job_chunk,
        " data_version=", data_version, " fit_version=", fit_version, " running on ", Sys.info()["nodename"], ")")

#source("~/R/config.R")
source("/projects/R/config.R")
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

rdata_filepath <- file.path(scratch_path, paste0(project_id, '_msglm_data_', ms_folder, '_', data_version, '.RData'))
message('Loading data from ', rdata_filepath)
load(rdata_filepath)

if (Sys.getenv('SLURM_CPUS_PER_TASK') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK')) # SLURM way
} else if (Sys.getenv('NSLOTS') != '') {
  mcmc_nchains <- as.integer(Sys.getenv('NSLOTS')) # SGE way
} else {
  mcmc_nchains <- 8
}

require(dplyr)
require(msglm)
require(rstan)
require(tidyr)
require(stringr)

#modelobj <- "protgroup"
#quantobj <- "protgroup"
modelobj <- "protregroup"
quantobj <- "pepmod"
source(file.path(project_scripts_path, 'setup_modelobj.R'))

sel_object_ids <- modelobjs_df[[modelobj_idcol]][[job_chunk]]
message(sel_object_ids, " ", modelobj, " ID(s): ",
        paste0(sort(unique(dplyr::pull(modelobjs_df[modelobjs_df[[modelobj_idcol]] %in% sel_object_ids, ], object_label))), collapse=' '))
if (modelobj == "protregroup") {
msdata.df <- dplyr::filter(msdata$protregroup2pepmod, protregroup_id %in% sel_object_ids & is_specific) %>%
  dplyr::inner_join(dplyr::select(msdata[[str_c(quantobj, "s")]], pepmod_id, any_of(quantobj_idcol))) %>%
  dplyr::inner_join(dplyr::select(msdata[[str_c(quantobj, "_intensities")]], !!quantobj_idcol, msrun, intensity)) %>%
  dplyr::mutate(object_id = protregroup_id)
} else if (modelobj == "protgroup") {
msdata.df <- dplyr::filter(msdata$protgroups, protgroup_id %in% sel_object_ids) %>%
  dplyr::inner_join(dplyr::select(msdata$protgroup_intensities, protgroup_id, msrun, intensity)) %>%
  dplyr::mutate(object_id = protgroup_id)
}
msdata.df <- msdata.df %>%
  dplyr::inner_join(dplyr::select(msdata$msruns, msrun, condition, treatment, timepoint))

message('Preparing MS GLM data...')
model_data <- list()
model_data$mschannels <- dplyr::select(msdata$msruns, condition, msrun, replicate) %>%
  dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::arrange(condition, replicate, msrun) %>% dplyr::distinct() %>%
  dplyr::mutate(mschannel_ix = row_number(),
                msrun_ix = as.integer(factor(msrun, levels=unique(msrun))),
                msproto_ix = 1L)
experiment_shift_col <- 'total_msrun_shift'
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

model_data$objects <- dplyr::select(model_data$interactions, glm_object_ix, object_id) %>%
  dplyr::distinct() %>% dplyr::arrange(glm_object_ix)

if (modelobj == "protregroup") {
model_data$objects <- dplyr::mutate(model_data$objects, protregroup_id = object_id) %>%
                      dplyr::inner_join(dplyr::select(modelobjs_df, protregroup_id, protregroup_label, object_label,
                                                      majority_protein_acs, protein_acs, gene_names, protein_names, contains("is_"))) %>%
  #dplyr::select(-is_fit) %>%
  dplyr::arrange(glm_object_ix)

# arrange pepmodstates by object and by the number of quantitations
model_data$subobjects <- msdata.df %>%
    dplyr::inner_join(msdata[[str_c(quantobj, "s")]]) %>%
    dplyr::inner_join(msdata[[str_c("protregroup2", quantobj)]])
if (has_name(model_data$subobjects, "is_specific")) {
    model_data$subobjects <- filter(model_data$subobjects, is_specific)
} else {
    model_data$subobjects <- mutate(model_data$subobjects, is_specific=TRUE)
}
if (!has_name(model_data$subobjects, "charge")) {
    model_data$subobjects <- mutate(model_data$subobjects, charge=NA_integer_)
}
model_data$subobjects <- model_data$subobjects %>%
    dplyr::group_by_at(intersect(c("protregroup_id", "pepmod_id", quantobj_idcol, "charge", "is_specific"),
                                 names(model_data$subobjects))) %>%
    dplyr::summarise(n_quant = sum(!is.na(intensity)),
                     intensity_med = median(intensity, na.rm=TRUE)) %>% ungroup() %>%
    dplyr::inner_join(dplyr::select(model_data$objects, protregroup_id, glm_object_ix)) %>%
    dplyr::arrange(glm_object_ix, desc(is_specific), desc(n_quant), desc(intensity_med),
                   pepmod_id, charge) %>%
    dplyr::mutate(glm_subobject_ix = row_number()) %>%
    dplyr::filter(glm_subobject_ix <= 50) # remove less abundant subobjects of rich objects
} else if (modelobj == "protgroup") {
model_data$objects <- dplyr::mutate(model_data$objects, protgroup_id = object_id) %>%
                      dplyr::inner_join(dplyr::select(modelobjs_df, protgroup_id, protgroup_label, object_label,
                                                      majority_protein_acs, protein_acs, gene_names, protein_names, contains("is_"))) %>%
  #dplyr::select(-is_fit) %>%
  dplyr::arrange(glm_object_ix)
}

# entries for an interaction in all replicate experiments
model_data$observations <- dplyr::inner_join(model_data$interactions, model_data$mschannels) %>%
  dplyr::arrange(glm_object_ix, glm_iaction_ix, mschannel_ix) %>%
  dplyr::mutate(glm_observation_ix = seq_len(n()))

if (quantobj != "protgroup") {
  model_data$msdata <- dplyr::inner_join(model_data$observations, model_data$subobjects) %>%
    dplyr::left_join(msdata.df) %>%
    dplyr::arrange(glm_observation_ix, glm_subobject_ix) %>%
    dplyr::mutate(glm_subobservation_ix = seq_len(n()))
} else if (quantobj == "protgroup") {
  model_data$msdata <- dplyr::left_join(model_data$observations, msdata.df) %>%
    dplyr::arrange(glm_observation_ix)
}

model_data$msdata <- mutate(model_data$msdata,
                qdata_ix = if_else(!is.na(intensity), cumsum(!is.na(intensity)), NA_integer_),
                mdata_ix = if_else(is.na(intensity), cumsum(is.na(intensity)), NA_integer_))

model_data <- prepare_effects(model_data, underdefined_iactions=FALSE)

dims_info <- msglm.prepare_dims_info(model_data, object_cols=c('object_id', modelobj_idcol, "object_label",
                                                               "majority_protein_acs", "gene_names",
                                                               "is_viral", "is_contaminant", "is_reverse",
                                                               "protein_names"
                                                               ))
# remove unneeded data to free some memory
msdata <- NULL
gc()

msglm.stan_data <- stan.prepare_data(instr_calib, model_data,
                                     global_labu_shift = global_labu_shift,
                                     batch_tau=0.8, effect_repl_shift_tau=0.4)

message('Running STAN in NUTS mode...')
options(mc.cores=mcmc_nchains)
msglm.stan_fit <- stan.sampling(msglm.stan_data, adapt_delta=0.85, max_treedepth=15L,
                                stepsize=1E-5, iter=6000L, chains=mcmc_nchains)
#require(shinystan)
#launch_shinystan(shinystan::as.shinystan(msglm.stan_fit))

min.iteration <- as.integer(1.25 * msglm.stan_fit@sim$warmup)

msglm_results <- process.stan_fit(msglm.stan_fit, dims_info)

res_prefix <- paste0(project_id, "_", ms_folder, "_msglm", modelobj_suffix)
if (!dir.exists(file.path(scratch_path, res_prefix))) {
  dir.create(file.path(scratch_path, res_prefix))
}
rdata_filepath <- file.path(scratch_path, res_prefix, paste0(res_prefix, '_', fit_version, '_', job_chunk, '.RData'))
message('Saving STAN results to ', rdata_filepath, '...')
results_info <- list(project_id = project_id, ms_folder = ms_folder,
                     data_version = data_version, fit_version = fit_version,
                     job_name = job_name, job_chunk = job_chunk, modelobj = modelobj, quantobj = quantobj)
save(data_info, results_info,
     model_data, msglm.stan_data, msglm_results,
     dims_info, file = rdata_filepath)
message('Done.')
on.exit(unlink(tempdir(), force = TRUE), add=TRUE)
