#job.args <- c("shamad_ivip", "shamad_ivip_msglm", "20190911", "20190911", "0", "1436")
#job.args <- c("shamad_ivip", "shamad_ivip_msglm", "20191022", "20191022", "0", "1436")
Sys.setenv(TZ='Etc/GMT+1') # issue#612@rstan
if (!exists('job.args')) {
  job.args <- commandArgs(trailingOnly = TRUE)
}

project_id <- job.args[[1]]
message('Project ID=', project_id)

job_name <- as.character(job.args[[2]])
batch_no <- as.character(job.args[[3]])
data_version <- job.args[[4]]
job_version <- job.args[[5]]
fit_version <- job_version
job_id <- as.integer(job.args[[6]])
iter <- as.integer(job.args[[7]])
proc_id <- as.integer(Sys.getenv('SLURM_PROCID'))
dist <- as.character(Sys.getenv('SLURM_DISTRIBUTION'))
job_chunk <- iter + proc_id
message('This is iter: ', iter, ' with the proc_id: ', proc_id, ' producing the chunk num: ', job_chunk, ' and the cluster has dist: ', dist )
message('The # tasks per node: ',as.character(Sys.getenv('SLURM_NTASKS_PER_NODE')), ' with RAM per node: ', as.character(Sys.getenv('SLURM_MEM_PER_NODE')), ' RAM per cpu: ',  as.character(Sys.getenv('SLURM_MEM_PER_CPU')), ' and cpus per task: ', as.character(Sys.getenv('SLURM_CPUS_PER_TASK')))
message('---------------------')
message('There are ', as.character(Sys.getenv('SLURM_NTASK')), ' ntasks in SLURM and ', as.character(Sys.getenv('SLURM_STEP_NUM_TASKS')),' task counts in SLURM')
message('---------------------')
message('Job ', job_name, '(id=', job_id, '_', job_chunk,
        " data_version=", data_version, " fit_version=", fit_version, " running on ",Sys.info()["nodename"],")")

source("/projects/R/config.R")
#source("~/R/config.R") #local
source(file.path(base_scripts_path, 'R/misc/setup_base_paths.R'))
source(file.path(base_scripts_path, 'R/misc/setup_project_paths.R'))

rdata_filepath <- file.path(scratch_path, paste0(project_id,'_B',batch_no, '_msglm_data_', data_version, '.RData'))
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

model_obj <- "protgroup"

model_obj_df <- paste0(model_obj, "s")
model_objid_col <- paste0(model_obj, "_id")
chunk_suffix <- case_when(model_obj == "protgroup" ~ "",
                          model_obj == "protregroup" ~ "_prg",
                          TRUE ~ NA_character_)
quant_obj <- case_when(model_obj == "protgroup" ~ "protgroup",
                       model_obj == "protregroup" ~ "pepmodstate",
                       TRUE ~ NA_character_)
global_labu_shift <- case_when(quant_obj == "protgroup" ~ global_protgroup_labu_shift,
                               #quant_obj == "pepmodstate" ~ global_pepmodstate_labu_shift,
                               TRUE ~ NA_real_)

sel_object_ids <- msdata[[model_obj_df]][[model_objid_col]][[job_chunk]]
message(sel_object_ids, " ", model_obj, " ID(s): ",
        paste0(sort(unique(dplyr::pull(msdata[[model_obj_df]][msdata[[model_obj_df]][[model_objid_col]] %in% sel_object_ids, ], gene_names))), collapse=' '))
msdata.df <- dplyr::filter(msdata$protgroups, protgroup_id %in% sel_object_ids) %>%
    dplyr::inner_join(dplyr::select(msdata$protgroup_intensities, protgroup_id, msrun, intensity, intensity_norm)) %>%
    dplyr::inner_join(msdata$msruns %>% dplyr::select(msrun, condition, virus, virus_full_id, timepoint) %>% dplyr::distinct()) %>%
    #dplyr::select(-msrun) %>%
    dplyr::mutate(object_id = protgroup_id)
#msdata.df <- dplyr::select(msdata.df, -one_of("intensity")) %>% dplyr::rename(intensity = intensity_corr.F)

message('Preparing MS GLM data...')
model_data <- list()
model_data$mschannels <- dplyr::select(dplyr::filter(msdata$mschannels, mstag != "Sum"),
                                       virus, virus_full_id, timepoint, replicate, condition,
                                       msrun) %>%
  dplyr::inner_join(dplyr::select(total_msrun_shifts.df, msrun, total_msrun_shift)) %>%
  dplyr::arrange(condition, replicate, msrun) %>% dplyr::distinct() %>%
  dplyr::mutate(mschannel_ix = row_number(),
                msrun_ix = row_number(),
                msproto_ix = 1L)
experiment_shift_col <- 'total_msrun_shift'
model_data$mschannels$model_mschannel_shift <- model_data$mschannels[[experiment_shift_col]]
model_data$conditions <- conditions.df

model_data$conditions <- dplyr::select(model_data$mschannels, condition, virus_full_id, virus, timepoint) %>%
  dplyr::distinct() %>%
  dplyr::arrange(condition) %>%
  dplyr::mutate(condition_ix = as.integer(factor(condition, levels = rownames(conditionXeffect.mtx)))) %>%
  dplyr::arrange(condition_ix)
model_data$mschannels <- dplyr::inner_join(model_data$mschannels, model_data$conditions) %>%
  dplyr::arrange(mschannel_ix)

model_data$interactions <- tidyr::crossing(object_id = sel_object_ids,
                                           condition_ix = model_data$conditions$condition_ix) %>%
  dplyr::inner_join(dplyr::select(model_data$conditions, condition_ix, condition) %>% dplyr::distinct()) %>%
  dplyr::left_join(msdata.df %>% dplyr::select(condition, object_id) %>% dplyr::distinct() %>%
                   dplyr::mutate(is_virtual = FALSE)) %>%
  dplyr::mutate(is_virtual = is.na(is_virtual),
                iaction_id = paste(condition, object_id)) %>%
  dplyr::arrange(condition_ix, object_id) %>%
  dplyr::mutate(glm_iaction_ix = row_number(),
                glm_object_ix = as.integer(factor(object_id)))

model_data$objects <- dplyr::select(model_data$interactions, glm_object_ix, object_id) %>%
    dplyr::distinct() %>% dplyr::arrange(glm_object_ix) %>%
    dplyr::mutate(protgroup_id = object_id) %>%
    dplyr::inner_join(dplyr::select(msdata$protgroups, protgroup_id,
                                    majority_protein_acs, protein_acs, gene_names, #protein_names, FIXME
                                    contains("is_"))) %>%
    #dplyr::select(-is_fit) %>%
    dplyr::arrange(glm_object_ix)

# entries for an interaction in all replicate experiments
model_data$observations <- dplyr::inner_join(model_data$interactions, model_data$mschannels) %>%
  dplyr::arrange(glm_object_ix, glm_iaction_ix, mschannel_ix) %>%
  dplyr::mutate(glm_observation_ix = seq_len(n()))

model_data$msdata <- dplyr::left_join(model_data$observations, msdata.df) %>%
    dplyr::arrange(glm_observation_ix) %>%
    mutate(qdata_ix = if_else(!is.na(intensity), cumsum(!is.na(intensity)), NA_integer_),
           mdata_ix = if_else(is.na(intensity), cumsum(is.na(intensity)), NA_integer_))

model_data <- prepare_effects(model_data, underdefined_iactions=FALSE)

msglm.stan_data <- stan.prepare_data(instr_calib, model_data,
                                     global_labu_shift = global_labu_shift,
                                     batch_tau=0.8, effect_repl_shift_tau=0.4)

message('Running STAN in NUTS mode...')
options(mc.cores=mcmc_nchains)
msglm.stan_fit <- stan.sampling(msglm.stan_data, adapt_delta=0.9, max_treedepth=12L, iter=4000L, chains=mcmc_nchains)

#require(shinystan)
#launch_shinystan(shinystan::as.shinystan(msglm.stan_fit))

min.iteration <- as.integer(1.5 * msglm.stan_fit@sim$warmup)
dims_info <- msglm.prepare_dims_info(model_data, object_cols=c('object_id', model_objid_col, "majority_protein_acs", "gene_names",
                                                               "is_viral", "is_contaminant", "is_reverse"#,
                                                               #"protein_names",
))

msglm_results <- process.stan_fit(msglm.stan_fit, dims_info)

res_prefix <- paste0(project_id, "_B",batch_no,"_msglm", chunk_suffix)
if (!dir.exists(file.path(scratch_path, res_prefix))) {
  dir.create(file.path(scratch_path, res_prefix))
}
rdata_filepath <- file.path(scratch_path, res_prefix, paste0(res_prefix, '_', fit_version, '_', job_chunk, '.RData'))
message('Saving STAN results to ', rdata_filepath, '...')
results_info <- list(project_id = project_id, data_version = data_version, fit_version = fit_version,
                     job_name = job_name, job_chunk = job_chunk, model_obj = model_obj, quant_obj = quant_obj)
save(data_info, results_info,
     model_data, msglm.stan_data, msglm_results,
     dims_info, file = rdata_filepath)
message('Done.')
on.exit(unlink(tempdir(), force = TRUE), add=TRUE)
