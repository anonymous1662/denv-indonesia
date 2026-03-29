# code/02_run_batch.R

# Script runs one batch 
# Needs to be called many times in parallel
# One batch is run for each municipality

library(here)
source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))
source(here("code", "equilibrate_sim.R"))
source(here("code", "scenario_sim.R"))

config <- yaml::read_yaml(here("code", "config.yaml"))
args <- commandArgs(trailingOnly = TRUE)
batch_id   <- as.numeric(args[[1]])
batch_size <- 1

output_path  <- here("output", config$run_name)
prep_path    <- here("output",config$run_name, "prep")
batches_path <- here(output_path, "batches")
dir.create(batches_path, showWarnings = FALSE, recursive = TRUE)

prepared <- qread(here(prep_path, "inputs_prepared.qs"))

config          <- prepared$config
all_gids        <- prepared$all_gids
indo_demog_full <- prepared$indo_demog_full
indo_foi_full   <- prepared$indo_foi_full
betas_full      <- prepared$betas_full
wol_regions_all <- prepared$wol_regions

n_municipalities <- length(all_gids)
n_batches <- ceiling(n_municipalities / batch_size)

if (batch_id > n_batches) {
  stop(sprintf("batch_id=%d exceeds n_batches=%d (n_municipalities=%d, batch_size=%d)",
               batch_id, n_batches, n_municipalities, batch_size))
}

start_idx <- (batch_id - 1) * batch_size + 1
end_idx   <- min(batch_id * batch_size, n_municipalities)
batch_indices <- start_idx:end_idx

batch_gids <- all_gids[batch_indices]
betas      <- betas_full[batch_indices]
wol_regions_batch <- wol_regions_all[batch_indices]

cat(sprintf("\n--- Batch %d/%d: municipalities %d-%d (n=%d) ---\n",
            batch_id, n_batches, start_idx, end_idx, length(batch_gids)))

# Load model
cat("Loading generative model\n")
denv_mod <- odin2::odin(here("code/denv-model.R"),  workdir = tempfile(pattern = paste0("dust_batch", batch_id, "_")))

# Subset data
indo_demog <- indo_demog_full |> filter(adm2_code %in% batch_gids)
indo_foi   <- indo_foi_full   |> filter(gid_2 %in% batch_gids)

# Wrangle demog
demog_input <- wrangle_demography(
  indo_demog,
  year_start = config$scenario_start,
  year_end   = (config$scenario_start + config$scenario_years),
  pad_left   = config$calibration_years
)
demog_input <- round(demog_input)
cat(sprintf("  demog_input dimensions: %s\n", paste(dim(demog_input), collapse = " x ")))

# Set seed
set.seed(123)

# Equilibrate (once per batch)
cat("  Equilibrating...\n")
denv_pars <- setup_denv_pars(demog_input, config, betas)

equilibration <- equilibrate_sim(
  demog_input, config, betas,
  output_path, denv_mod, denv_pars
)

baseline_states   <- equilibration$baseline_states
denv_out_reduced  <- equilibration$denv_out_reduced
denv_pars         <- equilibration$denv_pars

scenarios_to_run <- c("baseline", "wolbachia", "vaccine_national", "wolbachia_vaccine_national")
batch_reduced <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)

# Update parameters for scenario runs
set.seed(123)
seeding_times <- matrix(nrow = ncol(demog_input), ncol = config$scenario_years)
seeding_ages <- matrix(nrow = ncol(demog_input), ncol = config$scenario_years)

for(i in 1:ncol(demog_input)){
  for(j in 1:config$scenario_years){
    seeding_times[i,j] <- round(runif(1, 2, 365), digits = 0) # no seeding on day 1 of the year as this is when ageing occurs in the model
    seeding_ages[i,j] <- round(runif(1, 10, 45), digits = 0)
  }}

pop_seeding_weights <- apply(demog_input[, , config$calibration_years, drop = FALSE], 2, sum)
seeding_value <- pmax(1, floor(config$seeding_value * pop_seeding_weights/sum(pop_seeding_weights)))
denv_pars$vacc_on <- 1
denv_pars$wol_on <- 0
denv_pars$scenario_years <- config$scenario_years
denv_pars$demog  <- demog_input[, , (config$calibration_years + 1):(config$calibration_years + config$scenario_years), drop = FALSE]
denv_pars$N_init <- matrix(as.vector(demog_input[, , config$calibration_years + 1]), nrow = nrow(demog_input), ncol = ncol(demog_input))
denv_pars$seeding_times <- seeding_times
denv_pars$seeding_ages <- seeding_ages
denv_pars$seeding_value <- seeding_value
denv_pars$vacc_regions <- rep(1, ncol(demog_input)) # vector of length n_regions coded 1 for regions targeted with vaccination

# Initialise list to store results for each scenario
results <- list()
start_time <- Sys.time()
failure_log <- FailureLog$new(entries = list())

for (sc in scenarios_to_run) {
  cat(sprintf("  Running scenario: %s\n", sc))
  
  denv_pars_sc <- denv_pars
  
  scenario_results <- run_with_retry(
    scenario_sim,
    max_attempts = 3,
    failure_log  = failure_log,
    batch        = batch_id,
    log_scenario = sc,
    demog_input  = demog_input,
    config       = config,
    betas        = betas,
    output_path  = output_path,
    denv_mod     = denv_mod,
    baseline_states = baseline_states,
    denv_pars       = denv_pars_sc,
    scenario        = sc,
    wol_regions     = wol_regions_batch
  )
  
  
  batch_reduced[[sc]] <- scenario_results$scenario_out_reduced
}

# Save batch artefact
batch_obj <- list(
  batch_id        = batch_id,
  batch_size      = batch_size,
  n_batches       = n_batches,
  indices         = batch_indices,
  gids            = batch_gids,
  wol_regions     = wol_regions_batch,
  denv_out_reduced_equil = denv_out_reduced,
  reduced_by_scenario = batch_reduced
)

batch_file <- here(batches_path, sprintf("batch_%03d.qs", batch_id))
qsave(batch_obj, batch_file)

# Save batch failure log
failure_df <- failure_log$as_df()
rownames(failure_df) <- NULL
fail_file <- here(batches_path, sprintf("batch_%03d_failures.qs", batch_id))
qsave(failure_df, fail_file)

cat("Saved batch output: ", batch_file, "\n")
cat("Saved batch failures: ", fail_file, "\n")