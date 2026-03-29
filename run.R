# code/run.R
# Master script to test code locally. Runs batches sequentially

library(here)
library(qs)
here::here()

# Load packages
source("code/00_load-packages.R")

# Set and create output paths
config       <- yaml::read_yaml(here("code", "config.yaml"))
run_name     <- config$run_name
output_path  <- here("output", run_name)
prep_path    <- file.path(output_path, "prep")
batches_path <- file.path(output_path, "batches")
figure_path  <- here("figures", run_name)

dir.create(prep_path,    recursive = TRUE, showWarnings = FALSE)
dir.create(batches_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_path,  recursive = TRUE, showWarnings = FALSE)

# Prepare inputs (run once)
source("code/01_prepare_inputs.R")
stopifnot(file.exists(file.path(prep_path, "inputs_prepared.qs")))

# Determine number of batches
# batch_size=1 means one municipality per batch (set in 02_run_batch.R)
prep      <- qs::qread(file.path(prep_path, "inputs_prepared.qs"))
n_batches <- length(prep$all_gids)
cat("n_municipalities =", n_batches, "\n")

# Pre-compile model once so all workers use cached version
cat("Pre-compiling model...\n")
odin2::odin(here("code/denv-model.R"))
cat("Model compiled and cached\n")

# Build batch commands
# 02_run_batch.R takes a single positional argument: batch_id
rscript <- file.path(R.home("bin"), "Rscript")
cmds <- sprintf(
  "%s %s %d",
  shQuote(rscript),
  shQuote(here("code/02_run_batch.R")),
  seq_len(n_batches)
)

# Run batches sequentially
cat("Running", n_batches, "batches sequentially\n")
results <- lapply(cmds, system)
# To run in parallel instead, replace the two lines above with:
# n_cores <- max(parallel::detectCores() - 1, 1)
# results <- parallel::mclapply(cmds, system, mc.cores = n_cores)

# Check all batches completed
batch_files <- list.files(batches_path, pattern = "^batch_\\d{3}\\.qs$", full.names = TRUE)
cat("Expected:", n_batches, "  Found:", length(batch_files), "\n")
stopifnot(length(batch_files) == n_batches)

# Combine + post-process
cat("Combining batches...\n")
source("code/03_combine_batches.R")

# Check final outputs
final_files <- c(
  "equil_summary.qs",
  "averted_results.qs",
  "averted_national.qs",
  "cumulative_national.qs",
  "mean_age_national.qs",
  "mean_age_regional.qs",
  "prop_susc_national.qs",
  "prop_susc_regional.qs"
)
print(data.frame(
  file   = final_files,
  exists = file.exists(file.path(output_path, final_files))
))
