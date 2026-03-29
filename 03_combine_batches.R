# code/03_combine_batches.R

# Combine batch outputs + run post-processing (run once)

library(here)
source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))

config <- yaml::read_yaml(here("code", "config.yaml"))

output_path  <- here("output", config$run_name)
prep_path    <- here("output", config$run_name, "prep")
batches_path <- here(output_path, "batches")

prepared <- qread(here(prep_path, "inputs_prepared.qs"))
config   <- prepared$config
all_gids <- prepared$all_gids

scenarios_to_run <- c("baseline", "wolbachia", "vaccine_national", "wolbachia_vaccine_national")

# Load all batches
batch_files <- list.files(batches_path, pattern = "^batch_\\d{3}\\.qs$", full.names = TRUE)
if (length(batch_files) == 0) stop("No batch_XXX.qs files found in batches/")

batches <- lapply(batch_files, qread)

# Ensure ordered by batch_id
batch_ids <- vapply(batches, function(x) x$batch_id, double(1))
ord <- order(batch_ids)
batches <- batches[ord]
batch_ids <- batch_ids[ord]

cat(sprintf("Found %d batches: %s\n", length(batches), paste(batch_ids, collapse = ", ")))


# Assemble equilibration
national_equil_cases <- Reduce("+", lapply(batches, function(b) {
  apply(b$denv_out_reduced_equil[["cases_year"]], 3, sum, na.rm = TRUE)
}))
national_equil_inf <- Reduce("+", lapply(batches, function(b) {
  apply(b$denv_out_reduced_equil[["inf_year"]], 3, sum, na.rm = TRUE)
}))
national_equil_hosp <- Reduce("+", lapply(batches, function(b) {
  apply(b$denv_out_reduced_equil[["hosp_year"]], 3, sum, na.rm = TRUE)
}))

equil_summary <- bind_rows(
  make_equil_df(national_equil_cases, "cases_year"),
  make_equil_df(national_equil_inf,   "inf_year"),
  make_equil_df(national_equil_hosp,  "hosp_year")
)

# Combine batches — one scenario at a time, write to disk, then free memory
cat("\nCombining results across batches...\n")
for (sc in scenarios_to_run) {
  cat(sprintf("  Combining scenario: %s\n", sc))
  
  sc_batches <- lapply(batches, function(b) list(
    reduced = b$reduced_by_scenario[[sc]],
    gids    = b$gids,
    indices = b$indices
  ))
  
  combined_sc <- combine_batches(sc_batches)
  qsave(combined_sc, here(output_path, sprintf("combined_%s.qs", sc)))
  rm(combined_sc, sc_batches); gc()
}
rm(batches); gc()

# Helper: load a single scenario from disk on demand
get_scenario <- function(sc) qread(here(output_path, sprintf("combined_%s.qs", sc)))

# Post-processing: averted outcomes by municipality
intervention_scenarios <- c("wolbachia", "vaccine_national", "wolbachia_vaccine_national")
outcomes <- c("cases_cumulative", "hosp_cumulative",
              "vacc_cases_cumulative", "vacc_hosp_cumulative")

averted_results <- setNames(vector("list", length(outcomes)), outcomes)

for (outcome in outcomes) {
  cat(sprintf("\nCalculating averted %s...\n", outcome))
  averted_results[[outcome]] <- setNames(vector("list", length(intervention_scenarios)), intervention_scenarios)
  
  baseline_data <- get_scenario("baseline")
  
  for (sc in intervention_scenarios) {
    cat(sprintf("  Scenario: %s\n", sc))
    
    intervention_data <- get_scenario(sc)
    
    sc_results <- setNames(vector("list", length(all_gids)), all_gids)
    
    for (gid in all_gids) {
      baseline_arr     <- baseline_data[[gid]][[outcome]]
      intervention_arr <- intervention_data[[gid]][[outcome]]
      
      baseline_total     <- apply(baseline_arr, c(3, 4), sum, na.rm = TRUE)
      intervention_total <- apply(intervention_arr, c(3, 4), sum, na.rm = TRUE)
      
      averted    <- (baseline_total - intervention_total)
      baseline_t <- baseline_total
      n_times    <- ncol(averted)
      
      rows <- vector("list", n_times)
      for (t in seq_len(n_times)) {
        bl <- baseline_t[, t]
        av <- averted[, t]
        rows[[t]] <- data.frame(
          gid_2      = gid,
          year       = t,
          median     = median(av, na.rm = TRUE),
          q025       = quantile(av, 0.025, na.rm = TRUE),
          q975       = quantile(av, 0.975, na.rm = TRUE),
          pct_median = if (all(bl == 0)) NA_real_ else median(ifelse(bl == 0, NA_real_, av / bl * 100), na.rm = TRUE),
          pct_q025   = if (all(bl == 0)) NA_real_ else quantile(ifelse(bl == 0, NA_real_, av / bl * 100), 0.025, na.rm = TRUE),
          pct_q975   = if (all(bl == 0)) NA_real_ else quantile(ifelse(bl == 0, NA_real_, av / bl * 100), 0.975, na.rm = TRUE)
        )
      }
      sc_results[[gid]] <- do.call(rbind, rows)
    }
    
    averted_results[[outcome]][[sc]] <- do.call(rbind, sc_results)
    rownames(averted_results[[outcome]][[sc]]) <- NULL
    
    rm(intervention_data, sc_results); gc()
  }
  
  rm(baseline_data); gc()
}

##### Cumulative outcomes by municipality (all scenarios) #####
# Same structure as averted_results but reports absolute cumulative totals
# per municipality per scenario, not differences from baseline
outcomes_cumulative <- c("cases_cumulative", "hosp_cumulative")
cumulative_results <- setNames(vector("list", length(outcomes_cumulative)), outcomes_cumulative)

for (outcome in outcomes_cumulative) {
  cat(sprintf("\nCalculating municipal cumulative %s...\n", outcome))
  cumulative_results[[outcome]] <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  for (sc in scenarios_to_run) {
    cat(sprintf("  Scenario: %s\n", sc))
    
    sc_data    <- get_scenario(sc)
    sc_results <- setNames(vector("list", length(all_gids)), all_gids)
    
    for (gid in all_gids) {
      arr   <- sc_data[[gid]][[outcome]]
      total <- apply(arr, c(3, 4), sum, na.rm = TRUE)  # [particle, time]
      n_times <- ncol(total)
      
      rows <- vector("list", n_times)
      for (t in seq_len(n_times)) {
        v <- total[, t]
        rows[[t]] <- data.frame(
          gid_2  = gid,
          year   = t,
          median = median(v, na.rm = TRUE),
          q025   = quantile(v, 0.025, na.rm = TRUE),
          q975   = quantile(v, 0.975, na.rm = TRUE)
        )
      }
      sc_results[[gid]] <- do.call(rbind, rows)
    }
    
    cumulative_results[[outcome]][[sc]] <- do.call(rbind, sc_results)
    rownames(cumulative_results[[outcome]][[sc]]) <- NULL
    
    rm(sc_data, sc_results); gc()
  }
}


##### National averted outcomes #####
outcomes_all <- c(
  "cases_cumulative", "hosp_cumulative",
  "vacc_cases_cumulative", "vacc_hosp_cumulative",
  "vacc_seroneg_cases_cumulative", "vacc_seroneg_hosp_cumulative",
  "vacc_seropos_cases_cumulative", "vacc_seropos_hosp_cumulative",
  "inf_cumulative"
)
averted_national <- setNames(vector("list", length(outcomes_all)), outcomes_all)

for (outcome in outcomes_all) {
  cat(sprintf("\nCalculating national averted %s...\n", outcome))
  averted_national[[outcome]] <- setNames(vector("list", length(intervention_scenarios)), intervention_scenarios)
  
  baseline_data <- get_scenario("baseline")
  
  for (sc in intervention_scenarios) {
    cat(sprintf("  Scenario: %s\n", sc))
    
    intervention_data <- get_scenario(sc)
    
    # Build per-GID [particle, time] matrices
    baseline_list <- lapply(all_gids, function(gid) {
      apply(baseline_data[[gid]][[outcome]], c(3, 4), sum, na.rm = TRUE)
    })
    
    intervention_list <- lapply(all_gids, function(gid) {
      apply(intervention_data[[gid]][[outcome]], c(3, 4), sum, na.rm = TRUE)
    })
    
    # Stack into 3D arrays [particle, time, gid] and sum across GIDs
    baseline_array     <- simplify2array(baseline_list)
    intervention_array <- simplify2array(intervention_list)
    
    national_baseline     <- apply(baseline_array, c(1, 2), sum, na.rm = TRUE)
    national_intervention <- apply(intervention_array, c(1, 2), sum, na.rm = TRUE)
    
    rm(baseline_list, intervention_list, baseline_array, intervention_array,
       intervention_data); gc()
    
    # Compute averted nationally (sum first, then difference)
    national_averted <- national_baseline - national_intervention
    
    n_times <- ncol(national_averted)
    rows <- vector("list", n_times)
    for (t in seq_len(n_times)) {
      av <- national_averted[, t]
      bl <- national_baseline[, t]
      rows[[t]] <- data.frame(
        year       = t,
        median     = median(av, na.rm = TRUE),
        q025       = quantile(av, 0.025, na.rm = TRUE),
        q975       = quantile(av, 0.975, na.rm = TRUE),
        pct_median = if (all(bl == 0)) NA_real_ else median(ifelse(bl == 0, NA_real_, av / bl * 100), na.rm = TRUE),
        pct_q025   = if (all(bl == 0)) NA_real_ else quantile(ifelse(bl == 0, NA_real_, av / bl * 100), 0.025, na.rm = TRUE),
        pct_q975   = if (all(bl == 0)) NA_real_ else quantile(ifelse(bl == 0, NA_real_, av / bl * 100), 0.975, na.rm = TRUE)
      )
    }
    
    averted_national[[outcome]][[sc]] <- do.call(rbind, rows)
    rownames(averted_national[[outcome]][[sc]]) <- NULL
    
    rm(national_averted, national_baseline, national_intervention); gc()
  }
  
  rm(baseline_data); gc()
}

##### National cumulative totals by scenario #####
cumulative_national <- setNames(vector("list", length(outcomes_all)), outcomes_all)

for (outcome in outcomes_all) {
  cat(sprintf("\nCalculating national cumulative %s...\n", outcome))
  cumulative_national[[outcome]] <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  for (sc in scenarios_to_run) {
    sc_data <- get_scenario(sc)
    
    sc_list  <- lapply(all_gids, function(gid) {
      apply(sc_data[[gid]][[outcome]], c(3, 4), sum, na.rm = TRUE)
    })
    sc_total <- Reduce("+", sc_list)
    
    rm(sc_data, sc_list); gc()
    
    n_times <- ncol(sc_total)
    rows <- vector("list", n_times)
    for (t in seq_len(n_times)) {
      rows[[t]] <- data.frame(
        year   = t,
        median = median(sc_total[, t], na.rm = TRUE),
        q025   = quantile(sc_total[, t], 0.025, na.rm = TRUE),
        q975   = quantile(sc_total[, t], 0.975, na.rm = TRUE)
      )
    }
    
    cumulative_national[[outcome]][[sc]] <- do.call(rbind, rows)
    rownames(cumulative_national[[outcome]][[sc]]) <- NULL
    
    rm(sc_total); gc()
  }
}

##### Mean age of cases #####
cat("\nCalculating mean age of cases...\n")
age_vector <- 0:90

# Load all scenarios for mean age calculations
combined_reduced_all <- setNames(lapply(scenarios_to_run, get_scenario), scenarios_to_run)
mean_age_national <- calc_national_mean_age_of_cases(combined_reduced_all, all_gids, scenarios_to_run, age_vector)
mean_age_regional <- calc_regional_mean_age_of_cases(combined_reduced_all, all_gids, scenarios_to_run, age_vector)
rm(combined_reduced_all); gc()

##### Proportion susceptible #####
cat("\nCalculating proportion susceptible...\n")
demog_for_susc <- wrangle_demography(
  prepared$indo_demog_full,
  year_start = config$scenario_start,
  year_end   = config$scenario_start + config$scenario_years,
  pad_left   = 0
)

combined_reduced_all <- setNames(lapply(scenarios_to_run, get_scenario), scenarios_to_run)
prop_susc_national_raw <- calc_national_prop_susceptible(combined_reduced_all, all_gids, scenarios_to_run, demog_for_susc)
prop_susc_regional_raw <- calc_regional_prop_susceptible(combined_reduced_all, all_gids, scenarios_to_run, demog_for_susc)
rm(combined_reduced_all); gc()

# Summarise national proportion susceptible across particles
prop_susc_national <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
for (sc in scenarios_to_run) {
  mat <- prop_susc_national_raw[[sc]]  # [particle, time]
  n_times <- ncol(mat)
  rows <- vector("list", n_times)
  for (t in seq_len(n_times)) {
    rows[[t]] <- data.frame(
      year   = t,
      median = median(mat[, t], na.rm = TRUE),
      q025   = quantile(mat[, t], 0.025, na.rm = TRUE),
      q975   = quantile(mat[, t], 0.975, na.rm = TRUE)
    )
  }
  prop_susc_national[[sc]] <- do.call(rbind, rows)
  rownames(prop_susc_national[[sc]]) <- NULL
}
rm(prop_susc_national_raw); gc()

# Summarise regional proportion susceptible across particles
prop_susc_regional <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
for (sc in scenarios_to_run) {
  sc_rows <- vector("list", length(all_gids))
  for (gid_idx in seq_along(all_gids)) {
    gid <- all_gids[gid_idx]
    mat <- prop_susc_regional_raw[[sc]][[gid]]  # [particle, time]
    n_times <- ncol(mat)
    rows <- vector("list", n_times)
    for (t in seq_len(n_times)) {
      rows[[t]] <- data.frame(
        gid_2  = gid,
        year   = t,
        median = median(mat[, t], na.rm = TRUE),
        q025   = quantile(mat[, t], 0.025, na.rm = TRUE),
        q975   = quantile(mat[, t], 0.975, na.rm = TRUE)
      )
    }
    sc_rows[[gid_idx]] <- do.call(rbind, rows)
  }
  prop_susc_regional[[sc]] <- do.call(rbind, sc_rows)
  rownames(prop_susc_regional[[sc]]) <- NULL
}
rm(prop_susc_regional_raw); gc()

##### Save final results #####
qsave(equil_summary,      here(output_path, "equil_summary.qs"))
qsave(averted_results,    here(output_path, "averted_results.qs"))
qsave(averted_national,   here(output_path, "averted_national.qs"))
qsave(cumulative_national,here(output_path, "cumulative_national.qs"))
qsave(cumulative_results, here(output_path, "cumulative_results.qs"))
qsave(mean_age_national,  here(output_path, "mean_age_national.qs"))
qsave(mean_age_regional,  here(output_path, "mean_age_regional.qs"))
qsave(prop_susc_national, here(output_path, "prop_susc_national.qs"))
qsave(prop_susc_regional, here(output_path, "prop_susc_regional.qs"))