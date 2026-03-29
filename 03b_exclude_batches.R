# code/03b_exclude_batches.R

# Re-derive all national summary outputs excluding problematic municipalities

library(here)
source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))

config <- yaml::read_yaml(here("code", "config.yaml"))

output_path <- here("output", config$run_name)
prep_path   <- here("output", config$run_name, "prep")

prepared <- qread(here(prep_path, "inputs_prepared.qs"))
config   <- prepared$config
all_gids <- prepared$all_gids

scenarios_to_run       <- c("baseline", "wolbachia", "vaccine_national", "wolbachia_vaccine_national")
intervention_scenarios <- c("wolbachia", "vaccine_national", "wolbachia_vaccine_national")

outcomes <- c("cases_cumulative", "hosp_cumulative",
              "vacc_cases_cumulative", "vacc_hosp_cumulative")

outcomes_cumulative <- c("cases_cumulative", "hosp_cumulative")

outcomes_all <- c(
  "cases_cumulative", "hosp_cumulative",
  "vacc_cases_cumulative", "vacc_hosp_cumulative",
  "vacc_seroneg_cases_cumulative", "vacc_seroneg_hosp_cumulative",
  "vacc_seropos_cases_cumulative", "vacc_seropos_hosp_cumulative",
  "inf_cumulative"
)

OUTLIER_THRESHOLD <- -1000

get_scenario <- function(sc) qread(here(output_path, sprintf("combined_%s.qs", sc)))

# Helper: for a given scenario, return the GIDs to use (all minus outliers)
get_gids <- function(sc) all_gids[!all_gids %in% outlier_gids[[sc]]]


###### Identify outlier municipalities per intervention scenario #####
cat("=== Identifying outlier municipalities ===\n\n")

averted_results <- qread(here(output_path, "averted_results.qs"))

outlier_gids <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
outlier_gids[["baseline"]] <- character(0)

for (sc in intervention_scenarios) {
  bad <- averted_results[["cases_cumulative"]][[sc]] |>
    dplyr::group_by(gid_2) |>
    dplyr::summarise(worst = min(median), .groups = "drop") |>
    dplyr::filter(worst < OUTLIER_THRESHOLD) |>
    dplyr::pull(gid_2)
  
  outlier_gids[[sc]] <- bad
  
  if (length(bad) == 0) {
    cat(sprintf("  %-40s  no outliers\n", sc))
  } else {
    cat(sprintf("  %-40s  excluding: %s\n", sc, paste(bad, collapse = ", ")))
  }
}

rm(averted_results); gc()
qsave(outlier_gids, here(output_path, "outlier_gids.qs"))
cat("\nSaved: outlier_gids.qs\n\n")


##### Averted outcomes by municipality #####
cat("=== Calculating averted outcomes by municipality (excl) ===\n\n")

averted_results_excl <- setNames(vector("list", length(outcomes)), outcomes)

for (outcome in outcomes) {
  cat(sprintf("\nCalculating averted %s...\n", outcome))
  averted_results_excl[[outcome]] <- setNames(vector("list", length(intervention_scenarios)), intervention_scenarios)
  
  baseline_data <- get_scenario("baseline")
  
  for (sc in intervention_scenarios) {
    cat(sprintf("  Scenario: %s\n", sc))
    gids_to_use <- get_gids(sc)
    
    intervention_data <- get_scenario(sc)
    sc_results        <- setNames(vector("list", length(gids_to_use)), gids_to_use)
    
    for (gid in gids_to_use) {
      baseline_total     <- apply(baseline_data[[gid]][[outcome]],     c(3, 4), sum, na.rm = TRUE)
      intervention_total <- apply(intervention_data[[gid]][[outcome]], c(3, 4), sum, na.rm = TRUE)
      
      averted    <- baseline_total - intervention_total
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
    
    averted_results_excl[[outcome]][[sc]] <- do.call(rbind, sc_results)
    rownames(averted_results_excl[[outcome]][[sc]]) <- NULL
    
    rm(intervention_data, sc_results); gc()
  }
  
  rm(baseline_data); gc()
}

qsave(averted_results_excl, here(output_path, "averted_results_excl.qs"))
cat("\nSaved: averted_results_excl.qs\n\n")

##### Cumulative outcomes by municipality #####
cat("=== Calculating cumulative outcomes by municipality (excl) ===\n\n")

cumulative_results_excl <- setNames(vector("list", length(outcomes_cumulative)), outcomes_cumulative)

for (outcome in outcomes_cumulative) {
  cat(sprintf("\nCalculating municipal cumulative %s...\n", outcome))
  cumulative_results_excl[[outcome]] <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  for (sc in scenarios_to_run) {
    cat(sprintf("  Scenario: %s\n", sc))
    gids_to_use <- get_gids(sc)
    
    sc_data    <- get_scenario(sc)
    sc_results <- setNames(vector("list", length(gids_to_use)), gids_to_use)
    
    for (gid in gids_to_use) {
      total   <- apply(sc_data[[gid]][[outcome]], c(3, 4), sum, na.rm = TRUE)
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
    
    cumulative_results_excl[[outcome]][[sc]] <- do.call(rbind, sc_results)
    rownames(cumulative_results_excl[[outcome]][[sc]]) <- NULL
    
    rm(sc_data, sc_results); gc()
  }
}

qsave(cumulative_results_excl, here(output_path, "cumulative_results_excl.qs"))
cat("\nSaved: cumulative_results_excl.qs\n\n")

##### National averted outcomes #####
cat("=== Calculating national averted outcomes (excl) ===\n\n")

averted_national_excl <- setNames(vector("list", length(outcomes_all)), outcomes_all)

for (outcome in outcomes_all) {
  cat(sprintf("\nCalculating national averted %s...\n", outcome))
  averted_national_excl[[outcome]] <- setNames(vector("list", length(intervention_scenarios)), intervention_scenarios)
  
  baseline_data <- get_scenario("baseline")
  
  for (sc in intervention_scenarios) {
    cat(sprintf("  Scenario: %s\n", sc))
    gids_to_use <- get_gids(sc)
    
    intervention_data <- get_scenario(sc)
    
    baseline_list <- lapply(gids_to_use, function(gid) {
      apply(baseline_data[[gid]][[outcome]], c(3, 4), sum, na.rm = TRUE)
    })
    intervention_list <- lapply(gids_to_use, function(gid) {
      apply(intervention_data[[gid]][[outcome]], c(3, 4), sum, na.rm = TRUE)
    })
    
    baseline_array     <- simplify2array(baseline_list)
    intervention_array <- simplify2array(intervention_list)
    
    national_baseline     <- apply(baseline_array, c(1, 2), sum, na.rm = TRUE)
    national_intervention <- apply(intervention_array, c(1, 2), sum, na.rm = TRUE)
    
    rm(baseline_list, intervention_list, baseline_array, intervention_array,
       intervention_data); gc()
    
    national_averted <- national_baseline - national_intervention
    n_times          <- ncol(national_averted)
    
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
    
    averted_national_excl[[outcome]][[sc]] <- do.call(rbind, rows)
    rownames(averted_national_excl[[outcome]][[sc]]) <- NULL
    
    rm(national_averted, national_baseline, national_intervention); gc()
  }
  
  rm(baseline_data); gc()
}

qsave(averted_national_excl, here(output_path, "averted_national_excl.qs"))
cat("\nSaved: averted_national_excl.qs\n\n")

##### National cumulative totals #####
cat("=== Calculating national cumulative totals (excl) ===\n\n")

cumulative_national_excl <- setNames(vector("list", length(outcomes_all)), outcomes_all)

for (outcome in outcomes_all) {
  cat(sprintf("\nCalculating national cumulative %s...\n", outcome))
  cumulative_national_excl[[outcome]] <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  for (sc in scenarios_to_run) {
    gids_to_use <- get_gids(sc)
    
    sc_data  <- get_scenario(sc)
    sc_list  <- lapply(gids_to_use, function(gid) {
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
    
    cumulative_national_excl[[outcome]][[sc]] <- do.call(rbind, rows)
    rownames(cumulative_national_excl[[outcome]][[sc]]) <- NULL
    
    rm(sc_total); gc()
  }
}

qsave(cumulative_national_excl, here(output_path, "cumulative_national_excl.qs"))
cat("Saved: cumulative_national_excl.qs\n\n")

##### Mean age of cases #####
cat("=== Calculating mean age of cases (excl) ===\n\n")

age_vector <- 0:90

mean_age_national_excl <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
mean_age_regional_excl <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)

for (sc in scenarios_to_run) {
  gids_to_use <- get_gids(sc)
  cat(sprintf("  %s: %d GIDs\n", sc, length(gids_to_use)))
  
  sc_data     <- get_scenario(sc)
  arr         <- sc_data[[gids_to_use[1]]][["cases_year"]]
  n_times     <- dim(arr)[4]
  n_particles <- dim(arr)[3]
  
  national_age_weighted <- matrix(0, nrow = n_particles, ncol = n_times)
  national_cases_total  <- matrix(0, nrow = n_particles, ncol = n_times)
  sc_results            <- setNames(vector("list", length(gids_to_use)), gids_to_use)
  
  for (gid in gids_to_use) {
    arr  <- sc_data[[gid]][["cases_year"]]
    arr3 <- array(arr[, 1, , ], dim = c(dim(arr)[1], n_particles, n_times))
    
    # Regional
    rows <- vector("list", n_times)
    for (t in seq_len(n_times)) {
      cases_by_age <- arr3[, , t]
      total_cases  <- colSums(cases_by_age)
      mean_ages    <- colSums(age_vector * cases_by_age) / total_cases
      mean_ages[is.nan(mean_ages)] <- NA
      rows[[t]] <- data.frame(
        gid_2  = gid,
        year   = t,
        median = median(mean_ages, na.rm = TRUE),
        q025   = quantile(mean_ages, 0.025, na.rm = TRUE),
        q975   = quantile(mean_ages, 0.975, na.rm = TRUE)
      )
      # Accumulate national
      national_cases_total[, t]  <- national_cases_total[, t]  + total_cases
      national_age_weighted[, t] <- national_age_weighted[, t] + colSums(age_vector * cases_by_age)
    }
    sc_results[[gid]] <- do.call(rbind, rows)
  }
  rm(sc_data); gc()
  
  # National summary
  rows <- vector("list", n_times)
  for (t in seq_len(n_times)) {
    mean_ages <- national_age_weighted[, t] / national_cases_total[, t]
    mean_ages[is.nan(mean_ages)] <- NA
    rows[[t]] <- data.frame(
      year   = t,
      median = median(mean_ages, na.rm = TRUE),
      q025   = quantile(mean_ages, 0.025, na.rm = TRUE),
      q975   = quantile(mean_ages, 0.975, na.rm = TRUE)
    )
  }
  mean_age_national_excl[[sc]] <- do.call(rbind, rows)
  rownames(mean_age_national_excl[[sc]]) <- NULL
  
  # Regional summary
  mean_age_regional_excl[[sc]] <- do.call(rbind, sc_results)
  rownames(mean_age_regional_excl[[sc]]) <- NULL
}

qsave(mean_age_national_excl, here(output_path, "mean_age_national_excl.qs"))
qsave(mean_age_regional_excl, here(output_path, "mean_age_regional_excl.qs"))
cat("Saved: mean_age_national_excl.qs\n")
cat("Saved: mean_age_regional_excl.qs\n\n")

##### Proportion susceptible #####
# Note - this output failed to save but I did not have enough time to look into why
cat("=== Calculating proportion susceptible (excl) ===\n\n")

demog_for_susc <- wrangle_demography(
  prepared$indo_demog_full,
  year_start = config$scenario_start,
  year_end   = config$scenario_start + config$scenario_years,
  pad_left   = 0
)

prop_susc_national_excl <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
prop_susc_regional_excl <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)

for (sc in scenarios_to_run) {
  gids_to_use <- get_gids(sc)
  idx_to_use  <- which(all_gids %in% gids_to_use)
  cat(sprintf("  %s: %d GIDs\n", sc, length(gids_to_use)))
  
  sc_data     <- get_scenario(sc)
  arr1        <- sc_data[[gids_to_use[1]]][["susceptibles"]]
  n_times     <- dim(arr1)[4]
  n_particles <- dim(arr1)[3]
  
  national_susc <- matrix(0, nrow = n_particles, ncol = n_times)
  national_pop  <- numeric(n_times)
  sc_rows       <- vector("list", length(all_gids))
  
  for (i in seq_along(gids_to_use)) {
    gid     <- gids_to_use[i]
    gid_idx <- idx_to_use[i]
    susc    <- sc_data[[gid]][["susceptibles"]]
    susc3   <- array(susc[, 1, , ], dim = c(dim(susc)[1], n_particles, n_times))
    
    pop_by_time <- apply(demog_for_susc[, gid_idx, , drop = FALSE], 3, sum)
    prop        <- matrix(0, nrow = n_particles, ncol = n_times)
    
    for (t in seq_len(n_times)) {
      national_susc[, t] <- national_susc[, t] + colSums(susc3[, , t])
      national_pop[t]    <- national_pop[t]    + pop_by_time[t]
      prop[, t]          <- colSums(susc3[, , t]) / pop_by_time[t]
    }
    
    rows <- vector("list", n_times)
    for (t in seq_len(n_times)) {
      rows[[t]] <- data.frame(
        gid_2  = gid,
        year   = t,
        median = median(prop[, t], na.rm = TRUE),
        q025   = quantile(prop[, t], 0.025, na.rm = TRUE),
        q975   = quantile(prop[, t], 0.975, na.rm = TRUE)
      )
    }
    sc_rows[[i]] <- do.call(rbind, rows)
  }
  rm(sc_data); gc()
  
  # National summary
  prop_mat <- national_susc / matrix(national_pop, nrow = n_particles, ncol = n_times, byrow = TRUE)
  rows <- vector("list", n_times)
  for (t in seq_len(n_times)) {
    rows[[t]] <- data.frame(
      year   = t,
      median = median(prop_mat[, t], na.rm = TRUE),
      q025   = quantile(prop_mat[, t], 0.025, na.rm = TRUE),
      q975   = quantile(prop_mat[, t], 0.975, na.rm = TRUE)
    )
  }
  prop_susc_national_excl[[sc]] <- do.call(rbind, rows)
  rownames(prop_susc_national_excl[[sc]]) <- NULL
  
  # Regional summary
  prop_susc_regional_excl[[sc]] <- do.call(rbind, sc_rows[seq_along(gids_to_use)])
  rownames(prop_susc_regional_excl[[sc]]) <- NULL
}

qsave(prop_susc_national_excl, here(output_path, "prop_susc_national_excl.qs"))
qsave(prop_susc_regional_excl, here(output_path, "prop_susc_regional_excl.qs"))
cat("Saved: prop_susc_national_excl.qs\n")
cat("Saved: prop_susc_regional_excl.qs\n\n")
