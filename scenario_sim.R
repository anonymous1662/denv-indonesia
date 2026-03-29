# code/scenario_sim.R

# Defines scenario_sim() function used in 02_run_batch.R

scenario_sim <- function(demog_input, config, betas, output_path, denv_mod, baseline_states, denv_pars, scenario, wol_regions) {
  
  n_cores <- max(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1)) - 4,1)
  
  if(!grepl("vaccine", scenario)){ # if the scenario name does not contain "vaccine", set vaccine effects to 0
    denv_pars$ve_inf <- c(0.0, 0.0)
    denv_pars$ve_symp <- c(0.0, 0.0)
    denv_pars$ve_hosp <- c(0.0, 0.0)
    denv_pars$ve_death <- c(0.0, 0.0)
  }
  
  if(!grepl("wolbachia", scenario)){
    denv_pars$wol_on <- 1
    denv_pars$wol_reduction <- 0
    denv_pars$wol_regions <- wol_regions
  }
  
  if (grepl("wolbachia", scenario)) { # if the scenario name contains "wolbachia"
    denv_pars$wol_on <- 1
    denv_pars$wol_regions <- wol_regions # vector of length(n_regions coded 1 for regions targeted with wMel)
    denv_pars$wol_reduction <- config$wol_reduction
  }
  
  scenario_sys <- dust_system_create(denv_mod,
                                     pars = denv_pars,
                                     n_particles = config$n_particles,
                                     n_threads = n_cores,
                                     time = 1,
                                     seed = 42)
  
  # Set state from equilibration end-point
  dust_system_set_state(scenario_sys, baseline_states)
  
  # Zero all output trackers so they reflect scenario-period only.
  # Compartments (S, I, C, R, V, N) are left intact so epidemic state
  # carries over correctly from equilibration.
  idx <- dust_unpack_index(scenario_sys)
  
  fields_to_zero <- c(
    # Cumulative trackers
    "inf_cumulative",
    "cases_cumulative", "hosp_cumulative", "deaths_cumulative",
    "vacc_cumulative", "vacc_cumulative_seronegative", "vacc_cumulative_seropositive",
    "vacc_cases_cumulative", "vacc_hosp_cumulative", "vacc_deaths_cumulative",
    "vacc_seroneg_cases_cumulative", "vacc_seroneg_hosp_cumulative", "vacc_seroneg_deaths_cumulative",
    "vacc_seropos_cases_cumulative", "vacc_seropos_hosp_cumulative", "vacc_seropos_deaths_cumulative",
    # Within-year trackers
    "inf_year",
    "cases_year", "hosp_year", "deaths_year",
    "vacc_cases_year", "vacc_hosp_year", "vacc_deaths_year",
    "vacc_seroneg_cases_year", "vacc_seroneg_hosp_year", "vacc_seroneg_deaths_year",
    "vacc_seropos_cases_year", "vacc_seropos_hosp_year", "vacc_seropos_deaths_year",
    # Daily incidence trackers
    "cases_i", "hosp_i", "deaths_i",
    "vacc_i",
    # Annual FOI accumulator
    "annual_foi_region"
  )
  
  state <- dust_system_state(scenario_sys)
  for (field in fields_to_zero) {
    if (!is.null(idx[[field]])) {
      state[idx[[field]], ] <- 0
    }
  }
  dust_system_set_state(scenario_sys, state)
  cat(sc, "RNG hash:", digest::digest(dust_system_rng_state(scenario_sys)), "\n")
  
  # Output at yearly intervals
  t <- c(seq(365, 365*config$scenario_years, by = 365)) # define output time points
  sim_raw <- dust_system_simulate(scenario_sys, t)
  print(Sys.time() - start_time)
  
  
  clean_values <- function(x, name = "", threshold = 1e15) {
    
    count_bad <- function(col) sum(!is.finite(col) | abs(col) > threshold, na.rm = TRUE)
    replace_bad <- function(col) replace(col, !is.finite(col) | abs(col) > threshold, NA)
    
    if (is.data.frame(x)) {
      numeric_cols <- sapply(x, is.numeric)
      n_bad <- sum(sapply(x[numeric_cols], count_bad))
      x[numeric_cols] <- lapply(x[numeric_cols], replace_bad)
    } else if (is.numeric(x)) {
      n_bad <- count_bad(x)
      x <- replace_bad(x)
    } else {
      n_bad <- 0
    }
    
    if (n_bad > 0) {
      message(sprintf("  [%s]: %d value(s) replaced with NA", name, n_bad))
    }
    
    return(x)
  }
  
  scenario_out <- dust_unpack_state(scenario_sys, sim_raw)
  scenario_out_reduced <- list(
    # Total population
    cases_year = scenario_out$cases_year,
    hosp_year = scenario_out$hosp_year,
    cases_cumulative = scenario_out$cases_cumulative,
    hosp_cumulative = scenario_out$hosp_cumulative,
    
    # Infections
    inf_year = scenario_out$inf_year,
    inf_cumulative = scenario_out$inf_cumulative,
    
    # Vaccinated individuals
    vacc_cases_year = scenario_out$vacc_cases_year,
    vacc_hosp_year = scenario_out$vacc_hosp_year,
    vacc_cases_cumulative = scenario_out$vacc_cases_cumulative,
    vacc_hosp_cumulative = scenario_out$vacc_hosp_cumulative,

    # Seronegative recipients
    vacc_seroneg_cases_year = scenario_out$vacc_seroneg_cases_year,
    vacc_seroneg_hosp_year = scenario_out$vacc_seroneg_hosp_year,
    vacc_seroneg_cases_cumulative = scenario_out$vacc_seroneg_cases_cumulative,
    vacc_seroneg_hosp_cumulative = scenario_out$vacc_seroneg_hosp_cumulative,
    
    # Seropositive recipients
    vacc_seropos_cases_year = scenario_out$vacc_seropos_cases_year,
    vacc_seropos_hosp_year = scenario_out$vacc_seropos_hosp_year,
    vacc_seropos_cases_cumulative = scenario_out$vacc_seropos_cases_cumulative,
    vacc_seropos_hosp_cumulative = scenario_out$vacc_seropos_hosp_cumulative,
    
    # Number of susceptibles and age of cases
    susceptibles = scenario_out$S
  )
  message("Cleaning values in scenario_out_reduced:")
  scenario_out_reduced <- mapply(
    clean_values,
    x    = scenario_out_reduced,
    name = names(scenario_out_reduced),
    SIMPLIFY = FALSE
  )
# Store results for this scenario (no separate infections object)
  results <- list(scenario_out_reduced = scenario_out_reduced)
  #cat(paste0("Finished simulations for: ", scenario, "\n"))
  
  # Return all scenario results
  return(results)
}
  