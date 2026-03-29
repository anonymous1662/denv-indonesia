# code/utils.R
# Defines a set of functions

# Function splits population in each adm2 into year groups from age groups
expand_age_groups <- function(demog, max_age){
  
  demog_single <- demog |>
    dplyr::group_by(adm2_code, year) |>
    dplyr::arrange(age_lower, .by_group = TRUE) |>
    dplyr::mutate(
      age_upper = dplyr::lead(age_lower, default = max_age + 1) - 1,
      age_upper = pmin(age_upper, max_age),
      age_width = age_upper - age_lower + 1,
      population_single = population / age_width
    ) |>
    dplyr::ungroup()
  
  expanded_demog <- demog_single |>
    dplyr::rowwise() |>
    dplyr::reframe(
      adm2_code = adm2_code,
      year = year,
      age = seq(age_lower, age_upper),
      population = population_single
    )
  
  expanded_demog
}


wrangle_demography <- function(demog, year_start, year_end, pad_left = 0){
  demog <- demog |> 
    filter(year >= year_start & year <= year_end) |> 
    arrange(adm2_code)
  
  demog <- xtabs(population ~ age + adm2_code + year, data = demog)
  
  # Ensure demog is always 3D (age x adm2 x year) even with a single adm2
  if(length(dim(demog)) == 2){
    demog <- array(demog, dim = c(dim(demog)[1], 1, dim(demog)[2]),
                   dimnames = list(dimnames(demog)[[1]], 
                                   unique(demog_orig$adm2_code),  
                                   dimnames(demog)[[2]]))
  }
  
  if(pad_left != 0){ 
    demog <- abind(array(demog[,,1], dim = c(dim(demog)[1:2], pad_left)), demog, along = 3)
  }
  
  return(demog)
}

setup_denv_pars <- function(demog_input, config, betas) {
  
  # Define seeding parameters
  set.seed(123)
  seeding_times <- matrix(nrow = ncol(demog_input), ncol = config$calibration_years)
  seeding_ages <- matrix(nrow = ncol(demog_input), ncol = config$calibration_years)
  
  for(i in 1:ncol(demog_input)){
    for(j in 1:config$calibration_years){
      seeding_times[i,j] <- round(runif(1, 2, 365), digits = 0)
      seeding_ages[i,j] <- round(runif(1, 10, 45), digits = 0)
    }
  }
  
  pop_seeding_weights <- apply(demog_input[, , 1, drop = FALSE], 2, sum)
  seeding_value <- pmax(1, floor(config$seeding_value * pop_seeding_weights/sum(pop_seeding_weights)))
  
  # Create parameter list
  denv_pars <- list(
    dt = as.double(1),
    scenario_years = config$calibration_years,
    N_init = matrix(as.vector(demog_input[, , 1]), nrow = nrow(demog_input), ncol = ncol(demog_input)),
    demog  = demog_input[, , 1:config$calibration_years, drop = FALSE],
    beta = betas,
    gamma = config$gamma,
    amp_seas = config$amp_seas,
    phase_seas = config$phase_seas,
    seeding_times = seeding_times,
    seeding_ages = seeding_ages,
    seeding_value = seeding_value,
    n_age = nrow(demog_input),
    n_regions = ncol(demog_input),
    n_serotypes = 4,
    symp_prop = config$symp_prop,
    hosp_prop = config$hosp_prop,
    death_prop = rep(0, nrow(demog_input)),
    vacc_on = 0,
    vacc_start = config$vacc_start,
    catch_up_end = as.numeric(config$vacc_start) + as.numeric(config$catch_up_duration),
    vacc_age = config$vacc_age,
    vacc_coverage = config$vacc_coverage,
    catch_up_age = numeric(nrow(demog_input)),
    vacc_regions = numeric(ncol(demog_input)),
    ve_inf = as.double(config$ve_inf),
    ve_symp = config$ve_symp,
    ve_hosp = config$ve_hosp,
    ve_death = config$ve_death,
    ve_wane = config$ve_wane,
    wol_on = 0,
    wol_start = config$wol_start,
    wol_reduction = config$wol_reduction,
    wol_regions = numeric(ncol(demog_input))
  )
  
  return(denv_pars)
}



##### Run_with_retry funtion #####
run_with_retry <- function(fn, max_attempts = 3, failure_log = NULL, batch = NULL, log_scenario = NULL, ...) {
  attempt <- 1
  last_error <- NULL
  
  while (attempt <= max_attempts) {
    tryCatch({
      result <- fn(...)
      
      if (attempt > 1 & !is.null(failure_log)) {
        failure_log$add(data.frame(
          batch     = if (!is.null(batch)) batch else NA,
          scenario  = if (!is.null(log_scenario)) log_scenario else NA,
          attempt   = attempt,
          outcome   = "retry_success",
          error_msg = NA,
          timestamp = as.character(Sys.time())
        ))
      }
      return(result)
      
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      last_error <<- e
      
      is_binomial_error <- grepl("binomial|binom|particle.*error", error_msg, ignore.case = TRUE) ||
        grepl("Invalid call to binomial", error_msg, fixed = TRUE)
      
      if (is_binomial_error) {
        
        if (!is.null(failure_log)) {
          failure_log$add(data.frame(
            batch     = if (!is.null(batch)) batch else NA,
            scenario  = if (!is.null(log_scenario)) log_scenario else NA,
            attempt   = attempt,
            outcome   = ifelse(attempt < max_attempts, "binomial_retry", "binomial_failed"),
            error_msg = substr(error_msg, 1, 250),
            timestamp = as.character(Sys.time())
          ))
        }
        
        if (attempt < max_attempts) {
          attempt <<- attempt + 1
          Sys.sleep(0.1)
        } else {
          stop(e)
        }
        
      } else {
        
        if (!is.null(failure_log)) {
          failure_log$add(data.frame(
            batch     = if (!is.null(batch)) batch else NA,
            scenario  = if (!is.null(log_scenario)) log_scenario else NA,
            attempt   = attempt,
            outcome   = "non_binomial_failed",
            error_msg = substr(error_msg, 1, 250),
            timestamp = as.character(Sys.time())
          ))
        }
        stop(e)
      }
    })
  }
  
  if (!is.null(last_error)) stop(last_error)
}


##### Combine batches function #####
combine_batches <- function(batch_list) {
  
  all_gids_ordered <- unlist(lapply(batch_list, function(b) b$gids))
  field_names <- names(batch_list[[1]]$reduced)  # $reduced is the named list of arrays
  
  combined <- setNames(vector("list", length(all_gids_ordered)), all_gids_ordered)
  
  for (batch in batch_list) {
    n_batch_regions <- length(batch$gids)
    for (i in seq_along(batch$gids)) {
      gid <- batch$gids[i]
      combined[[gid]] <- setNames(vector("list", length(field_names)), field_names)
      
      for (field in field_names) {
        combined[[gid]][[field]] <- batch$reduced[[field]][, i, , , drop = FALSE]
      }
    }
  }
  
  return(combined)
}

##### Get args ######
get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) stop(paste("Missing value for", flag))
  args[[hit + 1]]
}

##### Failure log #####
FailureLog <- setRefClass("FailureLog",
                          fields = list(entries = "list"),
                          methods = list(
                            add = function(row) {
                              entries[[length(entries) + 1]] <<- row
                            },
                            as_df = function() {
                              if (length(entries) == 0) {
                                return(data.frame(
                                  batch = integer(), scenario = character(), attempt = integer(),
                                  outcome = character(), error_msg = character(), timestamp = character()
                                ))
                              }
                              do.call(rbind, entries)
                            }
                          )
)

##### Make equilibration df #####
make_equil_df <- function(vec, outcome_name) {
  data.frame(
    outcome = outcome_name,
    time    = seq_along(vec),
    value   = vec
  )
}

##### Mean age of cases (national) #####
calc_national_mean_age_of_cases <- function(combined_reduced, all_gids, scenarios_to_run, age_vector) {
  
  mean_age_national <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  for (sc in scenarios_to_run) {
    arr         <- combined_reduced[[sc]][[all_gids[1]]][["cases_year"]]
    n_times     <- dim(arr)[4]
    n_particles <- dim(arr)[3]
    
    national_age_weighted <- matrix(0, nrow = n_particles, ncol = n_times)
    national_cases_total  <- matrix(0, nrow = n_particles, ncol = n_times)
    
    for (gid in all_gids) {
      arr  <- combined_reduced[[sc]][[gid]][["cases_year"]]
      arr3 <- array(arr[, 1, , ], dim = c(dim(arr)[1], n_particles, n_times))
      
      for (t in seq_len(n_times)) {
        national_cases_total[, t]  <- national_cases_total[, t]  + colSums(arr3[, , t])
        national_age_weighted[, t] <- national_age_weighted[, t] + colSums(age_vector * arr3[, , t])
      }
    }
    
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
    mean_age_national[[sc]] <- do.call(rbind, rows)
    rownames(mean_age_national[[sc]]) <- NULL
  }
  return(mean_age_national)
}

##### Mean age of cases (municipality) #####
calc_regional_mean_age_of_cases <- function(combined_reduced, all_gids, scenarios_to_run, age_vector) {
  
  mean_age_regional <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  for (sc in scenarios_to_run) {
    arr         <- combined_reduced[[sc]][[all_gids[1]]][["cases_year"]]
    n_times     <- dim(arr)[4]
    n_particles <- dim(arr)[3]
    
    sc_results <- setNames(vector("list", length(all_gids)), all_gids)
    
    for (gid in all_gids) {
      arr  <- combined_reduced[[sc]][[gid]][["cases_year"]]
      arr3 <- array(arr[, 1, , ], dim = c(dim(arr)[1], n_particles, n_times))
      
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
      }
      sc_results[[gid]] <- do.call(rbind, rows)
    }
    mean_age_regional[[sc]] <- do.call(rbind, sc_results)
    rownames(mean_age_regional[[sc]]) <- NULL
  }
  return(mean_age_regional)
}

##### Proportion susceptible (national) — returns per-particle values #####
calc_national_prop_susceptible <- function(combined_reduced, all_gids, scenarios_to_run, demog_for_susc) {
  
  prop_susc_national <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  arr1        <- combined_reduced[[scenarios_to_run[1]]][[all_gids[1]]][["susceptibles"]]
  n_times     <- dim(arr1)[4]
  n_particles <- dim(arr1)[3]
  
  for (sc in scenarios_to_run) {
    national_susc <- matrix(0, nrow = n_particles, ncol = n_times)
    national_pop  <- numeric(n_times)
    
    for (gid_idx in seq_along(all_gids)) {
      gid   <- all_gids[gid_idx]
      susc  <- combined_reduced[[sc]][[gid]][["susceptibles"]]
      susc3 <- array(susc[, 1, , ], dim = c(dim(susc)[1], n_particles, n_times))
      
      for (t in seq_len(n_times)) {
        national_susc[, t] <- national_susc[, t] + colSums(susc3[, , t])
        national_pop[t]    <- national_pop[t]    + sum(demog_for_susc[, gid_idx, t])
      }
    }
    
    # Return per-particle proportions — summarised in 03_combine_batches.R
    prop_susc_national[[sc]] <- national_susc / matrix(national_pop, nrow = n_particles, ncol = n_times, byrow = TRUE)
  }
  return(prop_susc_national)
}

##### Proportion susceptible (municipality) — returns per-particle values #####
calc_regional_prop_susceptible <- function(combined_reduced, all_gids, scenarios_to_run, demog_for_susc) {
  
  prop_susc_regional <- setNames(vector("list", length(scenarios_to_run)), scenarios_to_run)
  
  arr1        <- combined_reduced[[scenarios_to_run[1]]][[all_gids[1]]][["susceptibles"]]
  n_times     <- dim(arr1)[4]
  n_particles <- dim(arr1)[3]
  
  for (sc in scenarios_to_run) {
    sc_results <- setNames(vector("list", length(all_gids)), all_gids)
    
    for (gid_idx in seq_along(all_gids)) {
      gid   <- all_gids[gid_idx]
      susc  <- combined_reduced[[sc]][[gid]][["susceptibles"]]
      susc3 <- array(susc[, 1, , ], dim = c(dim(susc)[1], n_particles, n_times))
      
      # Per-particle proportion for this municipality
      pop_by_time <- apply(demog_for_susc[, gid_idx, , drop = FALSE], 3, sum)
      prop <- matrix(0, nrow = n_particles, ncol = n_times)
      for (t in seq_len(n_times)) {
        prop[, t] <- colSums(susc3[, , t]) / pop_by_time[t]
      }
      
      # Store as per-particle matrix for summarising in 03_combine_batches.R
      sc_results[[gid]] <- prop
    }
    prop_susc_regional[[sc]] <- sc_results
  }
  return(prop_susc_regional)
}

run_with_retry <- function(fn, max_attempts = 3, failure_log = NULL, batch = NULL, log_scenario = NULL, ...) {
  attempt <- 1
  last_error <- NULL
  
  while (attempt <= max_attempts) {
    tryCatch({
      result <- fn(...)
      
      if (attempt > 1 & !is.null(failure_log)) {
        failure_log$add(data.frame(
          batch     = if (!is.null(batch)) batch else NA,
          scenario  = if (!is.null(log_scenario)) log_scenario else NA,
          attempt   = attempt,
          outcome   = "retry_success",
          error_msg = NA,
          timestamp = as.character(Sys.time())
        ))
      }
      return(result)
      
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      last_error <<- e
      
      is_binomial_error <- grepl("binomial|binom|particle.*error", error_msg, ignore.case = TRUE) ||
        grepl("Invalid call to binomial", error_msg, fixed = TRUE)
      
      if (is_binomial_error) {
        
        if (!is.null(failure_log)) {
          failure_log$add(data.frame(
            batch     = if (!is.null(batch)) batch else NA,
            scenario  = if (!is.null(log_scenario)) log_scenario else NA,
            attempt   = attempt,
            outcome   = ifelse(attempt < max_attempts, "binomial_retry", "binomial_failed"),
            error_msg = substr(error_msg, 1, 250),
            timestamp = as.character(Sys.time())
          ))
        }
        
        if (attempt < max_attempts) {
          attempt <<- attempt + 1
          Sys.sleep(0.1)
        } else {
          stop(e)
        }
        
      } else {
        
        if (!is.null(failure_log)) {
          failure_log$add(data.frame(
            batch     = if (!is.null(batch)) batch else NA,
            scenario  = if (!is.null(log_scenario)) log_scenario else NA,
            attempt   = attempt,
            outcome   = "non_binomial_failed",
            error_msg = substr(error_msg, 1, 250),
            timestamp = as.character(Sys.time())
          ))
        }
        stop(e)
      }
    })
  }
  
  if (!is.null(last_error)) stop(last_error)
}
