# code/equilibrate_sim.R

# Defines equilibtate_sim() function used in 02_run_batch.R

equilibrate_sim <- function(demog_input, config, betas, output_path, denv_mod, denv_pars) {
  
  start_time <- Sys.time()
  cat(paste0("Running simulations for: equilibration \n"))
  
  # create dust system to run model
  denv_sys <- dust_system_create(denv_mod,
                                 pars = denv_pars,
                                 n_particles = 1,
                                 n_threads = 1) # n_threads is the number of parallel threads to run, useful when n_particles>1
  
  dust_system_set_state_initial(denv_sys) # use initial conditions defined in code
  
  
  t <- c(1, seq(365, 365*config$calibration_years, by = 365)) # define output time points for equilibration
  denv_out <- dust_system_simulate(denv_sys, t) # runs model and returns packed state array at times t
  baseline_states <- denv_out[,length(t)]
  qsave(denv_out[, length(t)], here(output_path, paste0("baseline_states.qs"))) # save baseline states
  
  denv_out <- dust_unpack_state(denv_sys, denv_out) # unpack output
  denv_out_reduced <- list(
#    cases_cumulative = denv_out$cases_cumulative,
#    hosp_cumulative = denv_out$hosp_cumulative,
    cases_year = denv_out$cases_year,
    hosp_year = denv_out$hosp_year,
    inf_year = denv_out$inf_year,
    S = denv_out$S
#    inf_cumulative = denv_out$inf_cumulative
    )  # reduce number of outputs saved
#  qsave(denv_out_reduced, here(output_path, paste0("reduced_equilibration.qs"))) # to plot dynamics over equilibration period
  
  print(Sys.time() - start_time)
  cat(paste0("Finished equilibration\n"))
  
  # Return useful objects
  return(list(
    baseline_states = baseline_states,
    denv_out_reduced = denv_out_reduced,
    denv_pars = denv_pars
  ))
}

