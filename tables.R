# code/tables.R

# Tables of cumulative cases and hospitalisations by scenario

library(dplyr)
library(tidyr)
library(qs)
library(yaml)
library(readr)


# Paths
path        <- ""
config      <- yaml::read_yaml(file.path(path, "code", "config.yaml"))
output_path <- file.path(path, "output", "filtered_500")
table_path  <- file.path(path, "tables", "filtered_500")
dir.create(table_path, recursive = TRUE, showWarnings = FALSE)


# Load results
cumulative_national <- qread(file.path(output_path, "cumulative_national_excl.qs"))


# Labels
scenarios       <- c("baseline", "wolbachia", "vaccine_national", "wolbachia_vaccine_national")
scenario_labels <- c("Baseline", "Wolbachia", "Vaccine National", "Wolbachia + Vaccine")
scenario_map    <- setNames(scenario_labels, scenarios)


# Helpers
fmt_num <- function(x, digits = 0) {
  format(round(x, digits), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_qi <- function(median, q025, q975, digits = 0) {
  paste0(
    fmt_num(median, digits), " (",
    fmt_num(q025, digits), "\u2013", fmt_num(q975, digits), ")"
  )
}

make_summary_table <- function(results_list, outcomes, outcome_labels) {
  bind_rows(lapply(names(outcomes), function(outcome_name) {
    bind_rows(lapply(scenarios, function(sc) {
      df <- results_list[[outcome_name]][[sc]]
      if (is.null(df)) return(NULL)

      df |>
        filter(year == max(year, na.rm = TRUE)) |>
        mutate(
          scenario = scenario_map[sc],
          outcome  = outcome_labels[[outcome_name]],
          summary  = fmt_qi(median, q025, q975)
        ) |>
        select(scenario, outcome, summary)
    }))
  })) |>
    mutate(scenario = factor(scenario, levels = scenario_labels)) |>
    arrange(scenario) |>
    pivot_wider(names_from = outcome, values_from = summary) |>
    rename(Scenario = scenario)
}


##### Population-level table #####
population_outcomes <- c(
  "cases_cumulative" = "Cases",
  "hosp_cumulative"  = "Hospitalisations"
)

population_table <- make_summary_table(
  results_list    = cumulative_national,
  outcomes        = population_outcomes,
  outcome_labels  = population_outcomes
)

write_csv(
  population_table,
  file.path(table_path, "table_population_actual_cases_hospitalisations.csv")
)

cat("\nPopulation-level table:\n")
print(population_table, n = Inf)


##### Individual-level table #####
individual_outcomes <- c(
  "vacc_cases_cumulative" = "Cases",
  "vacc_hosp_cumulative"  = "Hospitalisations"
)

individual_table <- make_summary_table(
  results_list    = cumulative_national,
  outcomes        = individual_outcomes,
  outcome_labels  = individual_outcomes
)

write_csv(
  individual_table,
  file.path(table_path, "table_supplement_individual_actual_cases_hospitalisations.csv")
)

cat("\nIndividual-level table (supplement):\n")
print(individual_table, n = Inf)


# Export HTML tables if gt is installed
if (requireNamespace("gt", quietly = TRUE)) {
  library(gt)

  gt_population <- population_table |>
    gt() |>
    tab_header(
      title = md("**Table X. Population-level cumulative cases and hospitalisations by scenario**"),
      subtitle = "Values are median (95% quantile interval) over 500 stochastic simulations for the 15-year cumulative total."
    )

  gt_individual <- individual_table |>
    gt() |>
    tab_header(
      title = md("**Table S1. Individual-level cumulative cases and hospitalisations by scenario**"),
      subtitle = "Values are median (95% quantile interval) over 500 stochastic simulations for the 15-year cumulative total in vaccinated individuals, or the equivalent subgroup in non-vaccination scenarios."
    )

  gtsave(gt_population, file.path(table_path, "table_population_actual_cases_hospitalisations.html"))
  gtsave(gt_individual, file.path(table_path, "table_supplement_individual_actual_cases_hospitalisations.html"))
}

cat("\nSaved:\n")
cat(" - table_population_actual_cases_hospitalisations.csv\n")
cat(" - table_supplement_individual_actual_cases_hospitalisations.csv\n")
cat("and, if gt is installed:\n")
cat(" - table_population_actual_cases_hospitalisations.html\n")
cat(" - table_supplement_individual_actual_cases_hospitalisations.html\n")