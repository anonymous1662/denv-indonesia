# code/report_figures.R

# Plot non-map figures
# 1) Final-year percentage averted bar charts
# 2) National cumulative burden across 15 years
# 3) Mean age of cases
# 4) Cumulative burden in vaccinated / equivalent population across 15 years

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(qs)
library(yaml)
library(grid)

path <- ""
config <- yaml::read_yaml(file.path(path, "code", "config.yaml"))

output_path <- file.path(path, "output", "filtered_500")
figure_path <- file.path(path, "figures", "report")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)


##### Load data #####
cumulative_national <- qread(here(output_path, "cumulative_national_excl.qs"))
averted_national    <- qread(here(output_path, "averted_national_excl.qs"))
mean_age_national   <- qread(here(output_path, "mean_age_national_excl.qs"))

##### Shared labels and colours: intervention-only figures #####
scenarios <- c("wolbachia", "vaccine_national", "wolbachia_vaccine_national")

intervention_labels <- c(
  "wolbachia"                  = "wMel",
  "vaccine_national"           = "Vacc.",
  "wolbachia_vaccine_national" = "wMel &\nVacc."
)

scenario_order <- c("wMel", "Vacc.", "wMel &\nVacc.")

scenario_colors <- c(
  "wMel"          = "#D55E00",
  "Vacc."         = "#009E73",
  "wMel &\nVacc." = "#0072B2"
)

##### Shared labels and colours: all-scenario figures #####
scenarios_all <- c("baseline", "wolbachia", "vaccine_national", "wolbachia_vaccine_national")

intervention_labels_all <- c(
  "baseline"                   = "Baseline",
  "wolbachia"                  = "wMel",
  "vaccine_national"           = "Vacc.",
  "wolbachia_vaccine_national" = "wMel & Vacc."
)

scenario_order_all <- c("Baseline", "wMel", "Vacc.", "wMel & Vacc.")

scenario_colors_all <- c(
  "Baseline"      = "#999999",
  "wMel"          = "#D55E00",
  "Vacc."         = "#009E73",
  "wMel & Vacc."  = "#0072B2"
)

scenario_linetypes_all <- c(
  "Baseline"      = "solid",
  "wMel"          = "dashed",
  "Vacc."         = "dashed",
  "wMel & Vacc."  = "dashed"
)


##### Figure 1: 2x2 bar chart of % averted at final year #####
extract_final_pct_averted <- function(outcome_name, population_group_label) {
  lapply(scenarios, function(sc) {
    df <- averted_national[[outcome_name]][[sc]]
    if (is.null(df)) return(NULL)
    
    df |>
      filter(year == max(year, na.rm = TRUE)) |>
      transmute(
        scenario = intervention_labels[sc],
        outcome = case_when(
          outcome_name %in% c("cases_cumulative", "vacc_cases_cumulative") ~ "Cases",
          outcome_name %in% c("hosp_cumulative", "vacc_hosp_cumulative")   ~ "Hospitalisations"
        ),
        population_group = population_group_label,
        median = pct_median,
        q025   = pct_q025,
        q975   = pct_q975
      )
  }) |>
    bind_rows() |>
    mutate(
      scenario = factor(scenario, levels = scenario_order),
      outcome = factor(outcome, levels = c("Cases", "Hospitalisations")),
      population_group = factor(
        population_group,
        levels = c("Population-level impact", "Individual-level impact")
      )
    )
}

bar_df <- bind_rows(
  extract_final_pct_averted("cases_cumulative",      "Population-level impact"),
  extract_final_pct_averted("hosp_cumulative",       "Population-level impact"),
  extract_final_pct_averted("vacc_cases_cumulative", "Individual-level impact"),
  extract_final_pct_averted("vacc_hosp_cumulative",  "Individual-level impact")
)

p_bar_main <- ggplot(bar_df, aes(x = scenario, y = median, fill = scenario)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.15, linewidth = 0.5) +
  facet_grid(population_group ~ outcome, scales = "fixed") +
  scale_fill_manual(
    values = scenario_colors,
    labels = c(
      "wMel"          = expression(italic(w) * "Mel"),
      "Vacc."         = "Vacc.",
      "wMel &\nVacc." = expression(italic(w) * "Mel & Vacc.")
    )
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = seq(0, 100, by = 20),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Proportion averted (%)") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position      = "bottom",
    legend.title         = element_blank(),
    legend.text          = element_text(size = 11),
    legend.key.size      = unit(0.8, "cm"),
    strip.background     = element_rect(fill = "#4a4a4a", colour = "grey30", linewidth = 0.6),
    strip.text           = element_text(face = "plain", size = 14, colour = "white"),
    panel.border         = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
    panel.grid.major.x   = element_line(colour = "grey80", linewidth = 0.4),
    panel.spacing        = unit(0.3, "lines"),
    panel.spacing.y      = unit(2, "lines"),
    axis.text.x          = element_text(face = "plain"),
    axis.ticks.x         = element_line(colour = "grey40", linewidth = 0.4),
    axis.ticks.length.x  = unit(3, "pt"),
    plot.margin          = margin(20, 20, 20, 20)
  )

ggsave(
  file.path(figure_path, "figure_main_pct_averted_bar.png"),
  p_bar_main,
  width = 10,
  height = 8,
  dpi = 300
)

cat("Saved: figure_main_pct_averted_bar.png\n")

##### Shared line-plot helpers #####
line_plot_cumulative <- function(df, y_label, y_format = comma) {
  ggplot(
    df,
    aes(
      x = year,
      y = median,
      colour = scenario,
      fill = scenario,
      linetype = scenario
    )
  ) +
    geom_ribbon(
      aes(ymin = q025, ymax = q975),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    
    scale_colour_manual(values = scenario_colors_all, drop = FALSE) +
    scale_fill_manual(values = scenario_colors_all, drop = FALSE) +
    scale_linetype_manual(values = scenario_linetypes_all, drop = FALSE) +
    
    scale_x_continuous(
      breaks = c(2025, 2030, 2035, 2040),
      limits = c(2025, 2040)
    ) +
    scale_y_continuous(labels = y_format) +
    
    labs(
      x = "Year",
      y = y_label,
      colour = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 13),
      legend.key.size = unit(0.8, "cm"),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
      panel.spacing = unit(0.1, "lines"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
      
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y = element_text(margin = margin(r = 12))
    )
}

line_plot_report <- function(df, y_label, y_format = comma) {
  ggplot(
    df,
    aes(
      x = year,
      y = median,
      colour = scenario,
      fill = scenario,
      linetype = scenario
    )
  ) +
    geom_ribbon(
      aes(ymin = q025, ymax = q975),
      alpha = 0.15,
      colour = NA
    ) +
    geom_line(linewidth = 0.9) +
    
    scale_colour_manual(values = scenario_colors_all, drop = FALSE) +
    scale_fill_manual(values = scenario_colors_all, drop = FALSE) +
    scale_linetype_manual(values = scenario_linetypes_all, drop = FALSE) +
    
    scale_x_continuous(
      breaks = c(2025, 2030, 2035, 2040),
      limits = c(2025, 2040)
    ) +
    scale_y_continuous(
      limits = c(18, 25),
      breaks = seq(18, 25, by = 2),
      labels = y_format
    ) +
    
    labs(
      x = "Year",
      y = y_label,
      colour = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.8, "cm"),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
      panel.spacing = unit(0.1, "lines"),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
      
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.title.y = element_text(margin = margin(r = 12))
    )
}

##### Figure 2: National cumulative burden across 15 years #####
cumulative_df <- bind_rows(lapply(scenarios_all, function(sc) {
  bind_rows(lapply(c("cases_cumulative", "hosp_cumulative"), function(outcome) {
    df <- cumulative_national[[outcome]][[sc]]
    if (is.null(df)) return(NULL)
    
    df |>
      mutate(
        scenario = intervention_labels_all[sc],
        outcome  = outcome,
        year     = config$scenario_start + year - 1
      )
  }))
})) |>
  mutate(scenario = factor(scenario, levels = scenario_order_all))

p_cumulative_cases <- line_plot_cumulative(
  cumulative_df |> filter(outcome == "cases_cumulative"),
  "Cumulative cases"
)

p_cumulative_hosp <- line_plot_cumulative(
  cumulative_df |> filter(outcome == "hosp_cumulative"),
  "Cumulative hospitalisations"
)

combined_cumulative <- (p_cumulative_cases | p_cumulative_hosp) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_path, "figure_national_cumulative_burden.png"),
  combined_cumulative,
  width = 12,
  height = 5,
  dpi = 300
)

cat("Saved: figure_national_cumulative_burden.png\n")


##### Build mean_age_df #####
mean_age_df <- bind_rows(lapply(scenarios_all, function(sc) {
  mean_age_national[[sc]] |>
    mutate(
      scenario = intervention_labels_all[sc],
      year     = config$scenario_start + year - 1
    )
})) |>
  mutate(scenario = factor(scenario, levels = scenario_order_all))


##### Figure 3: Mean age of cases #####
p_mean_age <- line_plot_report(
  mean_age_df,
  y_label = "Mean age of cases"
)

ggsave(
  file.path(figure_path, "figure_age.png"),
  p_mean_age,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Saved: figure_age.png\n")


##### Figure 4: Cumulative burden in vaccinated / equivalent population across 15 years #####
individual_cumulative_df <- bind_rows(lapply(scenarios_all, function(sc) {
  bind_rows(lapply(c("vacc_cases_cumulative", "vacc_hosp_cumulative"), function(outcome) {
    df <- cumulative_national[[outcome]][[sc]]
    if (is.null(df)) return(NULL)
    
    df |>
      mutate(
        scenario = intervention_labels_all[sc],
        outcome  = outcome,
        year     = config$scenario_start + year - 1
      )
  }))
})) |>
  mutate(scenario = factor(scenario, levels = scenario_order_all))

p_individual_cumulative_cases <- line_plot_cumulative(
  individual_cumulative_df |> filter(outcome == "vacc_cases_cumulative"),
  "Cumulative cases"
)

p_individual_cumulative_hosp <- line_plot_cumulative(
  individual_cumulative_df |> filter(outcome == "vacc_hosp_cumulative"),
  "Cumulative hospitalisations"
)

combined_individual_cumulative <- (p_individual_cumulative_cases | p_individual_cumulative_hosp) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_path, "figure_individual_cumulative_burden.png"),
  combined_individual_cumulative,
  width = 12,
  height = 5,
  dpi = 300
)

cat("Saved: figure_individual_cumulative_burden.png\n")


##### Figure 5: 2x2 cumulative burden figure #####
# Row a = Population-level burden
# Row b = Individual-level burden

library(ggh4x)

burden_2x2_df <- bind_rows(
  cumulative_df |>
    filter(outcome %in% c("cases_cumulative", "hosp_cumulative")) |>
    mutate(
      burden_level = "Population-level",
      outcome_label = case_when(
        outcome == "cases_cumulative" ~ "Cases",
        outcome == "hosp_cumulative"  ~ "Hospitalisations"
      )
    ),
  individual_cumulative_df |>
    filter(outcome %in% c("vacc_cases_cumulative", "vacc_hosp_cumulative")) |>
    mutate(
      burden_level = "Individual-level",
      outcome_label = case_when(
        outcome == "vacc_cases_cumulative" ~ "Cases",
        outcome == "vacc_hosp_cumulative"  ~ "Hospitalisations"
      )
    )
) |>
  mutate(
    burden_level = factor(
      burden_level,
      levels = c("Population-level", "Individual-level")
    ),
    outcome_label = factor(
      outcome_label,
      levels = c("Cases", "Hospitalisations")
    ),
    scenario = factor(scenario, levels = scenario_order_all)
  )

p_burden_2x2 <- ggplot(
  burden_2x2_df,
  aes(
    x = year,
    y = median,
    colour = scenario,
    fill = scenario,
    linetype = scenario
  )
) +
  geom_ribbon(
    aes(ymin = q025, ymax = q975),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  ggh4x::facet_grid2(
    burden_level ~ outcome_label,
    scales = "free_y",
    independent = "y"
  ) +
  scale_colour_manual(values = scenario_colors_all, drop = FALSE) +
  scale_fill_manual(values = scenario_colors_all, drop = FALSE) +
  scale_linetype_manual(values = scenario_linetypes_all, drop = FALSE) +
  scale_x_continuous(
    breaks = c(2025, 2030, 2035, 2040),
    limits = c(2025, 2040)
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Year",
    y = "Cumulative burden",
    colour = NULL,
    fill = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "bottom",
    legend.title       = element_blank(),
    legend.text        = element_text(size = 11),
    legend.key.size    = unit(0.8, "cm"),
    strip.background   = element_rect(fill = "#4a4a4a", colour = "grey30", linewidth = 0.6),
    strip.text         = element_text(face = "plain", size = 13, colour = "white"),
    panel.border       = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
    panel.spacing      = unit(0.6, "lines"),
    plot.margin        = margin(20, 20, 20, 20),
    axis.title.x       = element_text(margin = margin(t = 10)),
    axis.title.y       = element_text(margin = margin(r = 10))
  )

ggsave(
  file.path(figure_path, "figure_burden_2x2.png"),
  p_burden_2x2,
  width = 12,
  height = 8,
  dpi = 300
)

cat("Saved: figure_burden_2x2.png\n")
