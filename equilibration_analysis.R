# code/equilibration_analysis.R
# Produces plots of the equilibration period

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(qs)


path        <- "/Users/dave/Documents/GitHub/denv-indo"
config      <- yaml::read_yaml(file.path(path, "code", "config.yaml"))
output_path <- file.path(path, "output", "main_500")
figure_path <- file.path(path, "figures", "main_500")
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

denv_out_reduced <- qread(file.path(output_path, "equil_summary.qs"))

##### Report theme #####
theme_report <- theme_bw(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.title = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    
    axis.title.x = element_text(size = 15, margin = margin(t = 12)),
    axis.title.y = element_text(size = 15, margin = margin(r = 12)),
    axis.text = element_text(size = 12, colour = "black")
  )


##### Plot cases and hospitalisations #####
p_case <- denv_out_reduced |>
  filter(outcome %in% c("cases_year", "hosp_year")) |>
  mutate(
    outcome = recode(outcome, cases_year = "Cases", hosp_year = "Hospitalisations")
  ) |>
  ggplot(aes(x = time, y = value, colour = outcome)) +
  geom_line(linewidth = 0.6) +
  scale_y_continuous(
    labels = label_number(),
    limits = c(0, 60000000)
  ) +
  labs(
    x      = "Years",
    y      = "Count",
    colour = "Outcome"
  ) +
  theme_report


ggsave(file.path(figure_path, "equilibration_cases.png"),
       p_case, width = 8, height = 5, dpi = 300)

cat("Saved: equilibration_cases.png\n")

##### Plot infections #####
p_inf <- denv_out_reduced |>
  filter(outcome == "inf_year") |>
  ggplot(aes(x = time, y = value / 1e6)) +
  geom_line(linewidth = 0.6, colour = "steelblue") +
  scale_y_continuous(labels = label_number()) +
  labs(
    x = "Years",
    y = "Infections (millions)"
  ) +
  theme_report

ggsave(file.path(figure_path, "equilibration_infections.png"),
       p_inf, width = 8, height = 5, dpi = 300)

cat("Saved: equilibration_infections.png\n")
