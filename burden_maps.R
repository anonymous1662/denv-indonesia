# code/burden_maps.R

# Plot maps
# 1) Maps of cases averted per 10,000
# 2) Maps of hospitalisations averted per 10,000


library(here)
library(dplyr)
library(ggplot2)
library(scales)
library(qs)
library(sf)
library(janitor)
library(yaml)

path <- ""
config <- yaml::read_yaml(file.path(path, "code", "config.yaml"))

output_path <- file.path(path, "output", "filtered_500")
figure_path <- file.path(path, "figures", "report")
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)


##### Load data #####
averted_results <- qread(here(output_path, "averted_results_excl.qs"))

indo_adm2_shp <- st_read(
  here(path, "data", "geography", "gadm41_IDN_shp", "gadm41_IDN_2.shp"),
  quiet = TRUE
) |>
  clean_names()

pop_summary <- qread(here(path, "data", "demography3", "demog_adm2.qs")) |>
  filter(year == 2025) |>
  group_by(adm2_code) |>
  summarise(population = sum(population), .groups = "drop")


##### Shared labels #####
scenarios <- c("wolbachia", "vaccine_national", "wolbachia_vaccine_national")

intervention_labels_map <- c(
  "wolbachia"                  = "italic(w)*'Mel'",
  "vaccine_national"           = "'Vacc.'",
  "wolbachia_vaccine_national" = "italic(w)*'Mel'~'&'~'Vacc.'"
)



scenario_order_map <- c("italic(w)*'Mel'", "'Vacc.'", "italic(w)*'Mel'~'&'~'Vacc.'")


##### Map helpers #####
build_map_sf <- function(outcome_key) {
  bind_rows(lapply(scenarios, function(sc) {
    averted_results[[outcome_key]][[sc]] |>
      filter(year == max(year, na.rm = TRUE)) |>
      mutate(scenario = intervention_labels_map[sc])
  })) |>
    left_join(pop_summary, by = c("gid_2" = "adm2_code")) |>
    mutate(
      averted_per_10000 = 10000 * median / population,
      scenario = factor(scenario, levels = scenario_order_map)
    )
}

join_shapefile <- function(df) {
  indo_adm2_shp |>
    mutate(gid_2 = as.character(gid_2)) |>
    left_join(df, by = "gid_2") |>
    filter(!is.na(scenario))
}

map_theme <- theme_void(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#4a4a4a", colour = "grey30", linewidth = 0.6),
    strip.text       = element_text(face = "plain", size = 22, colour = "white"),
    panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
    legend.position  = "right",
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 16),
    legend.margin    = margin(0, 0, 0, 15),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(20, 20, 20, 20),
    panel.spacing    = unit(1, "lines")
  )

make_map_plot <- function(sf_data, cap_val, legend_title, low_val) {
  ggplot(sf_data) +
    geom_sf(aes(fill = averted_capped), colour = "grey30", linewidth = 0.05) +
    facet_wrap(~scenario, ncol = 1, labeller = label_parsed)+
    scale_fill_viridis_c(
      option    = "viridis",
      direction = -1,
      begin     = 0.1,
      end       = 0.95,
      limits    = c(low_val, cap_val),
      name      = legend_title,
      na.value  = "grey85",
      labels    = scales::comma,
      guide     = guide_colorbar(barwidth = 0.8, barheight = 12, ticks = FALSE)
    ) +
    map_theme
}


##### Figure 1: Cases averted map #####
map_cases <- build_map_sf("cases_cumulative")
cap_cases <- quantile(map_cases$averted_per_10000, 0.95, na.rm = TRUE)

map_cases <- map_cases |>
  mutate(averted_capped = pmin(averted_per_10000, cap_cases, na.rm = FALSE))

p_cases <- make_map_plot(
  join_shapefile(map_cases),
  cap_cases,
  "Cases averted\nper 10,000",
  -10
)

ggsave(
  file.path(figure_path, "figure_map_cases_averted_3x1.png"),
  p_cases,
  width = 12,
  height = 14,
  dpi = 1200
)

cat("Saved: figure_map_cases_averted_3x1.png\n")


##### Figure 2: Hospitalisations averted map #####
map_hosp <- build_map_sf("hosp_cumulative")
cap_hosp <- quantile(map_hosp$averted_per_10000, 0.95, na.rm = TRUE)

map_hosp <- map_hosp |>
  mutate(averted_capped = pmin(averted_per_10000, cap_hosp, na.rm = FALSE))

p_hosp <- make_map_plot(
  join_shapefile(map_hosp),
  cap_hosp,
  "Hospitalisations\naverted per 10,000",
  -1
)

ggsave(
  file.path(figure_path, "figure_map_hosp_averted_3x1.png"),
  p_hosp,
  width = 12,
  height = 14,
  dpi = 1200
)

cat("Saved: figure_map_hosp_averted_3x1.png\n")
