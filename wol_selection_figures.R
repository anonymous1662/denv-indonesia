# code/wol_selection_figures.R

# Top 100 municipalities map + scatter plots

library(here)
library(sf)
library(ggplot2)
library(viridis)
library(dplyr)
library(patchwork)
library(scales)
library(qs)
library(readxl)
library(janitor)
library(ggrepel)
library(stringr)

source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))
source(here("code", "calculate-R0_fn.R"))

figure_path <- here("figures", "choose_regions")
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

##### Load data #####
indo_adm2_shp <- st_read(
  here("data", "geography", "gadm41_IDN_shp", "gadm41_IDN_2.shp"),
  quiet = TRUE
) |>
  clean_names()

indo_foi_full <- read_xlsx(
  here("data", "foi", "imperial-estimates", "Estimates_for_Gavi_20_Aug_25.xlsx"),
  sheet = 3
) |>
  clean_names() |>
  filter(country == "Indonesia") |>
  rename(name_1 = region, name_2 = district, gid_2 = id_2, gid_1 = id_1) |>
  arrange(gid_2) |>
  mutate(
    foi       = total_foi / 4,
    foi_plot  = total_foi / 4
  )

indo_demog_full <- qread(here("data", "demography3", "demog_adm2.qs")) |>
  arrange(adm2_code)

indo_demog_full <- expand_age_groups(indo_demog_full, max_age = 90)

pwd_adm2 <- qread(here("data", "pwd_adm2.qs"))


##### Calculate R0 for all municipalities #####
all_gids <- unique(indo_foi_full$gid_2)
r0s_full <- numeric(length(all_gids))

cat(sprintf("Calculating R0 for %d municipalities...\n", length(all_gids)))

for (i in seq_along(all_gids)) {
  id <- all_gids[i]
  
  foi_value <- indo_foi_full |>
    filter(gid_2 == id) |>
    pull(foi)
  
  age_proportions <- indo_demog_full |>
    filter(adm2_code == id, year == 2025) |>
    mutate(prop = population / sum(population)) |>
    pull(prop)
  
  age_lower <- unique(indo_demog_full$age)
  age_upper <- c(age_lower[-1], 91)
  
  r0s_full[i] <- calculate_R0(
    FOI   = foi_value,
    n_j   = age_proportions,
    u_lim = age_upper,
    l_lim = age_lower,
    phis  = c(1, 1, 1, 1)
  )
  
  if (i %% 50 == 0) {
    cat(sprintf("  Processed %d/%d municipalities\n", i, length(all_gids)))
  }
}

##### Build municipality-level summary #####
indo_adm2_shp <- indo_adm2_shp |>
  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6)

pop_summary <- indo_demog_full |>
  filter(year == 2025) |>
  group_by(adm2_code) |>
  summarise(population = sum(population), .groups = "drop")

map_data <- indo_adm2_shp |>
  left_join(pop_summary, by = c("gid_2" = "adm2_code")) |>
  left_join(data.frame(gid_2 = all_gids, R0 = r0s_full), by = "gid_2") |>
  left_join(indo_foi_full |> select(gid_2, foi), by = "gid_2") |>
  left_join(pwd_adm2, by = "gid_2") |>
  mutate(
    pop_density          = population / area_km2,
    pop_weighted_density = pwd
  )


##### Island groups and colours #####
island_map <- c(
  "IDN.1"  = "Sumatra", "IDN.3"  = "Sumatra", "IDN.5"  = "Sumatra",
  "IDN.8"  = "Sumatra", "IDN.16" = "Sumatra", "IDN.17" = "Sumatra",
  "IDN.24" = "Sumatra", "IDN.30" = "Sumatra", "IDN.31" = "Sumatra",
  "IDN.32" = "Sumatra",
  "IDN.4"  = "Java",    "IDN.7"  = "Java",    "IDN.9"  = "Java",
  "IDN.10" = "Java",    "IDN.11" = "Java",    "IDN.33" = "Java",
  "IDN.12" = "Borneo",  "IDN.13" = "Borneo",  "IDN.14" = "Borneo",
  "IDN.34" = "Borneo",  "IDN.35" = "Borneo",
  "IDN.6"  = "Sulawesi","IDN.25" = "Sulawesi","IDN.26" = "Sulawesi",
  "IDN.27" = "Sulawesi","IDN.28" = "Sulawesi","IDN.29" = "Sulawesi",
  "IDN.2"  = "Bali & Nusa Tenggara", "IDN.20" = "Bali & Nusa Tenggara",
  "IDN.21" = "Bali & Nusa Tenggara",
  "IDN.18" = "Maluku",  "IDN.19" = "Maluku",
  "IDN.22" = "Papua",   "IDN.23" = "Papua"
)

island_colors <- c(
  "Java"                 = "#E63946",
  "Sumatra"              = "#2196F3",
  "Borneo"               = "#2D6A4F",
  "Sulawesi"             = "#FF9F1C",
  "Bali & Nusa Tenggara" = "#9B5DE5",
  "Maluku"               = "#00B4D8",
  "Papua"                = "#606C38"
)


##### Weighted ranking — top 100 only #####
shp_data <- indo_adm2_shp |>
  st_drop_geometry() |>
  select(gid_2, area_km2) |>
  left_join(pop_summary, by = c("gid_2" = "adm2_code")) |>
  left_join(pwd_adm2, by = "gid_2") |>
  mutate(pop_density = population / area_km2)

top100_data <- indo_foi_full |>
  select(gid_2, total_foi, foi_plot) |>
  left_join(shp_data, by = "gid_2") |>
  filter(!is.na(pwd), !is.na(total_foi)) |>
  mutate(
    foi_scaled = (total_foi - min(total_foi, na.rm = TRUE)) /
      (max(total_foi, na.rm = TRUE) - min(total_foi, na.rm = TRUE)),
    pwd_scaled = (log(pwd + 0.001) - min(log(pwd + 0.001), na.rm = TRUE)) /
      (max(log(pwd + 0.001), na.rm = TRUE) - min(log(pwd + 0.001), na.rm = TRUE))
  ) |>
  mutate(
    foi_pwd      = foi_scaled * pwd_scaled,
    rank         = row_number(desc(foi_pwd)),
    top100       = rank <= 100,
    gid_1_prefix = str_extract(gid_2, "IDN\\.\\d+"),
    island       = island_map[gid_1_prefix]
  ) |>
  arrange(rank)


##### Map: top 100 municipalities highlighted #####
top100_map_data <- map_data |>
  left_join(top100_data |> select(gid_2, rank, island, foi_pwd), by = "gid_2") |>
  mutate(top100 = !is.na(rank) & rank <= 100)

p_top100_map <- ggplot() +
  geom_sf(data = top100_map_data, fill = "grey80", colour = "grey95", linewidth = 0.05) +
  geom_sf(
    data      = top100_map_data |> filter(top100),
    aes(fill  = island),
    colour    = "grey90",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = island_colors, name = "Island", na.value = "grey60") +
  annotate(
    "text", x = -Inf, y = -Inf,
    label = "Top 100 municipalities in bright colours",
    hjust = -0.05, vjust = -1, size = 7, fontface = "italic", colour = "grey30"
  ) +
  theme_void(base_size = 20) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.8),
    legend.position  = "none",
    plot.margin      = margin(50, 20, 2, 20)
  )

ggsave(
  here(figure_path, "map_top100_highlighted.png"),
  p_top100_map,
  width = 14, height = 7, dpi = 800
)


##### Scatter: top 100 with highlighted cities #####
highlight_cities <- c(
  "Jakarta Barat", "Kota Bandung", "Kota Semarang",
  "Kota Yogyakarta", "Bontang", "Kota Kupang"
)

top100_labelled <- top100_data |>
  left_join(indo_foi_full |> select(gid_2, name_2), by = "gid_2") |>
  mutate(highlight = name_2 %in% highlight_cities)

cat("Highlighted cities found:\n")
print(top100_labelled |> filter(highlight) |> select(name_2, rank, gid_2))

p_rank_highlighted <- ggplot(top100_labelled, aes(x = pwd, y = foi_plot)) +
  geom_point(aes(colour = island), alpha = 0.5, size = 2, shape = 16) +
  geom_point(
    data  = top100_labelled |> filter(top100 & !highlight),
    aes(colour = island),
    size = 4, shape = 21, fill = NA, stroke = 1.5, alpha = 0.8
  ) +
  geom_point(
    data  = top100_labelled |> filter(highlight),
    aes(fill = island),
    colour = "black", size = 5, shape = 23, stroke = 1
  ) +
  scale_colour_manual(
    values = island_colors, name = "Island",
    guide  = guide_legend(nrow = 3, override.aes = list(shape = 15, size = 5, alpha = 1))
  ) +
  scale_fill_manual(values = island_colors, guide = "none") +
  geom_point(
    data = data.frame(
      x    = NA_real_,
      y    = NA_real_,
      type = factor(
        c("All municipalities", "Top 100", "Highlighted cities"),
        levels = c("All municipalities", "Top 100", "Highlighted cities")
      )
    ),
    aes(x = x, y = y, shape = type),
    size = 4, colour = "grey40"
  ) +
  scale_shape_manual(
    name   = NULL,
    values = c("All municipalities" = 16, "Top 100" = 21, "Highlighted cities" = 23),
    guide  = guide_legend(
      override.aes = list(
        colour = c("grey40", "grey40", "black"),
        fill   = c("grey40", NA,       "grey40"),
        size   = c(2,        5,        4),
        stroke = c(0,        1.5,      1),
        alpha  = c(0.5,      0.8,      1)
      )
    )
  ) +
  geom_label_repel(
    data               = top100_labelled |> filter(highlight & name_2 != "Kota Semarang"),
    aes(label = name_2, colour = island),
    size               = 5,
    fontface           = "bold",
    fill               = alpha("white", 0.8),
    box.padding        = 4,
    point.padding      = 0,
    segment.linewidth  = 0.5,
    segment.curvature  = 0,
    min.segment.length = 0,
    max.overlaps       = Inf,
    force              = 2,
    force_pull         = 0.5,
    show.legend        = FALSE
  ) +
  geom_label_repel(
    data               = top100_labelled |> filter(name_2 == "Kota Semarang"),
    aes(label = name_2, colour = island),
    size               = 5,
    fontface           = "bold",
    fill               = alpha("white", 0.8),
    nudge_x            = 0,
    nudge_y            = 0.2,
    box.padding        = 4,
    point.padding      = 2,
    segment.linewidth  = 0.5,
    segment.curvature  = -0.1,
    min.segment.length = 0,
    max.overlaps       = Inf,
    show.legend        = FALSE
  ) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  labs(x = "Population-weighted density", y = "FOI per serotype") +
  guides(fill = "none") +
  theme_minimal(base_size = 18) +
  theme(
    panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.8),
    panel.background  = element_rect(fill = "white", colour = NA),
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.text       = element_text(size = 20),
    legend.title      = element_text(size = 22),
    axis.title.x      = element_text(size = 26, margin = margin(t = 30)),
    axis.title.y      = element_text(size = 26, margin = margin(r = 30)),
    axis.text.x       = element_text(size = 22, margin = margin(t = 8)),
    axis.text.y       = element_text(size = 22, margin = margin(r = 8)),
    axis.ticks        = element_line(colour = "grey40", linewidth = 0.5),
    axis.ticks.length = unit(4, "pt"),
    plot.margin       = margin(50, 20, 2, 20)
  )

ggsave(
  here(figure_path, "scatter_top100_highlighted_cities.png"),
  p_rank_highlighted,
  width = 12, height = 8, dpi = 800
)


##### Scatter without labels #####
p_rank_no_labels <- ggplot(top100_labelled, aes(x = pwd, y = foi_plot)) +
  geom_point(aes(colour = island), alpha = 0.5, size = 2, shape = 16) +
  geom_point(
    data  = top100_labelled |> filter(top100 & !highlight),
    aes(colour = island),
    size = 4, shape = 21, fill = NA, stroke = 1.5, alpha = 0.8
  ) +
  geom_point(
    data  = top100_labelled |> filter(highlight),
    aes(fill = island),
    colour = "black", size = 5, shape = 23, stroke = 1
  ) +
  scale_colour_manual(
    values = island_colors, name = "Island",
    guide  = guide_legend(nrow = 3, override.aes = list(shape = 15, size = 5, alpha = 1))
  ) +
  scale_fill_manual(values = island_colors, guide = "none") +
  geom_point(
    data = data.frame(
      x    = NA_real_,
      y    = NA_real_,
      type = factor(
        c("All municipalities", "Top 100", "Highlighted cities"),
        levels = c("All municipalities", "Top 100", "Highlighted cities")
      )
    ),
    aes(x = x, y = y, shape = type),
    size = 4, colour = "grey40"
  ) +
  scale_shape_manual(
    name   = NULL,
    values = c("All municipalities" = 16, "Top 100" = 21, "Highlighted cities" = 23),
    guide  = guide_legend(
      override.aes = list(
        colour = c("grey40", "grey40", "black"),
        fill   = c("grey40", NA,       "grey40"),
        size   = c(2,        5,        4),
        stroke = c(0,        1.5,      1),
        alpha  = c(0.5,      0.8,      1)
      )
    )
  ) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  labs(x = "Population-weighted density", y = "FOI per serotype") +
  guides(fill = "none") +
  theme_minimal(base_size = 18) +
  theme(
    panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.8),
    panel.background  = element_rect(fill = "white", colour = NA),
    legend.position   = "bottom",
    legend.box        = "horizontal",
    legend.text       = element_text(size = 20),
    legend.title      = element_text(size = 22),
    axis.title.x      = element_text(size = 26, margin = margin(t = 30)),
    axis.title.y      = element_text(size = 26, margin = margin(r = 30)),
    axis.text.x       = element_text(size = 22, margin = margin(t = 8)),
    axis.text.y       = element_text(size = 22, margin = margin(r = 8)),
    axis.ticks        = element_line(colour = "grey40", linewidth = 0.5),
    axis.ticks.length = unit(4, "pt"),
    plot.margin       = margin(50, 20, 2, 20)
  )

ggsave(
  here(figure_path, "scatter_top100_no_labels.png"),
  p_rank_no_labels,
  width = 12, height = 8, dpi = 800
)


##### Add panel tags #####
p_rank_highlighted <- p_rank_highlighted +
  labs(tag = "a") +
  theme(
    plot.tag = element_text(face = "bold", size = 28),
    plot.tag.position = c(0.02, 1.03)
  )

p_top100_map <- p_top100_map +
  labs(tag = "b") +
  theme(
    plot.tag = element_text(face = "bold", size = 28),
    plot.tag.position = c(0.02, 0.98)
  )


##### Combined with labels — scatter on top, map on bottom #####
p_combined_labels <- (p_rank_highlighted) / (p_top100_map) +
  plot_layout(heights = c(1.2, 1), guides = "keep") &
  theme(
    legend.box    = "horizontal",
    legend.text   = element_text(size = 16),
    legend.title  = element_text(size = 18),
    legend.margin = margin(t = 20),
    plot.margin   = margin(5, 20, 5, 20)
  )

p_combined_labels_padded <- wrap_elements(p_combined_labels) +
  theme(plot.margin = margin(40, 30, 60, 30))

ggsave(
  here(figure_path, "combined_top100_with_labels.png"),
  p_combined_labels_padded,
  width = 16, height = 17, dpi = 1200
)


##### Combined without labels — scatter on top, map on bottom #####
p_combined_no_labels <- (p_rank_no_labels) / (p_top100_map) +
  plot_layout(heights = c(1.2, 1), guides = "keep") &
  theme(
    legend.box    = "horizontal",
    legend.text   = element_text(size = 16),
    legend.title  = element_text(size = 18),
    legend.margin = margin(t = 20),
    plot.margin   = margin(5, 20, 5, 20)
  )

p_combined_no_labels_padded <- wrap_elements(p_combined_no_labels) +
  theme(plot.margin = margin(40, 30, 60, 30))

ggsave(
  here(figure_path, "combined_top100_no_labels.png"),
  p_combined_no_labels_padded,
  width = 16, height = 17, dpi = 800
)

cat("Saved: all four plot variants\n")


##### Population coverage #####
total_pop <- pop_summary |>
  summarise(total = sum(population, na.rm = TRUE)) |>
  pull(total)

pop_coverage_joined <- top100_data |>
  select(gid_2, rank) |>
  left_join(pop_summary, by = c("gid_2" = "adm2_code"))

pop_top100 <- pop_coverage_joined |>
  filter(rank <= 100) |>
  pull(population) |>
  sum(na.rm = TRUE)

cat(sprintf(
  "\nTop 100 municipalities: %.1f%% (%s people)\n",
  pop_top100 / total_pop * 100,
  comma(pop_top100)
))

cumulative_coverage <- top100_data |>
  select(gid_2, rank) |>
  left_join(pop_summary, by = c("gid_2" = "adm2_code")) |>
  arrange(rank) |>
  mutate(
    cumulative_pop = cumsum(population),
    cumulative_pct = cumulative_pop / total_pop * 100
  )

p_cumulative <- ggplot(cumulative_coverage, aes(x = rank, y = cumulative_pct)) +
  geom_line(linewidth = 1, colour = "#2196F3") +
  geom_vline(xintercept = 100, linetype = "dashed", colour = "#FF9F1C", linewidth = 1) +
  geom_point(
    data = cumulative_coverage |> filter(rank == 100),
    colour = "#FF9F1C",
    size = 4,
    shape = 4,
    stroke = 1.5
  ) +
  annotate(
    "text",
    x = 110,
    y = cumulative_coverage$cumulative_pct[cumulative_coverage$rank == 100] - 6,  # more space from line
    label = sprintf(
      "Top 100\n%.1f%%",
      cumulative_coverage$cumulative_pct[cumulative_coverage$rank == 100]
    ),
    colour = "#FF9F1C",
    size = 4.5,
    hjust = -0.1
  ) +
  scale_x_continuous(breaks = seq(0, 500, 50)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, NA)) +
  labs(
    x = "Number of municipalities (ordered by rank)",
    y = "Proportion of population covered (%)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    
    # space between axis titles and axes
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    
    axis.ticks = element_line(colour = "grey30", linewidth = 0.6),   # add ticks
    axis.ticks.length = unit(4, "pt"),
    
    panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    plot.margin = margin(40, 40, 40, 40)
  )

ggsave(
  here(figure_path, "cumulative_population_coverage.png"),
  p_cumulative,
  width = 10, height = 6, dpi = 800
)
