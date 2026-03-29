# code/wol_selection_heatmaps.R

# Produces heatmaps of the metrics used in Wolbachia targeting
# pwd and FOI heatmaps
# scaled log(pwd) and scaled FOI
# joint metric

library(here)
library(dplyr)
library(ggplot2)
library(sf)
library(viridis)
library(scales)
library(qs)
library(readxl)
library(janitor)
library(patchwork)

figure_path <- here("figures", "heatmaps")
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)

# Load municipality shapefile
indo_adm2_shp <- st_read(
  here("data", "geography", "gadm41_IDN_shp", "gadm41_IDN_2.shp"),
  quiet = TRUE
) |>
  clean_names() |>
  mutate(gid_2 = as.character(gid_2))

# Load FOI
indo_foi_full <- read_xlsx(
  here("data", "foi", "imperial-estimates", "Estimates_for_Gavi_20_Aug_25.xlsx"),
  sheet = 3
) |>
  clean_names() |>
  filter(country == "Indonesia") |>
  transmute(
    gid_2 = as.character(id_2),
    foi = total_foi / 4
  ) |>
  distinct(gid_2, .keep_all = TRUE)

# Load municipality-level population-weighted density
pwd_adm2 <- qread(here("data", "pwd_adm2.qs")) |>
  transmute(
    gid_2 = as.character(gid_2),
    pop_weighted_density = as.numeric(pwd)
  ) |>
  distinct(gid_2, .keep_all = TRUE)

# Join data
map_data <- indo_adm2_shp |>
  left_join(indo_foi_full, by = "gid_2") |>
  left_join(pwd_adm2, by = "gid_2")

cat(sprintf("FOI missing for %d municipalities\n", sum(is.na(map_data$foi))))
cat(sprintf("PWD missing for %d municipalities\n", sum(is.na(map_data$pop_weighted_density))))

# Create scaled FOI, pwd and joint metric
map_data <- map_data |>
  mutate(
    foi_scaled = if_else(
      is.na(foi),
      NA_real_,
      (foi - min(foi, na.rm = TRUE)) /
        (max(foi, na.rm = TRUE) - min(foi, na.rm = TRUE))
    ),
    pwd_scaled = if_else(
      is.na(pop_weighted_density),
      NA_real_,
      (log(pop_weighted_density + 0.001) - min(log(pop_weighted_density + 0.001), na.rm = TRUE)) /
        (max(log(pop_weighted_density + 0.001), na.rm = TRUE) - min(log(pop_weighted_density + 0.001), na.rm = TRUE))
    ),
    joint_metric = foi_scaled * pwd_scaled
  )

map_theme <- theme_void(base_size = 20) +
  theme(
    panel.background     = element_rect(fill = "white", colour = NA),
    plot.background      = element_rect(fill = "white", colour = NA),
    panel.border         = element_rect(fill = NA, colour = "black", linewidth = 0.8),
    
    legend.position      = "right",
    legend.direction     = "vertical",
    legend.title         = element_text(size = 17, lineheight = 0.9, hjust = 0),
    legend.text          = element_text(size = 15),
    
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.margin    = margin(0, 0, 0, 0),
    legend.box.spacing   = unit(14, "mm"),
    legend.spacing.y     = unit(2, "pt"),
    legend.justification = c(0, 0.5),
    
    plot.margin          = margin(18, 24, 12, 30)
  )

# Original heatmap function for raw variables
make_heatmap <- function(sf_data, value_col, legend_text, upper_quantile = 0.95) {
  
  vals <- sf_data[[value_col]]
  vals <- vals[is.finite(vals) & !is.na(vals) & vals > 0]
  
  lower_lim <- min(vals)
  upper_lim <- as.numeric(quantile(vals, upper_quantile, na.rm = TRUE))
  
  plot_data <- sf_data |>
    mutate(fill_capped = if_else(
      is.na(.data[[value_col]]) | .data[[value_col]] <= 0,
      NA_real_,
      pmin(.data[[value_col]], upper_lim)
    ))
  
  ggplot(plot_data) +
    geom_sf(aes(fill = fill_capped), colour = "grey90", linewidth = 0.04) +
    scale_fill_viridis_c(
      option    = "viridis",
      direction = -1,
      trans     = "log10",
      limits    = c(lower_lim, upper_lim),
      oob       = squish,
      na.value  = "grey90",
      labels    = comma,
      name      = legend_text,
      guide     = guide_colorbar(
        title.position = "top",
        title.hjust    = 0,
        barwidth       = 0.8,
        barheight      = 7.2
      )
    ) +
    coord_sf(clip = "off") +
    map_theme
}

# New function for scaled variables in [0, 1]
make_scaled_heatmap <- function(sf_data, value_col, legend_text) {
  ggplot(sf_data) +
    geom_sf(aes(fill = .data[[value_col]]), colour = "grey90", linewidth = 0.04) +
    scale_fill_viridis_c(
      option    = "viridis",
      direction = -1,
      limits    = c(0, 1),
      oob       = squish,
      na.value  = "grey90",
      labels    = label_number(accuracy = 0.1),
      name      = legend_text,
      guide     = guide_colorbar(
        title.position = "top",
        title.hjust    = 0,
        barwidth       = 0.8,
        barheight      = 7.2
      )
    ) +
    coord_sf(clip = "off") +
    map_theme +
    theme(
      legend.title = element_text(
        size = 17,
        lineheight = 0.9,
        hjust = 0,
        margin = margin(b = 12)   # increase gap between title and colour bar
      )
    )
}


##### FOI and pwd heatmaps #####

p_foi <- make_heatmap(map_data, "foi", "FOI") +
  labs(tag = "a") +
  theme(
    plot.tag          = element_text(face = "bold", size = 24),
    plot.tag.position = c(-0.04, 1.02),
    plot.margin       = margin(30, 24, 4, 30)
  )

p_pwd <- make_heatmap(map_data, "pop_weighted_density", "PWD") +
  labs(tag = "b") +
  theme(
    plot.tag          = element_text(face = "bold", size = 24),
    plot.tag.position = c(-0.04, 1.02),
    plot.margin       = margin(4, 24, 18, 30)
  )

p_combined <- p_foi / p_pwd +
  plot_layout(
    heights = c(1, 1),
    guides  = "keep"
  ) &
  theme(
    plot.margin = margin(20, 40, 15, 30)
  )

ggsave(
  here(figure_path, "heatmap_foi_pwd_indonesia_combined.png"),
  p_combined,
  width  = 12.6,
  height = 9.2,
  dpi    = 1200
)

##### Scaled maps: a = scaled FOI, b = scaled PWD, c = joint metric #####

p_foi_scaled <- make_scaled_heatmap(map_data, "foi_scaled", "Scaled FOI") +
  labs(tag = "a") +
  theme(
    plot.tag          = element_text(face = "bold", size = 24),
    plot.tag.position = c(-0.04, 1.02),
    plot.margin       = margin(30, 24, 4, 30)
  )

p_pwd_scaled <- make_scaled_heatmap(map_data, "pwd_scaled", "Scaled log(PWD)") +
  labs(tag = "b") +
  theme(
    plot.tag          = element_text(face = "bold", size = 24),
    plot.tag.position = c(-0.04, 1.02),
    plot.margin       = margin(4, 24, 4, 30)
  )

p_joint <- make_scaled_heatmap(map_data, "joint_metric", "Joint metric") +
  labs(tag = "c") +
  theme(
    plot.tag          = element_text(face = "bold", size = 24),
    plot.tag.position = c(-0.04, 1.02),
    plot.margin       = margin(4, 24, 18, 30)
  )

p_scaled_combined <- p_foi_scaled / p_pwd_scaled / p_joint +
  plot_layout(
    heights = c(1, 1, 1),
    guides  = "keep"
  ) &
  theme(
    plot.margin = margin(20, 40, 15, 30)
  )

ggsave(
  here(figure_path, "heatmap_scaled_foi_pwd_joint_indonesia_combined.png"),
  p_scaled_combined,
  width  = 12.6,
  height = 13.2,
  dpi    = 1200
)

#ggsave(
#  here(figure_path, "heatmap_joint_metric_indonesia.png"),
#  p_joint,
#  width  = 12.6,
#  height = 4.4,
#  dpi    = 800
#)
