# code/00_process_worldpop.R

# Pre-processing: compute population-weighted density
# from WorldPop 1km age-sex rasters for Indonesia (2025)
# Run once locally before submitting cluster jobs.
# Output: data/demography/pwd_adm2.qs

library(here)
library(terra)
library(sf)
library(exactextractr)
library(dplyr)
library(janitor)
library(qs)
source(here("code", "00_load-packages.R"))

# Load adm2 shapefile
indo_adm2_shp <- st_read(here("data", "geography", "gadm41_IDN_shp", "gadm41_IDN_2.shp")) |>
  clean_names()

# Load and sum WorldPop 1km age-sex rasters
# Sum across all age-sex groups for total population per 1km cell
worldpop_dir <- here("data", "demography", "idn_agesex_structures_2025_CN_1km_R2025A_UA_v1")
tif_files    <- list.files(worldpop_dir, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("Found %d WorldPop age-sex raster files\n", length(tif_files)))

r_stack  <- rast(tif_files)
r_totpop <- app(r_stack, fun = sum, na.rm = TRUE)
names(r_totpop) <- "population"
rm(r_stack)

# Reproject shapefile to match WorldPop CRS (WGS84)
indo_adm2_proj <- st_transform(indo_adm2_shp, crs(r_totpop))

# Compute per-cell density
r_area    <- cellSize(r_totpop, unit = "km")
r_density <- r_totpop / r_area
names(r_density) <- "density"

# Compute population-weighted mean cell density per municipality
# pwd = sum(pop_i * density_i) / sum(pop_i)
# where i indexes 1km cells within each municipality
cat("Extracting population-weighted density per municipality...\n")
pwd_values <- exact_extract(
  x = c(r_density, r_totpop),
  y = indo_adm2_proj,
  fun = function(df) {
    pop   <- df$population * df$coverage_fraction
    dens  <- df$density
    valid <- !is.na(pop) & !is.na(dens) & pop > 0
    if (sum(valid) == 0) return(NA_real_)
    sum(pop[valid] * dens[valid]) / sum(pop[valid])
  },
  summarize_df = TRUE
)

pwd_df <- data.frame(
  gid_2 = indo_adm2_proj$gid_2,
  pwd   = pwd_values
)

cat(sprintf("Computed pwd for %d municipalities (%d NA)\n",
            nrow(pwd_df), sum(is.na(pwd_df$pwd))))

# Save output
qsave(pwd_df, here("data", "demography", "pwd_adm2.qs"))
cat("Saved to: data/demography/pwd_adm2.qs\n")

