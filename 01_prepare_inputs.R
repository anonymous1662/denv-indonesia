# code/01_prepare_inputs.R

# Prepare inputs for cluster batches (run once)
# Requires pwd_adm2.qs produced locally by 00_process_worldpop.R

library(here)
source(here("code", "00_load-packages.R"))
source(here("code", "utils.R"))
source(here("code", "calculate-R0_fn.R"))
config <- yaml::read_yaml(here("code", "config.yaml"))
prep_path <- here("output", config$run_name, "prep")
dir.create(prep_path, recursive = T)

# Load pre-computed population-weighted density
pwd_df <- qread(here("data", "pwd_adm2.qs"))

# Load FOI
indo_foi_full <- read_xlsx(
  here("data", "Estimates_for_Gavi_20_Aug_25.xlsx"),
  sheet = 3
) |>
  clean_names() |>
  filter(country == "Indonesia") |>
  rename(name_1 = region, name_2 = district, gid_2 = id_2, gid_1 = id_1) |>
  arrange(gid_2)

indo_demog_full <- qread(here("data", "demog_adm2.qs")) |>
  arrange(adm2_code) |>
  filter(adm2_code %in% unique(indo_foi_full$gid_2))

# Expand age groups
indo_demog_full <- expand_age_groups(indo_demog_full, max_age = 90)

# Build Wolbachia ranking (scaled FOI x scaled log(population-weighted density))
pop_summary <- indo_demog_full |>
  filter(year == 2025) |>
  group_by(adm2_code) |>
  summarise(population = sum(population), .groups = "drop")

wol_ranking <- indo_foi_full |>
  select(gid_2, total_foi) |>
  left_join(pop_summary, by = c("gid_2" = "adm2_code")) |>
  left_join(pwd_df,      by = "gid_2") |>
  filter(!is.na(population), !is.na(total_foi), !is.na(pwd)) |>
  mutate(
    foi_scaled = (total_foi - min(total_foi)) / (max(total_foi) - min(total_foi)),
    pwd_scaled = (log(pwd) - min(log(pwd))) / (max(log(pwd)) - min(log(pwd)))
  ) |>
  mutate(foi_pwd = foi_scaled * pwd_scaled) |>
  arrange(desc(foi_pwd)) |>
  mutate(rank = row_number())

top_gids <- wol_ranking |>
  filter(rank <= config$n_wol_regions) |>
  pull(gid_2)

indo_foi_full <- indo_foi_full |>
  mutate(wol_region = as.integer(gid_2 %in% top_gids))

cat(sprintf("Wolbachia targeted in %d of %d municipalities\n",
            sum(indo_foi_full$wol_region), nrow(indo_foi_full)))

# Define all_gids using demography sort order to ensure betas, wol_regions, 
# and demog_input columns are all in the same order throughout the pipeline.
all_gids <- dimnames(wrangle_demography(
  indo_demog_full,
  year_start = config$scenario_start,
  year_end   = config$scenario_start + config$scenario_years,
  pad_left   = 0
))[[2]]

stopifnot(
  "all_gids not fully covered by FOI data" = all(all_gids %in% indo_foi_full$gid_2),
  "FOI data has gids not in demography"    = all(indo_foi_full$gid_2 %in% all_gids),
  "length mismatch"                        = length(all_gids) == length(unique(indo_foi_full$gid_2))
)
cat(sprintf("all_gids defined from demography order: %d municipalities\n", length(all_gids)))

# Calculate R0 + betas for all municipalities
r0s_full   <- numeric(length(all_gids))
betas_full <- numeric(length(all_gids))

cat(sprintf("Calculating R0 for %d municipalities...\n", length(all_gids)))

# Precompute age bounds
age_lower <- sort(unique(indo_demog_full$age))
age_upper <- dplyr::lead(age_lower)
age_upper[length(age_upper)] <- 91

for (i in seq_along(all_gids)) {
  id <- all_gids[i]
  
  foi_value <- indo_foi_full |>
    filter(gid_2 == id) |>
    mutate(foi = total_foi / 4) |>
    pull(foi)
  
  age_proportions <- indo_demog_full |>
    filter(adm2_code == id & year == 2025) |>
    mutate(prop = population / sum(population)) |>
    pull(prop)
  
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

betas_full <- r0s_full * 0.2
cat("R0 calculation complete.\n")

# Wolbachia vector in the same order as all_gids
wol_regions <- indo_foi_full |>
  filter(gid_2 %in% all_gids) |>
  arrange(match(gid_2, all_gids)) |>
  pull(wol_region)

prepared <- list(
  config          = config,
  all_gids        = all_gids,
  indo_demog_full = indo_demog_full,
  indo_foi_full   = indo_foi_full,
  r0s_full        = r0s_full,
  betas_full      = betas_full,
  wol_regions     = wol_regions
)

# Save prepared inputs
qsave(prepared, here(prep_path, "inputs_prepared.qs"))
cat("Saved prepared inputs to:", here(prep_path, "inputs_prepared.qs"), "\n")
