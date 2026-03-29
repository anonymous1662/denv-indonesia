# code/process-demography.R
# Provided by Emilie Finch

library(terra)
library(exactextractr)

# Load in Indonesia shape files for adm1 and adm2
indo_adm1_shp <- st_read(here("data", "geography", "gadm41_IDN_shp", "gadm41_IDN_1.shp")) |> 
  clean_names()
indo_adm2_shp <- st_read(here("data", "geography", "gadm41_IDN_shp", "gadm41_IDN_2.shp")) |> 
  clean_names()

# Process ISIMP population raster under ssp2 climate scenario
# from: https://data.isimip.org/search/tree/ISIMIP2b/SecondaryInputData/socioeconomic/pop/ssp2soc/

r <- rast(here("data", "demography", "population_ssp2soc_2.5min_annual_2006-2100.nc4"))
names(r) <- as.character(2006:2100)
indonesia_ext <- ext(94, 142, -12, 8)
r_idn <- crop(r, indonesia_ext)

indo_adm1_shp <- st_transform(indo_adm1_shp, crs(r))
indo_adm2_shp <- st_transform(indo_adm2_shp, crs(r))
rm(r)

# Admin 2
adm2_pop <- exact_extract(r_idn, indo_adm2_shp, 'sum')
adm2_pop$adm2_code <- indo_adm2_shp$gid_2
adm2_pop$adm2_name <- indo_adm2_shp$name_2

adm2_pop <- adm2_pop |> 
  pivot_longer(-c(adm2_name, adm2_code), names_to = "year", values_to = "population") |>
  mutate(year = as.integer(str_remove(year, "sum\\.")))

# Admin 1
adm1_pop <- exact_extract(r_idn, indo_adm1_shp, 'sum')
adm1_pop$adm1_code <- indo_adm1_shp$gid_1
adm1_pop$adm1_name <- indo_adm1_shp$name_1

adm1_pop <- adm1_pop |> 
  pivot_longer(-c(adm1_name, adm1_code), names_to = "year", values_to = "population") |>
  mutate(year = as.integer(str_remove(year, "sum\\.")))


# Add national-level age structure projections

age_str <- read.csv(here("data", "demography", "PROJresult_AGE_SSP2_V13.csv")) |> 
  filter(region == "reg360") |> # filter for Indonesia
  rename(year = Time) |> 
  group_by(year, agest) |> 
  summarise(pop = sum(pop)) |> 
  mutate(age_lower = case_when(agest %in% c(-5, 0) ~ "0",
                               agest > 90 ~ "90", 
                               T ~ as.character(agest))) |> 
  group_by(age_lower, year) |>  # add any other grouping vars e.g. province
  summarise(pop = sum(pop, na.rm = TRUE), .groups = "drop") |> 
  mutate(age_lower = as.integer(age_lower)) |> 
  group_by(age_lower) |> 
  complete(year = seq(min(year), max(year), by = 1)) |>
  mutate(pop = approx(year[!is.na(pop)], # linearly interpolate between time points
                             pop[!is.na(pop)],
                             xout = year)$y) |> 
  group_by(year) |> 
  mutate(age_weight = pop/sum(pop)) |> 
  ungroup()

demog_adm2 <- expand_grid(year = unique(age_str$year),
                     adm2_code = unique(adm2_pop$adm2_code),
                     age_lower = unique(age_str$age_lower)) |> 
  left_join(adm2_pop |> rename(total_population = population)) |> 
  left_join(age_str |> select(age_lower, year, age_weight)) |> 
  mutate(population = total_population * age_weight) |> 
  select(adm2_code, year, age_lower, population) 


demog_adm1 <- expand_grid(year = unique(age_str$year),
                          adm1_code = unique(adm1_pop$adm1_code),
                          age_lower = unique(age_str$age_lower)) |> 
  left_join(adm1_pop |> rename(total_population = population)) |> 
  left_join(age_str |> select(age_lower, year, age_weight)) |> 
  mutate(population = total_population * age_weight) |> 
  select(adm1_code, year, age_lower, population) 

qsave(demog_adm1, here("data", "demography", "demog_adm1.qs"))
qsave(demog_adm2, here("data", "demography", "demog_adm2.qs"))
