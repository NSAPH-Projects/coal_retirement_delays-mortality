# ============================================================
# Project configuration
# ============================================================

# ----------------------------
# Years
# ----------------------------
years_historical <- 1999:2020
years_future <- 2026:2035
so2_baseline_years <- 2023:2025

# ----------------------------
# Coordinate systems
# ----------------------------
# Grid coordinates are in NAD83 / Conus Albers in the source files.
grid_crs <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# ----------------------------
# Directories
# ----------------------------
dir_raw <- "data"
dir_raw_annual_byunit <- file.path(dir_raw, "annual", "byunit")
dir_raw_annual_total <- file.path(dir_raw, "annual", "total")
dir_intermediate <- "data_intermediate"
dir_output <- "output"
dir_tables <- file.path(dir_output, "tables")
dir_figures <- file.path(dir_output, "figures")
dir_diagnostics <- file.path(dir_output, "diagnostics")

for (d in c(dir_intermediate, dir_output, dir_tables, dir_figures, dir_diagnostics)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ----------------------------
# Raw input files
# ----------------------------
byunit_fst_pattern <- "grids_pm25_byunit_%d.fst"
total_fst_pattern <- "grids_pm25_total_%d.fst"
zcta_shp <- "tl_2020_us_zcta520/tl_2020_us_zcta520.shp"
so2_fp <- file.path(dir_raw, "annual_emission_plant.csv")
scenario_fp <- file.path(dir_raw, "coal_generator_retirement_scenarios.csv")
unit_so2_fp <- file.path(dir_raw, "annual_emission_unit.csv")
demo_fp <- file.path(dir_raw, "demographics.csv")
medicare_fp <- file.path(dir_raw, "medicare_mortality_rate.csv")
plant_location_fp <- file.path(dir_raw, "pm25_facility_state.csv")

# ----------------------------
# Intermediate files
# ----------------------------
grid_zcta_weights_fp <- file.path(dir_intermediate, "grid_to_zcta_weights_area_normalized.csv")
target_plant_zcta_pm25_fp <- file.path(dir_intermediate, "target_plant_zcta_pm25_1999_2020.fst")
model_gamma_fp <- file.path(dir_intermediate, "log_model_gamma_estimate.csv")
smearing_fp <- file.path(dir_intermediate, "smearing_factors_by_plant.csv")
model_alpha_fp <- file.path(dir_intermediate, "log_model_alpha_plant_zcta_fixed_effects.fst")
future_so2_fp <- file.path(dir_intermediate, "future_so2_original_modified_2026_2035.csv")
future_pm25_fp <- file.path(dir_intermediate, "future_plant_zcta_coal_pm25_point_estimate_2026_2035.fst")
future_zcta_pm25_diff_fp <- file.path(dir_intermediate, "future_zcta_coal_pm25_difference_modified_minus_original_2026_2035.csv")
zcta_pm25_metrics_fp <- file.path(dir_intermediate, "zcta_cumulative_additional_and_avoidable_share_2026_2035.csv")


state_shp <- "tl_2020_us_state/tl_2020_us_state.shp"
plant_map_fp <- file.path(dir_raw, "coal_plants_map_plot_df_so2.csv")

# ----------------------------
# Health-impact coefficient
# ----------------------------
# Coal PM2.5 mortality response from the peer-reviewed coal PM2.5 CRF.
# Values are relative risks per 1 ug/m3 coal PM2.5. The model uses log(RR).
beta_tbl <- tibble::tibble(
  beta_type = c("lower", "central", "upper"),
  beta_E = log(c(1.0113, 1.0125, 1.0137))
)

# ----------------------------
# Target delayed-retirement plants
# ----------------------------
target_plants <- tibble::tibble(
  plant_id = c(
    298, 470, 564, 889, 1012, 1710, 3845, 4041, 4050, 4158,
    4162, 6021, 6085, 6090, 6139, 6178, 6249, 8023, 8224, 50835
  ),
  plant_name_target = c(
    "Limestone",
    "Comanche (CO)",
    "Stanton Energy Center",
    "Baldwin Energy Complex",
    "F B Culley",
    "J H Campbell",
    "Transalta Centralia Generation",
    "South Oak Creek",
    "Edgewater",
    "Dave Johnston",
    "Naughton",
    "Craig (CO)",
    "R M Schahfer",
    "Sherburne County",
    "Welsh",
    "Coleto Creek",
    "Winyah",
    "Columbia (WI)",
    "North Valmy",
    "TES Filer City Station"
  )
)

# ----------------------------
# Plot settings
# ----------------------------
color_original <- "#6BAED6"
color_modified <- "#8B1A1A"
color_reference <- "#2C7FB8"
color_subgroup <- "#F28E2B"
