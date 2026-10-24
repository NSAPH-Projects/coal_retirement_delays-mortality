# ============================================================
# 02_fit_model_and_predict_future_pm25.R
# Fit historical SO2-PM2.5 model and predict future scenario PM2.5
# ============================================================

source("R/00_config.R")
source("R/00_functions.R")

load_required_packages(c(
  "fst",
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "lubridate",
  "fixest",
  "tibble",
  "ggplot2",
  "scales",
  "sf",
  "patchwork"
))

# ------------------------------------------------------------
# 1. Read historical plant-ZCTA PM2.5 panel
# ------------------------------------------------------------

pm25_panel <- read_fst(target_plant_zcta_pm25_fp) |>
  mutate(
    year = as.integer(year),
    plant_id = as.integer(plant_id),
    plant_name = as.character(plant_name),
    zcta = standardize_zcta(zcta),
    coal_pm25_plant = as.numeric(coal_pm25_plant)
  )

# ------------------------------------------------------------
# 2. Read plant-year SO2 emissions
# ------------------------------------------------------------

so2_raw <- read_csv(so2_fp, show_col_types = FALSE) |>
  clean_names_simple()

so2_col <- names(so2_raw)[
  str_detect(
    names(so2_raw),
    "^so2_mass|so2.*mass.*short.*tons|so2.*tons"
  )
][1]

if (is.na(so2_col)) {
  stop("No SO2 emissions column found in so2_raw.")
}

if (!"facility_name" %in% names(so2_raw)) {
  so2_raw$facility_name <- NA_character_
}

so2_panel <- so2_raw |>
  transmute(
    plant_id = as.integer(facility_id),
    year = as.integer(year),
    plant_name_emission = as.character(facility_name),
    so2_tons = as.numeric(.data[[so2_col]])
  )

# ------------------------------------------------------------
# 3. Fit log-log model with plant-ZCTA fixed effects
# ------------------------------------------------------------

model_df <- pm25_panel |>
  left_join(so2_panel, by = c("plant_id", "year")) |>
  mutate(
    plant_zcta_id = paste(plant_id, zcta, sep = "_")
  )

reg_df_log <- model_df |>
  filter(
    !is.na(coal_pm25_plant),
    !is.na(so2_tons),
    coal_pm25_plant > 0,
    so2_tons > 0
  ) |>
  mutate(
    log_pm25 = log(coal_pm25_plant),
    log_so2 = log(so2_tons)
  )

m_log <- feols(
  log_pm25 ~ log_so2 | plant_zcta_id,
  data = reg_df_log,
  cluster = ~ plant_id
)

print(summary(m_log))
print(confint(m_log, parm = "log_so2", level = 0.95))

# Retransformation (smearing) correction. exp(fitted) is the conditional
# median, not the mean. Duan smearing corrects this with mean(exp(residual))
# and makes no distributional assumption. Residual variance differs across
# plants (SD 0.35-0.63), so the factor is computed per plant.
resid_vec <- as.numeric(resid(m_log))
keep_idx <- seq_len(nrow(reg_df_log))
if (!is.null(m_log$obs_selection)) {
  for (sel in m_log$obs_selection) keep_idx <- keep_idx[sel]
}
stopifnot(length(keep_idx) == length(resid_vec))

smear_by_plant <- tibble(
  plant_id  = reg_df_log$plant_id[keep_idx],
  resid_log = resid_vec
) |>
  group_by(plant_id) |>
  summarise(smear_duan = mean(exp(resid_log)), .groups = "drop")

write_csv(smear_by_plant, smearing_fp)
print(range(smear_by_plant$smear_duan))
print(nrow(smear_by_plant))

message(sprintf(
  "Duan smearing by plant: %.3f - %.3f (pooled %.5f)",
  min(smear_by_plant$smear_duan), max(smear_by_plant$smear_duan),
  mean(exp(resid_vec))
))

gamma_hat <- as.numeric(coef(m_log)["log_so2"])
gamma_ci <- confint(m_log, parm = "log_so2", level = 0.95)

alpha_table <- tibble(
  plant_zcta_id = names(fixef(m_log)$plant_zcta_id),
  alpha_hat = as.numeric(fixef(m_log)$plant_zcta_id)
) |>
  mutate(
    plant_id = as.integer(str_extract(plant_zcta_id, "^[0-9]+")),
    zcta = standardize_zcta(str_extract(plant_zcta_id, "(?<=_).*"))
  ) |>
  select(
    plant_id,
    zcta,
    plant_zcta_id,
    alpha_hat
  )

write_csv(
  tibble(
    term = "log_so2",
    estimate = gamma_hat,
    conf_low_95 = gamma_ci[1, 1],
    conf_high_95 = gamma_ci[1, 2],
    n_obs = nobs(m_log),
    n_plants = n_distinct(reg_df_log$plant_id),
    n_plant_zcta_fe = n_distinct(reg_df_log$plant_zcta_id)
  ),
  model_gamma_fp
)

write_fst(alpha_table, model_alpha_fp, compress = 80)

# ------------------------------------------------------------
# 4. Build future SO2 scenarios from recent SO2 and retirement dates
# ------------------------------------------------------------

gen_tbl <- read_csv(scenario_fp, show_col_types = FALSE) |>
  clean_names_simple() |>
  mutate(
    plant_id = as.integer(plant_id),
    generator_id = as.character(generator_id),
    original_retirement_date = as.Date(original_retirement_date),
    modified_retirement_date = as.Date(modified_retirement_date),
    indefinitely_delayed = as.integer(indefinitely_delayed),
    across(ends_with("_2023_2025"), as.numeric)
  )

if (!setequal(unique(gen_tbl$plant_id), target_plants$plant_id)) {
  stop("Generator scenario file does not cover exactly the 20 target plants.")
}

# Plant baseline emissions are the sum over coal generators. Unlike the
# plant-level CAMPD total this excludes gas-fired units and units that
# retired before 2023, so the baseline matches the generator set that the
# future scenarios operate on.
so2_recent <- gen_tbl |>
  group_by(plant_id) |>
  summarise(so2_mean_baseline = sum(so2_mean_2023_2025), .groups = "drop")

# EIA-860M reports retirement dates per generator. Collapsing them to a single
# plant date either erases a single-unit delay or applies it to the whole
# plant's emissions. Instead each generator gets its own operating fraction and
# these are averaged with unit-level SO2 as the weight. Emissions, not
# capacity, because scrubber status and utilisation differ within a plant:
# Sherburne County unit 1 is 44% of plant capacity but 5% of plant SO2.
# A missing date means the generator keeps operating through 2035.
future_so2 <- expand_grid(
  gen_tbl,
  year = years_future
) |>
  rowwise() |>
  mutate(
    frac_original = annual_operating_fraction(
      year = year,
      retirement_date = original_retirement_date,
      indefinitely_delayed = as.integer(is.na(original_retirement_date))
    ),
    frac_modified = annual_operating_fraction(
      year = year,
      retirement_date = modified_retirement_date,
      indefinitely_delayed = indefinitely_delayed
    )
  ) |>
  ungroup() |>
  group_by(plant_id, plant_name, year) |>
  summarise(
    operating_fraction_original =
      sum(so2_mean_2023_2025 * frac_original) / sum(so2_mean_2023_2025),
    operating_fraction_modified =
      sum(so2_mean_2023_2025 * frac_modified) / sum(so2_mean_2023_2025),
    .groups = "drop"
  ) |>
  left_join(so2_recent, by = "plant_id") |>
  mutate(
    so2_original = so2_mean_baseline * operating_fraction_original,
    so2_modified = so2_mean_baseline * operating_fraction_modified
  )

write_csv(future_so2, future_so2_fp)

future_so2_long <- future_so2 |>
  select(
    plant_id,
    plant_name,
    year,
    so2_original,
    so2_modified
  ) |>
  pivot_longer(
    cols = c(so2_original, so2_modified),
    names_to = "scenario",
    values_to = "so2_future"
  ) |>
  mutate(
    scenario = recode(
      scenario,
      so2_original = "original",
      so2_modified = "modified"
    )
  )

# ------------------------------------------------------------
# 5. Predict future plant-ZCTA coal PM2.5
# ------------------------------------------------------------

future_coal_pm25_point <- alpha_table |>
  inner_join(future_so2_long, by = "plant_id") |>
  left_join(smear_by_plant, by = "plant_id") |>
  mutate(
    gamma_hat = gamma_hat,
    smear_duan = tidyr::replace_na(smear_duan, 1),
    pred_log_coal_pm25 = if_else(
      !is.na(alpha_hat) & !is.na(so2_future) & so2_future > 0,
      alpha_hat + gamma_hat * log(so2_future),
      NA_real_
    ),
    # Smearing multiplies the concentration and then passes through the CRF;
    # it is not applied to the death count.
    pred_coal_pm25 = if_else(
      !is.na(pred_log_coal_pm25),
      exp(pred_log_coal_pm25) * smear_duan,
      0
    )
  ) |>
  select(
    scenario,
    year,
    plant_id,
    plant_name,
    zcta,
    so2_future,
    gamma_hat,
    alpha_hat,
    smear_duan,
    pred_log_coal_pm25,
    pred_coal_pm25
  ) |>
  arrange(
    scenario,
    year,
    plant_id,
    zcta
  )

write_fst(future_coal_pm25_point, future_pm25_fp, compress = 80)

# ------------------------------------------------------------
# 6. Aggregate to ZCTA-year scenario totals and modified-minus-original difference
# ------------------------------------------------------------

future_zcta_pm25_by_scenario <- future_coal_pm25_point |>
  group_by(scenario, year, zcta) |>
  summarise(
    coal_pm25_zcta_total = sum(pred_coal_pm25, na.rm = TRUE),
    .groups = "drop"
  )

future_zcta_pm25_diff <- future_zcta_pm25_by_scenario |>
  pivot_wider(
    names_from = scenario,
    values_from = coal_pm25_zcta_total,
    values_fill = 0
  ) |>
  mutate(
    coal_pm25_original = original,
    coal_pm25_modified = modified,
    additional_coal_pm25 = modified - original
  ) |>
  select(
    year,
    zcta,
    coal_pm25_original,
    coal_pm25_modified,
    additional_coal_pm25
  ) |>
  arrange(year, zcta)

write_csv(future_zcta_pm25_diff, future_zcta_pm25_diff_fp)

zcta_metrics <- future_zcta_pm25_diff |>
  group_by(zcta) |>
  summarise(
    cumulative_original = sum(coal_pm25_original, na.rm = TRUE),
    cumulative_modified = sum(coal_pm25_modified, na.rm = TRUE),
    cumulative_additional = sum(additional_coal_pm25, na.rm = TRUE),
    mean_annual_additional = mean(additional_coal_pm25, na.rm = TRUE),
    max_annual_additional = max(additional_coal_pm25, na.rm = TRUE),
    n_years_positive = sum(additional_coal_pm25 > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    avoidable_share_modified = case_when(
      cumulative_modified > eps ~ cumulative_additional / cumulative_modified,
      cumulative_modified <= eps & abs(cumulative_additional) <= eps ~ 0,
      TRUE ~ NA_real_
    ),
    avoidable_share_modified_pct = 100 * avoidable_share_modified
  )

write_csv(zcta_metrics, zcta_pm25_metrics_fp)

future_so2 |>
  group_by(year) |>
  summarise(
    so2_original_total = sum(so2_original, na.rm = TRUE),
    so2_modified_total = sum(so2_modified, na.rm = TRUE),
    additional_so2 = so2_modified_total - so2_original_total,
    .groups = "drop"
  ) |>
  write_csv(file.path(dir_diagnostics, "future_so2_yearly_totals.csv"))

# ------------------------------------------------------------
# 7. National cumulative additional coal-PM2.5 map
# ------------------------------------------------------------

# --------------------------------------------------
# 7.1 Aggregate annual ZCTA differences to 2026-2035 cumulative exposure
# --------------------------------------------------

zcta_shp <- file.path(
  "tl_2020_us_zcta520",
  "tl_2020_us_zcta520.shp"
)

state_shp <- file.path(
  "tl_2020_us_state",
  "tl_2020_us_state.shp"
)

zcta_cumulative_additional_pm25 <- future_zcta_pm25_diff %>%
  mutate(
    zcta = str_pad(as.character(zcta), width = 5, side = "left", pad = "0"),
    year = as.integer(year)
  ) %>%
  filter(year %in% 2026:2035) %>%
  group_by(zcta) %>%
  summarise(
    cumulative_original_coal_pm25 =
      sum(coal_pm25_original, na.rm = TRUE),
    cumulative_revised_coal_pm25 =
      sum(coal_pm25_modified, na.rm = TRUE),
    cumulative_additional_coal_pm25 =
      sum(additional_coal_pm25, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  zcta_cumulative_additional_pm25,
  file.path(
    dir_tables,
    "zcta_cumulative_additional_coal_pm25_2026_2035.csv"
  )
)

# --------------------------------------------------
# 7.2 Read ZCTA geometry and attach cumulative additional PM2.5
# --------------------------------------------------

zcta_map_additional <- st_read(zcta_shp, quiet = TRUE) %>%
  transmute(
    zcta = str_pad(
      as.character(ZCTA5CE20),
      width = 5,
      side = "left",
      pad = "0"
    ),
    geometry = geometry
  ) %>%
  st_transform(grid_crs) %>%
  left_join(
    zcta_cumulative_additional_pm25,
    by = "zcta"
  )

zcta_geometry <- zcta_map_additional %>%
  select(zcta, geometry)

# --------------------------------------------------
# 7.3 Read CONUS state boundaries
# --------------------------------------------------

state_map_raw <- st_read(state_shp, quiet = TRUE)

if ("STUSPS" %in% names(state_map_raw)) {
  state_map <- state_map_raw %>%
    filter(
      !STUSPS %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS")
    )
} else if ("STATEFP" %in% names(state_map_raw)) {
  state_map <- state_map_raw %>%
    filter(
      !STATEFP %in% c("02", "15", "72", "78", "66", "69", "60")
    )
} else {
  stop("State shapefile must contain either STUSPS or STATEFP.")
}

state_map <- state_map %>%
  st_transform(grid_crs) %>%
  st_make_valid()

conus_bbox <- st_bbox(state_map)

# --------------------------------------------------
# 7.4 Delayed-plant locations
# --------------------------------------------------
# Plant locations are taken from pm25_facility_state.csv.
# This file contains repeated facility-year rows, so we keep one
# unique latitude/longitude record per facility.

facility_location_fp <- file.path(
  dir_raw,
  "pm25_facility_state.csv"
)

if (file.exists(facility_location_fp)) {
  
  facility_locations <- read_csv(
    facility_location_fp,
    show_col_types = FALSE
  ) %>%
    transmute(
      plant_id = as.character(FacID),
      plant_name = as.character(Facility.Name),
      state_facility = as.character(state_facility),
      lon = as.numeric(Facility.Longitude),
      lat = as.numeric(Facility.Latitude)
    ) %>%
    filter(
      !is.na(plant_id),
      !is.na(plant_name),
      !is.na(lon),
      !is.na(lat)
    ) %>%
    distinct(
      plant_id,
      plant_name,
      state_facility,
      lon,
      lat
    )
  
  if (exists("target_plants")) {
    
    target_plants_clean <- target_plants %>%
      mutate(
        plant_id = as.character(plant_id)
      )
    
    plant_points_df <- facility_locations %>%
      semi_join(
        target_plants_clean,
        by = "plant_id"
      )
    
    if (nrow(plant_points_df) == 0 && "plant_id" %in% names(target_plants_clean)) {
      plant_points_df <- facility_locations %>%
        semi_join(
          target_plants_clean %>% distinct(plant_id),
          by = "plant_id"
        )
    }
    
    if (nrow(plant_points_df) == 0 && "plant_name" %in% names(target_plants_clean)) {
      plant_points_df <- facility_locations %>%
        semi_join(
          target_plants_clean %>% distinct(plant_name),
          by = "plant_name"
        )
    }
    
  } else {
    plant_points_df <- facility_locations
  }
  
  plant_points_sf <- plant_points_df %>%
    st_as_sf(
      coords = c("lon", "lat"),
      crs = 4326,
      remove = FALSE
    ) %>%
    st_transform(grid_crs)
  
  print(
    plant_points_df %>%
      select(plant_id, plant_name, state_facility, lon, lat) %>%
      arrange(plant_name)
  )
  
} else {
  plant_points_sf <- NULL
  message("pm25_facility_state.csv not found. Drawing map without plant points.")
}

# --------------------------------------------------
# 7.5 Plot cumulative additional coal-PM2.5 map
# --------------------------------------------------

p_map_cumulative_additional_pm25 <- ggplot() +
  geom_sf(
    data = zcta_map_additional,
    aes(fill = cumulative_additional_coal_pm25),
    color = NA
  ) +
  geom_sf(
    data = state_map,
    fill = NA,
    color = "grey80",
    linewidth = 0.18
  ) +
  {
    if (!is.null(plant_points_sf)) {
      geom_sf(
        data = plant_points_sf,
        shape = 21,
        size = 1.5,
        fill = "#C51B3A",
        color = "white",
        stroke = 0.35
      )
    }
  } +
  scale_fill_gradientn(
    colours = c(
      "#FFFFFF",
      "#FFF7BC",
      "#FEE391",
      "#FEC44F",
      "#FE9929",
      "#EC7014",
      "#CC4C02",
      "#8C2D04",
      "#252525"
    ),
    values = scales::rescale(
      c(0, 0.005, 0.02, 0.04, 0.10, 0.20, 0.40, 0.60, 1.00),
      from = c(0, 1)
    ),
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = label_number(accuracy = 0.01),
    limits = c(0, 1.1),
    oob = scales::squish,
    na.value = "grey92",
    name = expression(
      "Coal PM"[2.5] *
        "(" * mu * "g/m"^3 * "-years)"
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(4.8, "in"),
      barheight = unit(0.16, "in"),
      ticks = TRUE
    )
  ) +
  coord_sf(
    crs = grid_crs,
    xlim = c(conus_bbox["xmin"], conus_bbox["xmax"]),
    ylim = c(conus_bbox["ymin"], conus_bbox["ymax"]),
    datum = NA,
    expand = FALSE
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7.5, hjust = 0.5),
    legend.text = element_text(size = 6),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(8, 8, 8, 8)
  )

# --------------------------------------------------
# 7.6 Save cumulative additional coal-PM2.5 map
# --------------------------------------------------

ggsave(
  file.path(
    dir_figures,
    "map_cumulative_additional_coal_pm25_2026_2035.png"
  ),
  p_map_cumulative_additional_pm25,
  width = 7.25,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 8. Plant-specific exposure footprint maps: top 3 plants
# ------------------------------------------------------------

selected_plant_names <- c(
  "Welsh",
  "Coleto Creek",
  "Limestone"
)

plant_label_lookup_top3 <- c(
  "Welsh" = "Welsh (TX)",
  "Coleto Creek" = "Coleto Creek (TX)",
  "Limestone" = "Limestone (TX)"
)

# --------------------------------------------------
# 8.1 Calculate plant-specific cumulative additional PM2.5
# --------------------------------------------------

plant_zcta_cumulative_top3 <- future_coal_pm25_point %>%
  mutate(
    zcta = str_pad(as.character(zcta), width = 5, side = "left", pad = "0"),
    year = as.integer(year),
    scenario = str_to_lower(as.character(scenario))
  ) %>%
  filter(
    year %in% 2026:2035,
    plant_name %in% selected_plant_names
  ) %>%
  select(
    plant_id,
    plant_name,
    year,
    zcta,
    scenario,
    pred_coal_pm25
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = pred_coal_pm25,
    values_fill = 0
  ) %>%
  mutate(
    coal_pm25_original = original,
    coal_pm25_modified = modified,
    additional_coal_pm25 = modified - original
  ) %>%
  group_by(plant_id, plant_name, zcta) %>%
  summarise(
    cumulative_original_coal_pm25 = sum(coal_pm25_original, na.rm = TRUE),
    cumulative_revised_coal_pm25 = sum(coal_pm25_modified, na.rm = TRUE),
    cumulative_additional_coal_pm25 = sum(additional_coal_pm25, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    plant_label = recode(
      plant_name,
      !!!plant_label_lookup_top3
    )
  )

write_csv(
  plant_zcta_cumulative_top3,
  file.path(
    dir_tables,
    "plant_specific_cumulative_additional_coal_pm25_top3_2026_2035.csv"
  )
)

# --------------------------------------------------
# 8.2 Attach ZCTA geometry
# --------------------------------------------------

zcta_map_top3 <- zcta_geometry %>%
  inner_join(
    plant_zcta_cumulative_top3,
    by = "zcta"
  ) %>%
  filter(
    !is.na(cumulative_additional_coal_pm25),
    cumulative_additional_coal_pm25 > 0
  )

# --------------------------------------------------
# 8.2 Attach ZCTA geometry with small-gap interpolation
# --------------------------------------------------
# This interpolation is used only for visualization.
# It fills small internal gaps by assigning the value of the nearest
# observed positive-exposure ZCTA within max_interp_dist_km_top3.

max_interp_dist_km_top3 <- 75

plant_zcta_cumulative_top3_obs <- plant_zcta_cumulative_top3 %>%
  filter(
    !is.na(cumulative_additional_coal_pm25),
    cumulative_additional_coal_pm25 > 0
  ) %>%
  select(
    plant_id,
    plant_name,
    plant_label,
    zcta,
    cumulative_original_coal_pm25,
    cumulative_revised_coal_pm25,
    cumulative_additional_coal_pm25
  ) %>%
  mutate(
    interpolated = FALSE
  )

zcta_centroids_top3 <- zcta_geometry %>%
  st_transform(grid_crs) %>%
  st_point_on_surface() %>%
  select(zcta, geometry)

fill_small_gaps_one_top3_plant <- function(plant_nm) {
  
  obs <- plant_zcta_cumulative_top3_obs %>%
    filter(plant_name == plant_nm)
  
  if (nrow(obs) == 0) {
    return(NULL)
  }
  
  obs_pts <- zcta_centroids_top3 %>%
    inner_join(
      obs %>%
        select(
          plant_id,
          plant_name,
          plant_label,
          zcta,
          cumulative_original_coal_pm25,
          cumulative_revised_coal_pm25,
          cumulative_additional_coal_pm25
        ),
      by = "zcta"
    )
  
  if (nrow(obs_pts) == 0) {
    return(NULL)
  }
  
  missing_pts <- zcta_centroids_top3 %>%
    anti_join(
      obs %>% select(zcta),
      by = "zcta"
    )
  
  nearest_idx <- st_nearest_feature(
    missing_pts,
    obs_pts
  )
  
  nearest_dist_km <- as.numeric(
    st_distance(
      missing_pts,
      obs_pts[nearest_idx, ],
      by_element = TRUE
    )
  ) / 1000
  
  fill_tbl <- missing_pts %>%
    st_drop_geometry() %>%
    mutate(
      nearest_dist_km = nearest_dist_km,
      cumulative_original_coal_pm25 =
        obs_pts$cumulative_original_coal_pm25[nearest_idx],
      cumulative_revised_coal_pm25 =
        obs_pts$cumulative_revised_coal_pm25[nearest_idx],
      cumulative_additional_coal_pm25 =
        obs_pts$cumulative_additional_coal_pm25[nearest_idx]
    ) %>%
    filter(nearest_dist_km <= max_interp_dist_km_top3) %>%
    transmute(
      plant_id = obs$plant_id[1],
      plant_name = obs$plant_name[1],
      plant_label = obs$plant_label[1],
      zcta,
      cumulative_original_coal_pm25,
      cumulative_revised_coal_pm25,
      cumulative_additional_coal_pm25,
      interpolated = TRUE
    )
  
  fill_tbl
}

plant_zcta_cumulative_top3_interp <- bind_rows(
  plant_zcta_cumulative_top3_obs,
  bind_rows(
    lapply(
      unique(plant_zcta_cumulative_top3_obs$plant_name),
      fill_small_gaps_one_top3_plant
    )
  )
)

zcta_map_top3 <- zcta_geometry %>%
  inner_join(
    plant_zcta_cumulative_top3_interp,
    by = "zcta"
  )

# --------------------------------------------------
# 8.3 Plant points for selected plants
# --------------------------------------------------

plant_ids_top3 <- plant_zcta_cumulative_top3 %>%
  distinct(plant_id, plant_label) %>%
  mutate(plant_id = as.character(plant_id))

plant_points_top3_sf <- plant_points_sf %>%
  mutate(plant_id = as.character(plant_id)) %>%
  inner_join(
    plant_ids_top3,
    by = "plant_id"
  )

# --------------------------------------------------
# 8.4 Map context: Texas and surrounding states
# --------------------------------------------------

if ("STUSPS" %in% names(state_map)) {
  state_map_context <- state_map %>%
    filter(STUSPS %in% c("TX", "OK", "AR", "LA", "NM", "KS", "MO", "MS"))
} else {
  state_map_context <- state_map
}

context_bbox <- st_bbox(state_map_context)

fill_max <- quantile(
  zcta_map_top3$cumulative_additional_coal_pm25,
  probs = 0.999,
  na.rm = TRUE
)

# --------------------------------------------------
# 8.5 Plot plant-specific exposure footprints
# --------------------------------------------------

p_map_top3_plant_footprints <- ggplot() +
  geom_sf(
    data = zcta_map_top3,
    aes(fill = cumulative_additional_coal_pm25),
    color = NA
  ) +
  geom_sf(
    data = state_map_context,
    fill = NA,
    color = "grey75",
    linewidth = 0.20
  ) +
  geom_sf(
    data = plant_points_top3_sf,
    shape = 21,
    size = 1.5,
    fill = "#C51B3A",
    color = "white",
    stroke = 0.35
  ) +
  facet_wrap(~ plant_label, nrow = 1) +
  scale_fill_gradientn(
    colours = c(
      "#FFFDE7",
      "#FFF7BC",
      "#FEE391",
      "#FEC44F",
      "#FE9929",
      "#EC7014",
      "#CC4C02",
      "#8C2D04",
      "#252525"
    ),
    values = scales::rescale(c(
      0.000,
      0.005,
      0.015,
      0.030,
      0.045,
      0.060,
      0.075,
      0.090,
      fill_max
    )),
    limits = c(0, fill_max),
    oob = scales::squish,
    na.value = "grey95",
    labels = label_number(accuracy = 0.001),
    name = expression(
      "Coal-" * PM[2.5] *
        " (" * mu * "g/m"^3 * "-years)"
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(4.8, "in"),
      barheight = unit(0.1, "in"),
      ticks = TRUE
    )
  ) +
  coord_sf(
    crs = grid_crs,
    xlim = c(context_bbox["xmin"], context_bbox["xmax"]),
    ylim = c(context_bbox["ymin"], context_bbox["ymax"]),
    datum = NA,
    expand = FALSE
  ) +
  theme_void(base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 8,
      face = "bold",
      margin = margin(b = 10)
    ),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7.5, hjust = 0.5),
    legend.text = element_text(size = 6),
    panel.spacing = unit(0.4, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

# --------------------------------------------------
# 8.6 Save top-3 plant-specific exposure footprint map
# --------------------------------------------------

ggsave(
  file.path(
    dir_figures,
    "map_top3_plant_specific_additional_coal_pm25_2026_2035.png"
  ),
  p_map_top3_plant_footprints,
  width = 7.25,
  height = 3.2,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 9. Plant-specific exposure footprint maps: all 20 plants
# ------------------------------------------------------------
# ------------------------------------------------------------
# 9. Plant-specific exposure footprint maps: all 20 plants
# ------------------------------------------------------------
# Supplementary figure.
# All panels use the same CONUS map extent.
# Colors use a common absolute scale from 0 to 0.250 ug/m3-years.
# Small internal gaps are filled for visualization using nearest-neighbor
# interpolation within 75 km.

# --------------------------------------------------
# 9.1 Calculate plant-specific cumulative additional PM2.5
# --------------------------------------------------

plant_zcta_cumulative_all20 <- future_coal_pm25_point %>%
  mutate(
    zcta = str_pad(as.character(zcta), width = 5, side = "left", pad = "0"),
    year = as.integer(year),
    scenario = str_to_lower(as.character(scenario))
  ) %>%
  filter(year %in% 2026:2035) %>%
  select(
    plant_id,
    plant_name,
    year,
    zcta,
    scenario,
    pred_coal_pm25
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = pred_coal_pm25,
    values_fill = 0
  ) %>%
  mutate(
    coal_pm25_original = original,
    coal_pm25_modified = modified,
    additional_coal_pm25 = modified - original
  ) %>%
  group_by(plant_id, plant_name, zcta) %>%
  summarise(
    cumulative_original_coal_pm25 = sum(coal_pm25_original, na.rm = TRUE),
    cumulative_revised_coal_pm25 = sum(coal_pm25_modified, na.rm = TRUE),
    cumulative_additional_coal_pm25 = sum(additional_coal_pm25, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  plant_zcta_cumulative_all20,
  file.path(
    dir_tables,
    "plant_specific_cumulative_additional_coal_pm25_all20_2026_2035.csv"
  )
)

# --------------------------------------------------
# 9.2 Plant labels and manual panel order
# --------------------------------------------------

plant_label_lookup_all20 <- c(
  "Welsh" = "Welsh (TX)",
  "Coleto Creek" = "Coleto Creek (TX)",
  "Limestone" = "Limestone (TX)",
  "Dave Johnston" = "Dave Johnston (WY)",
  "North Valmy" = "North Valmy (NV)",
  "Stanton Energy Center" = "Stanton Energy Center (FL)",
  "Baldwin Energy Complex" = "Baldwin Energy Complex (IL)",
  "Naughton" = "Naughton (WY)",
  "Columbia (WI)" = "Columbia (WI)",
  "Transalta Centralia Generation" = "Centralia (WA)",
  "Winyah" = "Winyah (SC)",
  "Edgewater" = "Edgewater (WI)",
  "Comanche (CO)" = "Comanche (CO)",
  "TES Filer City Station" = "TES Filer City Station (MI)",
  "J H Campbell" = "J H Campbell (MI)",
  "F B Culley" = "F B Culley (IN)",
  "R M Schahfer" = "R M Schahfer (IN)",
  "South Oak Creek" = "South Oak Creek (WI)",
  "Sherburne County" = "Sherburne County (MN)",
  "Craig (CO)" = "Craig (CO)"
)

plant_order_all20 <- c(
  "Welsh (TX)",
  "Coleto Creek (TX)",
  "Limestone (TX)",
  "Dave Johnston (WY)",
  "North Valmy (NV)",
  "Stanton Energy Center (FL)",
  "Baldwin Energy Complex (IL)",
  "Naughton (WY)",
  "Columbia (WI)",
  "Centralia (WA)",
  "Winyah (SC)",
  "Edgewater (WI)",
  "Comanche (CO)",
  "TES Filer City Station (MI)",
  "J H Campbell (MI)",
  "F B Culley (IN)",
  "R M Schahfer (IN)",
  "South Oak Creek (WI)",
  "Sherburne County (MN)",
  "Craig (CO)"
)

plant_label_lookup_all20_tbl <- tibble(
  plant_name = names(plant_label_lookup_all20),
  plant_label = unname(plant_label_lookup_all20)
)

plant_order_tbl_all20 <- tibble(
  plant_label = plant_order_all20,
  panel_order = seq_along(plant_order_all20)
)

# --------------------------------------------------
# 9.3 Observed plant-ZCTA data for plotting
# --------------------------------------------------

plant_zcta_plot_all20_obs <- plant_zcta_cumulative_all20 %>%
  filter(cumulative_additional_coal_pm25 > 0) %>%
  select(
    plant_id,
    plant_name,
    zcta,
    cumulative_additional_coal_pm25
  ) %>%
  left_join(
    plant_label_lookup_all20_tbl,
    by = "plant_name"
  ) %>%
  mutate(
    plant_label = coalesce(plant_label, plant_name),
    interpolated = FALSE
  )

# Check label consistency before plotting.
missing_from_manual_order <- plant_zcta_plot_all20_obs %>%
  distinct(plant_label) %>%
  anti_join(
    plant_order_tbl_all20,
    by = "plant_label"
  )

if (nrow(missing_from_manual_order) > 0) {
  print(missing_from_manual_order)
  stop("Some plant labels are not included in plant_order_all20.")
}

# --------------------------------------------------
# 9.4 Fill small internal gaps for visualization
# --------------------------------------------------

max_interp_dist_km <- 75

zcta_centroids_all20 <- zcta_geometry %>%
  st_transform(grid_crs) %>%
  st_point_on_surface() %>%
  select(zcta, geometry)

fill_small_gaps_one_plant <- function(plant_nm) {
  
  obs <- plant_zcta_plot_all20_obs %>%
    filter(plant_name == plant_nm)
  
  if (nrow(obs) == 0) {
    return(NULL)
  }
  
  obs_pts <- zcta_centroids_all20 %>%
    inner_join(
      obs %>%
        select(
          plant_id,
          plant_name,
          plant_label,
          zcta,
          cumulative_additional_coal_pm25
        ),
      by = "zcta"
    )
  
  if (nrow(obs_pts) == 0) {
    return(NULL)
  }
  
  missing_pts <- zcta_centroids_all20 %>%
    anti_join(
      obs %>% select(zcta),
      by = "zcta"
    )
  
  nearest_idx <- st_nearest_feature(
    missing_pts,
    obs_pts
  )
  
  nearest_dist_km <- as.numeric(
    st_distance(
      missing_pts,
      obs_pts[nearest_idx, ],
      by_element = TRUE
    )
  ) / 1000
  
  nearest_values <- obs_pts$cumulative_additional_coal_pm25[nearest_idx]
  
  fill_tbl <- missing_pts %>%
    st_drop_geometry() %>%
    mutate(
      nearest_dist_km = nearest_dist_km,
      cumulative_additional_coal_pm25 = nearest_values
    ) %>%
    filter(nearest_dist_km <= max_interp_dist_km) %>%
    transmute(
      plant_id = obs$plant_id[1],
      plant_name = obs$plant_name[1],
      plant_label = obs$plant_label[1],
      zcta,
      cumulative_additional_coal_pm25,
      interpolated = TRUE
    )
  
  fill_tbl
}

plant_zcta_plot_all20_interp <- bind_rows(
  plant_zcta_plot_all20_obs,
  bind_rows(
    lapply(
      unique(plant_zcta_plot_all20_obs$plant_name),
      fill_small_gaps_one_plant
    )
  )
)

# --------------------------------------------------
# 9.5 Common absolute display scale
# --------------------------------------------------

fill_max_abs_all20 <- 0.250

plant_zcta_plot_all20_interp <- plant_zcta_plot_all20_interp %>%
  mutate(
    cumulative_additional_coal_pm25_display =
      pmin(cumulative_additional_coal_pm25, fill_max_abs_all20),
    plant_label_facet = factor(
      plant_label,
      levels = plant_order_all20
    )
  )

# --------------------------------------------------
# 9.6 Simplify geometry for all-20 supplementary map
# --------------------------------------------------

zcta_needed_all20 <- unique(plant_zcta_plot_all20_interp$zcta)

zcta_geometry_all20_simple <- zcta_geometry %>%
  filter(zcta %in% zcta_needed_all20) %>%
  select(zcta, geometry) %>%
  st_transform(grid_crs) %>%
  st_simplify(
    dTolerance = 5000,
    preserveTopology = TRUE
  )

state_map_simple <- state_map %>%
  st_transform(grid_crs) %>%
  st_simplify(
    dTolerance = 10000,
    preserveTopology = TRUE
  )

zcta_map_all20 <- zcta_geometry_all20_simple %>%
  inner_join(
    plant_zcta_plot_all20_interp,
    by = "zcta"
  )

gc()

# --------------------------------------------------
# 9.7 Plant points for all 20 plants
# --------------------------------------------------

if (is.null(plant_points_sf)) {
  stop("plant_points_sf is NULL. Plant locations are required for this figure.")
}

plant_ids_all20 <- plant_zcta_cumulative_all20 %>%
  distinct(plant_id, plant_name) %>%
  mutate(plant_id = as.character(plant_id))

plant_points_all20 <- plant_points_sf %>%
  mutate(plant_id = as.character(plant_id)) %>%
  inner_join(
    plant_ids_all20,
    by = "plant_id",
    suffix = c("_location", "_scenario")
  ) %>%
  mutate(
    plant_name = plant_name_scenario
  ) %>%
  left_join(
    plant_label_lookup_all20_tbl,
    by = "plant_name"
  ) %>%
  mutate(
    plant_label = coalesce(plant_label, plant_name),
    plant_label_facet = factor(
      plant_label,
      levels = plant_order_all20
    )
  )

missing_points_from_manual_order <- plant_points_all20 %>%
  st_drop_geometry() %>%
  distinct(plant_label) %>%
  anti_join(
    plant_order_tbl_all20,
    by = "plant_label"
  )

if (nrow(missing_points_from_manual_order) > 0) {
  print(missing_points_from_manual_order)
  stop("Some plant point labels are not included in plant_order_all20.")
}

# --------------------------------------------------
# 9.8 Plot all 20 plant-specific maps
# --------------------------------------------------

p_map_all20_plant_footprints <- ggplot() +
  geom_sf(
    data = zcta_map_all20,
    aes(fill = cumulative_additional_coal_pm25_display),
    color = NA
  ) +
  geom_sf(
    data = state_map_simple,
    fill = NA,
    color = "grey80",
    linewidth = 0.12
  ) +
  geom_sf(
    data = plant_points_all20,
    shape = 21,
    size = 1.1,
    fill = "#C51B3A",
    color = "white",
    stroke = 0.30
  ) +
  facet_wrap(
    ~ plant_label_facet,
    ncol = 4
  ) +
  scale_fill_gradientn(
    colours = c(
      "white",
      "#FFFDE7",
      "#FFF7BC",
      "#FEE391",
      "#FEC44F",
      "#FE9929",
      "#EC7014",
      "#CC4C02",
      "#8C2D04",
      "#252525"
    ),
    values = scales::rescale(c(
      0.000,
      0.0005,
      0.001,
      0.005,
      0.015,
      0.030,
      0.055,
      0.080,
      0.110,
      0.150
    )),
    limits = c(0, 0.150),
    oob = scales::squish,
    breaks = c(0.0000, 0.0500, 0.1000, 0.1500),
    labels = function(x) sprintf("%.4f", x),
    na.value = "white",
    name = expression(
      "Coal-" * PM[2.5] *
        " (" * mu * "g/m"^3 * "-years)"
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      barwidth = unit(3.6, "in"),
      barheight = unit(0.10, "in"),
      ticks = TRUE
    )
  ) +
  coord_sf(
    crs = grid_crs,
    xlim = c(conus_bbox["xmin"], conus_bbox["xmax"]),
    ylim = c(conus_bbox["ymin"], conus_bbox["ymax"]),
    datum = NA,
    expand = FALSE
  ) +
  theme_void(base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_text(
      size = 7,
      face = "bold",
      margin = margin(b = 4)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7, hjust = 0.5),
    legend.text = element_text(size = 6),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.spacing = unit(0.35, "lines"),
    plot.margin = margin(6, 6, 6, 6)
  )

# --------------------------------------------------
# 9.9 Save 4 x 5 supplementary figure
# --------------------------------------------------

ggsave(
  file.path(
    dir_figures,
    "map_all20_plant_specific_additional_coal_pm25_2026_2035.png"
  ),
  p_map_all20_plant_footprints,
  width = 8.5,
  height = 10,
  dpi = 600,
  bg = "white"
)
