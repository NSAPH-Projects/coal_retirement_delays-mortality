# ============================================================
# 07_make_attributable_share_pm25_map.R
# Fig. S5
#
# Projection period: 2026-2035
# ============================================================

source("R/00_config.R")
source("R/00_functions.R")

load_required_packages(c(
  "dplyr",
  "readr",
  "sf",
  "ggplot2",
  "scales",
  "grid"
))

# ------------------------------------------------------------
# 0. Settings and output paths
# ------------------------------------------------------------

numerical_tolerance <- 1e-10

plant_point_color <- if (exists("color_modified")) {
  color_modified
} else {
  "#C51B3A"
}

attributable_share_table_fp <- file.path(
  dir_tables,
  "zcta_coal_pm25_attributable_share_2026_2035.csv"
)

attributable_share_diagnostics_fp <- file.path(
  dir_diagnostics,
  "fig_s3_attributable_share_diagnostics.csv"
)

fig_s3_png_fp <- file.path(
  dir_figures,
  "fig_s3_attributable_share_coal_pm25_2026_2035.png"
)

fig_s3_pdf_fp <- file.path(
  dir_figures,
  "fig_s3_attributable_share_coal_pm25_2026_2035.pdf"
)

# ------------------------------------------------------------
# 1. Check required inputs
# ------------------------------------------------------------

if (!file.exists(future_zcta_pm25_diff_fp)) {
  stop(
    "Missing future ZCTA PM2.5 file: ",
    future_zcta_pm25_diff_fp,
    "\nRun the future PM2.5 prediction script first."
  )
}

if (!file.exists(zcta_shp)) {
  stop(
    "Could not find the ZCTA shapefile: ",
    zcta_shp
  )
}

state_shp_candidates <- c(
  file.path(
    "tl_2020_us_state",
    "tl_2020_us_state.shp"
  ),
  file.path(
    dir_raw,
    "tl_2020_us_state",
    "tl_2020_us_state.shp"
  )
)

state_shp <- state_shp_candidates[
  file.exists(state_shp_candidates)
][1]

if (is.na(state_shp)) {
  stop(
    "Could not find the state shapefile. Checked:\n",
    paste(
      state_shp_candidates,
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# 2. Read future ZCTA-level coal PM2.5 results
# ------------------------------------------------------------

future_pm25_raw <- read_csv(
  future_zcta_pm25_diff_fp,
  show_col_types = FALSE
) |>
  clean_names_simple()

required_columns <- c(
  "year",
  "zcta",
  "coal_pm25_original",
  "coal_pm25_modified"
)

missing_columns <- setdiff(
  required_columns,
  names(future_pm25_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "The future PM2.5 file is missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

future_pm25 <- future_pm25_raw |>
  transmute(
    year = as.integer(year),
    
    zcta = standardize_zcta(zcta),
    
    coal_pm25_original =
      as.numeric(coal_pm25_original),
    
    coal_pm25_revised =
      as.numeric(coal_pm25_modified)
  ) |>
  filter(
    year %in% years_future,
    !is.na(zcta),
    is.finite(coal_pm25_original),
    is.finite(coal_pm25_revised)
  ) |>
  mutate(
    additional_coal_pm25 =
      coal_pm25_revised -
      coal_pm25_original
  )

if (nrow(future_pm25) == 0) {
  stop(
    "No valid future coal PM2.5 observations remained for ",
    min(years_future),
    "-",
    max(years_future),
    "."
  )
}

# ------------------------------------------------------------
# 3. Check revised-minus-original exposure differences
# ------------------------------------------------------------

material_negative_rows <- future_pm25 |>
  filter(
    additional_coal_pm25 <
      -numerical_tolerance
  )

if (nrow(material_negative_rows) > 0) {
  warning(
    nrow(material_negative_rows),
    " ZCTA-year observations have revised exposure below ",
    "original exposure by more than the numerical tolerance. ",
    "Review the scenario results before interpreting Fig. S5."
  )
}

# Set only negligible floating-point negative values to zero.

future_pm25 <- future_pm25 |>
  mutate(
    additional_coal_pm25 = case_when(
      additional_coal_pm25 < 0 &
        additional_coal_pm25 >=
        -numerical_tolerance ~ 0,
      
      TRUE ~ additional_coal_pm25
    )
  )

# ------------------------------------------------------------
# 4. Calculate cumulative ZCTA-level exposure
# ------------------------------------------------------------

zcta_attributable_share <- future_pm25 |>
  group_by(zcta) |>
  summarise(
    cumulative_original_coal_pm25 =
      sum(
        coal_pm25_original,
        na.rm = TRUE
      ),
    
    cumulative_revised_coal_pm25 =
      sum(
        coal_pm25_revised,
        na.rm = TRUE
      ),
    
    cumulative_additional_coal_pm25 =
      sum(
        additional_coal_pm25,
        na.rm = TRUE
      ),
    
    n_projection_years =
      n_distinct(year),
    
    .groups = "drop"
  ) |>
  mutate(
    attributable_share_raw = if_else(
      is.finite(
        cumulative_revised_coal_pm25
      ) &
        cumulative_revised_coal_pm25 > eps,
      
      cumulative_additional_coal_pm25 /
        cumulative_revised_coal_pm25,
      
      NA_real_
    )
  )

# ------------------------------------------------------------
# 5. Validate and clean attributable shares
# ------------------------------------------------------------

materially_invalid_shares <- zcta_attributable_share |>
  filter(
    is.finite(attributable_share_raw),
    attributable_share_raw <
      -numerical_tolerance |
      attributable_share_raw >
      1 + numerical_tolerance
  )

if (nrow(materially_invalid_shares) > 0) {
  warning(
    nrow(materially_invalid_shares),
    " ZCTAs have attributable shares materially outside 0-1. ",
    "These values will be set to missing in Fig. S5. ",
    "Review revised and original exposure estimates."
  )
}

zcta_attributable_share <- zcta_attributable_share |>
  mutate(
    attributable_share = case_when(
      !is.finite(attributable_share_raw) ~
        NA_real_,
      
      attributable_share_raw >=
        -numerical_tolerance &
        attributable_share_raw < 0 ~
        0,
      
      attributable_share_raw >
        1 &
        attributable_share_raw <=
        1 + numerical_tolerance ~
        1,
      
      attributable_share_raw >= 0 &
        attributable_share_raw <= 1 ~
        attributable_share_raw,
      
      TRUE ~
        NA_real_
    ),
    
    attributable_share_pct =
      100 * attributable_share
  )

# ------------------------------------------------------------
# 6. Save ZCTA-level results
# ------------------------------------------------------------

write_csv(
  zcta_attributable_share,
  attributable_share_table_fp
)

# ------------------------------------------------------------
# 7. Diagnostics
# ------------------------------------------------------------

defined_share_values <- zcta_attributable_share |>
  filter(
    is.finite(attributable_share_pct)
  ) |>
  pull(attributable_share_pct)

if (length(defined_share_values) == 0) {
  stop(
    "No finite attributable-share values were calculated."
  )
}

attributable_share_diagnostics <- tibble(
  n_zcta_in_prediction_data =
    nrow(zcta_attributable_share),
  
  n_revised_exposure_positive =
    sum(
      zcta_attributable_share$
        cumulative_revised_coal_pm25 > eps,
      na.rm = TRUE
    ),
  
  n_revised_exposure_zero =
    sum(
      zcta_attributable_share$
        cumulative_revised_coal_pm25 <= eps,
      na.rm = TRUE
    ),
  
  n_share_defined =
    length(defined_share_values),
  
  n_materially_invalid_share =
    nrow(materially_invalid_shares),
  
  minimum_pct =
    min(
      defined_share_values,
      na.rm = TRUE
    ),
  
  p10_pct =
    as.numeric(
      quantile(
        defined_share_values,
        0.10,
        na.rm = TRUE
      )
    ),
  
  p25_pct =
    as.numeric(
      quantile(
        defined_share_values,
        0.25,
        na.rm = TRUE
      )
    ),
  
  median_pct =
    median(
      defined_share_values,
      na.rm = TRUE
    ),
  
  p75_pct =
    as.numeric(
      quantile(
        defined_share_values,
        0.75,
        na.rm = TRUE
      )
    ),
  
  p90_pct =
    as.numeric(
      quantile(
        defined_share_values,
        0.90,
        na.rm = TRUE
      )
    ),
  
  p95_pct =
    as.numeric(
      quantile(
        defined_share_values,
        0.95,
        na.rm = TRUE
      )
    ),
  
  p99_pct =
    as.numeric(
      quantile(
        defined_share_values,
        0.99,
        na.rm = TRUE
      )
    ),
  
  maximum_pct =
    max(
      defined_share_values,
      na.rm = TRUE
    )
)

print(attributable_share_diagnostics)

write_csv(
  attributable_share_diagnostics,
  attributable_share_diagnostics_fp
)

# ------------------------------------------------------------
# 8. Read and prepare ZCTA geometry
# ------------------------------------------------------------

zcta_geometry <- st_read(
  zcta_shp,
  quiet = TRUE
)

zcta_geometry <- zcta_geometry %>%
  st_simplify(dTolerance = 500, preserveTopology = TRUE)

if (!"ZCTA5CE20" %in% names(zcta_geometry)) {
  stop(
    "The ZCTA shapefile does not contain ZCTA5CE20."
  )
}

zcta_map_attributable <- zcta_geometry |>
  transmute(
    zcta = standardize_zcta(
      ZCTA5CE20
    ),
    geometry = geometry
  ) |>
  st_make_valid() |>
  st_transform(grid_crs) |>
  left_join(
    zcta_attributable_share,
    by = "zcta"
  )

# ------------------------------------------------------------
# 9. Read and prepare CONUS state boundaries
# ------------------------------------------------------------

state_map_raw <- st_read(
  state_shp,
  quiet = TRUE
)

if ("STUSPS" %in% names(state_map_raw)) {
  
  state_map <- state_map_raw |>
    filter(
      !STUSPS %in% c(
        "AK",
        "HI",
        "PR",
        "VI",
        "GU",
        "MP",
        "AS"
      )
    )
  
} else if ("STATEFP" %in% names(state_map_raw)) {
  
  state_map <- state_map_raw |>
    filter(
      !STATEFP %in% c(
        "02",
        "15",
        "72",
        "78",
        "66",
        "69",
        "60"
      )
    )
  
} else {
  
  stop(
    "State shapefile must contain ",
    "either STUSPS or STATEFP."
  )
}

state_map <- state_map |>
  st_make_valid() |>
  st_transform(grid_crs)

conus_bbox <- st_bbox(state_map)

# ------------------------------------------------------------
# 10. Read plant locations
# ------------------------------------------------------------

plant_points_sf <- NULL

if (file.exists(plant_location_fp)) {
  
  facility_locations <- read_csv(
    plant_location_fp,
    show_col_types = FALSE
  ) |>
    transmute(
      plant_id =
        as.integer(FacID),
      
      plant_name =
        as.character(Facility.Name),
      
      longitude =
        as.numeric(Facility.Longitude),
      
      latitude =
        as.numeric(Facility.Latitude)
    ) |>
    filter(
      plant_id %in%
        target_plants$plant_id,
      is.finite(longitude),
      is.finite(latitude)
    ) |>
    distinct(
      plant_id,
      .keep_all = TRUE
    )
  
  if (nrow(facility_locations) > 0) {
    
    plant_points_sf <- facility_locations |>
      st_as_sf(
        coords = c(
          "longitude",
          "latitude"
        ),
        crs = 4326,
        remove = FALSE
      ) |>
      st_transform(grid_crs)
  }
}

# ------------------------------------------------------------
# 11. Define continuous color scale
# ------------------------------------------------------------
#
# The mapped variable remains continuous from 0% to 100%.
# Color anchors are concentrated near 100% because the
# attributable-share distribution is concentrated toward
# the upper end of the range.

attributable_share_breaks <- c(
  0,
  25,
  50,
  75,
  90,
  95,
  99,
  100
)

attributable_share_colors <- c(
  "#FFFFCC",
  "#FFEDA0",
  "#FED976",
  "#FEB24C",
  "#FD8D3C",
  "#FC4E2A",
  "#BD0026",
  "#67000D"
)

# ------------------------------------------------------------
# 12. Create Fig. S5
# ------------------------------------------------------------
#
# Colored ZCTAs:
#   share of revised-scenario coal PM2.5 exposure attributable
#   to retirement-schedule changes.
#
# Gray ZCTAs:
#   no modeled future value, zero cumulative revised exposure,
#   or an attributable share outside the valid 0-1 range.
#
# White areas:
#   outside Census ZCTA coverage.

p_fig_s3 <- ggplot() +
  
  geom_sf(
    data = zcta_map_attributable,
    aes(
      fill = attributable_share_pct
    ),
    color = NA
  ) +
  
  geom_sf(
    data = state_map,
    fill = NA,
    color = "grey75",
    linewidth = 0.18
  ) +
  
  {
    if (!is.null(plant_points_sf)) {
      
      geom_sf(
        data = plant_points_sf,
        shape = 21,
        size = 1.5,
        fill = plant_point_color,
        color = "white",
        stroke = 0.35
      )
    }
  } +
  
  scale_fill_gradientn(
    colours = attributable_share_colors,
    
    values = scales::rescale(
      attributable_share_breaks,
      from = c(0, 100)
    ),
    
    limits = c(0, 100),
    
    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),
    
    labels = c(
      "0%",
      "25%",
      "50%",
      "75%",
      "100%"
    ),
    
    oob = scales::squish,
    na.value = "grey92",
    
    name = expression(
      "Share of revised coal " *
        PM[2.5] *
        "-years attributable to schedule changes"
    )
  ) +
  
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      
      barwidth = grid::unit(
        5.2,
        "in"
      ),
      
      barheight = grid::unit(
        0.16,
        "in"
      ),
      
      ticks = TRUE
    )
  ) +
  
  coord_sf(
    crs = grid_crs,
    
    xlim = c(
      conus_bbox["xmin"],
      conus_bbox["xmax"]
    ),
    
    ylim = c(
      conus_bbox["ymin"],
      conus_bbox["ymax"]
    ),
    
    datum = NA,
    expand = FALSE
  ) +
  
  theme_void(
    base_family = "sans"
  ) +
  
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    
    legend.title = element_text(
      size = 7.5,
      hjust = 0.5,
      margin = margin(
        b = 3
      )
    ),
    
    legend.text = element_text(
      size = 6
    ),
    
    legend.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    legend.box.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )

# ------------------------------------------------------------
# 13. Save Fig. S5
# ------------------------------------------------------------

ggsave(
  filename = fig_s3_png_fp,
  plot = p_fig_s3,
  width = 7.25,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = fig_s3_pdf_fp,
  plot = p_fig_s3,
  width = 7.25,
  height = 4.8,
  device = cairo_pdf,
  bg = "white"
)

message(
  "Saved Fig. S5 PNG: ",
  fig_s3_png_fp
)

message(
  "Saved Fig. S5 PDF: ",
  fig_s3_pdf_fp
)

message(
  "Saved ZCTA-level attributable-share table: ",
  attributable_share_table_fp
)

message(
  "Saved diagnostics: ",
  attributable_share_diagnostics_fp
)