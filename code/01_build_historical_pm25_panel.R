# ============================================================
# 01_build_historical_pm25_panel.R
# Build historical plant-ZCTA coal PM2.5 panel for target plants
# ============================================================

source("R/00_config.R")
source("R/00_functions.R")

load_required_packages(c(
  "fst", "sf", "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble"
))

# ------------------------------------------------------------
# 1. Build grid metadata from the first annual by-unit file
# ------------------------------------------------------------
first_fst <- file.path(dir_raw_annual_byunit, sprintf(byunit_fst_pattern, years_historical[1]))

xy0 <- read_fst(first_fst, columns = c("x", "y")) |>
  distinct(x, y) |>
  arrange(y, x)

grid_meta <- xy0 |>
  mutate(grid_id = row_number())

x_unique <- sort(unique(grid_meta$x))
y_unique <- sort(unique(grid_meta$y))
dx <- median(diff(x_unique), na.rm = TRUE)
dy <- median(diff(y_unique), na.rm = TRUE)

grid_polys <- st_sf(
  grid_id = grid_meta$grid_id,
  geometry = st_sfc(
    lapply(seq_len(nrow(grid_meta)), function(i) {
      make_square(grid_meta$x[i], grid_meta$y[i], dx, dy)
    }),
    crs = grid_crs
  )
)

write_fst(grid_meta, file.path(dir_intermediate, "grid_metadata.fst"), compress = 80)

# ------------------------------------------------------------
# 2. Create or read grid-to-ZCTA area weights
# ------------------------------------------------------------
if (file.exists(grid_zcta_weights_fp)) {
  message("Reading existing grid-ZCTA weights: ", grid_zcta_weights_fp)
  weights_area <- read_csv(grid_zcta_weights_fp, show_col_types = FALSE) |>
    mutate(
      grid_id = as.integer(grid_id),
      zcta = standardize_zcta(zcta),
      weight_area = as.numeric(weight_area)
    )
} else {
  message("Computing grid-ZCTA intersections. This may take time.")

  zcta <- st_read(zcta_shp, quiet = TRUE) |>
    st_transform(grid_crs)

  zcta_id_col <- "ZCTA5CE20"
  if (!zcta_id_col %in% names(zcta)) {
    stop("Could not find ZCTA column ", zcta_id_col, " in ", zcta_shp)
  }

  zcta <- zcta |>
    select(all_of(zcta_id_col)) |>
    st_make_valid()

  intersections <- st_intersection(
    grid_polys |> select(grid_id),
    zcta |> select(all_of(zcta_id_col))
  ) |>
    mutate(int_area = as.numeric(st_area(geometry))) |>
    st_drop_geometry()

  weights_area <- intersections |>
    group_by(.data[[zcta_id_col]]) |>
    mutate(weight_area = int_area / sum(int_area, na.rm = TRUE)) |>
    ungroup() |>
    select(grid_id, all_of(zcta_id_col), weight_area) |>
    rename(zcta = all_of(zcta_id_col)) |>
    mutate(
      grid_id = as.integer(grid_id),
      zcta = standardize_zcta(zcta),
      weight_area = as.numeric(weight_area)
    )

  write_csv(weights_area, grid_zcta_weights_fp)
}

# ------------------------------------------------------------
# 3. Process one year of by-unit PM2.5 into plant-ZCTA PM2.5
# ------------------------------------------------------------
process_one_year_plant <- function(yy) {
  message("Processing year: ", yy)
  fst_file <- file.path(dir_raw_annual_byunit, sprintf(byunit_fst_pattern, yy))
  all_cols <- get_fst_colnames(fst_file)
  unit_cols_all <- setdiff(all_cols, c("x", "y"))
  unit_col_map <- tibble(
    unit_col = unit_cols_all,
    unit_id = str_remove(unit_cols_all, "^X"),
    plant_id = extract_plant_id_from_unit_col(unit_cols_all)
  ) |>
    filter(plant_id %in% target_plants$plant_id)
  
  dat <- read_fst(fst_file, columns = c("x", "y", unit_col_map$unit_col)) |>
    left_join(grid_meta, by = c("x", "y"))

  if (any(is.na(dat$grid_id))) {
    stop("Some x/y coordinates did not match grid_meta in year ", yy)
  }

  plant_grid <- lapply(sort(unique(unit_col_map$plant_id)), function(pid) {
    cols_pid <- unit_col_map |>
      filter(plant_id == pid) |>
      pull(unit_col)

    vals <- dat[, cols_pid, drop = FALSE]
    tibble(
      grid_id = dat$grid_id,
      plant_id = pid,
      coal_pm25_plant_grid = rowSums(as.matrix(vals), na.rm = TRUE)
    )
  }) |>
    bind_rows()

  plant_grid |>
    left_join(weights_area, by = "grid_id") |>
    mutate(pm25_weighted = coal_pm25_plant_grid * weight_area) |>
    group_by(year = yy, plant_id, zcta) |>
    summarise(coal_pm25_plant = sum(pm25_weighted, na.rm = TRUE), .groups = "drop") |>
    left_join(target_plants, by = "plant_id") |>
    transmute(
      year,
      plant_id,
      plant_name = plant_name_target,
      zcta = standardize_zcta(zcta),
      coal_pm25_plant = as.numeric(coal_pm25_plant)
    )
}


# ------------------------------------------------------------
# 4A. Build all-coal historical ZCTA-year PM2.5 panel
# ------------------------------------------------------------
# This panel is not restricted to target delayed-retirement plants.
# It sums all unit-level coal PM2.5 columns at each grid cell and
# then applies the same grid-to-ZCTA area weights used above.

all_coal_zcta_pm25_fp <- file.path(
  dir_intermediate,
  "all_coal_zcta_pm25_1999_2020.fst"
)

all_coal_zcta_pm25_csv_fp <- file.path(
  dir_intermediate,
  "all_coal_zcta_pm25_1999_2020.csv"
)

all_coal_zcta_mean_pm25_fp <- file.path(
  dir_intermediate,
  "all_coal_zcta_mean_pm25_1999_2020.csv"
)

# Reload grid metadata if this block is run independently.
if (!exists("grid_meta")) {
  message("grid_meta not found in environment. Reading from intermediate file.")
  
  grid_meta <- read_fst(
    file.path(dir_intermediate, "grid_metadata.fst")
  ) |>
    mutate(
      grid_id = as.integer(grid_id),
      x = as.numeric(x),
      y = as.numeric(y)
    )
}

# Reload grid-ZCTA weights if this block is run independently.
if (!exists("weights_area")) {
  message("weights_area not found in environment. Reading from grid_zcta_weights_fp.")
  
  weights_area <- read_csv(
    grid_zcta_weights_fp,
    show_col_types = FALSE
  ) |>
    mutate(
      grid_id = as.integer(grid_id),
      zcta = standardize_zcta(zcta),
      weight_area = as.numeric(weight_area)
    )
}

grid_meta_all <- grid_meta |>
  select(x, y, grid_id) |>
  mutate(
    x = as.numeric(x),
    y = as.numeric(y),
    grid_id = as.integer(grid_id)
  )

weights_area_all <- weights_area |>
  select(grid_id, zcta, weight_area) |>
  mutate(
    grid_id = as.integer(grid_id),
    zcta = standardize_zcta(zcta),
    weight_area = as.numeric(weight_area)
  )

process_one_year_all_coal <- function(yy, grid_meta_ref, weights_area_ref) {
  
  message("Processing all-coal year: ", yy)
  
  fst_file <- file.path(
    dir_raw_annual_byunit,
    sprintf(byunit_fst_pattern, yy)
  )
  
  if (!file.exists(fst_file)) {
    stop("FST file does not exist: ", fst_file)
  }
  
  all_cols <- get_fst_colnames(fst_file)
  
  if (!all(c("x", "y") %in% all_cols)) {
    stop("FST file must contain x and y columns: ", fst_file)
  }
  
  unit_cols_all <- setdiff(all_cols, c("x", "y"))
  
  if (length(unit_cols_all) == 0) {
    stop("No unit-level PM2.5 columns found in year ", yy)
  }
  
  message("  Unit columns found: ", length(unit_cols_all))
  
  dat <- read_fst(
    fst_file,
    columns = c("x", "y", unit_cols_all)
  ) |>
    mutate(
      x = as.numeric(x),
      y = as.numeric(y)
    ) |>
    left_join(
      grid_meta_ref,
      by = c("x", "y")
    )
  
  if (any(is.na(dat$grid_id))) {
    stop("Some x/y coordinates did not match grid_meta in year ", yy)
  }
  
  unit_cols_numeric <- unit_cols_all[
    vapply(dat[unit_cols_all], is.numeric, logical(1))
  ]
  
  if (length(unit_cols_numeric) == 0) {
    stop("No numeric unit-level PM2.5 columns found in year ", yy)
  }
  
  grid_total <- tibble(
    grid_id = dat$grid_id,
    coal_pm25_grid_total = rowSums(
      as.matrix(dat[, unit_cols_numeric, drop = FALSE]),
      na.rm = TRUE
    )
  )
  
  zcta_year <- grid_total |>
    left_join(
      weights_area_ref,
      by = "grid_id"
    ) |>
    filter(
      !is.na(zcta),
      !is.na(weight_area),
      is.finite(weight_area)
    ) |>
    mutate(
      pm25_weighted = coal_pm25_grid_total * weight_area
    ) |>
    group_by(zcta) |>
    summarise(
      year = yy,
      coal_pm25_zcta = sum(pm25_weighted, na.rm = TRUE),
      .groups = "drop"
    ) |>
    transmute(
      year = as.integer(year),
      zcta = standardize_zcta(zcta),
      coal_pm25_zcta = as.numeric(coal_pm25_zcta)
    )
  
  zcta_year
}

all_coal_zcta_pm25_panel <- purrr::map_dfr(
  years_historical,
  \(yy) process_one_year_all_coal(
    yy = yy,
    grid_meta_ref = grid_meta_all,
    weights_area_ref = weights_area_all
  )
) |>
  arrange(year, zcta)

write_fst(
  all_coal_zcta_pm25_panel,
  all_coal_zcta_pm25_fp,
  compress = 80
)

write_csv(
  all_coal_zcta_pm25_panel,
  all_coal_zcta_pm25_csv_fp
)

all_coal_zcta_mean_pm25_1999_2020 <- all_coal_zcta_pm25_panel |>
  group_by(zcta) |>
  summarise(
    coal_pm25_mean_1999_2020 = mean(coal_pm25_zcta, na.rm = TRUE),
    n_years_observed = n_distinct(year),
    .groups = "drop"
  ) |>
  arrange(zcta)

write_csv(
  all_coal_zcta_mean_pm25_1999_2020,
  all_coal_zcta_mean_pm25_fp
)

cat("All-coal ZCTA-year rows:", nrow(all_coal_zcta_pm25_panel), "\n")
cat("All-coal unique years:", n_distinct(all_coal_zcta_pm25_panel$year), "\n")
cat("All-coal unique ZCTAs:", n_distinct(all_coal_zcta_pm25_panel$zcta), "\n")

all_coal_zcta_mean_pm25_1999_2020 |>
  summarise(
    n_zcta = n_distinct(zcta),
    min_pm25 = min(coal_pm25_mean_1999_2020, na.rm = TRUE),
    p25_pm25 = quantile(coal_pm25_mean_1999_2020, 0.25, na.rm = TRUE),
    median_pm25 = median(coal_pm25_mean_1999_2020, na.rm = TRUE),
    mean_pm25 = mean(coal_pm25_mean_1999_2020, na.rm = TRUE),
    p75_pm25 = quantile(coal_pm25_mean_1999_2020, 0.75, na.rm = TRUE),
    max_pm25 = max(coal_pm25_mean_1999_2020, na.rm = TRUE)
  ) |>
  print()


# ------------------------------------------------------------
# 4. Build the full historical plant-ZCTA-year panel
# ------------------------------------------------------------
plant_zcta_pm25_panel <- map_dfr(years_historical, process_one_year_plant)

write_fst(plant_zcta_pm25_panel, target_plant_zcta_pm25_fp, compress = 80)
