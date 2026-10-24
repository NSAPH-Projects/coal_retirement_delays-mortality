# ==================================================
# Baseline heterogeneity and historical coal PM2.5 map
#
# ==================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(scales)
library(sf)

# --------------------------------------------------
# 0) Settings
# --------------------------------------------------

if (!exists("dir_tables")) {
  dir_tables <- file.path("output", "tables")
}

if (!exists("dir_figures")) {
  dir_figures <- file.path("output", "figures")
}

if (!exists("crs_grid")) {
  crs_grid <- 5070
}

if (!exists("zcta_shp")) {
  zcta_shp <- file.path(
    dir_raw,
    "tl_2020_us_zcta520",
    "tl_2020_us_zcta520.shp"
  )
}

dir.create(dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

to_share <- function(x) {
  x <- as.numeric(x)
  if_else(!is.na(x) & x > 1, x / 100, x)
}

# ==================================================
# 1) Read ZIP-year coal PM2.5 and treat ZIP as ZCTA
# ==================================================

zip_pm25_fp <- file.path(dir_raw, "annual/coalPM25_zip_year.csv")

zip_pm25_raw <- read_csv(
  zip_pm25_fp,
  show_col_types = FALSE
)

zip_col <- "zip"
year_col <- "year"
pm25_col <- "coal_PM25"

zip_year_pm25 <- zip_pm25_raw %>%
  transmute(
    zcta = str_pad(
      as.character(.data[[zip_col]]),
      width = 5,
      side = "left",
      pad = "0"
    ),
    year = as.integer(.data[[year_col]]),
    coal_pm25 = as.numeric(.data[[pm25_col]])
  ) %>%
  filter(
    year %in% years_historical,
    !is.na(zcta),
    !is.na(year),
    is.finite(coal_pm25)
  )

demo_zcta <- demo %>%
  transmute(
    zcta = str_pad(
      as.character(zcta),
      width = 5,
      side = "left",
      pad = "0"
    )
  ) %>%
  distinct()

zip_year_pm25_matched <- zip_year_pm25 %>%
  inner_join(demo_zcta, by = "zcta")

matched_zcta <- zip_year_pm25_matched %>%
  distinct(zcta)

coverage_check <- demo %>%
  mutate(
    zcta = str_pad(
      as.character(zcta),
      width = 5,
      side = "left",
      pad = "0"
    )
  ) %>%
  summarise(
    total_zcta = n_distinct(zcta),
    matched_zcta_n = n_distinct(zcta[zcta %in% matched_zcta$zcta]),
    matched_zcta_share = matched_zcta_n / total_zcta,
    
    total_population = sum(population, na.rm = TRUE),
    matched_population =
      sum(population[zcta %in% matched_zcta$zcta], na.rm = TRUE),
    matched_population_share = matched_population / total_population,
    
    total_medicare = sum(n_bene, na.rm = TRUE),
    matched_medicare =
      sum(n_bene[zcta %in% matched_zcta$zcta], na.rm = TRUE),
    matched_medicare_share = matched_medicare / total_medicare
  )

print(
  coverage_check %>%
    mutate(
      matched_zcta_share = percent(matched_zcta_share, accuracy = 0.1),
      matched_population_share = percent(matched_population_share, accuracy = 0.1),
      matched_medicare_share = percent(matched_medicare_share, accuracy = 0.1)
    )
)

coal_pm25_zcta_mean_1999_2020 <- zip_year_pm25_matched %>%
  group_by(zcta) %>%
  summarise(
    coal_pm25_mean_1999_2020 =
      mean(coal_pm25, na.rm = TRUE),
    n_years_observed =
      n_distinct(year),
    .groups = "drop"
  )

print(summary(coal_pm25_zcta_mean_1999_2020$coal_pm25_mean_1999_2020))
print(nrow(coal_pm25_zcta_mean_1999_2020))

# ==================================================
# 2) Restrict demographics to matched ZIP/ZCTA universe
# ==================================================

demo_2020 <- demo %>%
  mutate(
    zcta = str_pad(
      as.character(zcta),
      width = 5,
      side = "left",
      pad = "0"
    )
  ) %>%
  semi_join(
    coal_pm25_zcta_mean_1999_2020,
    by = "zcta"
  )

required_demo_cols <- c(
  "zcta",
  "population",
  "female_pct",
  "male_pct",
  "pop_white",
  "pop_black",
  "pop_asian",
  "pop_native",
  "pop_hispanic",
  "pop_under_20_years",
  "pop_35_49_years",
  "pop_50_64_years",
  "pop_over_65_years"
)

missing_demo_cols <- setdiff(required_demo_cols, names(demo_2020))

if (length(missing_demo_cols) > 0) {
  stop(
    paste0(
      "Missing required columns in demo_2020: ",
      paste(missing_demo_cols, collapse = ", ")
    )
  )
}

demo_2020_sub <- demo_2020 %>%
  mutate(
    zcta = str_pad(as.character(zcta), width = 5, side = "left", pad = "0"),
    
    population = as.numeric(population),
    
    female_share_raw = to_share(female_pct),
    male_share_raw = to_share(male_pct),
    
    pop_white = replace_na(as.numeric(pop_white), 0),
    pop_black = replace_na(as.numeric(pop_black), 0),
    pop_asian = replace_na(as.numeric(pop_asian), 0),
    pop_native = replace_na(as.numeric(pop_native), 0),
    pop_hispanic = replace_na(as.numeric(pop_hispanic), 0),
    
    pop_under_20_years = replace_na(as.numeric(pop_under_20_years), 0),
    pop_35_49_years = replace_na(as.numeric(pop_35_49_years), 0),
    pop_50_64_years = replace_na(as.numeric(pop_50_64_years), 0),
    pop_over_65_years = replace_na(as.numeric(pop_over_65_years), 0)
  ) %>%
  filter(
    !is.na(zcta),
    !is.na(population),
    population > 0
  ) %>%
  mutate(
    sex_share_sum = male_share_raw + female_share_raw,
    male_share = if_else(
      sex_share_sum > 0,
      male_share_raw / sex_share_sum,
      NA_real_
    ),
    female_share = if_else(
      sex_share_sum > 0,
      female_share_raw / sex_share_sum,
      NA_real_
    ),
    
    pop_20_34_years = pmax(
      population -
        pop_under_20_years -
        pop_35_49_years -
        pop_50_64_years -
        pop_over_65_years,
      0
    )
  )

# ==================================================
# 3) Build subgroup populations: Sex, Race, Age
# ==================================================

# ------------------------------
# Sex
# ------------------------------

sex_pop <- bind_rows(
  demo_2020_sub %>%
    transmute(
      zcta,
      group_type = "Sex",
      group = "Male",
      subgroup_population = population * male_share
    ),
  demo_2020_sub %>%
    transmute(
      zcta,
      group_type = "Sex",
      group = "Female",
      subgroup_population = population * female_share
    )
)

# ------------------------------
# Race
# ------------------------------
# Five race categories are shown.
# Other/residual is excluded from the figure.
# Counts are normalized within each ZCTA over the displayed categories.

race_pop <- demo_2020_sub %>%
  transmute(
    zcta,
    population,
    White = pop_white,
    Black = pop_black,
    Asian = pop_asian,
    Native = pop_native,
    Hispanic = pop_hispanic
  ) %>%
  pivot_longer(
    cols = c(White, Black, Asian, Native, Hispanic),
    names_to = "group",
    values_to = "raw_group_count"
  ) %>%
  group_by(zcta) %>%
  mutate(
    race_total = sum(raw_group_count, na.rm = TRUE),
    subgroup_share = if_else(
      race_total > 0,
      raw_group_count / race_total,
      NA_real_
    ),
    subgroup_population = population * subgroup_share
  ) %>%
  ungroup() %>%
  transmute(
    zcta,
    group_type = "Race",
    group,
    subgroup_population
  )

# ------------------------------
# Age
# ------------------------------

age_pop <- demo_2020_sub %>%
  transmute(
    zcta,
    `Under 20` = pop_under_20_years,
    `20-34` = pop_20_34_years,
    `35-49` = pop_35_49_years,
    `50-64` = pop_50_64_years,
    `65+` = pop_over_65_years
  ) %>%
  pivot_longer(
    cols = c(`Under 20`, `20-34`, `35-49`, `50-64`, `65+`),
    names_to = "group",
    values_to = "subgroup_population"
  ) %>%
  mutate(
    group_type = "Age"
  )

subgroup_pop_2020 <- bind_rows(
  sex_pop,
  race_pop,
  age_pop
) %>%
  filter(
    !is.na(subgroup_population),
    subgroup_population > 0
  )

# ==================================================
# 4) QC
# ==================================================

subgroup_pop_qc <- subgroup_pop_2020 %>%
  group_by(group_type, zcta) %>%
  summarise(
    subgroup_total = sum(subgroup_population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    demo_2020_sub %>% select(zcta, population),
    by = "zcta"
  ) %>%
  mutate(
    ratio_to_population = subgroup_total / population
  ) %>%
  group_by(group_type) %>%
  summarise(
    min_ratio = min(ratio_to_population, na.rm = TRUE),
    median_ratio = median(ratio_to_population, na.rm = TRUE),
    max_ratio = max(ratio_to_population, na.rm = TRUE),
    .groups = "drop"
  )

print(subgroup_pop_qc)

# ==================================================
# 5) Population-weighted overall reference
# ==================================================

overall_reference_all <- demo_2020_sub %>%
  inner_join(
    coal_pm25_zcta_mean_1999_2020,
    by = "zcta"
  ) %>%
  summarise(
    reference_population = sum(population, na.rm = TRUE),
    reference_exposure =
      sum(population * coal_pm25_mean_1999_2020, na.rm = TRUE) /
      sum(population, na.rm = TRUE)
  )

overall_pw_coal_pm25 <- overall_reference_all$reference_exposure

print(overall_reference_all)

# ==================================================
# 6) Subgroup population-weighted exposure
# ==================================================

subgroup_exposure_1999_2020 <- subgroup_pop_2020 %>%
  inner_join(
    coal_pm25_zcta_mean_1999_2020,
    by = "zcta"
  ) %>%
  mutate(
    weighted_coal_pm25 =
      subgroup_population * coal_pm25_mean_1999_2020
  ) %>%
  group_by(group_type, group) %>%
  summarise(
    subgroup_population_2020 =
      sum(subgroup_population, na.rm = TRUE),
    exposure_burden =
      sum(weighted_coal_pm25, na.rm = TRUE),
    pw_coal_pm25_mean_1999_2020 =
      exposure_burden / subgroup_population_2020,
    .groups = "drop"
  ) %>%
  mutate(
    reference_exposure = overall_pw_coal_pm25,
    difference_vs_overall =
      pw_coal_pm25_mean_1999_2020 - reference_exposure,
    ratio_vs_overall =
      pw_coal_pm25_mean_1999_2020 / reference_exposure,
    pct_difference_vs_overall =
      100 * (ratio_vs_overall - 1)
  )

domain_order <- c("Sex", "Race", "Age")

group_levels <- tibble(
  group_type = c(
    rep("Sex", 2),
    rep("Race", 5),
    rep("Age", 5)
  ),
  group = c(
    "Male", "Female",
    "White", "Black", "Asian", "Native", "Hispanic",
    "Under 20", "20-34", "35-49", "50-64", "65+"
  ),
  group_order = c(
    1, 2,
    1, 2, 3, 4, 5,
    1, 2, 3, 4, 5
  )
)

baseline_heterogeneity_1999_2020 <- subgroup_exposure_1999_2020 %>%
  left_join(
    group_levels,
    by = c("group_type", "group")
  ) %>%
  mutate(
    group_type = factor(group_type, levels = domain_order)
  ) %>%
  arrange(group_type, group_order)

print(baseline_heterogeneity_1999_2020, n = Inf)

write_csv(
  baseline_heterogeneity_1999_2020,
  file.path(
    dir_tables,
    "subgroup_coal_pm25_annual_mean_1999_2020_vs_overall_sex_race_age.csv"
  )
)

# ==================================================
# 7) Dumbbell plot: Sex, Race, Age
# ==================================================

plot_df <- baseline_heterogeneity_1999_2020 %>%
  arrange(group_type, group_order) %>%
  mutate(
    group_label = recode(
      group,
      "Native" = "AI/AN"
    ),
    label_text = case_when(
      pct_difference_vs_overall >= 0 ~
        paste0("+", number(pct_difference_vs_overall, accuracy = 0.1), "%"),
      TRUE ~
        paste0(number(pct_difference_vs_overall, accuracy = 0.1), "%")
    )
  )

gap <- 0.9

domain_meta <- plot_df %>%
  count(group_type, name = "n") %>%
  mutate(
    domain_order_id = match(as.character(group_type), domain_order)
  ) %>%
  arrange(domain_order_id) %>%
  mutate(
    start_y = lag(cumsum(n + gap), default = 0) + 1,
    end_y = start_y + n - 1,
    mid_y = (start_y + end_y) / 2,
    sep_y = end_y + gap / 2
  )

plot_df <- plot_df %>%
  left_join(
    domain_meta %>% select(group_type, start_y),
    by = "group_type"
  ) %>%
  group_by(group_type) %>%
  arrange(group_order, .by_group = TRUE) %>%
  mutate(
    y = start_y + row_number() - 1
  ) %>%
  ungroup()

separator_df <- domain_meta %>%
  filter(as.character(group_type) != tail(domain_order, 1))

x_min_data <- min(
  c(plot_df$pw_coal_pm25_mean_1999_2020, plot_df$reference_exposure),
  na.rm = TRUE
)

x_max_data <- max(
  c(plot_df$pw_coal_pm25_mean_1999_2020, plot_df$reference_exposure),
  na.rm = TRUE
)

x_range <- x_max_data - x_min_data

if (!is.finite(x_range) || x_range <= 0) {
  x_range <- max(abs(x_max_data), 1) * 0.20
}

# Tighter left-side label space
x_data_left <- x_min_data - 0.03 * x_range
x_right <- x_max_data + 0.16 * x_range
x_left <- x_min_data - 0.28 * x_range

# Move labels closer to the plotting region
x_domain_label <- x_min_data - 0.20 * x_range
x_group_label <- x_min_data - 0.07 * x_range

x_value_label_offset <- 0.022 * x_range

x_breaks <- breaks_extended(n = 5)(c(x_min_data, x_max_data))
x_breaks <- x_breaks[
  x_breaks >= x_data_left &
    x_breaks <= x_right
]

separator_segment_df <- separator_df %>%
  mutate(
    x = x_data_left,
    xend = x_right
  )

p_longterm_dumbbell <- ggplot(plot_df) +
  geom_vline(
    xintercept = x_breaks,
    color = "grey90",
    linewidth = 0.45
  ) +
  geom_vline(
    xintercept = overall_pw_coal_pm25,
    color = "grey35",
    linetype = "dashed",
    linewidth = 0.55
  ) +
  geom_segment(
    data = separator_segment_df,
    aes(
      x = x,
      xend = xend,
      y = sep_y,
      yend = sep_y
    ),
    inherit.aes = FALSE,
    color = "grey82",
    linewidth = 0.45
  ) +
  geom_segment(
    aes(
      x = reference_exposure,
      xend = pw_coal_pm25_mean_1999_2020,
      y = y,
      yend = y
    ),
    color = "grey60",
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_point(
    aes(
      x = reference_exposure,
      y = y
    ),
    shape = 16,
    size = 2.7,
    color = "#2C7FB8"
  ) +
  geom_point(
    aes(
      x = pw_coal_pm25_mean_1999_2020,
      y = y
    ),
    shape = 18,
    size = 3.3,
    color = "#F28E2B"
  ) +
  geom_text(
    aes(
      x = if_else(
        pw_coal_pm25_mean_1999_2020 >= reference_exposure,
        pw_coal_pm25_mean_1999_2020 + x_value_label_offset,
        reference_exposure + x_value_label_offset
      ),
      y = y,
      label = label_text
    ),
    hjust = 0,
    size = 3.1,
    color = "black"
  ) +
  geom_text(
    data = domain_meta,
    aes(
      x = x_domain_label,
      y = mid_y,
      label = as.character(group_type)
    ),
    inherit.aes = FALSE,
    hjust = 1,
    fontface = "bold",
    size = 3.5
  ) +
  geom_text(
    aes(
      x = x_group_label,
      y = y,
      label = group_label
    ),
    hjust = 1,
    size = 3.4,
    color = "grey20"
  ) +
  scale_y_reverse(
    breaks = NULL,
    expand = expansion(add = 0.45)
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = label_number(accuracy = 0.001)
  ) +
  coord_cartesian(
    xlim = c(x_left, x_right),
    clip = "off"
  ) +
  labs(
    x = expression("Weighted annual mean coal-" * PM[2.5] * "(" * mu * "g/m"^3 * ")"),
    y = NULL,
    ) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 20, 18, 35)
  )

p_longterm_dumbbell

ggsave(
  file.path(
    dir_figures,
    "figure_longterm_1999_2020_heterogeneity_vs_overall_sex_race_age.png"
  ),
  p_longterm_dumbbell,
  width = 7.25,
  height = 4.8,
  dpi = 600,
  bg = "white"
)


# ==================================================
# Build all-coal ZCTA-year PM2.5 from grid files
# Not restricted to the 20 delayed plants
# ==================================================
# Put this in 01 after grid_meta and weights_area are created.
#
# Required existing objects from 01:
#   weights_area
#     columns: grid_id, zcta, weight_area
#   historical_fst_files
#     named or identifiable by year
#   years_historical
#
# Each grid file should contain:
#   x, y, and coal PM2.5 columns for units/plants.
# ==================================================

library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(fst)

# --------------------------------------------------
# 1) Helper: extract year from file name
# --------------------------------------------------

extract_year_from_fp <- function(fp) {
  as.integer(str_extract(basename(fp), "19[0-9]{2}|20[0-9]{2}"))
}

# --------------------------------------------------
# 2) Build all-coal ZCTA-year PM2.5
# --------------------------------------------------

historical_fst_files <- file.path(
  dir_raw_annual_byunit,
  sprintf(byunit_fst_pattern, years_historical)
)

all_coal_zcta_year <- map_dfr(
  historical_fst_files,
  function(fp) {
    
    yy <- extract_year_from_fp(fp)
    
    if (!yy %in% years_historical) {
      return(NULL)
    }
    
    message("Aggregating all-coal grid PM2.5 to ZCTA for year: ", yy)
    
    dat <- read_fst(fp)
    
    # Make sure grid ordering matches grid_meta from 01.
    dat <- dat %>%
      arrange(y, x) %>%
      mutate(
        grid_id = row_number()
      )
    
    # Use every numeric PM2.5 column except coordinates and grid id.
    pm25_cols <- setdiff(
      names(dat),
      c("x", "y", "grid_id")
    )
    
    # Optional safety: keep numeric columns only.
    pm25_cols <- pm25_cols[
      vapply(dat[pm25_cols], is.numeric, logical(1))
    ]
    
    if (length(pm25_cols) == 0) {
      stop("No numeric PM2.5 columns found in: ", fp)
    }
    
    grid_total <- dat %>%
      transmute(
        grid_id,
        coal_pm25_grid_total =
          rowSums(
            as.matrix(across(all_of(pm25_cols))),
            na.rm = TRUE
          )
      )
    
    grid_total %>%
      inner_join(
        weights_area,
        by = "grid_id"
      ) %>%
      mutate(
        coal_pm25_weighted =
          coal_pm25_grid_total * weight_area
      ) %>%
      group_by(zcta) %>%
      summarise(
        year = yy,
        coal_pm25_zcta =
          sum(coal_pm25_weighted, na.rm = TRUE),
        .groups = "drop"
      )
  }
)

# --------------------------------------------------
# 3) Save ZCTA-year all-coal PM2.5
# --------------------------------------------------

all_coal_zcta_year <- all_coal_zcta_year %>%
  mutate(
    zcta = str_pad(as.character(zcta), width = 5, side = "left", pad = "0"),
    year = as.integer(year)
  ) %>%
  arrange(zcta, year)

write_csv(
  all_coal_zcta_year,
  file.path(
    dir_intermediate,
    "all_coal_pm25_zcta_year_1999_2020_from_grid.csv"
  )
)

# --------------------------------------------------
# 4) Build 1999-2020 ZCTA mean
# --------------------------------------------------

coal_pm25_zcta_mean_1999_2020_grid_all <- all_coal_zcta_year %>%
  filter(year %in% years_historical) %>%
  group_by(zcta) %>%
  summarise(
    coal_pm25_mean_1999_2020 =
      mean(coal_pm25_zcta, na.rm = TRUE),
    n_years_observed =
      n_distinct(year),
    .groups = "drop"
  )

write_csv(
  coal_pm25_zcta_mean_1999_2020_grid_all,
  file.path(
    dir_intermediate,
    "all_coal_pm25_zcta_mean_1999_2020_from_grid.csv"
  )
)

# --------------------------------------------------
# 5) Coverage check
# --------------------------------------------------

coal_pm25_zcta_mean_1999_2020_grid_all %>%
  summarise(
    n_zcta = n_distinct(zcta),
    min_pm25 = min(coal_pm25_mean_1999_2020, na.rm = TRUE),
    median_pm25 = median(coal_pm25_mean_1999_2020, na.rm = TRUE),
    mean_pm25 = mean(coal_pm25_mean_1999_2020, na.rm = TRUE),
    max_pm25 = max(coal_pm25_mean_1999_2020, na.rm = TRUE)
  ) %>%
  print()

# ==================================================
# 8) ZCTA-level map: 1999-2020 annual mean coal PM2.5
# ==================================================
# This map uses ZCTA-level coal PM2.5 directly.
# No population weighting is used in the map.

coal_pm25_zcta_mean_1999_2020 <- read_csv(
  file.path(
    dir_intermediate,
    "all_coal_zcta_mean_pm25_1999_2020.csv"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    zcta = standardize_zcta(zcta),
    coal_pm25_mean_1999_2020 = as.numeric(coal_pm25_mean_1999_2020),
    n_years_observed = as.integer(n_years_observed)
  )

coal_pm25_zcta_mean_1999_2020 %>%
  summarise(
    n_zcta = n_distinct(zcta),
    min_pm25 = min(coal_pm25_mean_1999_2020, na.rm = TRUE),
    median_pm25 = median(coal_pm25_mean_1999_2020, na.rm = TRUE),
    mean_pm25 = mean(coal_pm25_mean_1999_2020, na.rm = TRUE),
    max_pm25 = max(coal_pm25_mean_1999_2020, na.rm = TRUE)
  ) %>%
  print()

# --------------------------------------------------
# 2) Read ZCTA geometry and attach ZCTA-level coal PM2.5
# --------------------------------------------------

zcta_map_historical <- st_read(zcta_shp, quiet = TRUE) %>%
  transmute(
    zcta = standardize_zcta(ZCTA5CE20),
    geometry = geometry
  ) %>%
  st_transform(grid_crs) %>%
  left_join(
    coal_pm25_zcta_mean_1999_2020 %>%
      select(zcta, coal_pm25_mean_1999_2020, n_years_observed),
    by = "zcta"
  )

# --------------------------------------------------
# 3) Read CONUS state boundaries
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
# 4) QC map coverage
# --------------------------------------------------

zcta_map_historical %>%
  st_drop_geometry() %>%
  summarise(
    total_zcta = n(),
    zcta_with_pm25 = sum(!is.na(coal_pm25_mean_1999_2020)),
    zcta_missing_pm25 = sum(is.na(coal_pm25_mean_1999_2020)),
    missing_share = zcta_missing_pm25 / total_zcta
  ) %>%
  mutate(
    missing_share = percent(missing_share, accuracy = 0.1)
  ) %>%
  print()

# --------------------------------------------------
# 5) Plot
# --------------------------------------------------
# No cap is applied. The raw ZCTA-level 1999-2020 annual mean
# coal PM2.5 values are used directly.

p_map_historical_coal_pm25 <- ggplot() +
  geom_sf(
    data = zcta_map_historical,
    aes(fill = coal_pm25_mean_1999_2020),
    color = NA
  ) +
  geom_sf(
    data = state_map,
    fill = NA,
    color = "grey80",
    linewidth = 0.18
  ) +
  scale_fill_gradientn(
    colours = c(
      "#fff7bc",
      "#fec44f",
      "#f16913",
      "#b30000",
      "#4d0000"
    ),
    labels = label_number(accuracy = 0.01),
    na.value = "grey92",
    name = expression(
      "Coal PM"[2.5] *
        "(" * mu * "g/m"^3 * ")"
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
  labs(
    title = expression(bold("Long-term annual mean coal-" * PM[2.5] * ", 1999-2020"))
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
# 6) Save
# --------------------------------------------------

ggsave(
  file.path(
    dir_figures,
    "map_zcta_annual_mean_coal_pm25_1999_2020_section5.png"
  ),
  p_map_historical_coal_pm25,
  width = 7.25,
  height = 4.8,
  dpi = 600,
  bg = "white"
)







# ==================================================
# 8) ZCTA-level map: 1999-2020 annual mean coal PM2.5
# ==================================================
# This map uses the all-coal grid-to-ZCTA output generated in:
#   01_build_historical_pm25_panel.R
#
# It is not restricted to the 20 delayed-retirement plants.
# No population weighting is used in the map.
# No color-scale cap is applied.

# --------------------------------------------------
# 8.1 Settings and fallbacks
# --------------------------------------------------

if (!exists("grid_crs")) {
  if (exists("crs_grid")) {
    grid_crs <- crs_grid
  } else {
    grid_crs <- 5070
  }
}

if (!exists("state_shp")) {
  state_shp <- file.path(
    dir_raw,
    "tl_2020_us_state",
    "tl_2020_us_state.shp"
  )
}

if (!exists("zcta_shp")) {
  zcta_shp <- file.path(
    dir_raw,
    "tl_2020_us_zcta520",
    "tl_2020_us_zcta520.shp"
  )
}

if (!exists("standardize_zcta")) {
  standardize_zcta <- function(x) {
    stringr::str_pad(
      as.character(x),
      width = 5,
      side = "left",
      pad = "0"
    )
  }
}

# --------------------------------------------------
# 8.2 Read all-coal ZCTA-level historical mean PM2.5
# --------------------------------------------------

all_coal_zcta_mean_fp <- file.path(
  dir_intermediate,
  "all_coal_zcta_mean_pm25_1999_2020.csv"
)

if (!file.exists(all_coal_zcta_mean_fp)) {
  stop(
    paste0(
      "Missing file: ",
      all_coal_zcta_mean_fp,
      "\nRun the all-coal grid-to-ZCTA block in 01 first."
    )
  )
}

coal_pm25_zcta_mean_1999_2020_map <- read_csv(
  all_coal_zcta_mean_fp,
  show_col_types = FALSE
) %>%
  transmute(
    zcta = standardize_zcta(zcta),
    coal_pm25_mean_1999_2020 = as.numeric(coal_pm25_mean_1999_2020),
    n_years_observed = as.integer(n_years_observed)
  )

coal_pm25_zcta_mean_1999_2020_map %>%
  summarise(
    n_zcta = n_distinct(zcta),
    min_pm25 = min(coal_pm25_mean_1999_2020, na.rm = TRUE),
    p25_pm25 = quantile(coal_pm25_mean_1999_2020, 0.25, na.rm = TRUE),
    median_pm25 = median(coal_pm25_mean_1999_2020, na.rm = TRUE),
    mean_pm25 = mean(coal_pm25_mean_1999_2020, na.rm = TRUE),
    p75_pm25 = quantile(coal_pm25_mean_1999_2020, 0.75, na.rm = TRUE),
    max_pm25 = max(coal_pm25_mean_1999_2020, na.rm = TRUE)
  ) %>%
  print()

# --------------------------------------------------
# 8.3 Read ZCTA geometry and attach coal PM2.5
# --------------------------------------------------

zcta_map_historical <- st_read(zcta_shp, quiet = TRUE) %>%
  transmute(
    zcta = standardize_zcta(ZCTA5CE20),
    geometry = geometry
  ) %>%
  st_transform(grid_crs) %>%
  left_join(
    coal_pm25_zcta_mean_1999_2020_map,
    by = "zcta"
  )

# --------------------------------------------------
# 8.4 Read state boundaries
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
# 8.5 Restrict geometries to CONUS
# --------------------------------------------------
state_map_conus <- state_map

zcta_map_historical_conus <- zcta_map_historical %>%
  filter(
    !substr(zcta, 1, 3) %in% c("995", "996", "997", "998", "999",
                               "967", "968",
                               "006", "007", "008", "009")
  )

message(
  "CONUS ZCTAs retained: ",
  format(nrow(zcta_map_historical_conus), big.mark = ","),
  " of ", format(nrow(zcta_map_historical), big.mark = ",")
)

# --------------------------------------------------
# 8.7 Plot historical coal PM2.5 map
# --------------------------------------------------

p_map_historical_coal_pm25 <- ggplot() +
  geom_sf(
    data = zcta_map_historical_conus,
    aes(fill = coal_pm25_mean_1999_2020),
    color = NA
  ) +
  geom_sf(
    data = state_map_conus,
    fill = NA,
    color = "grey75",
    linewidth = 0.16
  ) +
  scale_fill_gradientn(
    colours = c(
      "#fff7bc",
      "#fec44f",
      "#f16913",
      "#b30000",
      "#4d0000"
    ),
    labels = label_number(accuracy = 0.01),
    na.value = "grey92",
    name = expression(
      "Coal " * PM[2.5] *
        " (" * mu * "g/m"^3 * ")"
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 0.5,
      title.vjust = 0.5,
      label.position = "bottom",
      barwidth = grid::unit(4.8, "in"),
      barheight = grid::unit(0.16, "in"),
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
  labs(
    title = expression(
      bold("Long-term annual mean coal-" * PM[2.5] * ", 1999-2020")
    )
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 10.5, hjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(
      size = 7.5,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(r = 7)
    ),
    legend.text = element_text(size = 6),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(fill = "white", color = NA),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(8, 8, 8, 8)
  )

# --------------------------------------------------
# 8.8 Save historical coal PM2.5 map
# --------------------------------------------------

ggsave(
  file.path(
    dir_figures,
    "map_zcta_annual_mean_coal_pm25_1999_2020_all_coal_grid.png"
  ),
  p_map_historical_coal_pm25,
  width = 7.25,
  height = 4.8,
  dpi = 600,
  bg = "white"
)
