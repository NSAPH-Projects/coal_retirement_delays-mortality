# ============================================================
# 08_run_emissions_retirement_sensitivities_balanced.R
# One-at-a-time sensitivity analyses for balanced future SO2 scaling assumptions
# and removed retirement dates, 2026-2035
#
# Sensitivity scenarios:
#   1. Main: 100% of plant-specific 2023-2025 mean SO2
#   2. Uniform 75% of plant-specific 2023-2025 mean SO2
#   3. Uniform 50% of plant-specific 2023-2025 mean SO2
#   4. Uniform 125% of plant-specific 2023-2025 mean SO2
#   5. Uniform 150% of plant-specific 2023-2025 mean SO2
#   6. Plant-specific recent-low SO2 (minimum observed in 2023-2025)
#   7. Plant-specific recent-high SO2 (maximum observed in 2023-2025)
#   8. Removed retirement date = original date + 2 years
#   9. Removed retirement date = original date + 5 years
# ============================================================

SENSITIVITY_SCRIPT_VERSION <- "balanced-scaling-minmax-v4-2026-07-20"
message("Running sensitivity script version: ", SENSITIVITY_SCRIPT_VERSION)

source("R/00_config.R")
source("R/00_functions.R")

load_required_packages(c(
  "fst",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "lubridate",
  "purrr",
  "tibble",
  "ggplot2",
  "scales"
))

# ------------------------------------------------------------
# 0. Output locations
# ------------------------------------------------------------

sensitivity_dir <- file.path(dir_output, "sensitivity")
sensitivity_tables_dir <- file.path(sensitivity_dir, "tables")
sensitivity_figures_dir <- file.path(sensitivity_dir, "figures")
sensitivity_diagnostics_dir <- file.path(sensitivity_dir, "diagnostics")

for (d in c(
  sensitivity_dir,
  sensitivity_tables_dir,
  sensitivity_figures_dir,
  sensitivity_diagnostics_dir
)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ------------------------------------------------------------
# 1. Check required inputs
# ------------------------------------------------------------

required_files <- c(
  so2_fp,
  scenario_fp,
  medicare_fp,
  demo_fp,
  model_gamma_fp,
  model_alpha_fp
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste0("- ", missing_files, collapse = "\n"),
    "\nRun the main analysis through script 03 before this sensitivity script."
  )
}

# ------------------------------------------------------------
# 2. Define sensitivity scenarios
# ------------------------------------------------------------
# Each scenario changes one assumption at a time.
#
# baseline_type:
#   mean = plant-specific 2023-2025 arithmetic mean SO2
#   min  = plant-specific minimum annual SO2 during 2023-2025
#   max  = plant-specific maximum annual SO2 during 2023-2025
#
# mean_multiplier applies equally to original and revised plans.
#
# removed_date_rule:
#   through_2035       = use the main assumption for removed dates
#   original_plus_years = replace removed dates with original + N years

sensitivity_specs <- tribble(
  ~sensitivity_id, ~sensitivity_label,                 ~baseline_type, ~mean_multiplier, ~removed_date_rule,    ~removed_date_shift_years,
  "main",            "Main analysis",                 "mean",         1.00,             "through_2035",       0L,
  "uniform_75",      "Uniform 75% emissions",         "mean",         0.75,             "through_2035",       0L,
  "uniform_50",      "Uniform 50% emissions",         "mean",         0.50,             "through_2035",       0L,
  "uniform_125",     "Uniform 125% emissions",        "mean",         1.25,             "through_2035",       0L,
  "uniform_150",     "Uniform 150% emissions",        "mean",         1.50,             "through_2035",       0L,
  "recent_low",      "Plant-specific recent-low",    "min",          1.00,             "through_2035",       0L,
  "recent_high",     "Plant-specific recent-high",   "max",          1.00,             "through_2035",       0L,
  "removed_plus_2",  "Removed date: original +2 y",  "mean",         1.00,             "original_plus_years", 2L,
  "removed_plus_5",  "Removed date: original +5 y",  "mean",         1.00,             "original_plus_years", 5L
)

write_csv(
  sensitivity_specs,
  file.path(sensitivity_tables_dir, "sensitivity_scenario_definitions.csv")
)

# ------------------------------------------------------------
# 3. Read saved historical-model estimates
# ------------------------------------------------------------

gamma_tbl_saved <- read_csv(model_gamma_fp, show_col_types = FALSE)

if (!all(c("term", "estimate") %in% names(gamma_tbl_saved))) {
  stop("The saved gamma file does not contain the expected term and estimate columns.")
}

gamma_hat <- gamma_tbl_saved |>
  filter(term == "log_so2") |>
  pull(estimate)

if (length(gamma_hat) != 1 || !is.finite(gamma_hat)) {
  stop("Could not identify one finite saved log_so2 estimate in ", model_gamma_fp)
}

alpha_table <- read_fst(model_alpha_fp) |>
  transmute(
    plant_id = as.integer(plant_id),
    zcta = standardize_zcta(zcta),
    alpha_hat = as.numeric(alpha_hat)
  ) |>
  filter(
    !is.na(plant_id),
    !is.na(zcta),
    is.finite(alpha_hat)
  )

smear_by_plant <- read_csv(smearing_fp, show_col_types = FALSE) |>
  mutate(plant_id = as.integer(plant_id), smear_duan = as.numeric(smear_duan))

message("Using saved SO2 elasticity gamma = ", signif(gamma_hat, 6))
message("Using ", format(nrow(alpha_table), big.mark = ","), " saved plant-ZCTA fixed effects.")

# ------------------------------------------------------------
# 4. Read SO2 data and construct plant-specific baselines
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
  stop("No SO2 emissions column found in ", so2_fp)
}

so2_panel <- so2_raw |>
  transmute(
    plant_id = as.integer(facility_id),
    year = as.integer(year),
    so2_tons = as.numeric(.data[[so2_col]])
  ) |>
  filter(plant_id %in% target_plants$plant_id)

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

# Plant baselines are sums over coal generators of the unit-level CAMPD
# values, so the recent-low and recent-high scenarios are now unit-specific.
so2_baselines <- gen_tbl |>
  group_by(plant_id) |>
  summarise(
    n_baseline_years = length(so2_baseline_years),
    so2_mean_baseline = sum(so2_mean_2023_2025),
    so2_min_baseline = sum(so2_min_2023_2025),
    so2_max_baseline = sum(so2_max_2023_2025),
    has_zero_baseline_year = any(so2_min_2023_2025 == 0),
    .groups = "drop"
  ) |>
  right_join(
    target_plants |>
      transmute(
        plant_id = as.integer(plant_id),
        plant_name_target = as.character(plant_name_target)
      ),
    by = "plant_id"
  ) |>
  mutate(
    recent_low_to_mean_ratio = if_else(
      is.finite(so2_mean_baseline) & so2_mean_baseline > 0,
      so2_min_baseline / so2_mean_baseline, NA_real_
    ),
    recent_high_to_mean_ratio = if_else(
      is.finite(so2_mean_baseline) & so2_mean_baseline > 0,
      so2_max_baseline / so2_mean_baseline, NA_real_
    )
  )

write_csv(
  so2_baselines,
  file.path(sensitivity_diagnostics_dir, "plant_so2_baseline_qc_2023_2025.csv")
)

if (any(is.na(so2_baselines$so2_mean_baseline))) {
  stop(
    "At least one target plant has no finite SO2 observation during ",
    min(so2_baseline_years), "-", max(so2_baseline_years),
    ". See plant_so2_baseline_qc_2023_2025.csv."
  )
}

if (any(so2_baselines$n_baseline_years < length(so2_baseline_years))) {
  warning(
    "At least one target plant has fewer than ", length(so2_baseline_years),
    " SO2 baseline years. Review plant_so2_baseline_qc_2023_2025.csv."
  )
}

if (any(so2_baselines$has_zero_baseline_year, na.rm = TRUE)) {
  warning(
    "At least one plant has a zero SO2 baseline year. The recent-low scenario will use zero for that plant. " ,
    "Confirm that zero is a valid reported value rather than missing or incomplete reporting."
  )
}

if (any(!is.finite(so2_baselines$so2_max_baseline))) {
  stop(
    "At least one target plant has no finite maximum SO2 baseline during ",
    min(so2_baseline_years), "-", max(so2_baseline_years),
    ". See plant_so2_baseline_qc_2023_2025.csv."
  )
}

# ------------------------------------------------------------
# 5. Read retirement scenarios
# ------------------------------------------------------------

# Plant-level representative dates, kept only for the diagnostic columns
# written out below. All emissions arithmetic is done per generator.
scenario_tbl <- gen_tbl |>
  mutate(
    affected = indefinitely_delayed == 1L |
      (!is.na(modified_retirement_date) & !is.na(original_retirement_date) &
         modified_retirement_date > original_retirement_date)
  ) |>
  filter(affected) |>
  arrange(plant_id, original_retirement_date) |>
  group_by(plant_id) |>
  summarise(
    plant_name = dplyr::first(plant_name),
    original_retirement_date = dplyr::first(original_retirement_date),
    modified_retirement_date = dplyr::first(modified_retirement_date),
    indefinitely_delayed = as.integer(dplyr::first(indefinitely_delayed)),
    .groups = "drop"
  )

if (n_distinct(scenario_tbl$plant_id) != nrow(target_plants)) {
  stop(
    "The retirement-scenario file does not contain exactly one row for each of the ",
    nrow(target_plants), " target plants."
  )
}

# ------------------------------------------------------------
# 6. Read Medicare population and baseline mortality inputs
# ------------------------------------------------------------

medicare_nat <- read_csv(medicare_fp, show_col_types = FALSE) |>
  transmute(
    year = as.integer(year),
    total_medicare_enrollment = as.numeric(total_medicare_enrollment)
  )

future_medicare_total <- medicare_nat |>
  filter(year %in% years_future)

if (nrow(future_medicare_total) != length(years_future)) {
  stop("Projected Medicare enrollment is missing at least one future year.")
}

demo <- read_csv(demo_fp, show_col_types = FALSE) |>
  transmute(
    zcta = standardize_zcta(zcta),
    n_bene = as.numeric(n_bene),
    mortality_rate_zcta = as.numeric(death_rate) * baseline_rate_multiplier
  ) |>
  filter(
    !is.na(zcta),
    is.finite(n_bene),
    n_bene > 0,
    is.finite(mortality_rate_zcta),
    mortality_rate_zcta >= 0
  )

zcta_inputs <- demo |>
  mutate(
    zcta_medicare_share = n_bene / sum(n_bene, na.rm = TRUE)
  )

future_zcta_medicare <- future_medicare_total |>
  crossing(zcta_inputs) |>
  mutate(
    medicare_population_future =
      total_medicare_enrollment * zcta_medicare_share
  ) |>
  select(
    year,
    zcta,
    medicare_population_future,
    mortality_rate_zcta
  )

# ------------------------------------------------------------
# 7. Functions to build and evaluate one sensitivity scenario
# ------------------------------------------------------------

build_future_so2 <- function(spec_row) {
  sensitivity_id <- spec_row$sensitivity_id[[1]]
  baseline_type <- spec_row$baseline_type[[1]]
  mean_multiplier <- spec_row$mean_multiplier[[1]]
  removed_date_rule <- spec_row$removed_date_rule[[1]]
  removed_date_shift_years <- spec_row$removed_date_shift_years[[1]]

  if (!removed_date_rule %in% c("original_plus_years", "through_2035")) {
    stop("Unknown removed_date_rule: ", removed_date_rule)
  }
  if (!baseline_type %in% c("mean", "min", "max")) {
    stop("Unknown baseline_type: ", baseline_type)
  }

  # Alternative retirement dates apply per generator, and only to generators
  # whose date was removed. Generators that kept a reported date keep it
  # (Craig units 2-3, Sherburne unit 3, Dave Johnston unit 3).
  gen_used <- gen_tbl |>
    mutate(
      so2_unit_used = dplyr::case_when(
        baseline_type == "mean" ~ so2_mean_2023_2025 * mean_multiplier,
        baseline_type == "min"  ~ so2_min_2023_2025,
        baseline_type == "max"  ~ so2_max_2023_2025
      ),
      modified_retirement_date_used = if_else(
        indefinitely_delayed == 1L & removed_date_rule == "original_plus_years",
        original_retirement_date %m+% lubridate::years(removed_date_shift_years),
        modified_retirement_date
      ),
      modified_indefinitely_delayed_used = if_else(
        indefinitely_delayed == 1L & removed_date_rule == "original_plus_years",
        0L, indefinitely_delayed
      )
    )

  expand_grid(
    gen_used,
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
        retirement_date = modified_retirement_date_used,
        indefinitely_delayed = modified_indefinitely_delayed_used
      )
    ) |>
    ungroup() |>
    group_by(plant_id, year) |>
    summarise(
      so2_baseline_used = sum(so2_unit_used),
      so2_original = sum(so2_unit_used * frac_original),
      so2_modified = sum(so2_unit_used * frac_modified),
      operating_fraction_original =
        sum(so2_unit_used * frac_original) / pmax(sum(so2_unit_used), eps),
      operating_fraction_modified =
        sum(so2_unit_used * frac_modified) / pmax(sum(so2_unit_used), eps),
      .groups = "drop"
    ) |>
    left_join(scenario_tbl, by = "plant_id") |>
    left_join(
      so2_baselines |>
        select(plant_id, n_baseline_years, so2_mean_baseline, so2_min_baseline,
               so2_max_baseline, has_zero_baseline_year,
               recent_low_to_mean_ratio, recent_high_to_mean_ratio),
      by = "plant_id"
    ) |>
    mutate(
      sensitivity_id = sensitivity_id,
      sensitivity_label = spec_row$sensitivity_label[[1]],
      baseline_type = baseline_type,
      mean_multiplier = mean_multiplier,
      removed_date_rule = removed_date_rule,
      removed_date_shift_years = removed_date_shift_years,
      modified_retirement_date_used = if_else(
        indefinitely_delayed == 1L & removed_date_rule == "original_plus_years",
        original_retirement_date %m+% lubridate::years(removed_date_shift_years),
        modified_retirement_date
      ),
      modified_indefinitely_delayed_used = if_else(
        indefinitely_delayed == 1L & removed_date_rule == "original_plus_years",
        0L, indefinitely_delayed
      )
    ) |>
    select(
      sensitivity_id, sensitivity_label, plant_id, plant_name, year,
      baseline_type, mean_multiplier, removed_date_rule,
      removed_date_shift_years, n_baseline_years, so2_mean_baseline,
      so2_min_baseline, so2_max_baseline, so2_baseline_used,
      recent_low_to_mean_ratio, recent_high_to_mean_ratio,
      has_zero_baseline_year, original_retirement_date,
      modified_retirement_date, modified_retirement_date_used,
      indefinitely_delayed, modified_indefinitely_delayed_used,
      operating_fraction_original, operating_fraction_modified,
      so2_original, so2_modified
    )
}

predict_zcta_pm25 <- function(future_so2_one) {
  sensitivity_id <- unique(future_so2_one$sensitivity_id)
  
  if (length(sensitivity_id) != 1) {
    stop("predict_zcta_pm25() requires exactly one sensitivity_id.")
  }
  
  out <- vector("list", length = 2 * length(years_future))
  counter <- 1L
  
  for (scenario_name in c("original", "modified")) {
    so2_var <- paste0("so2_", scenario_name)
    
    for (yy in years_future) {
      plant_year_so2 <- future_so2_one |>
        filter(year == yy) |>
        transmute(
          plant_id,
          so2_future = as.numeric(.data[[so2_var]])
        )
      
      zcta_one <- alpha_table |>
        inner_join(plant_year_so2, by = "plant_id") |>
        left_join(smear_by_plant, by = "plant_id") |>
        mutate(
          smear_duan = tidyr::replace_na(smear_duan, 1),
          pred_coal_pm25_positive = exp(
            alpha_hat + gamma_hat * log(
              pmax(so2_future, .Machine$double.xmin)
            )
          ) * smear_duan,
          pred_coal_pm25 = if_else(
            is.finite(so2_future) & so2_future > 0,
            pred_coal_pm25_positive,
            0
          )
        ) |>
        group_by(zcta) |>
        summarise(
          coal_pm25_zcta = sum(pred_coal_pm25, na.rm = TRUE),
          .groups = "drop"
        ) |>
        mutate(
          sensitivity_id = sensitivity_id,
          scenario = scenario_name,
          year = yy
        ) |>
        select(
          sensitivity_id,
          scenario,
          year,
          zcta,
          coal_pm25_zcta
        )
      
      out[[counter]] <- zcta_one
      counter <- counter + 1L
    }
  }
  
  bind_rows(out)
}

estimate_health_impacts <- function(zcta_pm25_one, future_so2_one, spec_row) {
  sensitivity_id <- spec_row$sensitivity_id[[1]]
  sensitivity_label <- spec_row$sensitivity_label[[1]]
  
  health_base <- future_zcta_medicare |>
    crossing(scenario = c("original", "modified")) |>
    left_join(
      zcta_pm25_one |>
        select(-sensitivity_id),
      by = c("scenario", "year", "zcta")
    ) |>
    mutate(
      # Assign the current sensitivity ID to every ZCTA, including
      # ZCTAs with no matched modeled PM2.5 value. Those ZCTAs receive
      # zero exposure below but must remain in the same population
      # denominator rather than forming a separate NA group.
      sensitivity_id = .env$sensitivity_id,
      sensitivity_label = .env$sensitivity_label,
      coal_pm25_zcta = replace_na(coal_pm25_zcta, 0)
    )
  
  annual_exposure <- health_base |>
    group_by(
      sensitivity_id,
      sensitivity_label,
      scenario,
      year
    ) |>
    summarise(
      medicare_population = sum(medicare_population_future, na.rm = TRUE),
      population_weighted_coal_pm25 =
        sum(
          medicare_population_future * coal_pm25_zcta,
          na.rm = TRUE
        ) /
        sum(medicare_population_future, na.rm = TRUE),
      .groups = "drop"
    )
  
  annual_deaths <- health_base |>
    crossing(beta_tbl) |>
    mutate(
      deaths =
        medicare_population_future *
        mortality_rate_zcta *
        (1 - exp(-beta_E * coal_pm25_zcta))
    ) |>
    group_by(
      sensitivity_id,
      sensitivity_label,
      beta_type,
      scenario,
      year
    ) |>
    summarise(
      total_deaths = sum(deaths, na.rm = TRUE),
      .groups = "drop"
    )
  
  cumulative_deaths <- annual_deaths |>
    group_by(
      sensitivity_id,
      sensitivity_label,
      beta_type,
      scenario
    ) |>
    summarise(
      cumulative_deaths = sum(total_deaths, na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = scenario,
      values_from = cumulative_deaths,
      values_fill = 0
    ) |>
    mutate(additional_deaths = modified - original)
  
  death_summary <- cumulative_deaths |>
    select(
      sensitivity_id,
      sensitivity_label,
      beta_type,
      original,
      modified,
      additional_deaths
    ) |>
    pivot_wider(
      names_from = beta_type,
      values_from = c(original, modified, additional_deaths)
    ) |>
    transmute(
      sensitivity_id,
      sensitivity_label,
      original_deaths_central = original_central,
      original_deaths_low = original_lower,
      original_deaths_high = original_upper,
      revised_deaths_central = modified_central,
      revised_deaths_low = modified_lower,
      revised_deaths_high = modified_upper,
      additional_deaths_central = additional_deaths_central,
      additional_deaths_low = additional_deaths_lower,
      additional_deaths_high = additional_deaths_upper
    )
  
  exposure_summary <- annual_exposure |>
    group_by(
      sensitivity_id,
      sensitivity_label,
      scenario
    ) |>
    summarise(
      cumulative_population_weighted_coal_pm25 =
        sum(population_weighted_coal_pm25, na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_wider(
      names_from = scenario,
      values_from = cumulative_population_weighted_coal_pm25,
      values_fill = 0
    ) |>
    transmute(
      sensitivity_id,
      sensitivity_label,
      cumulative_pw_coal_pm25_original = original,
      cumulative_pw_coal_pm25_revised = modified,
      additional_cumulative_pw_coal_pm25 = modified - original
    )
  
  so2_summary <- future_so2_one |>
    summarise(
      cumulative_so2_original = sum(so2_original, na.rm = TRUE),
      cumulative_so2_revised = sum(so2_modified, na.rm = TRUE),
      additional_cumulative_so2 =
        cumulative_so2_revised - cumulative_so2_original
    ) |>
    mutate(
      sensitivity_id = sensitivity_id,
      sensitivity_label = sensitivity_label,
      .before = 1
    )
  
  summary_one <- death_summary |>
    left_join(
      exposure_summary,
      by = c("sensitivity_id", "sensitivity_label")
    ) |>
    left_join(
      so2_summary,
      by = c("sensitivity_id", "sensitivity_label")
    ) |>
    left_join(
      sensitivity_specs,
      by = c("sensitivity_id", "sensitivity_label")
    )
  
  list(
    summary = summary_one,
    annual_deaths = annual_deaths,
    annual_exposure = annual_exposure
  )
}

run_one_sensitivity <- function(spec_row) {
  sensitivity_id <- spec_row$sensitivity_id[[1]]
  sensitivity_label <- spec_row$sensitivity_label[[1]]
  
  message("Running sensitivity: ", sensitivity_label, " [", sensitivity_id, "]")
  
  future_so2_one <- build_future_so2(spec_row)
  zcta_pm25_one <- predict_zcta_pm25(future_so2_one)
  health_one <- estimate_health_impacts(
    zcta_pm25_one = zcta_pm25_one,
    future_so2_one = future_so2_one,
    spec_row = spec_row
  )
  
  list(
    future_so2 = future_so2_one,
    zcta_pm25 = zcta_pm25_one,
    summary = health_one$summary,
    annual_deaths = health_one$annual_deaths,
    annual_exposure = health_one$annual_exposure
  )
}

# ------------------------------------------------------------
# 8. Run all sensitivity scenarios sequentially
# ------------------------------------------------------------
# Processing sequentially avoids creating one extremely large
# plant-ZCTA prediction table for all scenarios at once.

sensitivity_results <- vector("list", nrow(sensitivity_specs))

for (i in seq_len(nrow(sensitivity_specs))) {
  sensitivity_results[[i]] <- run_one_sensitivity(
    sensitivity_specs[i, , drop = FALSE]
  )
}

# ------------------------------------------------------------
# 9. Combine and save outputs
# ------------------------------------------------------------

sensitivity_summary <- map_dfr(sensitivity_results, "summary")
sensitivity_annual_deaths <- map_dfr(sensitivity_results, "annual_deaths")
sensitivity_annual_exposure <- map_dfr(sensitivity_results, "annual_exposure")
sensitivity_future_so2 <- map_dfr(sensitivity_results, "future_so2")

# Each sensitivity must yield exactly one summary row.
summary_qc <- sensitivity_summary |>
  count(sensitivity_id, sensitivity_label, name = "n_summary_rows")

if (
  nrow(sensitivity_summary) != nrow(sensitivity_specs) ||
  any(is.na(sensitivity_summary$sensitivity_id)) ||
  any(summary_qc$n_summary_rows != 1)
) {
  print(summary_qc)
  stop(
    "Sensitivity summary did not produce exactly one nonmissing row per scenario."
  )
}

main_additional_deaths <- sensitivity_summary |>
  filter(sensitivity_id == "main") |>
  pull(additional_deaths_central)

main_additional_exposure <- sensitivity_summary |>
  filter(sensitivity_id == "main") |>
  pull(additional_cumulative_pw_coal_pm25)

if (
  length(main_additional_deaths) != 1 ||
  !is.finite(main_additional_deaths) ||
  main_additional_deaths == 0
) {
  stop(
    "Could not identify one finite, nonzero main-analysis ",
    "additional-deaths estimate."
  )
}

if (
  length(main_additional_exposure) != 1 ||
  !is.finite(main_additional_exposure) ||
  main_additional_exposure == 0
) {
  stop(
    "Could not identify one finite, nonzero main-analysis ",
    "additional-exposure estimate."
  )
}

# The main denominators are validated above as finite, scalar, and nonzero.
# Use direct vectorized arithmetic here. Do not use dplyr::if_else() with
# a length-one condition and a multirow result vector.
sensitivity_summary <- sensitivity_summary |>
  mutate(
    scenario_order = match(
      sensitivity_id,
      sensitivity_specs$sensitivity_id
    ),
    absolute_change_from_main_deaths =
      additional_deaths_central - main_additional_deaths,
    percent_change_from_main_deaths =
      100 * absolute_change_from_main_deaths /
      main_additional_deaths,
    absolute_change_from_main_exposure =
      additional_cumulative_pw_coal_pm25 - main_additional_exposure,
    percent_change_from_main_exposure =
      100 * absolute_change_from_main_exposure /
      main_additional_exposure
  ) |>
  arrange(scenario_order)

if (nrow(sensitivity_summary) != 9L) {
  stop(
    "Expected 9 sensitivity-summary rows, but found ",
    nrow(sensitivity_summary), "."
  )
}

if (anyDuplicated(sensitivity_summary$sensitivity_id) > 0) {
  stop("Duplicate sensitivity_id values remain in the summary output.")
}

sensitivity_summary_formatted <- sensitivity_summary |>
  transmute(
    sensitivity_id,
    sensitivity_label,
    emissions_assumption = case_when(
      sensitivity_id == "main" ~ "100% of plant-specific 2023-2025 mean SO2",
      sensitivity_id == "uniform_75" ~ "75% of plant-specific 2023-2025 mean SO2",
      sensitivity_id == "uniform_50" ~ "50% of plant-specific 2023-2025 mean SO2",
      sensitivity_id == "uniform_125" ~ "125% of plant-specific 2023-2025 mean SO2",
      sensitivity_id == "uniform_150" ~ "150% of plant-specific 2023-2025 mean SO2",
      sensitivity_id == "recent_low" ~ "Plant-specific minimum annual SO2 during 2023-2025",
      sensitivity_id == "recent_high" ~ "Plant-specific maximum annual SO2 during 2023-2025",
      TRUE ~ "100% of plant-specific 2023-2025 mean SO2"
    ),
    removed_date_assumption = case_when(
      sensitivity_id %in% c(
        "main",
        "uniform_75",
        "uniform_50",
        "uniform_125",
        "uniform_150",
        "recent_low",
        "recent_high"
      ) ~ "Continued operation through 2035",
      sensitivity_id == "removed_plus_2" ~ "Original retirement date plus 2 years",
      sensitivity_id == "removed_plus_5" ~ "Original retirement date plus 5 years",
      TRUE ~ NA_character_
    ),
    additional_deaths_central,
    additional_deaths_low,
    additional_deaths_high,
    additional_deaths_display = paste0(
      number(additional_deaths_central, accuracy = 0.1),
      " (",
      number(additional_deaths_low, accuracy = 0.1),
      "-",
      number(additional_deaths_high, accuracy = 0.1),
      ")"
    ),
    percent_change_from_main_deaths,
    additional_cumulative_pw_coal_pm25,
    percent_change_from_main_exposure,
    cumulative_so2_original,
    cumulative_so2_revised,
    additional_cumulative_so2
  )

write_csv(
  sensitivity_summary,
  file.path(sensitivity_tables_dir, "sensitivity_summary_full.csv")
)

write_csv(
  sensitivity_summary_formatted,
  file.path(sensitivity_tables_dir, "table_sensitivity_summary.csv")
)

write_csv(
  sensitivity_annual_deaths,
  file.path(sensitivity_tables_dir, "sensitivity_annual_deaths_by_scenario.csv")
)

write_csv(
  sensitivity_annual_exposure,
  file.path(sensitivity_tables_dir, "sensitivity_annual_population_weighted_pm25.csv")
)

write_csv(
  sensitivity_future_so2,
  file.path(sensitivity_tables_dir, "sensitivity_future_so2_by_plant_year.csv")
)

plant_assumption_qc <- sensitivity_future_so2 |>
  distinct(
    sensitivity_id,
    sensitivity_label,
    plant_id,
    plant_name,
    baseline_type,
    mean_multiplier,
    n_baseline_years,
    so2_mean_baseline,
    so2_min_baseline,
    so2_max_baseline,
    so2_baseline_used,
    recent_low_to_mean_ratio,
    recent_high_to_mean_ratio,
    has_zero_baseline_year,
    original_retirement_date,
    modified_retirement_date,
    modified_retirement_date_used,
    indefinitely_delayed,
    modified_indefinitely_delayed_used
  ) |>
  arrange(
    match(sensitivity_id, sensitivity_specs$sensitivity_id),
    plant_id
  )

write_csv(
  plant_assumption_qc,
  file.path(sensitivity_diagnostics_dir, "sensitivity_assumptions_by_plant.csv")
)

# ------------------------------------------------------------
# 10. Check whether the recreated main result matches script 03
# ------------------------------------------------------------

main_reference_fp <- file.path(
  dir_tables,
  "total_cumulative_additional_deaths_ci_2026_2035.csv"
)

if (file.exists(main_reference_fp)) {
  main_reference <- read_csv(
    main_reference_fp,
    show_col_types = FALSE
  )
  
  if ("additional_deaths_central" %in% names(main_reference)) {
    reproducibility_check <- tibble(
      stored_main_additional_deaths =
        as.numeric(main_reference$additional_deaths_central[1]),
      recreated_main_additional_deaths = main_additional_deaths,
      absolute_difference =
        recreated_main_additional_deaths -
        stored_main_additional_deaths,
      relative_difference_pct = if_else(
        stored_main_additional_deaths != 0,
        100 * absolute_difference /
          stored_main_additional_deaths,
        NA_real_
      )
    )
    
    write_csv(
      reproducibility_check,
      file.path(
        sensitivity_diagnostics_dir,
        "main_result_reproducibility_check.csv"
      )
    )
    
    print(reproducibility_check)
    
    if (abs(reproducibility_check$relative_difference_pct) > 0.1) {
      warning(
        "The recreated main result differs from the stored main result by more than 0.1%. " ,
        "Review main_result_reproducibility_check.csv before interpreting sensitivities."
      )
    }
  }
}

# ------------------------------------------------------------
# 11. Make an SI-ready sensitivity figure
# ------------------------------------------------------------
# The lower and upper values reflect only the lower and upper
# coal-PM2.5 mortality-response coefficients, matching script 03.

plot_color <- if (exists("color_modified")) {
  color_modified
} else {
  "#C51B3A"
}

sensitivity_plot_data <- sensitivity_summary |>
  mutate(
    sensitivity_label = factor(
      sensitivity_label,
      levels = rev(sensitivity_specs$sensitivity_label)
    )
  )

p_sensitivity_deaths <- ggplot(
  sensitivity_plot_data,
  aes(
    x = additional_deaths_central,
    y = sensitivity_label
  )
) +
  geom_vline(
    xintercept = main_additional_deaths,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey45"
  ) +
  geom_errorbar(
    aes(
      xmin = additional_deaths_low,
      xmax = additional_deaths_high
    ),
    width = 0.16,
    linewidth = 0.55,
    orientation = "y",
    color = "black"
  ) +
  geom_point(
    size = 2.7,
    color = plot_color
  ) +
  scale_x_continuous(
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    title = "Sensitivity of cumulative additional Medicare deaths",
    subtitle = "Revised retirement scenario minus original scenario, 2026-2035",
    x = "Cumulative additional deaths",
    y = NULL,
    caption = paste0(
      "Points show central estimates; horizontal ranges use lower and upper coal-PM2.5 mortality-response coefficients. ",
      "The dashed line marks the main estimate."
    )
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, color = "grey35", hjust = 0),
    plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(
    sensitivity_figures_dir,
    "figure_sensitivity_additional_medicare_deaths.png"
  ),
  p_sensitivity_deaths,
  width = 7.25,
  height = 5.2,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 12. Print concise result table
# ------------------------------------------------------------

print(
  sensitivity_summary_formatted |>
    select(
      sensitivity_label,
      additional_deaths_display,
      percent_change_from_main_deaths,
      percent_change_from_main_exposure
    )
)

message(
  "Sensitivity analysis complete. Outputs saved under:\n",
  sensitivity_dir
)
