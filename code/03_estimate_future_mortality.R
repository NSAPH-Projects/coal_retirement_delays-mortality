# ============================================================
# 03_estimate_future_mortality.R
# Estimate future exposure burden and Medicare deaths by scenario
# ============================================================

source("R/00_config.R")
source("R/00_functions.R")

load_required_packages(c(
  "fst", "dplyr", "tidyr", "readr", "stringr", "tibble", 
  "ggplot2", "scales", "patchwork", "cowplot"
))

# ------------------------------------------------------------
# 1. Read projected national Medicare enrollment
# ------------------------------------------------------------
medicare_nat <- read_csv(medicare_fp, show_col_types = FALSE) |>
  mutate(
    year = as.integer(year),
    total_medicare_enrollment = as.numeric(total_medicare_enrollment)
  )

future_medicare_total <- medicare_nat |>
  filter(year %in% years_future) |>
  transmute(year, total_medicare_enrollment)

if (nrow(future_medicare_total) != length(years_future)) {
  stop("Projected Medicare enrollment is missing at least one future year.")
}

# ------------------------------------------------------------
# 2. Read demographics and construct future ZCTA Medicare population
# ------------------------------------------------------------
demo <- read_csv(demo_fp, show_col_types = FALSE) |>
  mutate(
    zcta = standardize_zcta(zcta),
    population = as.numeric(population),
    n_bene = as.numeric(n_bene),
    death_rate = as.numeric(death_rate),
    female_pct = as.numeric(female_pct),
    male_pct = as.numeric(male_pct),
    dual_pct = as.numeric(dual_pct),
    pop_under_20_years = as.numeric(pop_under_20_years),
    pop_35_49_years = as.numeric(pop_35_49_years),
    pop_50_64_years = as.numeric(pop_50_64_years),
    pop_over_65_years = as.numeric(pop_over_65_years),
    pop_white = as.numeric(pop_white),
    pop_black = as.numeric(pop_black),
    pop_asian = as.numeric(pop_asian),
    pop_native = as.numeric(pop_native),
    pop_hispanic = as.numeric(pop_hispanic)
  ) |>
  filter(
    !is.na(zcta),
    !is.na(n_bene), n_bene > 0,
    !is.na(death_rate), death_rate >= 0
  )

zcta_inputs <- demo |>
  mutate(
    zcta_medicare_share = n_bene / sum(n_bene, na.rm = TRUE),
    mortality_rate_zcta = death_rate
  ) |>
  select(zcta, n_bene, zcta_medicare_share, mortality_rate_zcta, everything())

future_zcta_medicare <- future_medicare_total |>
  crossing(zcta_inputs) |>
  mutate(medicare_population_future = total_medicare_enrollment * zcta_medicare_share)

write_csv(
  future_zcta_medicare |>
    select(year, zcta, medicare_population_future, mortality_rate_zcta),
  file.path(dir_tables, "future_zcta_medicare_population_and_mortality.csv")
)

# ------------------------------------------------------------
# 3. Optional subgroup exposure summaries under future scenarios
# ------------------------------------------------------------
subgroup_share <- build_medicare_subgroup_shares(zcta_inputs)

subgroup_share_check <- subgroup_share |>
  group_by(group_type, zcta) |>
  summarise(share_sum = sum(subgroup_share, na.rm = TRUE), .groups = "drop") |>
  group_by(group_type) |>
  summarise(
    min_share_sum = min(share_sum, na.rm = TRUE),
    median_share_sum = median(share_sum, na.rm = TRUE),
    max_share_sum = max(share_sum, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(subgroup_share_check, file.path(dir_diagnostics, "future_subgroup_share_qc.csv"))

future_subgroup_pop <- future_zcta_medicare |>
  select(year, zcta, medicare_population_future) |>
  left_join(subgroup_share, by = "zcta", relationship = "many-to-many") |>
  mutate(subgroup_population_future = medicare_population_future * subgroup_share) |>
  filter(!is.na(group_type), !is.na(group), !is.na(subgroup_population_future), subgroup_population_future > 0)

# ------------------------------------------------------------
# 4. Read future plant-ZCTA coal PM2.5 predictions
# ------------------------------------------------------------
plant_pm25 <- read_fst(future_pm25_fp) |>
  mutate(
    scenario = str_to_lower(as.character(scenario)),
    year = as.integer(year),
    plant_id = as.integer(plant_id),
    plant_name = as.character(plant_name),
    zcta = standardize_zcta(zcta),
    pred_coal_pm25 = as.numeric(pred_coal_pm25)
  ) |>
  filter(year %in% years_future, scenario %in% c("original", "modified")) |>
  select(scenario, year, plant_id, plant_name, zcta, pred_coal_pm25)

zcta_pm25_scenario <- plant_pm25 |>
  group_by(scenario, year, zcta) |>
  summarise(coal_pm25_zcta = sum(pred_coal_pm25, na.rm = TRUE), .groups = "drop")

write_csv(zcta_pm25_scenario, file.path(dir_tables, "zcta_year_coal_pm25_by_scenario.csv"))

# ------------------------------------------------------------
# 5. Subgroup population-weighted future exposure summaries
# ------------------------------------------------------------
subgroup_exposure_annual <- future_subgroup_pop |>
  left_join(zcta_pm25_scenario, by = c("year", "zcta")) |>
  mutate(
    coal_pm25_zcta = tidyr::replace_na(coal_pm25_zcta, 0),
    weighted_coal_pm25 = subgroup_population_future * coal_pm25_zcta
  ) |>
  group_by(group_type, group, scenario, year) |>
  summarise(
    subgroup_population = sum(subgroup_population_future, na.rm = TRUE),
    population_weighted_coal_pm25 =
      sum(weighted_coal_pm25, na.rm = TRUE) / sum(subgroup_population_future, na.rm = TRUE),
    .groups = "drop"
  )

subgroup_exposure_summary <- subgroup_exposure_annual |>
  group_by(group_type, group, scenario) |>
  summarise(
    beneficiary_years = sum(subgroup_population, na.rm = TRUE),
    mean_annual_pw_coal_pm25 = weighted.mean(
      population_weighted_coal_pm25,
      w = subgroup_population,
      na.rm = TRUE
    ),
    cumulative_pw_coal_pm25 = sum(population_weighted_coal_pm25, na.rm = TRUE),
    .groups = "drop"
  )

subgroup_effect_summary <- subgroup_exposure_summary |>
  select(group_type, group, scenario, cumulative_pw_coal_pm25) |>
  pivot_wider(names_from = scenario, values_from = cumulative_pw_coal_pm25, values_fill = 0) |>
  mutate(
    additional = modified - original,
    pct_increase = 100 * additional / original
  )

write_csv(subgroup_exposure_annual, file.path(dir_tables, "future_subgroup_exposure_annual_by_scenario.csv"))
write_csv(subgroup_exposure_summary, file.path(dir_tables, "future_subgroup_exposure_summary_by_scenario.csv"))
write_csv(subgroup_effect_summary, file.path(dir_tables, "future_subgroup_delayed_effect_summary.csv"))

# ------------------------------------------------------------
# 6. Estimate deaths at ZCTA-year-scenario level and allocate to plants
# ------------------------------------------------------------
zcta_total_pm25 <- plant_pm25 |>
  group_by(scenario, year, zcta) |>
  summarise(total_coal_pm25_zcta = sum(pred_coal_pm25, na.rm = TRUE), .groups = "drop")

plant_zcta_death_input <- plant_pm25 |>
  left_join(zcta_total_pm25, by = c("scenario", "year", "zcta")) |>
  left_join(
    future_zcta_medicare |>
      select(year, zcta, medicare_population_future, mortality_rate_zcta),
    by = c("year", "zcta")
  ) |>
  mutate(
    medicare_population_future = tidyr::replace_na(medicare_population_future, 0),
    mortality_rate_zcta = tidyr::replace_na(mortality_rate_zcta, 0),
    pred_coal_pm25 = tidyr::replace_na(pred_coal_pm25, 0),
    total_coal_pm25_zcta = tidyr::replace_na(total_coal_pm25_zcta, 0),
    plant_pm25_share = if_else(total_coal_pm25_zcta > 0, pred_coal_pm25 / total_coal_pm25_zcta, 0)
  ) |>
  crossing(beta_tbl) |>
  mutate(
    deaths_standalone =
      medicare_population_future * mortality_rate_zcta *
      (1 - exp(-beta_E * pred_coal_pm25)),
    deaths_total_zcta =
      medicare_population_future * mortality_rate_zcta *
      (1 - exp(-beta_E * total_coal_pm25_zcta)),
    deaths_allocated = deaths_total_zcta * plant_pm25_share
  )

# Scenario-year total deaths use distinct ZCTA totals to avoid plant-row duplication.
scenario_year_deaths_total <- plant_zcta_death_input |>
  distinct(beta_type, scenario, year, zcta, deaths_total_zcta) |>
  group_by(beta_type, scenario, year) |>
  summarise(total_deaths = sum(deaths_total_zcta, na.rm = TRUE), .groups = "drop")

scenario_year_deaths_diff <- scenario_year_deaths_total |>
  pivot_wider(names_from = scenario, values_from = total_deaths, values_fill = 0) |>
  mutate(additional_deaths = modified - original) |>
  arrange(beta_type, year)

write_csv(scenario_year_deaths_total, file.path(dir_tables, "scenario_year_total_deaths.csv"))
write_csv(scenario_year_deaths_diff, file.path(dir_tables, "scenario_year_total_deaths_modified_minus_original.csv"))

# Plant-year deaths. Allocated deaths sum back to the scenario-level total.
plant_year_deaths_by_scenario <- plant_zcta_death_input |>
  group_by(beta_type, scenario, year, plant_id, plant_name) |>
  summarise(
    plant_deaths_allocated = sum(deaths_allocated, na.rm = TRUE),
    plant_deaths_standalone = sum(deaths_standalone, na.rm = TRUE),
    plant_population_weighted_pm25 =
      sum(medicare_population_future * pred_coal_pm25, na.rm = TRUE) /
      sum(medicare_population_future, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(beta_type, scenario, year, desc(plant_deaths_allocated))

plant_year_deaths_diff <- plant_year_deaths_by_scenario |>
  select(
    beta_type, year, plant_id, plant_name, scenario,
    plant_deaths_allocated, plant_deaths_standalone, plant_population_weighted_pm25
  ) |>
  pivot_wider(
    names_from = scenario,
    values_from = c(plant_deaths_allocated, plant_deaths_standalone, plant_population_weighted_pm25),
    values_fill = 0
  ) |>
  mutate(
    additional_deaths_allocated = plant_deaths_allocated_modified - plant_deaths_allocated_original,
    additional_deaths_standalone = plant_deaths_standalone_modified - plant_deaths_standalone_original,
    additional_pw_pm25 = plant_population_weighted_pm25_modified - plant_population_weighted_pm25_original
  ) |>
  arrange(beta_type, year, desc(additional_deaths_allocated))

write_csv(plant_year_deaths_by_scenario, file.path(dir_tables, "plant_year_deaths_by_scenario.csv"))
write_csv(plant_year_deaths_diff, file.path(dir_tables, "plant_year_deaths_modified_minus_original.csv"))

# ------------------------------------------------------------
# 7. Cumulative plant-level and total summaries
# ------------------------------------------------------------
plant_cumulative_by_scenario <- plant_year_deaths_by_scenario |>
  group_by(beta_type, scenario, plant_id, plant_name) |>
  summarise(
    cumulative_deaths_allocated = sum(plant_deaths_allocated, na.rm = TRUE),
    cumulative_deaths_standalone = sum(plant_deaths_standalone, na.rm = TRUE),
    mean_annual_population_weighted_pm25 = mean(plant_population_weighted_pm25, na.rm = TRUE),
    cumulative_population_weighted_pm25 = sum(plant_population_weighted_pm25, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(beta_type, scenario, desc(cumulative_deaths_allocated))

plant_cumulative_diff <- plant_cumulative_by_scenario |>
  pivot_wider(
    names_from = scenario,
    values_from = c(
      cumulative_deaths_allocated,
      cumulative_deaths_standalone,
      mean_annual_population_weighted_pm25,
      cumulative_population_weighted_pm25
    ),
    values_fill = 0
  ) |>
  mutate(
    additional_cumulative_deaths_allocated =
      cumulative_deaths_allocated_modified - cumulative_deaths_allocated_original,
    additional_cumulative_deaths_standalone =
      cumulative_deaths_standalone_modified - cumulative_deaths_standalone_original,
    additional_mean_annual_population_weighted_pm25 =
      mean_annual_population_weighted_pm25_modified - mean_annual_population_weighted_pm25_original,
    additional_cumulative_population_weighted_pm25 =
      cumulative_population_weighted_pm25_modified - cumulative_population_weighted_pm25_original
  ) |>
  arrange(beta_type, desc(additional_cumulative_deaths_allocated))

write_csv(plant_cumulative_by_scenario, file.path(dir_tables, "plant_cumulative_deaths_by_scenario_2026_2035.csv"))
write_csv(plant_cumulative_diff, file.path(dir_tables, "plant_cumulative_deaths_modified_minus_original_2026_2035.csv"))

total_additional_deaths_ci <- scenario_year_deaths_total |>
  filter(beta_type %in% c("lower", "central", "upper")) |>
  group_by(beta_type, scenario) |>
  summarise(cumulative_deaths = sum(total_deaths, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = scenario, values_from = cumulative_deaths, values_fill = 0) |>
  mutate(additional_deaths = modified - original) |>
  select(beta_type, original, modified, additional_deaths) |>
  pivot_wider(names_from = beta_type, values_from = c(original, modified, additional_deaths)) |>
  transmute(
    original_central = original_central,
    original_low = original_lower,
    original_high = original_upper,
    modified_central = modified_central,
    modified_low = modified_lower,
    modified_high = modified_upper,
    additional_deaths_central = additional_deaths_central,
    additional_deaths_low = additional_deaths_lower,
    additional_deaths_high = additional_deaths_upper
  )

write_csv(total_additional_deaths_ci, file.path(dir_tables, "total_cumulative_additional_deaths_ci_2026_2035.csv"))
print(total_additional_deaths_ci)

# Plant allocated sum should match scenario-level additional deaths.
plant_allocated_sum <- plant_cumulative_diff |>
  group_by(beta_type) |>
  summarise(
    plant_sum_additional_deaths = sum(additional_cumulative_deaths_allocated, na.rm = TRUE),
    .groups = "drop"
  )

total_vs_plant_check <- scenario_year_deaths_diff |>
  group_by(beta_type) |>
  summarise(scenario_additional_deaths = sum(additional_deaths, na.rm = TRUE), .groups = "drop") |>
  left_join(plant_allocated_sum, by = "beta_type") |>
  mutate(
    difference = plant_sum_additional_deaths - scenario_additional_deaths,
    relative_difference = difference / scenario_additional_deaths
  )

write_csv(total_vs_plant_check, file.path(dir_diagnostics, "scenario_total_vs_plant_allocated_check.csv"))
print(total_vs_plant_check)


# ==================================================
# Optional figure: subgroup delayed-effect plot
# Sex, Race, and Age only
#
# Top panels:
#   cumulative population-weighted coal PM2.5
#   under original vs revised scenarios
#
# Bottom panels:
#   percent increase under revised plan relative to original plan
#
# Layout:
#   Sex | Race | Age
#   One shared legend at the top
# ==================================================

# --------------------------------------------------
# Standardize group labels for plotting
# --------------------------------------------------

subgroup_exposure_summary_plot <- subgroup_exposure_summary %>%
  mutate(
    group_type = recode(group_type, "Race/ethnicity" = "Race"),
    group = recode(group, "Native" = "AI/AN")
  ) %>%
  filter(
    group_type %in% c("Sex", "Race", "Age"),
    !group %in% c("Other")
  )

subgroup_effect_summary_plot <- subgroup_effect_summary %>%
  mutate(
    group_type = recode(group_type, "Race/ethnicity" = "Race"),
    group = recode(group, "Native" = "AI/AN")
  ) %>%
  filter(
    group_type %in% c("Sex", "Race", "Age"),
    !group %in% c("Other")
  )

domain_order <- c("Sex", "Race", "Age")

group_levels_plot <- tibble(
  group_type = c(
    rep("Sex", 2),
    rep("Race", 5),
    rep("Age", 5)
  ),
  group = c(
    "Male", "Female",
    "White", "Black", "Asian", "AI/AN", "Hispanic",
    "Under 20", "20-34", "35-49", "50-64", "65+"
  ),
  group_order = c(
    1, 2,
    1, 2, 3, 4, 5,
    1, 2, 3, 4, 5
  )
) %>%
  semi_join(
    bind_rows(
      subgroup_exposure_summary_plot %>% select(group_type, group),
      subgroup_effect_summary_plot %>% select(group_type, group)
    ) %>%
      distinct(),
    by = c("group_type", "group")
  )

# Shared y-axis limits across panels
top_y_max <- subgroup_exposure_summary_plot %>%
  filter(
    group_type %in% c("Sex", "Race", "Age"),
    scenario %in% c("original", "modified")
  ) %>%
  summarise(
    y = max(cumulative_pw_coal_pm25, na.rm = TRUE)
  ) %>%
  pull(y)

if (!is.finite(top_y_max) || top_y_max <= 0) {
  top_y_max <- 1
}

top_y_max <- top_y_max * 1.15

bottom_y_min <- subgroup_effect_summary_plot %>%
  filter(group_type %in% c("Sex", "Race", "Age")) %>%
  summarise(
    y = min(pct_increase, na.rm = TRUE)
  ) %>%
  pull(y)

bottom_y_max <- subgroup_effect_summary_plot %>%
  filter(group_type %in% c("Sex", "Race", "Age")) %>%
  summarise(
    y = max(pct_increase, na.rm = TRUE)
  ) %>%
  pull(y)

if (!is.finite(bottom_y_min)) bottom_y_min <- 240
if (!is.finite(bottom_y_max)) bottom_y_max <- 300

bottom_y_min <- min(235, floor(bottom_y_min / 10) * 10)
bottom_y_max <- max(300, ceiling(bottom_y_max / 10) * 10)

make_subgroup_panel_horizontal <- function(domain_name,
                                           exposure_summary,
                                           effect_summary,
                                           group_levels,
                                           top_y_max,
                                           bottom_y_min,
                                           bottom_y_max,
                                           show_y_title = FALSE) {
  
  group_order_domain <- group_levels %>%
    filter(group_type == domain_name) %>%
    arrange(group_order) %>%
    pull(group)
  
  bar_df <- exposure_summary %>%
    filter(
      group_type == domain_name,
      scenario %in% c("original", "modified")
    ) %>%
    mutate(
      group = factor(group, levels = group_order_domain),
      scenario = factor(
        scenario,
        levels = c("original", "modified"),
        labels = c("Original plan", "Revised plan")
      )
    ) %>%
    filter(!is.na(group), !is.na(scenario))
  
  eff_df <- effect_summary %>%
    filter(group_type == domain_name) %>%
    mutate(
      group = factor(group, levels = group_order_domain)
    ) %>%
    filter(!is.na(group))
  
  has_pct_ci <- all(c("pct_increase_low", "pct_increase_high") %in% names(eff_df))
  
  y_lab_top <- if (show_y_title) {
    expression("Cumulative weighted coal-" * PM[2.5] * "\n(" * mu * "g/m"^3 * "-years)")
  } else {
    NULL
  }
  
  y_lab_bottom <- if (show_y_title) {
    "Increase from original (%)"
  } else {
    NULL
  }
  
  p_top <- ggplot(
    bar_df,
    aes(
      x = group,
      y = cumulative_pw_coal_pm25,
      fill = scenario
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.72),
      width = 0.62
    ) +
    scale_fill_manual(
      values = c(
        "Original plan" = "#6BAED6",
        "Revised plan" = "#8B1A1A"
      ),
      drop = TRUE,
      na.translate = FALSE,
      name = NULL
    ) +
    scale_y_continuous(
      limits = c(0, top_y_max),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      x = NULL,
      y = y_lab_top
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 7, margin = margin(r = 5)),
      axis.text.y = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(4, 5, 1, 5)
    )
  
  p_bottom <- ggplot(
    eff_df,
    aes(
      x = group,
      y = pct_increase
    )
  )
  
  if (has_pct_ci) {
    p_bottom <- p_bottom +
      geom_errorbar(
        aes(
          ymin = pct_increase_low,
          ymax = pct_increase_high
        ),
        width = 0.12,
        linewidth = 0.4,
        color = "black"
      )
  }
  
  p_bottom <- p_bottom +
    geom_point(
      shape = 23,
      size = 2.8,
      fill = "white",
      color = "black",
      stroke = 0.85
    ) +
    geom_text(
      aes(
        label = paste0("+", number(pct_increase, accuracy = 0.1), "%")
      ),
      vjust = -1.0,
      size = 2.35
    ) +
    scale_y_continuous(
      limits = c(bottom_y_min, bottom_y_max),
      labels = label_percent(scale = 1, accuracy = 1),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = NULL,
      y = y_lab_bottom
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        size = 7.4
      ),
      axis.title.y = element_text(size = 7, margin = margin(r = 5)),
      axis.text.y = element_text(size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(1, 5, 4, 5)
    )
  
  p_domain_title <- ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = domain_name,
      fontface = "bold",
      size = 3.4
    ) +
    theme_void(base_family = "sans") +
    theme(
      plot.margin = margin(0, 5, 0, 5)
    )
  
  p_top / plot_spacer() / p_bottom / p_domain_title +
    plot_layout(
      heights = c(3.0, 0.12, 1.25, 0.22)
    )
}

# Shared legend, built from a dummy plot
legend_plot <- tibble(
  group = factor(c("A", "B")),
  value = c(1, 1),
  scenario = factor(
    c("Original plan", "Revised plan"),
    levels = c("Original plan", "Revised plan")
  )
) %>%
  ggplot(
    aes(
      x = group,
      y = value,
      fill = scenario
    )
  ) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "Original plan" = "#6BAED6",
      "Revised plan" = "#8B1A1A"
    ),
    name = NULL
  ) +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8.8),
    legend.key.size = unit(0.35, "cm")
  )

shared_legend <- cowplot::get_legend(legend_plot)

p_sex_subgroup <- make_subgroup_panel_horizontal(
  domain_name = "Sex",
  exposure_summary = subgroup_exposure_summary_plot,
  effect_summary = subgroup_effect_summary_plot,
  group_levels = group_levels_plot,
  top_y_max = top_y_max,
  bottom_y_min = bottom_y_min,
  bottom_y_max = bottom_y_max,
  show_y_title = TRUE
)

p_race_subgroup <- make_subgroup_panel_horizontal(
  domain_name = "Race",
  exposure_summary = subgroup_exposure_summary_plot,
  effect_summary = subgroup_effect_summary_plot,
  group_levels = group_levels_plot,
  top_y_max = top_y_max,
  bottom_y_min = bottom_y_min,
  bottom_y_max = bottom_y_max,
  show_y_title = FALSE
)

p_age_subgroup <- make_subgroup_panel_horizontal(
  domain_name = "Age",
  exposure_summary = subgroup_exposure_summary_plot,
  effect_summary = subgroup_effect_summary_plot,
  group_levels = group_levels_plot,
  top_y_max = top_y_max,
  bottom_y_min = bottom_y_min,
  bottom_y_max = bottom_y_max,
  show_y_title = FALSE
)

p_subgroup_delayed_effect_1x3 <- (
  p_sex_subgroup | p_race_subgroup | p_age_subgroup
) +
  plot_layout(
    widths = c(0.70, 1.35, 1.35),
    guides = "collect"
  ) &
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal"
  )

p_subgroup_delayed_effect_1x3 <- p_subgroup_delayed_effect_1x3 +
  plot_annotation(
    theme = theme(
      plot.margin = margin(8, 10, 8, 10)
    )
  )

p_subgroup_delayed_effect_1x3

ggsave(
  file.path(dir_figures, "figure_subgroup_delayed_effect_1x3.png"),
  p_subgroup_delayed_effect_1x3,
  width = 11.5,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

# ============================================================
# Supplementary figure: plant-specific cumulative deaths
# by retirement scenario, 2026-2035
# ============================================================
# Uses central beta estimates and allocated plant-level deaths.
# Annual deaths are accumulated within each plant-scenario pair.
# ------------------------------------------------------------
# 2. Set plotting years and plant order
# ------------------------------------------------------------

# ------------------------------------------------------------
# Plant labels for plotting
# ------------------------------------------------------------

plant_label_lookup_all20 <- c(
  "Welsh" = "Welsh (TX)",
  "Coleto Creek" = "Coleto Creek (TX)",
  "Limestone" = "Limestone (TX)",
  "Dave Johnston" = "Dave Johnston (WY)",
  "North Valmy" = "North Valmy (NV)",
  "Stanton Energy Center" = "Stanton Energy Center (FL)",
  "Baldwin Energy Complex" = "Baldwin Energy Complex (IL)",
  "Naughton" = "Naughton (WY)",
  "Columbia" = "Columbia (WI)",
  "Columbia (WI)" = "Columbia (WI)",
  "Transalta Centralia Generation" = "Centralia (WA)",
  "Centralia" = "Centralia (WA)",
  "Winyah" = "Winyah (SC)",
  "Edgewater" = "Edgewater (WI)",
  "Comanche" = "Comanche (CO)",
  "Comanche (CO)" = "Comanche (CO)",
  "TES Filer City Station" = "TES Filer City Station (MI)",
  "J H Campbell" = "J H Campbell (MI)",
  "F B Culley" = "F B Culley (IN)",
  "R M Schahfer" = "R M Schahfer (IN)",
  "South Oak Creek" = "South Oak Creek (WI)",
  "Sherburne County" = "Sherburne County (MN)",
  "Craig" = "Craig (CO)",
  "Craig (CO)" = "Craig (CO)"
)

plant_label_lookup_all20_tbl <- tibble(
  plant_name = names(plant_label_lookup_all20),
  plant_label = unname(plant_label_lookup_all20)
)

years_plot_deaths <- 2026:2035

plant_order_cumulative_deaths <- plant_cumulative_diff |>
  filter(beta_type == "central") |>
  arrange(desc(additional_cumulative_deaths_allocated)) |>
  distinct(plant_id, plant_name, .keep_all = TRUE) |>
  left_join(
    plant_label_lookup_all20_tbl,
    by = "plant_name"
  ) |>
  mutate(
    plant_label = coalesce(plant_label, plant_name)
  ) |>
  pull(plant_label)

plant_frame_deaths <- plant_cumulative_diff |>
  filter(beta_type == "central") |>
  arrange(desc(additional_cumulative_deaths_allocated)) |>
  distinct(plant_id, plant_name) |>
  left_join(
    plant_label_lookup_all20_tbl,
    by = "plant_name"
  ) |>
  mutate(
    plant_label = coalesce(plant_label, plant_name)
  )
# ------------------------------------------------------------
# 3. Prepare annual plant-scenario deaths
# ------------------------------------------------------------

plant_deaths_annual_plot <- plant_year_deaths_by_scenario |>
  filter(
    beta_type == "central",
    scenario %in% c("original", "modified"),
    year %in% years_plot_deaths
  ) |>
  group_by(
    plant_id,
    plant_name,
    scenario,
    year
  ) |>
  summarise(
    annual_deaths_allocated = sum(plant_deaths_allocated, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    plant_label_lookup_all20_tbl,
    by = "plant_name"
  ) |>
  mutate(
    plant_label = coalesce(plant_label, plant_name)
  )

# ------------------------------------------------------------
# 4. Complete missing plant-year-scenario combinations
# ------------------------------------------------------------
# Missing years/scenarios are filled with zero annual deaths before cumsum.

plant_deaths_cumulative_plot <- plant_frame_deaths |>
  crossing(
    scenario = c("original", "modified"),
    year = 2025:2035
  ) |>
  left_join(
    plant_deaths_annual_plot |>
      select(
        plant_id,
        plant_name,
        scenario,
        year,
        annual_deaths_allocated
      ),
    by = c("plant_id", "plant_name", "scenario", "year")
  ) |>
  mutate(
    annual_deaths_allocated = tidyr::replace_na(annual_deaths_allocated, 0)
  ) |>
  arrange(
    plant_id,
    scenario,
    year
  ) |>
  group_by(
    plant_id,
    plant_name,
    plant_label,
    scenario
  ) |>
  mutate(
    cumulative_deaths_allocated = cumsum(annual_deaths_allocated)
  ) |>
  ungroup() |>
  mutate(
    plant_label = factor(
      plant_label,
      levels = plant_order_cumulative_deaths
    ),
    scenario = factor(
      scenario,
      levels = c("original", "modified"),
      labels = c(
        "Original retirement scenario",
        "Revised retirement scenario"
      )
    )
  )

write_csv(
  plant_deaths_cumulative_plot,
  file.path(
    dir_tables,
    "plant_specific_cumulative_allocated_deaths_by_scenario_2026_2035.csv"
  )
)


# ------------------------------------------------------------
# 5. Plot 20 plant-specific cumulative death trajectories
# ------------------------------------------------------------

p_plant_cumulative_deaths_all20 <- ggplot(
  plant_deaths_cumulative_plot,
  aes(
    x = year,
    y = cumulative_deaths_allocated,
    color = scenario,
    group = scenario
  )
) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 0.1) +
  facet_wrap(
    ~ plant_label,
    ncol = 4,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = c(2025, 2030, 2035),
    limits = c(2025, 2035)
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.10, 0.30))
  ) +
  scale_color_manual(
    values = c(
      "Original retirement scenario" = "#6BAED6",
      "Revised retirement scenario" = "#8B1A1A"
    ),
    name = NULL
  ) +
  labs(
    x = NULL,
    y = "Cumulative deaths"
  ) +
  theme_bw(base_family = "sans") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 7,
      face = "bold"
    ),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.2
    ),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.key.width = unit(0.35, "in"),
    panel.spacing = unit(0.55, "lines"),
    plot.margin = margin(6, 6, 6, 6)
  )

p_plant_cumulative_deaths_all20

# ------------------------------------------------------------
# 6. Save figure
# ------------------------------------------------------------

ggsave(
  file.path(
    dir_figures,
    "figure_s_plant_specific_cumulative_allocated_deaths_by_scenario_2026_2035.png"
  ),
  p_plant_cumulative_deaths_all20,
  width = 8.5,
  height = 10,
  dpi = 600,
  bg = "white"
)

