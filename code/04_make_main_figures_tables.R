# ============================================================
# 04_make_main_figures_tables.R
# Make final figures and formatted summary tables
# ============================================================

source("R/00_config.R")
source("R/00_functions.R")

load_required_packages(c(
  "dplyr", "tidyr", "ggplot2", "scales", "readr", "stringr", "forcats", "tibble"
))

# ------------------------------------------------------------
# 1. Read mortality output tables from script 03
# ------------------------------------------------------------
scenario_year_deaths_total <- read_csv(
  file.path(dir_tables, "scenario_year_total_deaths.csv"),
  show_col_types = FALSE
) |>
  mutate(
    beta_type = as.character(beta_type),
    scenario = str_to_lower(as.character(scenario)),
    year = as.integer(year),
    total_deaths = as.numeric(total_deaths)
  ) |>
  filter(
    beta_type %in% c("lower", "central", "upper"),
    scenario %in% c("original", "modified")
  )

plant_cumulative_by_scenario <- read_csv(
  file.path(dir_tables, "plant_cumulative_deaths_by_scenario_2026_2035.csv"),
  show_col_types = FALSE
)

plant_cumulative_diff <- read_csv(
  file.path(dir_tables, "plant_cumulative_deaths_modified_minus_original_2026_2035.csv"),
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Cumulative scenario deaths with beta lower/upper range
# ------------------------------------------------------------
cumulative_deaths_long <- scenario_year_deaths_total |>
  arrange(beta_type, scenario, year) |>
  group_by(beta_type, scenario) |>
  mutate(cumulative_deaths = cumsum(total_deaths)) |>
  ungroup()


baseline_2025 <- cumulative_deaths_long |>
  distinct(beta_type, scenario) |>
  mutate(year = 2025L, total_deaths = 0, cumulative_deaths = 0)

cumulative_deaths_long <- bind_rows(baseline_2025, cumulative_deaths_long) |>
  arrange(beta_type, scenario, year)

cumulative_deaths_plot <- cumulative_deaths_long |>
  select(beta_type, scenario, year, cumulative_deaths) |>
  pivot_wider(names_from = beta_type, values_from = cumulative_deaths) |>
  mutate(
    scenario_label = case_when(
      scenario == "original" ~ "Original plan",
      scenario == "modified" ~ "Revised plan",
      TRUE ~ scenario
    ),
    scenario_label = factor(scenario_label, levels = c("Original plan", "Revised plan"))
  )

write_csv(cumulative_deaths_plot, file.path(dir_tables, "cumulative_total_deaths_by_scenario_with_ci.csv"))

p_cumulative_deaths <- ggplot(
  cumulative_deaths_plot,
  aes(x = year, y = central, color = scenario_label, fill = scenario_label)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, linewidth = 0) +
  geom_line(linewidth = 1.25) +
  scale_color_manual(values = c("Original plan" = color_original, "Revised plan" = color_modified), name = NULL) +
  scale_fill_manual(values = c("Original plan" = color_original, "Revised plan" = color_modified), name = NULL) +
  scale_x_continuous(breaks = 2025:2035) +
  scale_y_continuous(labels = label_number(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Cumulative coal-PM2.5-attributable deaths, 2026-2035",
    subtitle = "Point estimates with lower and upper bounds from the mortality-response estimate",
    x = NULL,
    y = "Cumulative attributable deaths",
    caption = "Cumulative deaths are calculated from scenario-year total attributable deaths. Shaded bands use lower and upper coal-PM2.5 mortality response estimates. Year 2025 is shown as a zero baseline; modeled deaths accrue over 2026-2035."
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "grey35", hjust = 0),
    plot.caption = element_text(size = 8.5, color = "grey40", hjust = 0),
    legend.position = "top",
    legend.justification = "center",
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10.5),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(dir_figures, "figure_cumulative_total_attributable_deaths_by_scenario.png"),
  p_cumulative_deaths,
  width = 10.8,
  height = 6.4,
  dpi = 320
)

# ------------------------------------------------------------
# 3. Cumulative additional deaths only: revised minus original
# ------------------------------------------------------------
annual_additional_deaths <- scenario_year_deaths_total |>
  select(beta_type, scenario, year, total_deaths) |>
  pivot_wider(names_from = scenario, values_from = total_deaths, values_fill = 0) |>
  mutate(annual_additional_deaths = modified - original) |>
  arrange(beta_type, year) |>
  group_by(beta_type) |>
  mutate(cumulative_additional_deaths = cumsum(annual_additional_deaths)) |>
  ungroup()

baseline_additional_2025 <- annual_additional_deaths |>
  distinct(beta_type) |>
  mutate(
    year = 2025L,
    original = 0,
    modified = 0,
    annual_additional_deaths = 0,
    cumulative_additional_deaths = 0
  )

cumulative_additional_plot <- bind_rows(baseline_additional_2025, annual_additional_deaths) |>
  arrange(beta_type, year) |>
  select(beta_type, year, cumulative_additional_deaths) |>
  pivot_wider(names_from = beta_type, values_from = cumulative_additional_deaths) |>
  transmute(
    year,
    additional_deaths_central = central,
    additional_deaths_low = lower,
    additional_deaths_high = upper
  )

write_csv(
  cumulative_additional_plot,
  file.path(dir_tables, "cumulative_additional_deaths_trajectory_with_ci.csv")
)

p_cumulative_additional <- ggplot(cumulative_additional_plot, aes(x = year, y = additional_deaths_central)) +
  geom_ribbon(
    aes(ymin = additional_deaths_low, ymax = additional_deaths_high),
    fill = color_modified,
    alpha = 0.16,
    linewidth = 0
  ) +
  geom_line(color = color_modified, linewidth = 1.25) +
  scale_x_continuous(breaks = 2025:2035) +
  scale_y_continuous(labels = label_number(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Cumulative additional deaths from delayed retirements",
    subtitle = "Revised plan minus original plan, 2026-2035",
    x = NULL,
    y = "Cumulative additional deaths",
    caption = "Line shows the central estimate; shaded band uses lower and upper coal-PM2.5 mortality response estimates. Year 2025 is shown as a zero baseline."
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "grey35", hjust = 0),
    plot.caption = element_text(size = 8.5, color = "grey40", hjust = 0),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10.5),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(dir_figures, "figure_cumulative_additional_deaths_modified_minus_original.png"),
  p_cumulative_additional,
  width = 10.8,
  height = 6.4,
  dpi = 320
)


# ==================================================
# 4) Plot: cumulative deaths by scenario + additional deaths
# ==================================================
# Lines/ribbons:
#   cumulative attributable deaths under each scenario
#
# Bars/error bars:
#   cumulative additional deaths = modified - original

# --------------------------------------------------
# 4.1 Cumulative additional deaths for bars
# --------------------------------------------------

cumulative_additional_deaths_plot <- cumulative_deaths_plot %>%
  select(
    scenario,
    year,
    central,
    lower,
    upper
  ) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(central, lower, upper),
    values_fill = 0
  ) %>%
  mutate(
    additional_central =
      central_modified - central_original,
    additional_low =
      lower_modified - lower_original,
    additional_high =
      upper_modified - upper_original
  ) %>%
  arrange(year)

write_csv(
  cumulative_additional_deaths_plot,
  file.path(dir_tables, "cumulative_additional_deaths_modified_minus_original_with_ci.csv")
)

# --------------------------------------------------
# 4.2 Plot
# --------------------------------------------------

col_original <- "#6BAED6"
col_revised <- "#8B1A1A"

fill_original <- "#6BAED6"
fill_revised <- "#8B1A1A"

bar_additional_fill <- "#5F8F5F"
bar_additional_line <- "#2F5E2F"

p_cumulative_deaths_combined <- ggplot() +
  geom_col(
    data = cumulative_additional_deaths_plot,
    aes(
      x = year,
      y = additional_central,
      fill = "Additional deaths"
    ),
    alpha = 0.5,
    width = 0.72
  ) +
  geom_errorbar(
    data = cumulative_additional_deaths_plot,
    aes(
      x = year,
      ymin = additional_low,
      ymax = additional_high
    ),
    width = 0.14,
    linewidth = 0.35,
    color = bar_additional_line,
    alpha = 0.75
  ) +
  geom_ribbon(
    data = cumulative_deaths_plot %>% filter(scenario_label == "Original plan"),
    aes(
      x = year,
      ymin = lower,
      ymax = upper
    ),
    fill = "#6BAED6",
    alpha = 0.16,
    linewidth = 0,
    show.legend = FALSE
  ) +
  geom_ribbon(
    data = cumulative_deaths_plot %>% filter(scenario_label == "Revised plan"),
    aes(
      x = year,
      ymin = lower,
      ymax = upper
    ),
    fill = "#8B1A1A",
    alpha = 0.16,
    linewidth = 0,
    show.legend = FALSE
  ) +
  geom_line(
    data = cumulative_deaths_plot,
    aes(
      x = year,
      y = central,
      color = scenario_label
    ),
    linewidth = 1.25
  ) +
  scale_color_manual(
    values = c(
      "Original plan" = col_original,
      "Revised plan" = col_revised
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Additional deaths" = bar_additional_fill,
      "Original plan" = fill_original,
      "Revised plan" = fill_revised
    ),
    name = NULL
  ) +
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(linewidth = 1.25)
    ),
    fill = guide_legend(
      order = 2,
      override.aes = list(
        alpha = 0.55,
        color = NA
      )
    )
  ) +
  scale_x_continuous(
    breaks = 2025:2035,
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 1, big.mark = " "),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = expression(bold("Cumulative coal-" * PM[2.5] * " attributable deaths, 2026-2035")),
    x = NULL,
    y = "Cumulative deaths"
  ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 6),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

p_cumulative_deaths_combined

# --------------------------------------------------
# 4.3 Save
# --------------------------------------------------

ggsave(
  file.path(dir_figures, "figure_cumulative_deaths_with_additional_bars.png"),
  p_cumulative_deaths_combined,
  width = 7.25,
  height = 4.8,
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 5. Top five plants by additional allocated deaths
# ------------------------------------------------------------
top5_plants <- plant_cumulative_diff |>
  filter(beta_type == "central") |>
  arrange(desc(additional_cumulative_deaths_allocated)) |>
  slice_head(n = 20) |>
  select(plant_id, plant_name, additional_central = additional_cumulative_deaths_allocated)

top5_plants <- plant_cumulative_by_scenario |>
  mutate(
    beta_type = as.character(beta_type),
    scenario = str_to_lower(as.character(scenario))
  ) |>
  filter(
    beta_type == "central",
    scenario == "modified"
  ) |>
  arrange(desc(cumulative_deaths_allocated)) |>
  slice_head(n = 20) |>
  transmute(
    plant_id,
    plant_name,
    revised_central = cumulative_deaths_allocated
  )

plant_scenario_ci_top5 <- plant_cumulative_by_scenario |>
  semi_join(top5_plants, by = c("plant_id", "plant_name")) |>
  select(beta_type, scenario, plant_id, plant_name, cumulative_deaths_allocated) |>
  pivot_wider(names_from = beta_type, values_from = cumulative_deaths_allocated) |>
  mutate(
    deaths_central = central,
    deaths_low = lower,
    deaths_high = upper,
    scenario_label = case_when(
      scenario == "original" ~ "Original plan",
      scenario == "modified" ~ "Revised plan",
      TRUE ~ scenario
    )
  ) |>
  left_join(top5_plants, by = c("plant_id", "plant_name"))

plant_additional_ci_top5 <- plant_cumulative_diff |>
  semi_join(top5_plants, by = c("plant_id", "plant_name")) |>
  select(plant_id, plant_name, beta_type, additional_cumulative_deaths_allocated) |>
  pivot_wider(names_from = beta_type, values_from = additional_cumulative_deaths_allocated) |>
  transmute(
    plant_id,
    plant_name,
    additional_deaths_central = central,
    additional_deaths_low = lower,
    additional_deaths_high = upper
  )
  

write_csv(plant_additional_ci_top5, file.path(dir_tables, "top5_plant_additional_deaths_ci_2026_2035.csv"))

plant_order <- top5_plants |>
  arrange(revised_central) |>
  pull(plant_name)

plant_label_lookup <- c(
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

plant_order_labeled <- recode(plant_order, !!!plant_label_lookup)

plant_scenario_ci_top5 <- plant_scenario_ci_top5 |>
  mutate(
    plant_label = recode(plant_name, !!!plant_label_lookup),
    plant_label = factor(plant_label, levels = plant_order_labeled),
    scenario_label = factor(scenario_label, levels = c("Original plan", "Revised plan"))
  )

plant_additional_ci_top5 <- plant_additional_ci_top5 |>
  mutate(
    plant_label = recode(plant_name, !!!plant_label_lookup),
    plant_label = factor(plant_label, levels = plant_order_labeled),
    label_x = max(plant_scenario_ci_top5$deaths_high, na.rm = TRUE) * 1.08,
    label = paste0(
      "+", number(additional_deaths_central, accuracy = 0.1),
      " deaths\n(", number(additional_deaths_low, accuracy = 0.1),
      ", ", number(additional_deaths_high, accuracy = 0.1), ")"
    )
  )

p_top5_plant_deaths <- ggplot(
  plant_scenario_ci_top5,
  aes(x = deaths_central, y = plant_label, fill = scenario_label)
) +
  geom_col(position = position_dodge(width = 0.70), width = 0.58) +
  geom_errorbar(
    aes(xmin = deaths_low, xmax = deaths_high),
    position = position_dodge(width = 0.70),
    width = 0.18,
    linewidth = 0.45,
    orientation = "y",
    color = "black"
  ) +
  geom_text(
    data = plant_additional_ci_top5,
    aes(x = label_x, y = plant_label, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.0,
    lineheight = 0.9
  ) +
  scale_fill_manual(values = c("Original plan" = color_original, "Revised plan" = color_modified), name = NULL) +
  scale_x_continuous(labels = label_number(accuracy = 1), expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Plant-level deaths attributable to delayed retirements",
    x = "Attributable deaths, 2026-2035",
    y = NULL,
    ) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 45, 20, 10)
  ) +
  coord_cartesian(clip = "off")

p_top5_plant_deaths

ggsave(
  file.path(dir_figures, "figure_top5_plant_attributable_deaths_with_ci.png"),
  p_top5_plant_deaths,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)

 # ------------------------------------------------------------
# 5. Summary tables
# ------------------------------------------------------------
total_additional_deaths_ci <- read_csv(
  file.path(dir_tables, "total_cumulative_additional_deaths_ci_2026_2035.csv"),
  show_col_types = FALSE
) |>
  mutate(
    result_label = paste0(
      number(additional_deaths_central, accuracy = 0.1),
      " additional deaths (",
      number(additional_deaths_low, accuracy = 0.1),
      ", ",
      number(additional_deaths_high, accuracy = 0.1),
      ")"
    )
  )

write_csv(total_additional_deaths_ci, file.path(dir_tables, "table_total_additional_deaths_ci_2026_2035.csv"))

plant_contribution_share <- plant_cumulative_diff |>
  filter(beta_type == "central") |>
  arrange(desc(additional_cumulative_deaths_allocated)) |>
  mutate(
    total_additional_deaths = sum(additional_cumulative_deaths_allocated, na.rm = TRUE),
    share_of_total = additional_cumulative_deaths_allocated / total_additional_deaths,
    cumulative_share_of_total = cumsum(share_of_total)
  )

write_csv(plant_contribution_share, file.path(dir_tables, "table_all_plants_additional_deaths_contribution.csv"))
