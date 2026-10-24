# ============================================================
# 06_plant_retirement_status_map.R
# Plant retirement-status map with 2024 SO2 emissions
# ============================================================

source("R/00_config.R")
source("R/00_functions.R")

load_required_packages(c(
  "ggplot2", "dplyr", "readr", "maps", "scales", "grid", "stringr"
))

# ------------------------------------------------------------
# 1. Input data
# ------------------------------------------------------------
plants <- read_csv(
  plant_map_fp,
  show_col_types = FALSE
) %>%
  mutate(
    category = str_squish(as.character(category)),
    category = recode(
      category,
      "Retired" = "Retired",
      "Plan unchanged" = "No change",
      "Unchanged" = "No change",
      "Explicitly delayed" = "Delayed retirement",
      "Date removed / indefinitely delayed" = "Retirement date removed"
    ),
    category = factor(
      category,
      levels = c(
        "Retired",
        "No change",
        "Delayed retirement",
        "Retirement date removed"
      )
    ),
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude),
    so2_plot = pmax(as.numeric(so2_2024_tons), 0)
  )

plants <- read_csv(
  file.path(dir_raw, "coal_plants_map_plot_df_so2.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    category_raw = stringr::str_squish(as.character(category)),
    
    category = dplyr::case_when(
      category_raw == "Retired" ~ "Retired",
      category_raw %in% c("Plan unchanged", "Unchanged") ~ "No change",
      category_raw == "Explicitly delayed" ~ "Delayed retirement",
      category_raw == "Date removed / indefinitely delayed" ~ "Retirement date removed",
      TRUE ~ NA_character_
    ),
    
    category = factor(
      category,
      levels = c(
        "Retired",
        "No change",
        "Delayed retirement",
        "Retirement date removed"
      )
    ),
    
    longitude = as.numeric(longitude),
    latitude = as.numeric(latitude),
    so2_plot = pmax(as.numeric(so2_2024_tons), 0)
  ) %>%
  filter(
    is.finite(longitude),
    is.finite(latitude),
    is.finite(so2_plot),
    !is.na(category)
  ) %>%
  as.data.frame()

# ------------------------------------------------------------
# 2. QC
# ------------------------------------------------------------

plants %>%
  count(category)

plants %>%
  summarise(
    n_plants = n(),
    min_so2 = min(so2_plot, na.rm = TRUE),
    median_so2 = median(so2_plot, na.rm = TRUE),
    max_so2 = max(so2_plot, na.rm = TRUE)
  ) %>%
  print()

# ------------------------------------------------------------
# 3. Map data and aesthetics
# ------------------------------------------------------------

us_states <- ggplot2::map_data("state")

pal <- c(
  "Retired" = "#6BAED6",
  "No change" = "#9E9E9E",
  "Delayed retirement" = "#F28E2B",
  "Retirement date removed" = "#9B1D20"
)

# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------

p_plant_status_map <- ggplot() +
  geom_polygon(
    data = us_states,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "#F5F5F5",
    color = "#D0D0D0",
    linewidth = 0.16
  ) +
  geom_point(
    data = plants,
    aes(
      x = longitude,
      y = latitude,
      size = so2_plot,
      colour = category
    ),
    shape = 16,
    alpha = 0.88
  ) +
  coord_quickmap(
    xlim = c(-125, -66),
    ylim = c(24, 50),
    expand = FALSE
  ) +
  scale_colour_manual(
    values = pal,
    breaks = c(
      "Retired",
      "No change",
      "Delayed retirement",
      "Retirement date removed"
    ),
    name = NULL,
    drop = FALSE
  ) +
  scale_size_area(
    max_size = 12,
    breaks = c(1000, 5000, 10000),
    labels = comma_format(accuracy = 1),
    name = expression(SO[2] * " emissions (tons)")
  ) +
  guides(
    size = guide_legend(
      order = 1,
      title.position = "left",
      title.hjust = 0.5,
      override.aes = list(
        colour = "#252525",
        alpha = 1
      )
    ),
    colour = guide_legend(
      order = 2,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 5,
        alpha = 1
      )
    )
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    plot.caption = element_text(
      size = 8.5,
      color = "#666666",
      hjust = 0,
      margin = margin(t = 8)
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(8, "pt"),
    legend.spacing.y = unit(2, "pt"),
    legend.title = element_text(size = 7, hjust = 0.5),
    legend.text = element_text(size = 7),
    
    plot.margin = margin(8, 8, 8, 8)
  )

p_plant_status_map

# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------

ggsave(
  file.path(
    dir_figures,
    "map_coal_plants_retirement_status_so2.png"
  ),
  p_plant_status_map,
  width = 7.25,
  height = 4.8,
  dpi = 600,
  bg = "white"
)
