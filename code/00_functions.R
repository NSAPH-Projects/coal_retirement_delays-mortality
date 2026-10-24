# ============================================================
# Shared functions
# ============================================================

load_required_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      ". Install them before running the scripts."
    )
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

standardize_zcta <- function(x) {
  stringr::str_pad(as.character(x), width = 5, side = "left", pad = "0")
}

clean_names_simple <- function(df) {
  names(df) <- names(df) |>
    tolower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
  df
}

share_to_unit <- function(x) {
  x <- as.numeric(x)
  dplyr::if_else(!is.na(x) & x > 1, x / 100, x)
}

make_square <- function(x, y, dx, dy) {
  sf::st_polygon(list(matrix(
    c(
      x - dx / 2, y - dy / 2,
      x + dx / 2, y - dy / 2,
      x + dx / 2, y + dy / 2,
      x - dx / 2, y + dy / 2,
      x - dx / 2, y - dy / 2
    ),
    ncol = 2,
    byrow = TRUE
  )))
}

extract_plant_id_from_unit_col <- function(x) {
  unit_id <- stringr::str_remove(x, "^X")
  as.integer(stringr::str_extract(unit_id, "^\\d+"))
}

get_fst_colnames <- function(fp) {
  names(fst::read_fst(fp, from = 1, to = 1))
}

annual_operating_fraction <- function(year, retirement_date, indefinitely_delayed = 0) {
  if (is.na(retirement_date)) {
    return(1)
  }
  if (!is.na(indefinitely_delayed) && indefinitely_delayed == 1) {
    return(1)
  }

  retirement_year <- lubridate::year(retirement_date)
  retirement_month <- lubridate::month(retirement_date)

  if (year < retirement_year) {
    return(1)
  }
  if (year > retirement_year) {
    return(0)
  }

  # If a plant retires in month m, approximate operation as m/12 of that year.
  pmax(pmin(retirement_month / 12, 1), 0)
}

weighted_cor <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[ok]
  y <- y[ok]
  w <- w[ok]
  if (length(x) < 3) return(NA_real_)

  mx <- sum(w * x) / sum(w)
  my <- sum(w * y) / sum(w)
  cov_xy <- sum(w * (x - mx) * (y - my)) / sum(w)
  var_x <- sum(w * (x - mx)^2) / sum(w)
  var_y <- sum(w * (y - my)^2) / sum(w)
  cov_xy / sqrt(var_x * var_y)
}

read_grid_total_file <- function(fp, yy) {
  fst::read_fst(fp) |>
    dplyr::transmute(
      year = as.integer(yy),
      x = as.numeric(x),
      y = as.numeric(y),
      coal_pm25 = as.numeric(vals.out)
    ) |>
    dplyr::filter(is.finite(x), is.finite(y), is.finite(coal_pm25))
}

build_total_population_subgroups <- function(demo) {
  # Sex and age use total population. Dual eligible people are defined as
  # n_bene * dual_pct; non-dual is all other residents. Race uses five main
  # mutually exclusive race categories and normalizes them within ZCTA.
  demo2 <- demo |>
    dplyr::mutate(
      male_share_raw = share_to_unit(male_pct),
      female_share_raw = share_to_unit(female_pct),
      sex_share_sum = male_share_raw + female_share_raw,
      male_share = dplyr::if_else(sex_share_sum > 0, male_share_raw / sex_share_sum, NA_real_),
      female_share = dplyr::if_else(sex_share_sum > 0, female_share_raw / sex_share_sum, NA_real_),
      dual_share = share_to_unit(dual_pct),
      pop_20_34_years = pmax(
        population - pop_under_20_years - pop_35_49_years -
          pop_50_64_years - pop_over_65_years,
        0
      )
    )

  sex_pop <- dplyr::bind_rows(
    demo2 |>
      dplyr::transmute(
        zcta, group_type = "Sex", group = "Male",
        subgroup_population = population * male_share
      ),
    demo2 |>
      dplyr::transmute(
        zcta, group_type = "Sex", group = "Female",
        subgroup_population = population * female_share
      )
  )

  race_pop <- demo2 |>
    dplyr::transmute(
      zcta, population,
      White = pop_white,
      Black = pop_black,
      Asian = pop_asian,
      Native = pop_native,
      Hispanic = pop_hispanic
    ) |>
    tidyr::pivot_longer(
      cols = c(White, Black, Asian, Native, Hispanic),
      names_to = "group",
      values_to = "raw_group_count"
    ) |>
    dplyr::group_by(zcta) |>
    dplyr::mutate(
      race_total = sum(raw_group_count, na.rm = TRUE),
      subgroup_share = dplyr::if_else(race_total > 0, raw_group_count / race_total, NA_real_),
      subgroup_population = population * subgroup_share
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(zcta, group_type = "Race", group, subgroup_population)

  age_pop <- demo2 |>
    dplyr::transmute(
      zcta,
      `Under 20` = pop_under_20_years,
      `20-34` = pop_20_34_years,
      `35-49` = pop_35_49_years,
      `50-64` = pop_50_64_years,
      `65+` = pop_over_65_years
    ) |>
    tidyr::pivot_longer(
      cols = c(`Under 20`, `20-34`, `35-49`, `50-64`, `65+`),
      names_to = "group",
      values_to = "subgroup_population"
    ) |>
    dplyr::mutate(group_type = "Age")

  dual_pop <- demo2 |>
    dplyr::filter(!is.na(n_bene), n_bene >= 0, !is.na(dual_share)) |>
    dplyr::mutate(
      dual_count = n_bene * dual_share,
      nondual_count = pmax(population - dual_count, 0)
    ) |>
    dplyr::transmute(zcta, `Non-dual` = nondual_count, Dual = dual_count) |>
    tidyr::pivot_longer(
      cols = c(`Non-dual`, Dual),
      names_to = "group",
      values_to = "subgroup_population"
    ) |>
    dplyr::mutate(group_type = "Dual eligibility")

  dplyr::bind_rows(sex_pop, race_pop, dual_pop, age_pop) |>
    dplyr::filter(!is.na(subgroup_population), subgroup_population > 0)
}

build_medicare_subgroup_shares <- function(demo) {
  # Shares used for future Medicare subgroup exposure summaries.
  demo2 <- demo |>
    dplyr::mutate(
      male_share_raw = share_to_unit(male_pct),
      female_share_raw = share_to_unit(female_pct),
      sex_share_sum = male_share_raw + female_share_raw,
      male_share = dplyr::if_else(sex_share_sum > 0, male_share_raw / sex_share_sum, NA_real_),
      female_share = dplyr::if_else(sex_share_sum > 0, female_share_raw / sex_share_sum, NA_real_),
      dual_share = share_to_unit(dual_pct),
      pop_20_34_years = pmax(
        population - pop_under_20_years - pop_35_49_years -
          pop_50_64_years - pop_over_65_years,
        0
      )
    )

  sex_share <- dplyr::bind_rows(
    demo2 |> dplyr::transmute(zcta, group_type = "Sex", group = "Male", subgroup_share = male_share),
    demo2 |> dplyr::transmute(zcta, group_type = "Sex", group = "Female", subgroup_share = female_share)
  )

  dual_share <- dplyr::bind_rows(
    demo2 |> dplyr::transmute(zcta, group_type = "Dual eligibility", group = "Dual", subgroup_share = dual_share),
    demo2 |> dplyr::transmute(zcta, group_type = "Dual eligibility", group = "Non-dual", subgroup_share = pmax(1 - dual_share, 0))
  )

  race_share <- demo2 |>
    dplyr::transmute(
      zcta,
      White = pop_white,
      Black = pop_black,
      Asian = pop_asian,
      Native = pop_native,
      Hispanic = pop_hispanic
    ) |>
    tidyr::pivot_longer(
      cols = c(White, Black, Asian, Native, Hispanic),
      names_to = "group",
      values_to = "raw_group_count"
    ) |>
    dplyr::group_by(zcta) |>
    dplyr::mutate(
      race_total = sum(raw_group_count, na.rm = TRUE),
      subgroup_share = dplyr::if_else(race_total > 0, raw_group_count / race_total, NA_real_)
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(zcta, group_type = "Race", group, subgroup_share)

  age_share <- demo2 |>
    dplyr::transmute(
      zcta,
      `Under 20` = pop_under_20_years,
      `20-34` = pop_20_34_years,
      `35-49` = pop_35_49_years,
      `50-64` = pop_50_64_years,
      `65+` = pop_over_65_years
    ) |>
    tidyr::pivot_longer(
      cols = c(`Under 20`, `20-34`, `35-49`, `50-64`, `65+`),
      names_to = "group",
      values_to = "age_count"
    ) |>
    dplyr::group_by(zcta) |>
    dplyr::mutate(
      age_total = sum(age_count, na.rm = TRUE),
      subgroup_share = dplyr::if_else(age_total > 0, age_count / age_total, NA_real_)
    ) |>
    dplyr::ungroup() |>
    dplyr::transmute(zcta, group_type = "Age", group, subgroup_share)

  dplyr::bind_rows(sex_share, dual_share, race_share, age_share) |>
    dplyr::filter(!is.na(subgroup_share), subgroup_share >= 0)
}
