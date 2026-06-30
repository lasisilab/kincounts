library(readr)
library(dplyr)
library(tidyr)

project_root <- if (dir.exists("data")) "." else ".."
raw_file <- file.path(project_root, "data", "raw", "usa_00001.csv")
processed_dir <- file.path(project_root, "data", "processed")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_file)) {
  stop(
    "Missing raw IPUMS extract: ", raw_file, "\n",
    "Download the extract from IPUMS USA, save it there, and rerun this script.",
    call. = FALSE
  )
}

ipums_data <- read_csv(raw_file, show_col_types = FALSE)

required_cols <- c("YEAR", "AGE", "CHBORN", "RACE")
missing_cols <- setdiff(required_cols, names(ipums_data))
if (length(missing_cols) > 0) {
  stop("Raw extract is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

races <- c(
  "White",
  "Black/African American",
  "American Indian or Alaska Native",
  "Chinese",
  "Japanese",
  "Other Asian or Pacific Islander",
  "Other race, nec",
  "Two major races",
  "Three or more major races"
)

mother_data <- ipums_data %>%
  filter(CHBORN > 0) %>%
  mutate(
    n_child = CHBORN - 1L,
    Race = races[RACE]
  ) %>%
  transmute(
    Year = YEAR,
    Age = AGE,
    Race,
    n_child
  )

fertility_data <- bind_rows(
  mother_data %>%
    filter(Year %in% c(1960, 1970, 1980, 1990), Age >= 50, Age <= 59),
  mother_data %>%
    filter(Year == 1960, Age >= 60, Age <= 69) %>%
    transmute(
      Year = 1950,
      Age = Age - 10L,
      Race,
      n_child
    )
) %>%
  filter(!is.na(n_child)) %>%
  arrange(Year, Age)

cohort_lookup <- tibble(
  Year = c(1950, 1960, 1970, 1980, 1990),
  cohort = c("1891-1900", "1901-1910", "1911-1920", "1921-1930", "1931-1940")
)

cohort_sample_sizes <- fertility_data %>%
  count(Year, name = "n_women") %>%
  left_join(cohort_lookup, by = "Year") %>%
  select(Year, cohort, n_women)

fertility_counts <- fertility_data %>%
  count(Year, n_child, name = "n_women") %>%
  left_join(cohort_lookup, by = "Year") %>%
  select(Year, cohort, n_child, n_women) %>%
  arrange(Year, n_child)

fertility_pmf <- fertility_counts %>%
  group_by(Year, cohort) %>%
  mutate(probability = n_women / sum(n_women)) %>%
  ungroup() %>%
  transmute(Year, cohort, children = n_child, probability)

write_csv(cohort_sample_sizes, file.path(processed_dir, "ipums_cohort_sample_sizes.csv"))
write_csv(fertility_counts, file.path(processed_dir, "ipums_fertility_counts_by_cohort.csv"))
write_csv(fertility_pmf, file.path(processed_dir, "ipums_fertility_pmf_by_cohort.csv"))

fertility_pmf %>%
  group_split(Year) %>%
  lapply(function(dat) {
    year <- unique(dat$Year)
    write_csv(
      dat %>% transmute(children, probability),
      file.path(processed_dir, paste0("fertility_pmf_", year, ".csv"))
    )
  })

