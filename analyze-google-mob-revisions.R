#!/usr/bin/env Rscript

library(tidyverse)

read_us_mob_report <- function(path, locpat = "US-CA"){
  readr::read_csv(
    path,
    col_types = readr::cols_only(
      country_region = col_character(),
      sub_region_1 = col_character(),
      sub_region_2 = col_character(),
      iso_3166_2_code = col_character(),
      date = col_date(format = "%Y-%m-%d"),
      residential_percent_change_from_baseline = col_double()
    )
  ) %>% filter(is.na(sub_region_2)) %>%
    select(-sub_region_2) %>%
    filter(str_detect(iso_3166_2_code, locpat))
}


dnames <- dir("google-mobility-reports-wayback") 
paths <- file.path("google-mobility-reports-wayback", 
                   dnames, "US-states-mobility.csv")

dfs <- purrr::map(paths, read_us_mob_report)
names(dfs) <- dnames
df <- bind_rows(dfs, .id = "asof_date")

dfsum <- df %>% group_by(date) %>% arrange(asof_date) %>%
  summarise(
    n = n(),
    varrcb = var(residential_percent_change_from_baseline),
    minrcb = min(residential_percent_change_from_baseline),
    maxrcb = max(residential_percent_change_from_baseline),
    lastrcb = residential_percent_change_from_baseline[n()]
  )

stopifnot(all(dfsum$n == 1 | dfsum$varrcb == 0))

p <-
  ggplot(dfsum, aes(
    x = date,
    y = lastrcb,
  ))  + geom_point() + theme_minimal() + 
  labs(x = "Date", y = "Percent increase in time spent in residential areas")

ggsave(p,filename = "residential-mobility-plot.pdf", width = 4, height = 4)
