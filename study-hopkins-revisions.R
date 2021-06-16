#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

forecast_loc <- "06"
forecast_dates <- seq.Date(lubridate::ymd("2020-06-29"), 
                           lubridate::ymd("2021-04-26"),
                           by = 7)

hopdir <- map(forecast_dates, ~file.path("hopkins", .x))

tdat <- map(hopdir, load_hopkins, weekly = FALSE)

load_and_filter <- function(dir){
  load_hopkins(dir, weekly = FALSE) %>% 
    filter(location == forecast_loc) %>%
    filter(target_type == "day ahead inc case" |
             target_type == "day ahead inc death")
}

hoploc <- map(hopdir, load_and_filter)

names(hoploc) <- forecast_dates

hopall <- bind_rows(hoploc, .id = "issue_date") %>%
  pivot_wider(names_from = target_type, values_from = value) %>%
  ungroup() %>%
  mutate(wday = lubridate::wday(target_end_date)) %>%
  rename(cases = `day ahead inc case`, deaths = `day ahead inc death`) %>%
  select(target_end_date, wday, cases, deaths, issue_date)

hfd <- forecast_dates[forecast_dates > "2020-11-15"]

hhs <- map(hfd, ~load_health_data(forecast_date = .x, forecast_loc = "06"))
names(hhs) <- hfd

all <- bind_rows(hhs, .id = "issue_date") %>% 
  right_join(hopall, by = c("target_end_date", "issue_date"))

allsum <- all %>%
  pivot_longer(c(cases:deaths, hospitalizations), names_to = "variable") %>%
  group_by(target_end_date, variable) %>%
  arrange(issue_date) %>%
  summarise(min = min(value),
            max = max(value),
            last = value[n()])

allsum %>% 
  ggplot(aes(x = target_end_date, ymin = min, ymax = max, y = last)) + 
  geom_linerange() + geom_point() + facet_grid(variable~., scales = "free_y") +
  theme_minimal() + 
  labs(x = "Date", y = "Count")

hopsum %>% 
  ggplot(aes(x = target_end_date, ymin = mindeaths, ymax = maxdeaths, y = lastdeaths)) + 
  geom_linerange() + geom_point()