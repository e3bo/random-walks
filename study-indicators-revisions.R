#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

eval_date <- "2021-05-24"
forecast_loc <- "06"
forecast_dates <- c(seq.Date(lubridate::ymd("2020-06-29"), 
                           lubridate::ymd("2021-04-26"),
                           by = 7), 
                    lubridate::ymd(eval_date))

hopdir <- map(forecast_dates, ~file.path("hopkins", .x))

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

p <- allsum %>%
  ggplot(aes(
    x = target_end_date,
    ymin = min,
    ymax = max,
    y = last
  )) +
  geom_linerange() + geom_point() + facet_grid(variable ~ ., scales = "free_y") +
  theme_minimal() +
  labs(x = "Date", y = "Count")

ggsave("indicators-revisions-cal.png", plot = p, width = 6.5, height = 6)

q('no')
winsize <- 14
outlier_limits <- all %>% 
  filter(issue_date == eval_date ) %>% 
  select(-issue_date, -wday) %>%
  pivot_longer(c(cases, deaths, hospitalizations), names_to = "variable") %>%
  group_by(target_end_date, variable) %>% 
  arrange(target_end_date) %>%
  group_by(variable) %>%
  mutate(naive_res = value - lag(value), 
         naive_res7 = value - lag(value, n = 7L)) %>%
  summarize(rollsd_naive_res = slider::slide_dbl(naive_res, sd, .before = winsize - 1, .complete = TRUE),
            rollsd_naive_res7 = slider::slide_dbl(naive_res7, sd, .before = winsize  - 1, .complete = TRUE),
            window_end = target_end_date,
            value = value,
            .groups = "drop") %>%
  mutate(lower = value - 2 * rollsd_naive_res,
         upper = value + 2 * rollsd_naive_res)
  
df <- allsum %>% left_join(outlier_limits, by = c("target_end_date" = "window_end", "variable"))

df %>% ggplot(aes(x = target_end_date, y = last)) + geom_point() + geom_line(aes(y = lower)) + geom_line(aes(y = upper)) + facet_grid(variable ~., scales = "free_y")
