#!/usr/bin/env Rscript

library(tidyverse)
source("covidhub-common.R")



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
