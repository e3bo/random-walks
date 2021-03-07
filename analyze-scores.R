#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")
lambda <- 1 / seq(0.001, 0.1, length.out = 10)

dirnames <- paste0("lambda", sprintf("%06.2f", lambda), "-CEID-InfectionKalman")
dirnames2 <- paste0("lambda", sprintf("%06.2f", lambda), 
                    "-CEID-InfectionKalmanEmp")
load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat1 <- map_dfr(dirnames, load_from_dir)
fdst11 <- map_dfr(dirnames2, load_from_dir)
fdat2 <- bind_rows(fdat0, fdat1, fdst11)

truth_data <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         locations = unique(fdat2$location))

locations_to_exclude <- c("78", "72", "69", "66", "60")

scores <- score_forecasts(fdat2, truth_data, return_format = "wide") %>%
  filter(!location %in% locations_to_exclude) %>%
  select(model, horizon, location, target_variable, target_end_date, 
         coverage_50, coverage_95, abs_error, wis)

s1 <-scores %>% group_by(horizon, location, model, target_end_date) %>% 
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop")

s2 <-scores %>% 
  group_by(horizon, location, model) %>% 
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop_last") %>%
  mutate(relative_wis = wis / wis[str_detect(model, "CEID-Walk")])
write_csv(s2, "horizon-location-model.csv")

s3 <-scores %>% 
  group_by(location, model) %>% 
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop_last") %>%
  mutate(relative_wis = wis / wis[str_detect(model, "CEID-Walk")])
write_csv(s3, "location-model.csv")

ascv <- s3 %>%
  filter(str_detect(model, "^lambda.*Emp$")) %>%
  group_by(location) %>%
  filter(abs_error == min(abs_error)) %>%
  mutate(model = "CV-CEID-InfectionKalmanEmp")

s4 <-
  bind_rows(ascv, s3) %>% group_by(model) %>%
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop_last")  

write_csv(s4, "model.csv")
