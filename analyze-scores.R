#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")

lambda <- 20
par2name <- function(lambda){
  paste0("lambda", sprintf("%06.2f", lambda), 
         "-status-quo",
         "-CEID-InfectionKalman")
}
dirnames <- purrr::map(lambda, par2name)

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat1 <- map_dfr(dirnames, load_from_dir)
fdat2 <- bind_rows(fdat0, fdat1)
locations_to_exclude <- c("78", "72", "69", "66", "60")

truth_data1 <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         locations = unique(fdat1$location))

truth_data2 <- load_truth(truth_source = "HealthData",
                          target_variable = "inc hosp",
                          locations = unique(fdat1$location))

truth_data3 <- load_truth(truth_source = "JHU",
                          target_variable = "inc death",
                          locations = unique(fdat1$location))

truth_data <- bind_rows(truth_data1, truth_data2, truth_data3)

scores <- 
  score_forecasts(fdat2, truth_data, return_format = "wide") %>%
  filter(!location %in% locations_to_exclude) %>%
  select(model, forecast_date, horizon, location, target_variable, target_end_date, 
         coverage_50, coverage_95, abs_error, wis)

scores %>% filter(target_end_date == "2020-10-10" & horizon == "1") %>% arrange(wis)

scores %>% filter(target_end_date == "2020-10-17" & horizon == "2") %>% arrange(wis)

scores %>% filter(target_end_date == "2020-10-24" & horizon == "3") %>% arrange(wis)

scores %>% filter(target_end_date == "2020-10-31" & horizon == "4") %>% arrange(wis)

q('no')

tscv <-train_scores %>% 
  group_by(location, model, target_variable) %>% 
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop_last") %>%
  group_by(location, target_variable) %>%
  filter(abs_error == min(abs_error))

## val data will be used to evaluate the model with lambda selected from the
## training data

val_data <- bind_rows(fdat0, fdat1, fdat11) %>%
  filter(forecast_date >= "2020-11-30")

cv_mod <- tscv %>% select(location, model, target_variable) %>% 
  left_join(val_data, by = c("location", "model", "target_variable")) %>%
  mutate(model = "CV-CEID-InfectionKalman")

val_data <- bind_rows(val_data, cv_mod)

scores <- score_forecasts(val_data, truth_data, return_format = "wide") %>%
  filter(!location %in% locations_to_exclude) %>%
  select(model, horizon, location, target_variable, target_end_date, 
         coverage_50, coverage_95, abs_error, wis)

s1 <-scores %>% group_by(horizon, location, model, 
                         target_variable, target_end_date) %>% 
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop")

s2 <-scores %>% 
  group_by(horizon, location, model, target_variable) %>% 
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop_last")

s2 %>% filter(target_variable == "inc case") %>% 
  write_csv("horizon-location-model-cases.csv")
s2 %>% filter(target_variable == "inc hosp") %>% 
  write_csv("horizon-location-model-hosp.csv")
s2 %>% filter(target_variable == "inc death") %>% 
  write_csv("horizon-location-model-death.csv")

s3 <-scores %>% 
  group_by(location, model, target_variable) %>% 
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop_last")
write_csv(s3, "location-model.csv")

s4 <-
  scores %>% group_by(model, target_variable) %>%
  summarize(coverage_50 = mean(coverage_50), 
            coverage_95 = mean(coverage_95),
            abs_error = mean(abs_error),
            wis = mean(wis),
            .groups = "drop_last")  
write_csv(s4, "model.csv")
