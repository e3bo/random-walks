#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")

lambda <- 1 / seq(0.001, 0.1, length.out = 10)
agrid <- c(0.94, 0.95)
par2name <- function(lambda, a){
  paste0("lambda", sprintf("%06.2f", lambda), 
         "-a", sprintf("%02.2f", a),
         "-CEID-InfectionKalman")
}
dirnames <- outer(lambda, agrid, par2name)
dirnames2 <- paste0(dirnames, "Emp")

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat1 <- map_dfr(dirnames, load_from_dir)
fdat11 <- map_dfr(dirnames2, load_from_dir)
locations_to_exclude <- c("78", "72", "69", "66", "60")

## train_data will be used to select hyperparameters for each location, and variable
train_data <- fdat11 %>% filter(forecast_date <= "2020-11-30")

truth_data1 <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         locations = unique(train_data$location))

truth_data2 <- load_truth(truth_source = "HealthData",
                          target_variable = "inc hosp",
                          locations = unique(fdat1$location))

truth_data3 <- load_truth(truth_source = "JHU",
                          target_variable = "inc death",
                          locations = unique(train_data$location))

truth_data <- bind_rows(truth_data1, truth_data2, truth_data3)

train_scores <- 
  score_forecasts(train_data, truth_data, return_format = "wide") %>%
  filter(!location %in% locations_to_exclude) %>%
  select(model, horizon, location, target_variable, target_end_date, 
         coverage_50, coverage_95, abs_error, wis)

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
