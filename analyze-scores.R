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
  select(model, forecast_date, horizon, location, target_variable, 
         target_end_date, coverage_50, coverage_95, abs_error, wis)

ssums <- scores %>% 
  filter(forecast_date > "2020-10-01") %>%
  filter(target_variable %in% c( "inc case", "inc death")) %>% 
  group_by(horizon, model, target_variable) %>% 
  summarize(meanw = mean(wis))

ssumsh <- scores %>% 
  filter(forecast_date > "2020-12-07") %>%
  filter(target_variable == "inc hosp") %>% 
  group_by(horizon, model) %>% 
  summarize(meanw = mean(wis))

ssums %>% ggplot(aes(x = as.integer(horizon), y = meanw, color = model)) + 
  geom_point() + geom_line() +
  ggthemes::scale_colour_colorblind(name = "Model") +
  facet_grid(target_variable~., scales = "free_y") +
  labs(x = "Forecast horizon (weeks)", y = "Mean weighted interval score") +
  theme_minimal() + 
  theme(legend.position="top")

ssumsh %>% 
  ggplot(aes(x = as.integer(horizon), y = meanw, color = model)) + 
  geom_point() + 
  geom_line() + 
  labs(x = "Forecast horizon (days)", y = "Mean weighted interval score") +
  theme_minimal() + 
  theme(legend.position="top") + 
  ggthemes::scale_colour_colorblind(name = "Model")