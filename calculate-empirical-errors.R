#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

lambda <- 1 / seq(0.001, 0.1, length.out = 10)
dirnames <- paste0("lambda", sprintf("%06.2f", lambda), "-CEID-InfectionKalman")

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat1 <- map_dfr(dirnames, load_from_dir)

truth_data <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         locations = unique(fdat1$location))

locations_to_exclude <- c("78", "72", "69", "66", "60")

resids <- dplyr::left_join(
  x = fdat1,
  y = truth_data,
  by = c("location",
         "target_variable", "target_end_date")
) %>% dplyr::select(-c("model.y")) %>%
  dplyr::rename(model = model.x,
                prediction = value.x,
                true_value = value.y) %>% dplyr::filter(!is.na(true_value)) %>%
  filter(type == "point") %>%
  mutate(error = true_value - prediction) %>%
  filter(!location %in% locations_to_exclude) %>%
  select(model, location, horizon, target_variable, target_end_date, true_value, 
         error) %>%
  mutate(horizon = as.integer(horizon))

m <- lm(log(error^2 + 1) ~ model + location + horizon + true_value, 
        data = resids)

update_pred_intervals <- function(mod_name, sd_model){
  files <- dir(mod_name, full.names = TRUE)
  new_mod_name <- paste0(mod_name, "Emp")
  for (file in files){
    fcst <- read_forecast(file) %>%
      filter(!location %in% locations_to_exclude)
    fcst2 <- fcst %>% filter(type == "point") %>% 
      mutate(horizon = substring(target, 1, 1) %>% as.integer())
    pred <- fcst2 %>% select(location, value, horizon, target) %>% 
      rename(true_value=value)
    pred$model <- mod_name
    pred$sd <- sqrt(exp(predict(m, newdata = pred)))
    fcst3 <- left_join(fcst, select(pred, location, target, true_value, sd), 
                       by = c("target", "location"))
    fcst4 <- fcst3 %>% 
      mutate(value_new = ifelse(type == "quantile", 
                                qnorm(p = quantile, mean = true_value, sd = sd),
                                value),
             value_new = ifelse(value_new < 0, 0, value_new)) %>%
      select(-value, -true_value, -sd) %>%
      rename(value=value_new)
    if (!dir.exists(new_mod_name))
      dir.create(new_mod_name)
    new_file_name <- str_replace_all(file, "CEID-InfectionKalman", 
                                     "CEID-InfectionKalmanEmp")
    write_csv(fcst4, new_file_name)
  }
}

map(dirnames, update_pred_intervals, sd_model = m)