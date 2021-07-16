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
fdat2 <- bind_rows(fdat0, fdat1) %>% 
  mutate(model = ifelse(model == "lambda020.00-status-quo-CEID-InfectionKalman", 
                        "GISST", model))
locations_to_exclude <- c("78", "72", "69", "66", "60")

truth_data1 <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         data_location = "local_hub_repo",
                         local_repo_path = "~/hub",
                         locations = unique(fdat1$location))

truth_data2 <- load_truth(truth_source = "HealthData",
                          target_variable = "inc hosp",
                          data_location = "local_hub_repo",
                          local_repo_path = "~/hub",
                          locations = unique(fdat1$location))

truth_data3 <- load_truth(truth_source = "JHU",
                          target_variable = "inc death",
                          data_location = "local_hub_repo",
                          local_repo_path = "~/hub",
                          locations = unique(fdat1$location))

truth_data <- bind_rows(truth_data1, truth_data2, truth_data3)

scores <- 
  score_forecasts(fdat2, truth_data, return_format = "wide") %>%
  filter(!location %in% locations_to_exclude) %>%
  select(model, forecast_date, horizon, location, target_variable, 
         target_end_date, coverage_50, coverage_95, abs_error, wis) %>%
  mutate(facet_var = fct_recode(target_variable, 
                                "Cases" = "inc case",
                                "Deaths" = "inc death",
                                "Hospital admissions" = "inc hosp")) %>%
  filter(model != "CEID-Walk") %>% 
  group_by(forecast_date, horizon, location, target_variable) %>%
  mutate(nmodels = n()) %>%
  filter(nmodels == 3 | (target_variable == "inc hosp" & nmodels == 2))

add_theme_mods <- function(plt){
  plt +
  ggthemes::scale_colour_colorblind(name = "Model") +
  theme_minimal() + 
  theme(legend.position = "top")
}

pscores_cd_by_date <- (
  scores %>%
    filter(target_variable %in% c("inc case", "inc death")) %>%
    ggplot(aes(
      x = target_end_date, y = wis, color = model
    )) +
    geom_point() +
    scale_y_log10() +
    facet_grid(horizon ~ facet_var, scales = "free_y") +
    labs(x = "Observation date", y = "Weighted interval score")
) %>%
  add_theme_mods()

ggsave("cases-deaths-wis-by-date.png", pscores_cd_by_date, width = 5.2, 
       height = 7.5, dpi = 600)

pscores_h_by_date <-
  (
    scores %>%
      filter(target_variable %in% c("inc hosp")) %>%
      ggplot(aes(
        x = target_end_date, y = wis, color = model
      )) +
      geom_point() +
      scale_y_log10(breaks = c(10, 1000)) +
      coord_cartesian(clip = 'off') +
      facet_grid(as.integer(horizon) ~ .) +
      labs(x = "Observation date", y = "Weighted interval score")
  ) %>%
  add_theme_mods()
ggsave("hosp-wis-by-date.png", pscores_h_by_date, width = 5.2, dpi = 600)

ssums <- scores %>%
  filter(facet_var %in% c("Cases", "Deaths")) %>% 
  group_by(horizon, model, facet_var) %>% 
  summarize(meanw = mean(wis), .groups = "drop")

ssumsh <- scores %>% 
  filter(target_variable == "inc hosp") %>% 
  group_by(horizon, model) %>% 
  summarize(meanw = mean(wis), .groups = "drop")

pscores_cd <-
  (
    ssums %>% ggplot(aes(
      x = as.integer(horizon),
      y = meanw,
      color = model
    )) +
      geom_point() + geom_line() +
      facet_grid(facet_var ~ ., scales = "free_y") +
      labs(x = "Forecast horizon (weeks)", y = "Mean weighted interval score")
  ) %>%
  add_theme_mods()
ggsave("cases-deaths-wis.png", pscores_cd, width = 5.2, height=4, dpi = 600)

pscores_h <- (
  ssumsh %>%
    ggplot(aes(
      x = as.integer(horizon),
      y = meanw,
      color = model
    )) +
    geom_point() +
    geom_line() +
    labs(x = "Forecast horizon (days)", y = "Mean weighted interval score")
) %>%
  add_theme_mods()

ggsave("hosp-wis.png", pscores_h, width = 5.2, dpi = 600)