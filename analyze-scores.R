#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")

lambda <- c(158.49, 107.98, 73.56, 50.12, 34.15, 23.26, 15.85, 10.8)

dirnames <- paste0("lambda", sprintf("%06.2f", lambda), "-CEID-InfectionKalman")

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat1 <- map_dfr(dirnames, load_from_dir)
fdat2 <- bind_rows(fdat0, fdat1)

truth_data <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         locations = unique(fdat2$location))

scores <- score_forecasts(fdat2, truth_data, return_format = "long")
ws <- scores %>% filter(score_name == "wis")
as <- scores %>% filter(score_name == "abs_error")

ws1 <-
  ws %>% group_by(horizon, location, model, target_end_date) %>% 
  summarize(meanscore = mean(score_value), .groups = "drop")

ws1
ws1 %>% ggplot(aes(x = target_end_date, y = meanscore, color = model)) + 
  geom_point() + facet_grid(location ~ horizon, scales = "free_y") + 
  scale_x_date(guide = guide_axis(angle = 90)) + 
  labs(x = "target date", y = "wis")

ws2 <- ws %>% group_by(horizon, location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop")
ws2
ws2 %>% ggplot(aes(x = horizon, y = meanscore, fill = model)) + 
  geom_col(position = position_dodge()) + 
  facet_grid(location ~ ., scales = "free_y") + 
  labs(y = "WIS")

as2 <- as %>% group_by(horizon, location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop")
as2
as2 %>% ggplot(aes(x = horizon, y = meanscore, fill = model)) + 
  geom_col(position = position_dodge()) + 
  facet_grid(location ~ ., scales = "free_y") + 
  labs(y = "MAE")



ws3 <- ws %>% group_by(location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop")
ws3
ws3 %>% ggplot(aes(x = location, y = meanscore, color = model)) + 
  geom_point() + coord_flip() + labs(y = "WIS")

ws4 <- ws %>% group_by(model) %>% 
  summarize(meanscore = mean(score_value), .groups = "drop")
ws4
ws4 %>% ggplot(aes(x = model, y = meanscore)) + geom_col() + labs(y = "wis") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ws$lambda <- ws$model %>%
  str_extract("lambda\\d{3}\\.\\d{2}") %>% str_remove("lambda") %>% as.numeric()

ws5 <- ws %>% filter(str_detect(model, "InfectionKalman")) %>%
  group_by(model, lambda, forecast_date, location) %>%
  summarize(meanscore = mean(score_value), .groups = "drop") %>%
  left_join(hub_locations, by = c("location" = "fips"))

ws5 %>%
  ggplot(aes(
    x = forecast_date,
    color = as.factor(lambda),
    y = meanscore,
    group = model
  )) +
  scale_y_continuous(trans = "log1p") + 
  geom_point() + geom_line() + facet_wrap( ~ location_name, scales = "free_y") + 
  labs(x = "forecast date", y = "WIS", color = expression(lambda))
