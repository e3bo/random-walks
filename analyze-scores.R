#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")
fdat1 <-
  dir("lambda125.89-CEID-InfectionKalman", full.names = TRUE) %>% 
  load_forecast_files_repo()
fdat11 <-
  dir("lambda158.49-CEID-InfectionKalman", full.names = TRUE) %>% 
  load_forecast_files_repo()  
fdat2 <- bind_rows(fdat0, fdat1, fdat11)

truth_data <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         locations = unique(fdat2$location))

scores <- score_forecasts(fdat2, truth_data, return_format = "long")
ws <- scores %>% filter(score_name == "wis")

ws1 <-
  ws %>% group_by(horizon, location, model, target_end_date) %>% 
  summarize(meanscore = mean(score_value), .groups = "drop")

ws1
ws1 %>% ggplot(aes(x = target_end_date, y = meanscore, color = model)) + 
  geom_point() + facet_grid(location ~ horizon, scales = "free_y") + 
  labs(x = "target date", y = "wis")

ws2 <- ws %>% group_by(horizon, location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop")
ws2
ws2 %>% ggplot(aes(x = horizon, y = meanscore, fill = model)) + 
  geom_col(position = position_dodge()) + 
  facet_grid(location ~ ., scales = "free_y") + 
  labs(y = "wis")

ws3 <- ws %>% group_by(location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop")
ws3
ws3 %>% ggplot(aes(x = location, y = meanscore, color = model)) + 
  geom_point() + coord_flip() + labs(y = "wis")

ws4 <- ws %>% group_by(model) %>% 
  summarize(meanscore = mean(score_value), .groups = "drop")
ws4
ws4 %>% ggplot(aes(x = model, y = meanscore)) + geom_col() + labs(y = "wis")