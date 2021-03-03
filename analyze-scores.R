#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")
lambda <- 1 / seq(0.001, 0.1, length.out = 10)

dirnames <- paste0("lambda", sprintf("%06.2f", lambda), "-CEID-InfectionKalman")

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat1 <- map_dfr(dirnames, load_from_dir)
fdat2 <- bind_rows(fdat0, fdat1)

truth_data <- load_truth(truth_source = "JHU",
                         target_variable = "inc case",
                         locations = unique(fdat2$location))

locations_to_exclude <- c("78", "72", "69", "66", "60")

scores <- score_forecasts(fdat2, truth_data, return_format = "long") %>%
  filter(!location %in% locations_to_exclude)
ws <- scores %>% filter(score_name == "wis")
as <- scores %>% filter(score_name == "abs_error")


ws1 <-
  ws %>% group_by(horizon, location, model, target_end_date) %>% 
  summarize(meanscore = mean(score_value), .groups = "drop")

ws1

ws2 <- ws %>% group_by(horizon, location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop") %>%
  group_by(location, horizon) %>% 
  mutate(relscore = meanscore / meanscore[str_detect(model, "CEID-Walk")])

ws2
write_csv(ws2, "wis-horizon-location-model.csv")

as2 <- as %>% group_by(horizon, location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop")
as2


ws3 <- ws %>% group_by(location, model) %>%
  summarize(meanscore = mean(score_value), .groups = "drop") %>%
  group_by(location) %>% 
  mutate(relscore = meanscore / meanscore[str_detect(model, "CEID-Walk")])
ws3
#ws3 %>% ggplot(aes(x = model, y = meanscore)) + 
write_csv(ws3, "wis-model-location.csv")

wscv <- ws3 %>%
  filter(str_detect(model, "^lambda")) %>%
  group_by(location) %>%
  summarize(
    cv_lambda = model[which.min(meanscore)],
    model = "CEID-InfectionKalman",
    meanscore = min(meanscore), .groups = "drop"
  )

ws4 <- wscv %>% 
  select(-cv_lambda) %>% 
  bind_rows(ws3) %>%
  group_by(model) %>% 
  summarize(meanscore = mean(meanscore), .groups = "drop")
write_csv(ws4, "wis-model.csv")

ws4 %>% ggplot(aes(x = model, y = meanscore)) + geom_col() + labs(y = "wis") +
  scale_x_discrete(guide = guide_axis(angle = 45))