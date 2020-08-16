#!/usr/bin/env Rscript

# Setup environment

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")
ddt <- Sys.getenv("ddt")
targ_dir <- file.path("hopkins", ddt)
output_dir <- "forecasts"
datf <- load_hopkins(targ_dir) %>% rename(true_value = "value")
forecast_paths <- dir(output_dir, full.names = TRUE)

# Calculate errors

score <- function(path, tdf2) {
  fcst <- read_forecast(path) %>% mutate(target_type = str_remove(target, "^\\d+ "))
  
  ptdf <- fcst %>% filter(type == "point") %>%
    left_join(tdf2, by = c("target_end_date", "target_type", "location")) %>%
    filter(!is.na(true_value)) %>%
    mutate(error = value - true_value)

  probdf <- fcst %>% filter(type == "quantile") %>%
    left_join(tdf2, by = c("target_end_date", "target_type", "location")) %>%
    filter(!is.na(true_value)) %>%
    mutate(is_below = true_value < value,
           brier_loss = (is_below - quantile) ^ 2) %>%
    mutate(pinball_loss = purrr::map2_dbl(true_value - value, quantile, 
                                          verification::check.func))
  
  probsumdf <- probdf %>% group_by(forecast_date, target) %>% 
    summarise(mean_quantile_score = mean(pinball_loss),
           mean_brier_score = mean(brier_loss), .groups = "drop")
  
  list(point = ptdf, prob = probdf, probsum = probsumdf)
}

scores <- map(forecast_paths, score, tdf2 = datf)

residuals <- 
  map(scores, "prob") %>% 
  bind_rows() %>%
  mutate(loc_type = case_when(location == "US" ~ "national", 
                              nchar(location) == 2 ~ "state", 
                              TRUE ~ "county" ))

# Summarize errors

summary <- list()
summary$by_loc_type <-
  residuals %>%
  group_by(loc_type) %>%
  summarise(mean_qs = mean(pinball_loss), .groups = "drop")

summary$by_loc_type_targ_type <- 
  residuals %>% 
  group_by(loc_type, target_type) %>%
  summarise(mean_qs = mean(pinball_loss), .groups = "drop")

summary$by_loc_targ_fdt <-
  residuals %>%
  group_by(loc_type, target_type, forecast_date) %>%
  summarise(mean_qs = mean(pinball_loss), .groups = "drop")

# Create output

dir.create("metrics")
resids_path <- file.path("metrics", paste0(ddt, "-residuals.rds"))
saveRDS(residuals, resids_path)

summary_plot_path <- file.path("metrics", paste0(ddt, "-score-by-loc-type.csv"))
write_csv(summary$by_loc_type, path = summary_plot_path)

summary_plot_path <- file.path("metrics", paste0(ddt, "-score-by-loc-type-targ-type.csv"))
write_csv(summary$by_loc_type_targ_type, path = summary_plot_path)

summary_plot_path <- file.path("metrics", paste0(ddt, "-score-by-loc-type-targ-type-forecast-date.csv"))
write_csv(summary$by_loc_targ_fdt, path = summary_plot_path)
