#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

data_date <- Sys.getenv("ddt") %>% as.Date()
resids_name <- paste0(data_date, "-residuals.rds")
model_dirs <- c("drift", "no-drift")
score_paths <- file.path(model_dirs, "metrics", resids_name)

scores <- 
  map(score_paths, readRDS)
names(scores) <-  model_dirs

scores_weights <- 
  bind_rows(scores, .id = "model_id") %>%
  filter(data_date - forecast_date < 30) %>%
  filter(str_detect(target, "^1 wk")) %>%
  group_by(forecast_date, target, target_type, location, quantile, loc_type) %>% 
  summarise(n_models = n(),
            best_model = model_id[which.min(pinball_loss)],
            .groups = "drop")
  
weights <- 
  scores_weights %>% 
  group_by(target_type, location, quantile) %>% 
  count(best_model, .groups = "keep") %>%
  summarise(weight = n / sum(n), 
            model = best_model, .groups = "drop")
  
wts_name <- paste0(data_date, "-weights.csv")
weight_dir <- "model-weights"
wts_path <- file.path("model-weights", wts_name)

if (!dir.exists(weight_dir)) {
  dir.create(weight_dir)
}
write_csv(weights, path = wts_path)
