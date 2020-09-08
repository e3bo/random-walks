#! /usr/bin/env Rscript

# converts scored submissions to mixing weights which minimize pinball loss

tictoc::tic()

# Load packages and define local functions---------------------------------

suppressPackageStartupMessages(library(tidyverse))

# Read configuration ------------------------------------------------------

model_dirs <- c("drift", "no-drift")
score_paths <- file.path(model_dirs, "metrics", "2020-08-03-residuals.rds")

# Run analysis ------------------------------------------------------------

scores <- 
  map(score_paths, readRDS)
names(scores) <-  model_dirs

scores_weights <- bind_rows(scores, .id = "model_id") %>%
  group_by(forecast_date, target, target_type, location, quantile, loc_type) %>% 
  summarise(n_models = n(),
            best_model = model_id[which.min(pinball_loss)],
            .groups = "drop")
  
weights <- scores_weights %>% group_by(target_type, loc_type) %>% 
  count(best_model, .groups = "keep") %>%
  summarise(weight = n / sum(n), 
            model = best_model, .groups = "drop")
  
# Produce outputs ---------------------------------------------------------

write_csv(weights, path = "weights.csv")
time <- tictoc::toc()
walltime <- list(wall = time$toc - time$tic)
jsonlite::write_json(walltime, "optimal-weight-calc-time.json")
