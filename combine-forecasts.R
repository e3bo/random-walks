#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

source("covidhub-common.R")

forecast_date <- Sys.getenv("fdt") %>% as.Date()

wts_name <- paste0(forecast_date, "-weights.csv")
weights <- read_csv(
  file.path("model-weights", wts_name),
  col_types = cols(
    target_type = col_character(),
    loc_type = col_character(),
    quantile = col_double(),
    weight = col_double(),
    model = col_character()
  )
)
output_dir <- "forecasts"
model_dirs <- c("drift", "no-drift")

forecast_dirs <- file.path(model_dirs, output_dir)
fcst_name <- paste0(forecast_date, "-CEID-Walk.csv")
paths <- file.path(forecast_dirs, fcst_name)
names(paths) <- model_dirs

fcsts <- 
  map(paths, read_forecast) %>% 
  bind_rows(.id = "model") %>% 
  mutate(target_type = str_remove(target, "^[0-9]+ ")) %>%
  mutate(loc_type = case_when(location == "US" ~ "national", 
                              nchar(location) == 2 ~ "state", 
                              TRUE ~ "county")) %>%
  left_join(weights, by = c("target_type", "loc_type", "model", "quantile")) %>%
  mutate(weight = ifelse(is.na(weight), 0, weight)) %>%
  group_by(forecast_date, target, target_end_date, location, type, quantile) %>%
  summarise(value = sum(value * weight), .groups = "drop") %>%
  mutate(value = round(value, digits = 2))

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
opath <- file.path(output_dir, fcst_name)
write_csv(fcsts, opath)