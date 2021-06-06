#!/usr/bin/env Rscript

library(covidHubUtils)

models <- c("COVIDhub-ensemble", "COVIDhub-baseline", "CEID-Walk")
locs <-
  c(
    "06"
  )

fdts_mon <- seq.Date(as.Date("2020-06-08"), as.Date("2021-04-26"), by = "7 days")
fdts_tue <- fdts_mon + 1
fdts_sun <- fdts_mon - 1

fdat <- load_forecasts(models = models,
                       forecast_dates = c(fdts_sun, fdts_mon, fdts_tue),
                       locations = locs, 
                       targets = c(paste(1:4, "wk ahead inc case"),
                                   paste(1:4, "wk ahead inc death"),
                                   paste(c(1:28), "day ahead inc hosp"))) 

saveRDS(fdat, file = "other-model-forecasts.rds")
