#!/usr/bin/env Rscript

library(covidHubUtils)

models <- c("COVIDhub-ensemble", "CEID-Walk")
locs <-
  c(
    "01",
    "04",
    "06",
    "08",
    "12",
    "13",
    "17",
    "18",
    "21",
    "22",
    "24",
    "25",
    "26",
    "27",
    "29",
    "34",
    "36",
    "37",
    "42",
    "45",
    "47",
    "48",
    "51",
    "53",
    "55"
  )

fdts_mon <- seq.Date(as.Date("2020-07-20"), as.Date("2020-11-30"), by = "7 days")
fdts_tue <- fdts_mon + 1
fdts_sun <- fdts_mon - 1

fdat <- load_forecasts(models = models,
                       forecast_dates = c(fdts_sun, fdts_mon, fdts_tue),
                       locations = locs, 
                       targets = paste(1:4, "wk ahead inc case"))

saveRDS(fdat, file = "other-model-forecasts.rds")
