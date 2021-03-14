#!/usr/bin/env Rscript

library(epidatr)
library(dplyr)

forecast_date <- Sys.getenv("fdt", unset = "2020-11-16")
forecast_loc <- Sys.getenv("loc", unset = "36")

abb <- covidHubUtils::hub_locations %>% filter(fips == forecast_loc) %>% 
  pull(abbreviation)

issues <- forecast_date %>% stringr::str_remove_all("-")
dates <- paste0("20200101-", issues)

req <- healthdata(dates = dates, states = abb, issues = issues)

outdir <- file.path("healthdata", forecast_date, forecast_loc) 
if (!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
}
outname <- file.path(outdir, "epidata.csv")
readr::write_csv(req$epidata, outname)
