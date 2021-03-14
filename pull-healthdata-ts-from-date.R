#!/usr/bin/env Rscript

library(epidatr)
library(tidyverse)

forecast_date <- Sys.getenv("fdt", unset = "2020-11-16")
forecast_loc <- Sys.getenv("loc", unset = "36")

abb <- covidHubUtils::hub_locations %>% 
  filter(fips == forecast_loc) %>% 
  pull(abbreviation)
issue_date <- lubridate::ymd(forecast_date)
out <- tibble()

while (issue_date > "2020-11-02"){ # release date on healthdata.gov is 2020-11-03, although earliest issue currently available seems to be 11/16
  issue <- issue_date %>% stringr::str_remove_all("-")
  dates <- paste0("20200101-", issue)
  req <- try(healthdata(dates = dates, states = abb, issues = issue))
  if (inherits(req, "try-error")){
    issue_date <- issue_date - lubridate::ddays(1)
  } else {
    out <- req$epidata
    break
  }
}

outdir <- file.path("healthdata", forecast_date, forecast_loc) 
if (!dir.exists(outdir)){
  dir.create(outdir, recursive = TRUE)
}
outname <- file.path(outdir, "epidata.csv")
readr::write_csv(out, outname)
