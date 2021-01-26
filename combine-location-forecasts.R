#!/usr/bin/env Rscript

library(tidyverse)
source("covidhub-common.R")

covid_hub_forecaster_name <- "CEID-InfectionKalman"
out <-
  dir("forecasts") %>% stringr::str_match_all(sprintf(
    "(20\\d{2}-\\d{2}-\\d{2})-\\d{2}-%s.csv",
    covid_hub_forecaster_name
  ))
fdt <- purrr::map_chr(out, 2) %>% lubridate::as_date()
fname <- purrr::map_chr(out, 1)

splt <- split(fname, fdt)

if (!dir.exists(covid_hub_forecaster_name)){
  dir.create(covid_hub_forecaster_name)
}

combine <- function(date, fnames, chubname){
  outname <- paste0(date, "-", chubname, ".csv")
  srcs <- file.path("forecasts", fnames)
  dest <- file.path(chubname, outname)
  comb <- map(srcs, read_forecast) %>% bind_rows()
  write_csv(comb, dest)
}

map2(names(splt), splt, combine, chubname = covid_hub_forecaster_name)