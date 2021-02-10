#!/usr/bin/env Rscript

library(tidyverse)
source("covidhub-common.R")

combine <- function(date, fnames, chubname){
  outname <- paste0(date, "-", chubname, ".csv")
  srcs <- file.path("forecasts", fnames)
  dest <- file.path(chubname, outname)
  comb <- map(srcs, read_forecast) %>% bind_rows()
  write_csv(comb, dest)
}

lambda <- c(158.49, 107.98, 73.56, 50.12, 34.15, 23.26, 15.85, 10.8)
covid_hub_forecaster_name <- paste0("lambda", sprintf("%06.2f", lambda), "-CEID-InfectionKalman")

agg_fcsts <- function(chname) {
  pat <- sprintf("(20\\d{2}-\\d{2}-\\d{2})-fips\\d{2}-%s.csv", chname)
  out <- dir("forecasts") %>% 
    str_subset(chname) %>% 
    stringr::str_match_all(pat)
    
  fdt <- purrr::map_chr(out, 2) %>% lubridate::as_date()
  fname <- purrr::map_chr(out, 1)
  
  splt <- split(fname, fdt)
  if (!dir.exists(chname)) {
    dir.create(chname)
  }
  map2(names(splt), splt, combine, chubname = chname)
}

map(covid_hub_forecaster_name, agg_fcsts)