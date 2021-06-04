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

lambda <- 20
par2name <- function(lambda){
  paste0("lambda", sprintf("%06.2f", lambda), 
         "-status-quo",
         "-CEID-InfectionKalman")
}
covid_hub_forecaster_name <- purrr::map(lambda, par2name)

agg_fcsts <- function(chname) {
  pat <- sprintf("(20\\d{2}-\\d{2}-\\d{2})-fips\\d{2}/%s.csv", chname)
  out <- dir("forecasts", recursive = TRUE) %>% 
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
