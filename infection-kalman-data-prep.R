#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

forecast_date <- Sys.getenv("fdt", unset = "2020-10-12")
forecast_loc <- Sys.getenv("loc", unset = "36")
hopdir <- file.path("hopkins", forecast_date)
tictoc::tic("data loading")
tdat <- load_hopkins(hopdir, weekly = FALSE) 
tictoc::toc()

ltdat <- tdat %>% filter(location == forecast_loc) %>% 
  filter(target_type == "day ahead inc case")

ltdat2 <- ltdat %>% mutate(time = lubridate::decimal_date(target_end_date))
case_data <- ltdat2 %>% ungroup() %>% 
  mutate(wday = lubridate::wday(target_end_date)) %>%
  select(target_end_date, time, wday, value) %>% 
  rename(reports = value)

moving_average <- function(x, n = 7) {
  stats::filter(x, rep(1 / n, n), sides = 1)
}

case_data$smooth <- moving_average(case_data$reports)


data_fname = paste0(forecast_date, "--", forecast_loc, ".csv")
if(!dir.exists("data")) dir.create("data")
write_csv(case_data, path = file.path("data", data_fname))

wsize <- 60
gamma <- 365.25/9
tau_init <- case_data$reports %>% tail(n = wsize) %>% var()

pvar_df <- tribble(
  ~par, ~init, ~lower, ~upper,
  "E_0", 1e4, 10, 1e5,
  "I_0", 1e4, 10, 1e5,
  "tau", tau_init, tau_init * 1e-2, tau_init * 10
) %>% 
  bind_rows(tibble(par = paste0("b", seq_len(wsize)),
                   init = gamma,
                   lower = 0.1 * gamma,
                   upper = 4 * gamma))

init_fname <- paste0(forecast_date, "--", forecast_loc, ".csv")
if(!dir.exists("initial-pars")) dir.create("initial-pars")
write_csv(pvar_df, path = file.path("initial-pars", init_fname))
