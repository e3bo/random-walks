#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

forecast_date <- Sys.getenv("fdt", unset = "2020-10-12")
forecast_loc <- "36"
hopdir <- file.path("hopkins", forecast_date)
tictoc::tic("data loading")
tdat <- load_hopkins(hopdir, weekly = FALSE) 
tictoc::toc()

nys <- tdat %>% filter(location == forecast_loc) %>% 
  filter(target_type == "day ahead inc case")

nys2 <- nys %>% mutate(time = lubridate::decimal_date(target_end_date))
case_data <- nys2 %>% ungroup() %>% 
  mutate(wday = lubridate::wday(target_end_date)) %>%
  select(target_end_date, time, wday, value) %>% 
  rename(reports = value)

data_fname = paste0("data--", forecast_date, "--", forecast_loc, ".csv")
write_csv(case_data, path = data_fname)

wsize <- 60
gamma <- 365.25/9

pvar_df <- tribble(
  ~par, ~init, ~lower, ~upper,
  "E_0", 1e4, 10, 1e5,
  "I_0", 1e4, 10, 1e5,
  "tau", 1e-2, 1e-4, 1e3,
  "betasd", 1, 1e-8, 5
) %>% 
  bind_rows(tibble(par = paste0("rho", seq(2, 2)),
                   init = 0.4,
                   lower = 0,
                   upper = 1)) %>%
  bind_rows(tibble(par = paste0("b", seq_len(wsize)),
                   init = gamma,
                   lower = 0.1 * gamma,
                   upper = 4 * gamma))

init_fname <- paste0("initial-pars--", forecast_date, "--", forecast_loc, ".csv")
write_csv(pvar_df, path = init_fname)
