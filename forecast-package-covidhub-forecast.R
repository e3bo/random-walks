#!/usr/bin/env Rscript

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

mconfig <- ini::read.ini("model.ini")
eval_string <- function(x){
  str2lang(x) %>% eval()
}
mconfig$model_pars <- map(mconfig$model_pars, eval_string)

forecast_dates <- Sys.getenv("fdt") %>% as.Date()
hopdir <- file.path("hopkins", forecast_dates)

tdat <- load_hopkins(hopdir) %>% filter(str_detect(location, "^13"))

loglik_naive <- function(y){
  res <- diff(y)
  sigma <- sd(res, na.rm = TRUE)
  densities <- dnorm(res, sd = sigma, log = TRUE)
  sum(densities, na.rm = TRUE)
}

rw_forecast <- function(df, fdt, npaths = 100, 
                        h = 4, tailn = Inf, include_drift = FALSE){
  stopifnot(npaths > 0)
  last_end_date <- max(df$target_end_date)
  stopifnot(fdt - last_end_date < lubridate::ddays(3))
  df1 <- df %>% filter(target_end_date <= fdt)
  y <- ts(df1$value) %>% tail(n = tailn)
  object <- forecast::naive(y, h = 1)$model
  sim <- matrix(NA, nrow = npaths, ncol = h)
  for (i in 1:npaths) {
    sim[i, ] <- simulate(object, nsim = h, bootstrap = FALSE)
  }
  sim1 <- (sim + abs(sim)) / 2 #to zero out any negative values
  sim_dates <- seq(last_end_date + 7, last_end_date + h *7, by = 7)
  colnames(sim1) <- as.character(sim_dates)
  sim12 <- cbind(Rep = seq_len(nrow(sim1)), sim1) %>% as_tibble()
  simdf <- sim12 %>% pivot_longer(-Rep, names_to = "Date", 
                                  values_to = "value")
  simdf %>% mutate(Date = lubridate::ymd(Date))
}

reshape_paths_like_pompout <- function(x){
  x %>% mutate(target_type = recode(target_type, 
                                    `wk ahead inc death`= "deaths",
                                    `wk ahead inc case` = "cases")) %>% 
    pivot_wider(id_cols = c(Rep, Date), names_from = "target_type", 
                values_from = "value")
}

sim_paths <- function(fdt, tdat, mname, odir, ...) {

  tdat1 <- tdat %>% 
    filter(target_type %in% c("wk ahead inc death", "wk ahead inc case")) %>%
    group_by(location, target_type) %>%
    arrange(target_end_date) %>%
    nest() %>% 
    mutate(paths = map(data, rw_forecast, fdt = fdt, ...)) %>%
    select(-data) %>%
    unnest(cols = paths) %>%
    group_by(location) %>%
    nest() %>%
    mutate(paths2 = map(data, reshape_paths_like_pompout)) %>%
    select(-data) %>%
    mutate(fcst_df = map2(paths2, location, paths_to_forecast, hop = tdat))
    
  full <- bind_rows(tdat1$fcst_df)

  fname <- paste0(fdt, "-CEID-", mname, ".csv")
  if (!dir.exists(odir)) {
    dir.create(odir)
  }
  opath <- file.path(odir, fname)
  write_csv(full, opath)
}

paths <- map(forecast_dates,
    sim_paths,
    tdat = tdat,
    npaths = mconfig$model_pars$npaths,
    tailn = mconfig$model_pars$tailn,
    include_drift = mconfig$model_pars$include_drift,
    odir = "forecasts",
    mname = "Walk")

## Produce metrics
time <- tictoc::toc()
walltime <- list(wall = time$toc - time$tic)
dir.create("metrics")
mpath <- file.path("metrics", paste0(forecast_dates, 
                                     "-forecast-calc-time.json"))
jsonlite::write_json(walltime, mpath)
