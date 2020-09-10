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

forecast_date <- Sys.getenv("fdt") %>% as.Date()
hopdir <- file.path("hopkins", forecast_date)

tdat <- load_hopkins(hopdir) #%>% filter(str_detect(location, "^13"))

rw_forecast <- function(df, target_type, fdt, 
                        h = 4, tailn = Inf, include_drift = FALSE){
  last_end_date <- max(df$target_end_date)
  stopifnot(fdt - last_end_date < lubridate::ddays(3))
  df1 <- df %>% filter(target_end_date <= fdt)
  y <- ts(df1$value) %>% tail(n = tailn)
  object <- forecast::rwf(y, drift = include_drift, h = h)$model
  if (target_type == "wk ahead inc case"){
    prob <- c(0.025, 0.100, 0.250, 0.500, 0.750, 0.900, 0.975)
  } else {
    prob <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  }
  z <- qnorm(p = prob)
  nquantiles <- length(z)
  
  # start copying from forecast:::forecast.lagwalk
  lag <- object$par$lag
  fullperiods <- (h - 1) / lag + 1
  steps <- rep(1:fullperiods, rep(lag, fullperiods))[1:h]
  fc <- rep(object$future, fullperiods)[1:h] + steps * object$par$drift
  mse <- mean(object$residuals ^ 2, na.rm = TRUE)
  se <- sqrt(mse * steps + (steps * object$par$drift.se) ^ 2)
  # end copying
  
  quantile_values <- matrix(NA, nrow = h, ncol = nquantiles)
  for (i in 1:nquantiles) {
    quantile_values[, i] <- fc + z[i] * se
  }
  quantile_values <- cbind(quantile_values, fc) # for point forecast
  colnames(quantile_values) <- c(prob, "point")
  qdw <- cbind(tibble("wks_ahead" = steps), as_tibble(quantile_values))
  qdl <- pivot_longer(qdw, -wks_ahead, names_to = "quantile", 
                      values_to = "value") %>%
    mutate(value = ifelse(value < 0, 0, value)) %>%
    mutate(type = ifelse(quantile == "point", "point", "quantile")) %>%
    mutate(quantile = ifelse(quantile == "point", NA, quantile)) %>%
    add_column(forecast_date = fdt) %>%
    mutate(target = paste(wks_ahead, target_type)) %>%
    mutate(target_end_date = wks_ahead * 7 + last_end_date) %>%
    mutate(value = round(value, digits = 2)) %>%
    select(-wks_ahead)
  qdl
}

gen_forecasts <- function(fdt, tdat, mname, odir, ...) {

  full <- tdat %>%
    filter(target_type %in% c("wk ahead inc death", "wk ahead inc case")) %>%
    filter(nchar(location) != 5 | target_type == "wk ahead inc case") %>% # only case forecasts for counties 
    group_by(location, target_type) %>%
    arrange(target_end_date) %>%
    nest() %>%
    mutate(fcsts = map2(data, target_type, rw_forecast, fdt = fdt, ...)) %>%
    select(-data) %>%
    unnest(cols = fcsts) %>%
    ungroup %>%
    select(-target_type) %>%
    relocate(forecast_date,
             target,
             target_end_date,
             location,
             type,
             quantile,
             value)
  
  fname <- paste0(fdt, "-CEID-", mname, ".csv")
  if (!dir.exists(odir)) {
    dir.create(odir)
  }
  opath <- file.path(odir, fname)
  write_csv(full, opath)
}

gen_forecasts(forecast_date,
    tdat = tdat,
    tailn = mconfig$model_pars$tailn,
    include_drift = mconfig$model_pars$include_drift,
    odir = "forecasts",
    mname = "Walk")

## Produce metrics
time <- tictoc::toc()
walltime <- list(wall = time$toc - time$tic)
dir.create("metrics")
mpath <- file.path("metrics", paste0(forecast_date, 
                                     "-forecast-calc-time.json"))
jsonlite::write_json(walltime, mpath, auto_unbox = TRUE)
