#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

JuliaCall::julia_setup("/opt/julia-1.5.3/bin")
JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
JuliaCall::julia_eval("using DataFrames")

## main functions


## main script

forecast_date <- Sys.getenv("fdt", unset = "2021-03-29")
forecast_loc <- Sys.getenv("loc", unset = "36")

hopdir <- file.path("hopkins", forecast_date)
tdat <- load_hopkins(hopdir, weekly = FALSE)
ltdat <- tdat %>% filter(location == forecast_loc) %>% 
  filter(target_type == "day ahead inc case" | target_type == "day ahead inc death")
ltdat2 <- ltdat %>% mutate(time = lubridate::decimal_date(target_end_date)) %>%
  pivot_wider(names_from = target_type, values_from = value)

jhu_data <- ltdat2 %>% ungroup() %>% 
  mutate(wday = lubridate::wday(target_end_date)) %>%
  rename(cases = `day ahead inc case`, deaths = `day ahead inc death`) %>%
  select(target_end_date, time, wday, cases, deaths)

healthd <- file.path("healthdata", forecast_date, forecast_loc, "epidata.csv")
cov_thresh <- .5
tdat2 <- read_csv(healthd,
         col_types = cols_only(date = col_date("%Y%m%d"),
         previous_day_admission_adult_covid_confirmed = col_integer(),
         previous_day_admission_pediatric_covid_confirmed = col_integer(),
         previous_day_admission_adult_covid_confirmed_coverage = col_integer())
)

most_recent_coverage <- tdat2 %>% arrange(date) %>%
  pull(previous_day_admission_adult_covid_confirmed_coverage) %>%
  tail(n = 1)

tdat3 <- tdat2 %>%
  filter(previous_day_admission_adult_covid_confirmed_coverage >
           most_recent_coverage * cov_thresh) %>%
  filter(previous_day_admission_adult_covid_confirmed_coverage <
           most_recent_coverage / cov_thresh) %>%
  mutate(
    hospitalizations = previous_day_admission_adult_covid_confirmed +
      previous_day_admission_pediatric_covid_confirmed,
    target_end_date = date - lubridate::ddays(1)
  ) %>%
  select(target_end_date, hospitalizations)

mob_path <- file.path("covidcast-safegraph-home-prop-7dav",
                      forecast_date,
                      "epidata.csv")
abbr <- covidcast::fips_to_abbr(forecast_loc) %>% tolower()
mob_ts <- read_csv(mob_path, 
                   col_types = cols_only(geo_value = col_character(),
                                         time_value = col_date(),
                                         value = col_double())) %>% 
  filter(geo_value == abbr) %>% 
  select(time_value, value) %>%
  rename(target_end_date = time_value, 
         prophome = value)


vacc_path <- file.path("hopkins-vaccine",
                       forecast_date,
                       "vaccine_data_us_timeline.csv")
if (file.exists(vacc_path)) {
  vacc_ts <- read_csv(
    vacc_path,
    col_types = cols_only(
      FIPS = col_double(),
      Vaccine_Type = col_character(),
      Doses_admin = col_double(),
      Date = col_date()
    )
  ) %>%
    mutate(FIPS = sprintf("%02d", FIPS)) %>%
    filter(FIPS == forecast_loc & Vaccine_Type == "All") %>%
    select(Date, Doses_admin) %>%
    rename(target_end_date = Date, doses = Doses_admin)
  
  i <- 1 #replace leading NAs with 0
  while (is.na(vacc_ts$doses[i])) {
    vacc_ts$doses[i] <- 0
    i <- i + 1
  }
  
  obs_data <- left_join(jhu_data, tdat3, by = "target_end_date") %>%
    left_join(vacc_ts, by = "target_end_date") %>%
    left_join(mob_ts, by = "target_end_date") %>%
    mutate(prophome = lead(prophome, 3),
           prophome = zoo::na.fill(prophome, "extend"))
           
  obs_data$doses[lubridate::year(obs_data$target_end_date) == 2020] <-
    0
} else {
  stop("Fitting not implemented for dates without vaccine data")
}

#wind <- obs_data %>% slice(match(1, obs_data$cases > 0):n())
wind <- obs_data %>% slice(370:n())
wsize <- nrow(wind)
N <- covidHubUtils::hub_locations %>% filter(fips == forecast_loc) %>% 
  pull(population)

wfixed <- c(
  N = N,
  rho1 = 0.4,
  gamma = 365.25 / 4,
  gamma_h = 365.25 / 10,
  gamma_d = 365.25 / 10,
  eta = 365.25 / 4,
  t0 = min(wind$time) - 0.00273224
)

y <- wind %>% select(cases, hospitalizations, deaths)
x <- wind %>% select(time, doses, prophome)
x$dosesiqr <- x$doses /  diff(quantile(x$doses, c(.25, .75)))
x$prophomeiqr <- x$prophome / diff(quantile(x$prophome, c(.25, .75)))
#x$doses <- rnorm(x$doses)



winit <- initialize_estimates(x = x, y = y, wfixed = wfixed)

tictoc::tic("optimization")

fit <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x,
  betasd = 0.1,
  epsilon = 1e-1,
  max_iterations = 100,
  a = .9,
  y = y,
  pm = param_map,
  winit,
  invisible = 0
)



fit_over_betagrid <- function(a, betasdgrid) {
  fits <- list()
  fits[[1]] <- lbfgs::lbfgs(
    calc_kf_nll,
    calc_kf_grad,
    x = x,
    betasd = betasdgrid[1],
    epsilon = 1e-1,
    max_iterations = 2,
    a = a,
    y = y[, ],
    pm = param_map,
    winit,
    invisible = 0,
    orthantwise_c = 0,
    orthantwise_start = 10,
    orthantwise_end = 10
  )
  for (i in seq(2, length(betasdgrid))) {
    fits[[i]] <- lbfgs::lbfgs(
      calc_kf_nll,
      calc_kf_grad,
      x = x,
      betasd = betasdgrid[i],
      epsilon = 1e-4,
      max_iterations = 20,
      a = a,
      y = y[, ],
      pm = param_map,
      fits[[i - 1]]$par,
      invisible = 0,
      orthantwise_c = 0,
      orthantwise_start = 10,
      orthantwise_end = 10
    )
  }
  return(fits)
}

betasdgrid <- seq(0.001, 0.1, len = 10)
agrid <- c(0.94, 0.95)

fits <- map(agrid, fit_over_betagrid, betasdgrid = betasdgrid)
tictoc::toc()

ti <- max(wind$target_end_date) + lubridate::ddays(1)
fdtw <- lubridate::wday(forecast_date)
if (fdtw > 2){
  extra <- 7 - fdtw ## extra days to ensure we simulate up to the end of the 4 week ahead target
} else {
  extra <- 0
}
tf <- lubridate::ymd(forecast_date) + lubridate::ddays(28 + extra)
target_end_dates <- seq(from = ti, to = tf, by = "day")

target_end_times <- lubridate::decimal_date(target_end_dates)
target_wday <- lubridate::wday(target_end_dates)

fet <- tibble(target_end_times, target_wday, target_end_dates)

fit_dir <-
  file.path(
    "fits",
    paste0(
      forecast_date,
      "-fips",
      forecast_loc))

if (!dir.exists(fit_dir))
  dir.create(fit_dir, recursive = TRUE)

the_file <- file.path(fit_dir, "fit.RData")
save(x, y, wfixed, fits, fet, agrid, betasdgrid, file = the_file)

## 


dets <- kf_nll_details(w=fit$par, x=x, y=y, betasd = .1, a = 0.9, pm = param_map, fet = NULL)

par(mfrow = c(3, 1))
qqnorm(dets$ytilde_k[1, ] / sqrt(dets$S[1, 1, ]), sub = "Cases")
abline(0, 1)
qqnorm(dets$ytilde_k[2, ] / sqrt(dets$S[2, 2, ]), sub = "Hospitalizations")
abline(0, 1)
qqnorm(dets$ytilde_k[3, ] / sqrt(dets$S[3, 3, ]), sub = "Deaths")
abline(0, 1)

rho_t <- detect_frac(365.25 * (x$time - 2020.164))
plot(x$time, y$cases, xlab = "Time", ylab = "Cases")
pred_cases <- dets$xhat_kkmo["C", ] * rho_t + dets$xhat_kkmo["Hnew", ]
est_cases <- dets$xhat_kkmo["C", ] + dets$xhat_kkmo["Hnew", ]
se_cases <- sqrt(dets$S[1, 1, ])
lines(x$time, se_cases * 2 + pred_cases, col = "grey")
lines(x$time, pred_cases)
lines(x$time, est_cases, lty = 2)
lines(x$time,-se_cases * 2 + pred_cases, col = "grey")

plot(x$time, y$hospitalizations, xlab = "Time",
     ylab = "Hospitalizations")
pred_hosps <- dets$xhat_kkmo["Hnew", ]
se_hosps <- sqrt(dets$S[2, 2, ])
lines(x$time, se_hosps * 2 + pred_hosps, col = "grey")
lines(x$time, pred_hosps)
lines(x$time,-se_hosps * 2 + pred_hosps, col = "grey")

plot(x$time, y$deaths, xlab = "Time",
     ylab = "Deaths")
pred_deaths <- dets$xhat_kkmo["Drep", ]
se_deaths <- sqrt(dets$S[3, 3, ])
lines(x$time, se_deaths * 2 + pred_deaths, col = "grey")
lines(x$time, pred_deaths)
lines(x$time,-se_deaths * 2 + pred_deaths, col = "grey")

plot(exp(fit$par[-c(1:10)] + fit$par[9] * x$dosesiqr + fit$par[10] * x$prophomeiqr) / wfixed["gamma"], type = 'l', ylim = c(0, 2))
lines(exp(fit$par[-c(1:10)] + fit$par[9] * x$dosesiqr + 0 * x$prophomeiqr) / wfixed["gamma"], col = 2)
lines(exp(fit$par[-c(1:10)] + 0 * x$dosesiqr + fit$par[10] * x$prophomeiqr) / wfixed["gamma"], col = 3)

