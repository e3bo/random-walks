#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

JuliaCall::julia_setup("/opt/julia-1.5.3/bin")
JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
JuliaCall::julia_eval("using DataFrames")

## main functions

calc_kf_nll <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  if (ncol(y) == 3) {
    pvar <-
      c(
        p$logE0,
        p$logH0,
        p$logtauc,
        p$logtauh,
        p$logtaud,
        p$logchp,
        p$loghfp,
        p$loggammahd,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    nll <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.obj",
      "(pvar, z; N = N, η = η, γ = γ, a = a, betasd = betasd, just_nll = true)"
    ))
  } else if(ncol(y) == 1 && "cases" %in% names(y)) {
    pvar <-
      c(
        p$logE0,
        p$logtauc,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    JuliaCall::julia_assign("h0", exp(p$logH0))
    JuliaCall::julia_assign("τh", exp(p$logtauh))
    JuliaCall::julia_assign("τd", exp(p$logtaud))
    JuliaCall::julia_assign("chp", exp(p$logchp))
    JuliaCall::julia_assign("hfp", exp(p$loghfp))
    nll <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.obj",
      "(pvar, z; N = N, η = η, γ = γ, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τh, chp = chp, hfp = hfp, a = a, betasd = betasd, just_nll = true)"
    ))
  }

  nll
}

calc_kf_grad <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  if (ncol(y) == 3) {
    pvar <-
      c(
        p$logE0,
        p$logH0,
        p$logtauc,
        p$logtauh,
        p$logtaud,
        p$logchp,
        p$loghfp,
        p$loggammahd,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    g <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.grad",
      "(pvar, z; N = N, η = η, γ = γ, a = a, betasd = betasd, just_nll = true)"
    ))
  } else if(ncol(y) == 1 && "cases" %in% names(y)) {
    pvar <-
      c(
        p$logE0,
        p$logtauc,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    JuliaCall::julia_assign("h0", exp(p$logH0))
    JuliaCall::julia_assign("τh", exp(p$logtauh))
    JuliaCall::julia_assign("τd", exp(p$logtaud))
    JuliaCall::julia_assign("chp", exp(p$logchp))
    JuliaCall::julia_assign("hfp", exp(p$loghfp))
    g <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.grad",
      "(pvar, z; N = N, η = η, γ = γ, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τh, chp = chp, hfp = hfp, a = a, betasd = betasd, just_nll = true)"
    ))
  }
  g
}

calc_kf_hess <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  JuliaCall::julia_assign("logE0", p$logE0)
  JuliaCall::julia_assign("logtau", p$logtau)
  JuliaCall::julia_assign("bpars", p$bpars)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("N", p$N)
  JuliaCall::julia_assign("ρ", p$rho1)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("z", data.matrix(y))
  JuliaCall::julia_assign("a", a)
  JuliaCall::julia_assign("betasd", betasd)
  g <- JuliaCall::julia_eval(paste0(
    "InfectionKalman.hess",
    "(logE0, logtau, bpars, z; ρ = ρ, N = N, η = η, γ = γ, a = a, betasd = betasd)"
  ))
  g
}


## main script

forecast_date <- Sys.getenv("fdt", unset = "2020-11-16")
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
  filter(previous_day_admission_adult_covid_confirmed >
           most_recent_coverage * cov_thresh) %>%
  filter(previous_day_admission_adult_covid_confirmed <
           most_recent_coverage / cov_thresh) %>%
  mutate(
    hospitalizations = previous_day_admission_adult_covid_confirmed +
      previous_day_admission_pediatric_covid_confirmed,
    target_end_date = date - lubridate::ddays(1)
  ) %>%
  select(target_end_date, hospitalizations)

obs_data <- left_join(jhu_data, tdat3, by = "target_end_date")

wind <- obs_data %>% slice(match(1, obs_data$cases):n())
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
  t0 = rev(obs_data$time)[wsize + 1]
)

y <- wind %>% select(cases, hospitalizations, deaths)
x <- matrix(wind$time, ncol = 1)

tau_cases_init <- var(y$cases, na.rm = TRUE)
tau_hosp_init <- var(y$hospitalizations, na.rm = TRUE)
tau_deaths_init <- var(y$deaths, na.rm = TRUE)

E0init <- ((mean(y$cases) / wfixed["rho1"]) * (365.25 / wfixed["eta"])) %>% unname()
H0init <- (mean(y$hospitalizations, na.rm = TRUE) * 365.25 / wfixed["gamma_h"]) %>% unname()

chp_init <- sum(y$hospitalizations, na.rm = TRUE) / sum(y$cases, na.rm = TRUE)
hfp_init <- sum(y$deaths, na.rm = TRUE) / sum(y$hospitalizations, na.rm = TRUE)

binit <- unname(rep(log(wfixed["gamma"]), wsize))
names(binit) <- paste0("b", seq_len(wsize))
winit <- c(
  logE0 = log(E0init),
  logH0 = log(H0init),
  logtauc = log(tau_cases_init),
  logtauh = log(tau_hosp_init),
  logtaud = log(tau_deaths_init),
  logchp = log(chp_init),
  loghfp = log(hfp_init),
  loggammahd = log(365.25 / 2),
  binit
)

tictoc::tic("optimization")

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
    invisible = 0
  )
  for (i in seq(2, length(betasdgrid))) {
    fits[[i]] <- lbfgs::lbfgs(
      calc_kf_nll,
      calc_kf_grad,
      x = x,
      betasd = betasdgrid[i],
      epsilon = 1e-4,
      max_iterations = 2,
      a = a,
      y = y[, ],
      pm = param_map,
      fits[[i - 1]]$par,
      invisible = 0
    )
  }
  return(fits)
}

betasdgrid <- seq(0.001, 0.1, len = 10)
agrid <- c(0.94, 0.95)

fits <- map(agrid, fit_over_betagrid, betasdgrid = betasdgrid)
tictoc::toc()

ti <- max(wind$target_end_date) + lubridate::ddays(1)
tf <- lubridate::ymd(forecast_date) + lubridate::ddays(28)
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

saveRDS(x, file.path(fit_dir, "x.rds"))
saveRDS(y, file.path(fit_dir, "y.rds"))
saveRDS(wfixed, file.path(fit_dir, "wfixed.rds"))
saveRDS(fits, file.path(fit_dir, "fits.rds"))
saveRDS(fet, file.path(fit_dir, "fet.rds"))
saveRDS(agrid, file.path(fit_dir, "agrid.rds"))
saveRDS(betasdgrid, file.path(fit_dir, "betasdgrid.rds"))

q("no")
# Diagnostics, requires variables in above script and functions in write-formatted-forecasts.R
fits <- readRDS(file.path(fit_dir, "fits.rds"))
fet <- readRDS(file.path(fit_dir, "fet.rds"))
agrid <- readRDS(file.path(fit_dir, "agrid.rds"))
betasdgrid <- readRDS(file.path(fit_dir, "betasdgrid.rds"))

fitind1 <- 1
fitind2 <- 2
dets <- kf_nll_details(winit, x = x, y = y, param_map, 
                       betasd = betagrid[fitind2], a = agrid[fitind1],
                       fet)

nll <- calc_kf_nll(winit, x = x, y = y, param_map, 
                       betasd = betagrid[fitind2], a = agrid[fitind1])

fitind1 <- 1
fitind2 <- 1
dets <- kf_nll_details(fits[[fitind1]][[fitind2]]$par, x = x, y = y, param_map, 
                       betasd = betasdgrid[fitind2], a = agrid[fitind1],
                       fet)
par(mfrow = c(1,1))
qqnorm(dets$ytilde_k[1,] / sqrt(dets$S[1,1,]))
abline(0, 1)
qqnorm(dets$ytilde_k[2,] / sqrt(dets$S[2,2,]))
abline(0, 1)
qqnorm(dets$ytilde_k[3,] / sqrt(dets$S[3,3,]))
abline(0,1 )

rho_t <- detect_frac(seq_len(nrow(y)))
plot(x[,1], y$cases, xlab = "Time", ylab = "Cases")
pred_cases <- dets$xhat_kkmo["C",] * rho_t + dets$xhat_kkmo["Hnew",]
est_cases <- dets$xhat_kkmo["C",] + dets$xhat_kkmo["Hnew",]
se_cases <- sqrt(dets$S[1,1,])
lines(x[,1], se_cases * 2 + pred_cases, col = "grey")
lines(x[,1], pred_cases)
lines(x[,1], est_cases, lty = 2)
lines(x[,1], -se_cases * 2 + pred_cases, col = "grey")

plot(x[,1], y$hospitalizations, xlab = "Time", 
     ylab = "Hospitalizations")
pred_hosps <- dets$xhat_kkmo["Hnew",]
se_hosps <- sqrt(dets$S[2,2,])
lines(x[,1], se_hosps * 2 + pred_hosps, col = "grey")
lines(x[,1], pred_hosps)
lines(x[,1], -se_hosps * 2 + pred_hosps, col = "grey")

plot(x[,1], y$deaths, xlab = "Time", 
     ylab = "Deaths")
pred_deaths <- dets$xhat_kkmo["Drep",]
se_deaths <- sqrt(dets$S[3,3,])
lines(x[,1], se_deaths * 2 + pred_deaths, col = "grey")
lines(x[,1], pred_deaths)
lines(x[,1], -se_deaths * 2 + pred_deaths, col = "grey")

case_inds <- which(fet$target_wday == 7)
case_fcst <- create_forecast_df(means = dets$sim_means[1, case_inds,],
                                vars = dets$sim_cov[1, 1, case_inds,],
                                location = forecast_loc)

case_fcst %>% 
  ggplot(aes(x = target_end_date, y = as.numeric(value), color = quantile)) + 
  geom_line() + geom_point()

hosp_inds <- fet$target_end_dates %in% (lubridate::ymd(forecast_date) + 1:28)
stopifnot(sum(hosp_inds) == 28)

hosp_fcst <- create_forecast_df(means = dets$sim_means[2, hosp_inds,],
                                vars = dets$sim_cov[2, 2, hosp_inds,],
                                target_type = "hospitalizations",
                                location = forecast_loc)

hosp_fcst %>% 
  ggplot(aes(x = target_end_date, y = as.numeric(value), color = quantile)) + 
  geom_line() + geom_point()

death_fcst <- create_forecast_df(means = dets$sim_means[3, case_inds,],
                                vars = dets$sim_cov[3, 3, case_inds,],
                                target_type = "inc deaths",
                                location = forecast_loc)

death_fcst %>% 
  ggplot(aes(x = target_end_date, y = as.numeric(value), color = quantile)) + 
  geom_line() + geom_point()
