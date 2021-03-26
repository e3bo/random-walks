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
        p$logdoseeffect,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    nll <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.obj",
      "(pvar, cov, z; N = N, η = η, γ = γ, a = a, betasd = betasd, just_nll = true)"
    ))
  } else if(ncol(y) == 1 && "cases" %in% names(y)) {
    pvar <-
      c(
        p$logE0,
        p$logtauc,
        p$logdoseeffect,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
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
      "(pvar, cov, z; N = N, η = η, γ = γ, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τh, chp = chp, hfp = hfp, a = a, betasd = betasd, just_nll = true)"
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
        p$logdoseeffect,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    g <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.grad",
      "(pvar, cov, z; N = N, η = η, γ = γ, a = a, betasd = betasd, just_nll = true)"
    ))
  } else if(ncol(y) == 1 && "cases" %in% names(y)) {
    pvar <-
      c(
        p$logE0,
        p$logtauc,
        p$logdoseeffect,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
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
      "(pvar, cov, z; N = N, η = η, γ = γ, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τd, chp = chp, hfp = hfp, a = a, betasd = betasd, just_nll = true)"
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

forecast_date <- Sys.getenv("fdt", unset = "2021-03-20")
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

vacc_ts <- read_csv("vaccine_data_us_timeline.csv") %>%
  filter(FIPS == forecast_loc & Vaccine_Type == "All") %>%
  select(Date, Doses_admin) %>%
  rename(target_end_date=Date, doses = Doses_admin)

i <- 1 #replace leading NAs with 0
while(is.na(vacc_ts$doses[i])){
  vacc_ts$doses[i] <- 0
  i <- i + 1
}

obs_data <- left_join(jhu_data, tdat3, by = "target_end_date") %>% 
  left_join(vacc_ts, by = "target_end_date")
obs_data$doses[lubridate::year(obs_data$target_end_date) == 2020] <- 0

wind <- obs_data %>% slice(match(1, obs_data$cases > 0):n())
#wind <- obs_data %>% slice(match(1, obs_data$cases > 0):347)
#wind <- obs_data %>% slice(100:200)
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
x <- wind %>% select(time, doses)

tau_cases_init <- max(var(y$cases, na.rm = TRUE), 1)
tau_hosp_init <- max(var(y$hospitalizations, na.rm = TRUE), 1)
tau_deaths_init <- max(var(y$deaths, na.rm = TRUE), 1)

E0init <-
  ((mean(y$cases) / wfixed["rho1"]) * (365.25 / wfixed["eta"])) %>% unname()

hosp_obs <- !is.na(y$hospitalizations)
chp_init <- sum(y$hospitalizations, na.rm = TRUE) / sum(y$cases[hosp_obs], na.rm = TRUE)
hfp_init <- max(sum(y$deaths[hosp_obs], na.rm = TRUE) / sum(y$hospitalizations, na.rm = TRUE), 0.01)

H0init <- hfp_init / mean(y$deaths)
logdoseeffect_init = log(log(wfixed["N"]) - log(wfixed["N"] - 1)) # first dose reduces susceptibility by 1/N

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
  logdoseeffect = logdoseeffect_init,
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
fdtw <- lubridate::wday(forecast_date)
if (fdtw > 2){
  extra <- 7 - fdtw ## extra days to ensure we simulate up to the end of the 4 week ahead target
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
