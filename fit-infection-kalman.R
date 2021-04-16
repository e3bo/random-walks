#!/usr/bin/env R

## Environment prep

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

JuliaCall::julia_setup("/opt/julia-1.5.3/bin")
JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
JuliaCall::julia_eval("using DataFrames")

## Data prep

forecast_date <- Sys.getenv("fdt", unset = "2021-03-29")
forecast_loc <- Sys.getenv("loc", unset = "17")

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
wind <- obs_data %>% slice(41:n())
wsize <- nrow(wind)
N <- covidHubUtils::hub_locations %>% filter(fips == forecast_loc) %>% 
  pull(population)

case_to_death_lag <- 20

wfixed <- c(
  N = N,
  gamma = 365.25 / 4,
  gamma_h = 365.25 / (case_to_death_lag / 2),
  gamma_d = 365.25 / (case_to_death_lag / 2),
  eta = 365.25 / 4
)

y <- wind %>% select(cases, hospitalizations, deaths)
x0 <- wind %>% select(time, doses, prophome)
x0$dosesiqr <- x0$doses /  diff(quantile(x0$doses[x0$doses > 0], c(.25, .75)))
x0$prophomeiqr <- (x0$prophome - mean(x0$prophome)) / diff(quantile(x0$prophome, c(.25, .75)))
x0$bvecmap <- rep(1:14, each = 28)

set.seed(1)

# estimate the probability that a case is reported assuming that all deaths are reported
cvd <- data.frame(time = x0$time, r = y$cases / lead(y$deaths, case_to_death_lag)) %>% 
  mutate(logr = log(r)) %>% filter(is.finite(logr)) %>% filter(time < 2020.47)
mars <- earth::earth(logr~time, data = cvd)
#plot(logr ~ time, data = cvd)
#lines(cvd$time, predict(mars))
cvd$rhat <- exp(predict(mars)) %>% as.numeric()
right <- cvd %>% select("time", "rhat")

x <- left_join(x0, right, by = "time") %>% 
  mutate(rhat2 = zoo::na.fill(rhat, "extend"),
         rhot = rhat2 / (2 *rhat2[n()]))


#pdf <- data.frame(date = wind$target_end_date, p = x2$rhat3)
#ggplot(pdf, aes(x = date, y = p)) + geom_point() + ylab("Pr (case is reported)")

winit <- initialize_estimates(x = x, y = y, wfixed = wfixed)

## fitting
iter1 <- 1
iter2 <- 1
bsd <- 0.01

tictoc::tic("fit 1")
fit1 <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x,
  betasd = bsd,
  epsilon = 1e-3,
  max_iterations = iter1,
  y = y,
  pm = param_map,
  winit,
  invisible = 0
)
tt1 <- tictoc::toc()

tictoc::tic("hessian 1")
h1 <- calc_kf_hess(w=fit1$par, x=x, y=y, betasd=bsd, pm=param_map)
tictoc::toc()

#rbind(winit, fit$par, sqrt(diag(solve(h))))

tictoc::tic("fit 2")
fit2 <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x,
  betasd = bsd,
  epsilon = 1e-3,
  max_iterations = iter2,
  y = y,
  pm = param_map,
  fit1$par,
  invisible = 0
)
tt2 <- tictoc::toc()

tictoc::tic("hessian 2")
h2 <- calc_kf_hess(w=fit2$par, x=x, y=y, betasd=bsd, pm = param_map)
tictoc::toc()
#rbind(winit, fit$par, fit2$par, sqrt(diag(solve(h2))))

## Save outputs

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
save(x, y, wfixed, fit1, fit2, h1, h2, file = the_file)

## Save metrics

dets1 <- kf_nll_details(w=fit1$par, x=x, y=y, betasd = bsd, pm = param_map, fet = NULL)
mae1 <- rowMeans(abs(dets1$ytilde_k), na.rm = TRUE)
naive_error <- colMeans(abs(apply(y, 2, diff)), na.rm = TRUE)
mase1 <- mae1 / naive_error

dets2 <- kf_nll_details(w=fit2$par, x=x, y=y, betasd = bsd, pm = param_map, fet = NULL)
mae2 <- rowMeans(abs(dets2$ytilde_k), na.rm = TRUE)
mase2 <- mae2 / naive_error

mets <- list(
  fit1nll = fit1$value,
  fit1mae = mae1,
  fit1mase = mase1, 
  fit1iter = iter1,
  fit1walltime = unname(tt1$toc - tt1$tic),
  fit1hessposdef = all(eigen(h1)$values > 0),
  fit1convergence = fit1$convergence,
  fit2nll = fit2$value,
  fit2mae = mae2,
  fit2mase = mase2, 
  fit2iter = iter2,
  fit2hessposdef = all(eigen(h2)$values > 0),
  fit2walltime = unname(tt2$toc - tt2$tic),
  fit2convergence = fit2$convergence,
  delta = (fit2$value - fit1$value) / iter2
)

met_path <- file.path(fit_dir, "fit-metrics.json")
jsonlite::write_json(mets, path = met_path, auto_unbox = TRUE)

q('no')
## scraps
dets <- kf_nll_details(w=fit1$par, x=x, y=y, betasd = .01,
                       pm = param_map, fet = NULL)

par(mfrow = c(3, 1))

qqnorm(dets$ytilde_k[1, ] / sqrt(dets$S[1, 1, ]), sub = "Cases")
abline(0, 1)
qqnorm(dets$ytilde_k[2, ] / sqrt(dets$S[2, 2, ]), sub = "Hospitalizations")
abline(0, 1)
qqnorm(dets$ytilde_k[3, ] / sqrt(dets$S[3, 3, ]), sub = "Deaths")
abline(0, 1)

rho_t <- x2$rhot
plot(x2$time, y$cases, xlab = "Time", ylab = "Cases")
pred_cases <- dets$xhat_kkmo["C", ] * rho_t + dets$xhat_kkmo["Hnew", ]
est_cases <- dets$xhat_kkmo["C", ] + dets$xhat_kkmo["Hnew", ]
se_cases <- sqrt(dets$S[1, 1, ])
lines(x$time, se_cases * 2 + pred_cases, col = "grey")
lines(x$time, pred_cases)
lines(x$time, est_cases, lty = 2)
lines(x$time,-se_cases * 2 + pred_cases, col = "grey")

plot(x$time, y$hospitalizations, xlab = "Time",
     ylab = "Hospitalizations", ylim = c(0, 2200))
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

make_rt_plot <- function(ft, x) {
  par(mfrow = c(1, 1))
  plot(
    x$time,
    exp(
      ft$par[-c(1:10)][x$bvecmap] + ft$par[7] * x$dosesiqr + ft$par[8] * x$prophomeiqr
    ) / wfixed["gamma"],
    type = 'l',
    xlab = "Time",
    ylab = expression(R[t])
  )
  legend(
    "top",
    col = c(1, "orange", "blue", "grey"),
    lty = 1,
    legend = c(
      "MLE estimate",
      "MLE - (effect of mobility)",
      "MLE - (effect of vaccine)",
      "MLE random walk intercept"
    )
  )
  lines(x$time,
        exp(ft$par[-c(1:10)][x$bvecmap]+ ft$par[7] * x$dosesiqr + 0 * x$prophomeiqr) / wfixed["gamma"],
        col = "orange")
  lines(x$time,
        exp(ft$par[-c(1:10)][x$bvecmap] + 0 * x$dosesiqr + ft$par[8] * x$prophomeiqr) / wfixed["gamma"],
        col = "blue")
  lines(x$time,
        exp(ft$par[-c(1:10)][x$bvecmap] + 0 * x$dosesiqr + 0 * x$prophomeiqr) / wfixed["gamma"],
        col = "grey")
}

make_rt_plot(fit2, x)

plot.new()
est_tab <- rbind(initial = winit, MLE = fit2$par, sd = sqrt(diag(solve(h2)))) %>% signif(3)
gridExtra::grid.table(est_tab[, 1:10])

plot.new()
gridExtra::grid.table(est_tab[, -c(1:10)])



