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
wind <- obs_data %>% slice(41:n())
wsize <- nrow(wind)
N <- covidHubUtils::hub_locations %>% filter(fips == forecast_loc) %>% 
  pull(population)

wfixed <- c(
  N = N,
  gamma = 365.25 / 4,
  gamma_h = 365.25 / 5,
  gamma_d = 365.25 / 5,
  eta = 365.25 / 4
)

y <- wind %>% select(cases, hospitalizations, deaths)
x <- wind %>% select(time, doses, prophome)
x$dosesiqr <- x$doses /  diff(quantile(x$doses[x$doses > 0], c(.25, .75)))
x$prophomeiqr <- (x$prophome - mean(x$prophome)) / diff(quantile(x$prophome, c(.25, .75)))
set.seed(1)

# estimate the probability that a case is reported assuming that all deaths are reported
cvd <- data.frame(time = x$time, r = y$cases / lead(y$deaths, 10)) %>% 
  mutate(logr = log(r)) %>% filter(is.finite(logr)) %>% filter(time < 2020.47)
mars <- earth::earth(logr~time, data = cvd)
#plot(logr ~ time, data = cvd)
#lines(cvd$time, predict(mars))
cvd$rhat <- exp(predict(mars)) %>% as.numeric()
right <- cvd %>% select("time", "rhat")

x2 <- left_join(x, right, by = "time") %>% 
  mutate(rhat2 = zoo::na.fill(rhat, "extend"),
         rhot = rhat2 / (2 *rhat2[n()]))


#pdf <- data.frame(date = wind$target_end_date, p = x2$rhat3)
#ggplot(pdf, aes(x = date, y = p)) + geom_point() + ylab("Pr (case is reported)")

winit <- initialize_estimates(x = x2, y = y, wfixed = wfixed)

tictoc::tic("optimization")

system.time(fit <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x2,
  betasd = 0.01,
  epsilon = 1e-1,
  max_iterations = 10,
  a = .9,
  y = y,
  pm = param_map,
  winit,
  invisible = 0
))

system.time(h <- calc_kf_hess(w=fit$par, x=x2, y=y, betasd=0.01, a =.1, pm = param_map))
rbind(winit, fit$par, sqrt(diag(solve(h))))

fit2 <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x2,
  betasd = 0.01,
  epsilon = 1e-3,
  max_iterations = 90,
  a = 0.9,
  y = y,
  pm = param_map,
  fit$par,
  invisible = 0
)

system.time(h2 <- calc_kf_hess(w=fit2$par, x=x2, y=y, betasd=0.01, a =.1, pm = param_map))
rbind(winit, fit$par, fit2$par, sqrt(diag(solve(h2))))

fit3 <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x2,
  betasd = 0.01,
  epsilon = 1e-3,
  max_iterations = 600,
  a = 0.9,
  y = y,
  pm = param_map,
  fit2$par,
  invisible = 0
)

system.time(h3 <- calc_kf_hess(w=fit3$par, x=x2, y=y, betasd=0.01, a =.1, pm = param_map))
rbind(winit, fit$par, fit2$par, fit3$par, sqrt(diag(solve(h3))))


tictoc::toc()


## 


dets <- kf_nll_details(w=fit2$par, x=x2, y=y, betasd = .01, a = 0.9, pm = param_map, fet = NULL)

par(mfrow = c(3, 1))
qqnorm(dets$ytilde_k[1, ] / sqrt(dets$S[1, 1, ]), sub = "Cases")
abline(0, 1)
qqnorm(dets$ytilde_k[2, ] / sqrt(dets$S[2, 2, ]), sub = "Hospitalizations")
abline(0, 1)
qqnorm(dets$ytilde_k[3, ] / sqrt(dets$S[3, 3, ]), sub = "Deaths")
abline(0, 1)

rho_t <- x2$rhot
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

make_rt_plot <- function(ft, x) {
  par(mfrow = c(1, 1))
  plot(
    x$time,
    exp(
      rep(ft$par[-c(1:9)], each = 28) + ft$par[8] * x$dosesiqr + ft$par[9] * x$prophomeiqr
    ) / wfixed["gamma"],
    type = 'l',
    xlab = "Time",
    ylab = expression(R[t])
  )
  legend(
    "topright",
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
        exp(rep(ft$par[-c(1:9)], each = 28) + ft$par[8] * x$dosesiqr + 0 * x$prophomeiqr) / wfixed["gamma"],
        col = "orange")
  lines(x$time,
        exp(rep(ft$par[-c(1:9)], each = 28) + 0 * x$dosesiqr + ft$par[9] * x$prophomeiqr) / wfixed["gamma"],
        col = "blue")
  lines(x$time,
        exp(rep(ft$par[-c(1:9)], each = 28) + 0 * x$dosesiqr + 0 * x$prophomeiqr) / wfixed["gamma"],
        col = "grey")
}

make_rt_plot(fit3, x)

rho_t <- detect_frac(365.25 * (x$time - 2020.2))
plot(x$time, y$cases, xlab = "Time", ylab = "Cases", xlim = c(2020.6, 2020.8), ylim = c(0, 2200))
pred_cases <- dets$xhat_kkmo["C", ] * rho_t + dets$xhat_kkmo["Hnew", ]
est_cases <- dets$xhat_kkmo["C", ] + dets$xhat_kkmo["Hnew", ]
se_cases <- sqrt(dets$S[1, 1, ])
lines(x$time, se_cases * 2 + pred_cases, col = "grey")
lines(x$time, pred_cases)
lines(x$time, est_cases, lty = 2)
lines(x$time,-se_cases * 2 + pred_cases, col = "grey")

par(mfrow = c(2, 1))
plot(x$time, fit2$par[9] * x$prophomeiqr, xlab = "Time", 
     ylab = expression(paste("Effect of mobility on ", log(beta[t]))))
plot(x$time, fit2$par[8] * x$betanoise, 
     xlab = "Time",
     ylab = expression(paste("Effect of heterogeneity on ", log(beta[t]))))

plot.new()
est_tab <- rbind(initial = winit, MLE = fit2$par, sd = sqrt(diag(solve(h2)))) %>% signif(3)
gridExtra::grid.table(est_tab[, 1:9])

plot.new()
gridExtra::grid.table(est_tab[, -c(1:9)])
