#!/usr/bin/env R

## Environment prep

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

JuliaCall::julia_setup("/opt/julia-1.5.3/bin")
JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
JuliaCall::julia_eval("using DataFrames")

## Data prep

forecast_date_start <- Sys.getenv("fdtstart", unset = "2020-06-29")
forecast_date <- Sys.getenv("fdt", unset = "2020-06-29")
forecast_loc <- Sys.getenv("loc", unset = "06")

hopdir <- file.path("hopkins", forecast_date)
tdat <- load_hopkins(hopdir, weekly = FALSE)
ltdat <- tdat %>% filter(location == forecast_loc) %>%
  filter(target_type == "day ahead inc case" |
           target_type == "day ahead inc death")
ltdat2 <-
  ltdat %>% mutate(time = lubridate::decimal_date(target_end_date)) %>%
  pivot_wider(names_from = target_type, values_from = value)

jhu_data <- ltdat2 %>% ungroup() %>%
  mutate(wday = lubridate::wday(target_end_date)) %>%
  rename(cases = `day ahead inc case`, deaths = `day ahead inc death`) %>%
  select(target_end_date, time, wday, cases, deaths)

mob_path <- file.path("covidcast-safegraph-home-prop",
                      forecast_date,
                      "epidata.csv")
abbr <- covidcast::fips_to_abbr(paste0(forecast_loc, "000")) %>% tolower()
mob_ts <- read_csv(
  mob_path,
  col_types = cols_only(
    geo_value = col_character(),
    time_value = col_date(),
    value = col_double()
  )
) %>%
  filter(geo_value == abbr) %>%
  select(time_value, value) %>%
  rename(target_end_date = time_value,
         prophome = value)

if (forecast_date > "2020-11-15") {
  healthd <-
    file.path("healthdata", forecast_date, forecast_loc, "epidata.csv")
  cov_thresh <- .9
  tdat2 <- read_csv(
    healthd,
    col_types = cols_only(
      date = col_date("%Y%m%d"),
      previous_day_admission_adult_covid_confirmed = col_integer(),
      previous_day_admission_pediatric_covid_confirmed = col_integer(),
      previous_day_admission_adult_covid_confirmed_coverage = col_integer()
    )
  )
  most_recent_coverage <- tdat2 %>% arrange(date) %>%
    pull(previous_day_admission_adult_covid_confirmed_coverage) %>%
    tail(n = 1)
  tdat3 <- tdat2 %>%
    filter(
      previous_day_admission_adult_covid_confirmed_coverage >
        most_recent_coverage * cov_thresh
    ) %>%
    filter(
      previous_day_admission_adult_covid_confirmed_coverage <
        most_recent_coverage / cov_thresh
    ) %>%
    mutate(
      hospitalizations = previous_day_admission_adult_covid_confirmed +
        previous_day_admission_pediatric_covid_confirmed,
      target_end_date = date - lubridate::ddays(1)
    ) %>%
    select(target_end_date, hospitalizations)
  obs_data <-
    left_join(jhu_data, mob_ts, by = "target_end_date") %>%
    mutate(prophome = zoo::na.fill(prophome, "extend")) %>%
    left_join(tdat3, by = "target_end_date")
} else {
  obs_data <- left_join(jhu_data, mob_ts, by = "target_end_date") %>%
    mutate(prophome = zoo::na.fill(prophome, "extend"))
}


if (forecast_date > "2020-12-31"){
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
    mutate(
      prophome = lead(prophome, 3),
      prophome = zoo::na.fill(prophome, "extend")
    )
  
  obs_data$doses[lubridate::year(obs_data$target_end_date) == 2020] <-
    0
} else {
  stop("Fitting not implemented for dates without vaccine data")
}
}

#wind <- obs_data %>% slice(match(1, obs_data$cases > 0):n())



wind <- obs_data %>% slice(41:n())
wsize <- nrow(wind)
N <-
  covidHubUtils::hub_locations %>% filter(fips == forecast_loc) %>%
  pull(population)

case_to_death_lag <- 20

wfixed <- c(
  N = N,
  gamma = 365.25 / 4,
  gamma_h = 365.25 / (case_to_death_lag / 2),
  gamma_d = 365.25 / (case_to_death_lag / 2),
  chp = 0.03,
  eta = 365.25 / 4
)

if (forecast_date > "2020-11-15") {
  y <- wind %>% select(cases, hospitalizations, deaths)
} else {
  y <- wind %>% select(cases, deaths)
}

x0 <- wind %>% select(target_end_date, time, wday, prophome)
x0$prophomeiqr <-
  (x0$prophome - mean(x0$prophome)) / diff(quantile(x0$prophome, c(.25, .75)))
x0$bvecmap <- rep(seq_len(ceiling(wsize / 7)), each = 7) %>% head(wsize) %>% as.integer()

## removal of untrusted data points

if(forecast_loc == 33){
  isuntrusted <- x0$target_end_date >= "2020-08-16" & x0$target_end_date <= "2020-09-02"
  y$hospitalizations[isuntrusted] <- NA
  #rationale: inconsistent with deaths data from hopkins and hospitalization data from https://carlsonschool.umn.edu/mili-misrc-covid19-tracking-project/download-data
}

set.seed(1)

# estimate the probability that a case is reported assuming that all deaths are reported
cvd <-
  data.frame(time = x0$time,
             r = y$cases / lead(y$deaths, case_to_death_lag)) %>%
  mutate(logr = log(r)) %>% filter(is.finite(logr)) %>% filter(time < 2020.47)
mars <- earth::earth(logr ~ time, data = cvd)
#plot(logr ~ time, data = cvd)
#lines(cvd$time, predict(mars))
cvd$rhat <- exp(predict(mars)) %>% as.numeric() %>% cummax() # make monotonic
right <- cvd %>% select("time", "rhat")

x <- left_join(x0, right, by = "time") %>%
  mutate(rhat2 = zoo::na.fill(rhat, "extend"),
         rhot = rhat2 / (2 * rhat2[n()])) %>%
  select(-rhat,-rhat2)

stopifnot(all(x$rhot <= 1))
stopifnot(all(x$rhot >= 0))

#pdf <- data.frame(date = wind$target_end_date, p = x$rhot)
#ggplot(pdf, aes(x = date, y = p)) + geom_point() + ylab("Pr (case is reported)")

if (forecast_date_start == forecast_date){
  stopifnot(forecast_date < "2020-11-16") # initializer not implemented for data with hospitalizations
  winit <- initialize_estimates(x = x, y = y, wfixed = wfixed)
} else {
  fit_dir_start <-
    file.path("fits",
              paste0(forecast_date_start,
                     "-fips",
                     forecast_loc))
  fp <- file.path(fit_dir_start, "fit.RData")
  attach(fp, pos = 2)
  winit <- get("fit1", pos = 2)$par
  detach(pos = 2)
  datediff <- (lubridate::ymd(forecast_date) - lubridate::ymd(forecast_date_start)) / lubridate::ddays(1)
  sizediff <- x$bvecmap[wsize] - x$bvecmap[wsize - datediff]
  if (sizediff > 0){
    winit <- c(winit, rep(winit[length(winit)], sizediff)) # extend the bvec by repeating last value
  }
}

if (forecast_date >= "2020-11-16" && forecast_date_start < "2020-11-16"){
  tauh_init <- var(na.omit(diff(y$hospitalizations)))
  is_NA_hosps <- is.na(y$hospitalizations)
  chp_init <- sum(y$hospitalizations[!is_NA_hosps]) / sum(y$cases[!is_NA_hosps])
  winit0 <- winit
  winit <- c(winit0[1:3], log(tauh_init), winit0[4], qlogis(chp_init), winit0[5:length(winit0)])
}

if (forecast_date == "2020-06-29" && forecast_loc == "06") {
  winit <-
    c(
      8.84525442487002,
      -1.9410973262972,
      4.22, #12.4822240973189,
      5.62041865346974,
      -0.215231653815656,
      -0.985527447518198,
      2.8253984751295,
      3.59,
      4.68123234065946,
      4.68070637757278,
      4.67844012908105,
      4.67339401927969,
      4.66485980762512,
      4.6533866117208,
      4.63997699124128,
      4.62550077865006,
      4.6116270964086,
      4.59946956774516,
      4.5902455935586,
      4.58440251810939,
      4.58207964451701,
      4.58370322953008,
      4.58817592970655,
      4.59186894945956,
      4.59167559076008
    )
}

## fitting
iter1 <- 1000
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
  winit,
  wfixed = wfixed,
  invisible = 0
)
tt1 <- tictoc::toc()

tictoc::tic("hessian 1")
h1 <- calc_kf_hess(
  w = fit1$par,
  x = x,
  y = y,
  betasd = bsd,
  pm = param_map
)
tictoc::toc()


## Save outputs

fit_dir <-
  file.path("fits",
            paste0(forecast_date,
                   "-fips",
                   forecast_loc))

if (!dir.exists(fit_dir))
  dir.create(fit_dir, recursive = TRUE)

the_file <- file.path(fit_dir, "fit.RData")
save(x, y, winit, wfixed, fit1, h1, bsd, file = the_file)

## Save metrics

dets1 <-
  kf_nll_details(
    w = fit1$par,
    x = x,
    y = y,
    betasd = bsd,
    fixed = wfixed,
    fet = NULL
  )
mae1 <- rowMeans(abs(dets1$ytilde_k), na.rm = TRUE)
naive_error <- colMeans(abs(apply(y, 2, diff)), na.rm = TRUE)
naive_error_weekly <- colMeans(abs(y[-c(1:7),] - y[-((nrow(y) - 6):nrow(y)),]), na.rm = TRUE)

if (ncol(y) == 2) {
  mase1 <- mae1 / c(naive_error["cases"], NA, naive_error["deaths"])
  mase1w <-
    mae1 / c(naive_error_weekly["cases"], hosps = NA, naive_error_weekly["deaths"])
} else {
  mase1 <- mae1 / naive_error
  mase1w <- mae1 / naive_error_weekly
}

g1 <-
  calc_kf_grad(
    w = fit1$par,
    x = x,
    y = y,
    betasd = bsd,
    wfixed = wfixed
  )

mets <- list(
  fit1nll = fit1$value,
  fit1xnorm = sqrt(sum(fit1$par ^ 2)),
  fit1gnorm = sqrt(sum(g1 ^ 2)),
  fit1mae_cases = mae1[1],
  fit1mae_hosps = mae1[2],
  fit1mae_death = mae1[3],
  fit1mase_cases = mase1[1],
  fit1mase_hosps = mase1[2],
  fit1mase_death = mase1[3],
  fit1mase_cases_w = mase1w[1],
  fit1mase_hosps_w = mase1w[2],
  fit1mase_death_w = mase1w[3],  
  fit1iter = iter1,
  fit1walltime = unname(tt1$toc - tt1$tic),
  fit1hessposdef = all(eigen(h1)$values > 0),
  fit1convergence = fit1$convergence
)

met_path <- file.path(fit_dir, "fit-metrics.json")
jsonlite::write_json(mets, path = met_path, auto_unbox = TRUE)