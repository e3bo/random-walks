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
  
  obs_data0 <-
    left_join(jhu_data, mob_ts, by = "target_end_date") %>%
    mutate(prophome = zoo::na.fill(prophome, "extend")) %>%
    left_join(tdat3, by = "target_end_date")
  
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
    obs_data <-
      left_join(obs_data0, vacc_ts, by = "target_end_date")
    obs_data$doses[lubridate::year(obs_data$target_end_date) == 2020] <-
      0
  } else {
    obs_data <- obs_data0
  }
} else {
  obs_data <- left_join(jhu_data, mob_ts, by = "target_end_date") %>%
    mutate(prophome = zoo::na.fill(prophome, "extend"))
}

wind <- obs_data %>% filter(target_end_date >= "2020-03-02")
wsize <- nrow(wind)
N <-
  covidHubUtils::hub_locations %>% filter(fips == forecast_loc) %>%
  pull(population)

case_to_death_lag <- 20

wfixed <- c(
  N = N,
  γ = 365.25 / 4,
  γ_h = 365.25 / (case_to_death_lag / 2),
  γ_d = 365.25 / (case_to_death_lag / 2),
  γ_z = 365.25 / 1, 
  p_h = 0.03,
  τ_h = 10, # wfixed p_h and τ_h values only are used if there are no hospitalization data, otherwise they are estimated
  η = 365.25 / 4
)

if (forecast_date > "2020-11-15") {
  z <- wind %>% select(cases, hospitalizations, deaths)
} else {
  z <- wind %>% select(cases, deaths)
}

if("doses" %in% names(wind)) {
  x0 <- wind %>% select(target_end_date, time, wday, prophome, doses) %>%
    mutate(doses = zoo::na.approx(doses, rule = 2))
  sel <- x0$doses > 0 
  x0$dosesiqr <-
    (x0$doses - mean(x0$doses)) / diff(quantile(x0$doses[sel], c(.25, .75)))
} else {
  x0 <- wind %>% select(target_end_date, time, wday, prophome)
}

x0$prophomeiqr <-
  (x0$prophome - mean(x0$prophome)) / diff(quantile(x0$prophome, c(.25, .75)))
x0$β_0map <- rep(seq_len(ceiling(wsize / 7)), each = 7) %>% head(wsize) %>% as.integer()
x0$τ_cmap <- rep(seq_len(ceiling(wsize / 28)), each = 28) %>% head(wsize) %>% as.integer()

if (ncol(z) == 3){
 # hwsize <- wsize - sum(is.na(z$hospitalizations))
#  x0$p_hmap <- c(rep(1L, wsize - hwsize), rep(seq_len(ceiling(hwsize / 28)), each = 28)) %>% head(wsize) %>% as.integer()
  x0$p_hmap <- rep(seq_len(ceiling(wsize / 28)), each = 28) %>% head(wsize) %>% as.integer() 
} else {
  x0$p_hmap <- 1L
}

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
             r = z$cases / lead(z$deaths, case_to_death_lag)) %>%
  mutate(logr = log(r)) %>% filter(is.finite(logr)) %>% filter(time < 2020.47)
mars <- earth::earth(logr ~ time, data = cvd)
#plot(logr ~ time, data = cvd)
#lines(cvd$time, predict(mars))
cvd$rhat <- exp(predict(mars)) %>% as.numeric() %>% cummax() # make monotonic
right <- cvd %>% select("time", "rhat")

x <- left_join(x0, right, by = "time") %>%
  mutate(rhat2 = zoo::na.fill(rhat, "extend"),
         ρ = rhat2 / (2 * rhat2[n()])) %>%
  select(-rhat,-rhat2)

stopifnot(all(x$ρ <= 1))
stopifnot(all(x$ρ >= 0))

#pdf <- data.frame(date = wind$target_end_date, p = x$rhot)
#ggplot(pdf, aes(x = date, y = p)) + geom_point() + ylab("Pr (case is reported)")

if (forecast_date_start == forecast_date) {
  stopifnot(forecast_date < "2020-11-16") # initializer not implemented for data with hospitalizations
  winit <- initialize_estimates(x = x, y = z, wfixed = wfixed)
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
  datediff <-
    (lubridate::ymd(forecast_date) - lubridate::ymd(forecast_date_start)) / lubridate::ddays(1)
  nβ_0 <- x$β_0map[wsize]
  
  sizediff <- nβ_0 - x$β_0map[wsize - datediff]
  if (sizediff > 0) {
    winit <-
      c(winit, rep(winit[length(winit)], sizediff)) # extend β_0 by repeating last value
  }
  
  np_h <- x$p_hmap[wsize]
  if (forecast_date >= "2020-11-16" &&
      forecast_date_start < "2020-11-16") {
    # initialize hospitalization parameters
    τ_hinit <- var(na.omit(diff(z$hospitalizations)))
    is_NA_hosps <- is.na(z$hospitalizations)
    p_hinit <-
      sum(z$hospitalizations[!is_NA_hosps]) / sum(z$cases[!is_NA_hosps])
    winit0 <- winit
    nw0 <- length(winit0)
    winit <- c(winit0[1:2],
               log(τ_hinit),
               winit0[3:(nw0 - nβ_0)],
               rep(qlogis(p_hinit), times = np_h),
               winit0[(nw0 - nβ_0 + 1):nw0])
  } else {
    sizediff_p_h <- np_h - x$p_hmap[wsize - datediff]
    if (sizediff_p_h > 0) {
      #extend p_h
      nw <- length(winit)
      last_val <- winit[nw - nβ_0]
      winit <-
        c(winit[1:(nw - nβ_0)], rep(last_val, times = sizediff_p_h), 
          winit[(nw - nβ_0 + 1):nw])
    }
  }
  
  if (forecast_date > "2021-02-08" &&
      forecast_date_start <= "2021-02-08"){
    # add parameter for doseeffect
    winit <- c(winit[1:5], -1, winit[6:length(winit)])
  }
  
  nτ_c <- x$τ_cmap[wsize]
  sizediff_τ_c <- nτ_c - x$τ_cmap[wsize - datediff]
  if (sizediff_τ_c > 0) {
    #extend τ_c by repeating last value
    nw <- length(winit)
    last_val <- winit[nw - nβ_0 - np_h]
    winit <-
      c(winit[1:(nw - nβ_0 - np_h)], rep(last_val, times = sizediff_τ_c), winit[(nw - nβ_0 - np_h + 1):nw])
  }
}

## fitting
iter1 <- 20
β_0sd <- 0.1
τ_csd <- 0.05
p_hsd <- 0.5

#wnew <- c(winit[1:4], winit[6:(11 + 11 -1)], rep(winit[5], 11), winit[(11 + 11):(length(winit))])

tictoc::tic("fit 1")
fit1 <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  cov = x,
  p_hsd = p_hsd,
  β_0sd = β_0sd,
  τ_csd = τ_csd, 
  epsilon = 1e-3,
  max_iterations = iter1,
  z = z,
  winit,
  wfixed = wfixed,
  invisible = 0
)
tt1 <- tictoc::toc()

if(FALSE){
tictoc::tic("hessian 1")
h1 <- calc_kf_hess(
  w = fit1$par,
  cov = x,
  z = z,
  p_hsd = p_hsd,
  β_0sd = β_0sd,
  τ_csd = τ_csd,
  wfixed = wfixed
)
tictoc::toc()
} else{
  h1 <- diag(nrow = length(winit))
}
## Save outputs

fit_dir <-
  file.path("fits",
            paste0(forecast_date,
                   "-fips",
                   forecast_loc))

if (!dir.exists(fit_dir))
  dir.create(fit_dir, recursive = TRUE)

the_file <- file.path(fit_dir, "fit.RData")
save(x, z, winit, wfixed, fit1, h1, p_hsd, β_0sd, τ_csd, file = the_file)

## Save metrics

dets1 <-
  calc_kf_nll_r(
    w = fit1$par,
    cov = x,
    z = z,
    p_hsd = p_hsd,
    β_0sd = β_0sd,
    τ_csd = τ_csd,
    wfixed = wfixed,
    just_nll = FALSE
  )
mae1 <- rowMeans(abs(dets1$ytilde_k), na.rm = TRUE)
naive_error <- colMeans(abs(apply(z, 2, diff)), na.rm = TRUE)
naive_error_weekly <- colMeans(abs(z[-c(1:7),] - z[-((nrow(z) - 6):nrow(z)),]), na.rm = TRUE)

if (ncol(z) == 2) {
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
    cov = x,
    z = z,
    p_hsd = p_hsd,
    β_0sd = β_0sd,
    τ_csd = τ_csd,
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
