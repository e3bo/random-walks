#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

make_rt_plot <- function(ft, cov, no_hosps, wfixed) {
  par(mfrow = c(1, 1))
  nβ_0 <- tail(cov$β_0map, 1)
  np <- length(ft$par)
  inds <- seq(np - nβ_0 + 1, np)
  intercept <- ft$par[inds][cov$β_0map]
  
  if ("dosesiqr" %in% names(cov)){
    X <- cbind(cov$prophomeiqr, cov$dosesiqr)
    effects <- -exp(ft$par[c(5, 6)])
  } else {
    X <- cbind(cov$prophomeiqr)
    if (no_hosps){
      effects <- -exp(ft$par[4])
    } else {
      effects <- -exp(ft$par[5])
    }
  }

  num_all <- exp(intercept + X %*% effects)
  plot(
    cov$time,
    num_all / wfixed["γ"],
    type = 'l',
    xlab = "Time",
    ylab = expression(R[t])
  )
  if (length(effects) == 2) {
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
    effects_no_dose <- effects_no_mob <- effects
    effects_no_mob[1] <- 0
    num_nomob <- exp(intercept + X %*% effects_no_mob)
    lines(cov$time,
          num_nomob / wfixed["γ"],
          col = "orange")
    effects_no_dose[2] <- 0
    num_nodose <- exp(intercept + X %*% effects_no_dose)
    lines(cov$time,
          num_nodose / wfixed["γ"],
          col = "blue")
  } else {
    legend(
      "top",
      col = c(1, "orange", "grey"),
      lty = 1,
      legend = c(
        "MLE estimate",
        "MLE - (effect of mobility)",
        "MLE random walk intercept"
      )
    )
    effects_no_mob <- effects
    effects_no_mob[1] <- 0
    num_nomob <- exp(intercept + X %*% effects_no_mob)
    lines(cov$time,
          num_nomob / wfixed["γ"],
          col = "orange")
  }
  lines(cov$time,
        exp(intercept) / wfixed["γ"],
        col = "grey")
}

create_forecast_df <- function(means,
                               vars,
                               target_type = "cases",
                               fdt = forecast_date,
                               location) {
  # takes in a vector of means and variance for the normal predicted
  # distribution of the next n weeks, and outputs a representation of
  # this forecast in the standard covidhub format
  if (target_type == "cases") {
    probs <- c(0.025, 0.100, 0.250, 0.500, 0.750, 0.900, 0.975)
    maxn <- 8
    targ_string <- " wk ahead inc case"
    stride <- 7
  } else if (target_type == "hospitalizations") {
    probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
    maxn <- 129
    targ_string <- " day ahead inc hosp"
    stride <- 1
  } else if (target_type == "inc deaths") {
    probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
    maxn <- 8
    targ_string <- " wk ahead inc death"
    stride <- 7
  } else {
    stop("not implemented")
  }
  n <- length(means)
  stopifnot(length(vars) == n)
  stopifnot(n <= maxn)
  step_ahead <- seq_len(n)
  targets <- paste0(step_ahead, targ_string)
  dte <- lubridate::ymd(fdt)
  if (target_type %in% c("cases", "inc deaths")) {
    fdt_wd <- lubridate::wday(dte) # 1 = Sunday, 7 = Sat.
    fdt_sun <- lubridate::ymd(dte) - (fdt_wd - 1)
    take_back_step <- fdt_wd <= 2
    if (take_back_step) {
      week0_sun <- fdt_sun - 7
    } else {
      week0_sun <- fdt_sun
    }
    time_zero <- week0_sun + 6
  } else {
    time_zero <- dte
  }
  t1 <- expand_grid(quantile = probs, h = step_ahead) %>%
    mutate(target = paste0(h, targ_string),
           target_end_date = time_zero + h * stride) %>%
    add_column(
      location = as.character(location),
      forecast_date = fdt,
      type = "quantile"
    ) %>%
    mutate(q1 = map2_dbl(quantile,
                         h,
                         ~ qnorm(
                           p = .x,
                           mean = means[.y],
                           sd = sqrt(vars[.y])
                         )),
           value = ifelse(q1 < 0, 0, q1)) %>%
    mutate(value = format(value, digits = 4),
           quantile = as.character(quantile))
  t2 <- t1 %>% filter(quantile == "0.5") %>%
    mutate(type = "point",
           quantile = NA_character_)
  bind_rows(t1, t2) %>% dplyr::select(forecast_date,
                                      target,
                                      target_end_date,
                                      location,
                                      type,
                                      quantile,
                                      value)
}

write_forecasts <-
  function(fit, cov, z, winit, wfixed, hess, fets, p_hsd, β_0sd, τ_csd, fdt, forecast_loc) {
    
    dets <-
      calc_kf_nll_r(
        w = fit$par,
        cov = cov,
        z = z,
        p_hsd = p_hsd,
        β_0sd = β_0sd,
        τ_csd = τ_csd,
        wfixed = wfixed,
        just_nll = FALSE,
        fets = fets
      )
    
    make_fit_plots(
      dets,
      cov = cov,
      winit = winit, 
      h = hess,
      p_hsd = p_hsd,
      β_0sd = β_0sd,
      τ_csd = τ_csd,      
      fdt = fdt,
      forecast_loc = forecast_loc,
      fit = fit,
      wfixed = wfixed
    )
    case_inds <- which(fets$target_wday == 7)
    if (lubridate::wday(fdt) > 2) {
      case_inds <-
        case_inds[-1] # drop first epiweek, which has been observed too much 
    }
    case_fcst <-
      create_forecast_df(
        means = dets$sim_means[1, case_inds],
        vars = dets$sim_cov[1, 1, case_inds],
        location = forecast_loc,
        fdt = fdt
      )
    stopifnot(setequal(
      fets$target_end_dates[case_inds],
      case_fcst$target_end_date %>% unique()
    ))
    hosp_inds <- fets$target_end_dates %in%
      (lubridate::ymd(forecast_date) + 1:28)
    hosp_fcst <-
      create_forecast_df(
        means = dets$sim_means[2, hosp_inds],
        vars = dets$sim_cov[2, 2, hosp_inds],
        target_type = "hospitalizations",
        location = forecast_loc,
        fdt = fdt
      )
    stopifnot(setequal(
      fets$target_end_dates[hosp_inds],
      hosp_fcst$target_end_date %>% unique()
    ))
    
    death_inds <- case_inds
    death_fcst <-
      create_forecast_df(
        means = dets$sim_means[3, death_inds],
        vars = dets$sim_cov[3, 3, death_inds],
        target_type = "inc deaths",
        location = forecast_loc,
        fdt = fdt
      )
    fcst <- bind_rows(case_fcst, hosp_fcst, death_fcst)
    
    lambda <- 1 / β_0sd
    fcst_dir <-
      file.path("forecasts",
                paste0(forecast_date,
                       "-fips",
                       forecast_loc))
    
    fcst_name <- paste0("lambda",
                        sprintf("%06.2f", lambda),
                        "-status-quo",
                        "-CEID-InfectionKalman.csv")
    fcst_path <- file.path(fcst_dir, fcst_name)
    if (!dir.exists(fcst_dir))
      dir.create(fcst_dir, recursive = TRUE)
    write_csv(x = fcst, path = fcst_path)
  }

make_fit_plots <- function(dets, cov, winit, h, p_hsd, β_0sd, τ_csd, fdt, forecast_loc, fit, wfixed) {
  plot_dir <-
    file.path("plots",
              paste0(forecast_date,
                     "-fips",
                     forecast_loc),
              paste0("lambda",
                     sprintf("%06.2f", 1 / β_0sd)))
  
  if (!dir.exists(plot_dir))
    dir.create(plot_dir, recursive = TRUE)
  
  no_hosps <- all(is.na(dets$ytilde_k[2,]))
  
  plot_path1 <- file.path(plot_dir, "qqplots.png")
  
  png(
    plot_path1,
    width = 7.5,
    height = 10,
    units = "in",
    res = 90
  )
  
  par(mfrow = c(3, 1))
  qqnorm(dets$ytilde_k[1,] / sqrt(dets$S[1, 1,]), sub = "Cases")
  abline(0, 1)
  if(no_hosps){
    plot.new()
  } else {
    qqnorm(dets$ytilde_k[2,] / sqrt(dets$S[2, 2,]), sub = "Hospitalizations")
    abline(0, 1)
  }
  qqnorm(dets$ytilde_k[3,] / sqrt(dets$S[3, 3,]), sub = "Deaths")
  abline(0, 1)
  dev.off()
  
  plot_path2 <- file.path(plot_dir, "fitted-time-series.png")
  png(
    plot_path2,
    width = 7.5,
    height = 10,
    units = "in",
    res = 90
  )
  par(mfrow = c(3, 1))
  
  plot(x$time, z$cases, xlab = "Time", ylab = "Cases")
  pred_cases <- dets$xhat_kkmo["Z_r",]
  se_cases <- sqrt(dets$S[1, 1,])
  lines(x$time, se_cases * 2 + pred_cases, col = "grey")
  lines(x$time, pred_cases)
  lines(x$time, -se_cases * 2 + pred_cases, col = "grey")
  
  if (no_hosps) {
    plot.new()
  } else {
    plot(x$time,
         z$hospitalizations,
         xlab = "Time",
         ylab = "Hospitalizations")
    pred_hosps <- dets$xhat_kkmo["A", ]
    se_hosps <- sqrt(dets$S[2, 2, ])
    lines(x$time, se_hosps * 2 + pred_hosps, col = "grey")
    lines(x$time, pred_hosps)
    lines(x$time,-se_hosps * 2 + pred_hosps, col = "grey")
  }
  
  plot(x$time, z$deaths, xlab = "Time",
       ylab = "Deaths")
  pred_deaths <- dets$xhat_kkmo["D_r",]
  se_deaths <- sqrt(dets$S[3, 3,])
  lines(x$time, se_deaths * 2 + pred_deaths, col = "grey")
  lines(x$time, pred_deaths)
  lines(x$time, -se_deaths * 2 + pred_deaths, col = "grey")
  dev.off()
  
  plot_path3 <- file.path(plot_dir, "case-forecasts.png")
  case_inds <- which(fets$target_wday == 7)
  case_fcst <-
    create_forecast_df(
      means = dets$sim_means[1, case_inds],
      vars = dets$sim_cov[1, 1, case_inds],
      location = forecast_loc,
      fdt = fdt
    )
  p3 <- case_fcst %>%
    ggplot(aes(
      x = target_end_date,
      y = as.numeric(value),
      color = quantile
    )) +
    geom_line() + geom_point() + labs(y = "Weekly cases")
  ggsave(plot_path3, p3)
  
  plot_path4 <- file.path(plot_dir, "hosp-forecasts.png")
  hosp_inds <-
    fets$target_end_dates %in% (lubridate::ymd(forecast_date) + 1:28)
  stopifnot(sum(hosp_inds) == 28)
  
  hosp_fcst <-
    create_forecast_df(
      means = dets$sim_means[2, hosp_inds],
      vars = dets$sim_cov[2, 2, hosp_inds],
      target_type = "hospitalizations",
      location = forecast_loc,
      fdt = fdt
    )
  
  p4 <- hosp_fcst %>%
    ggplot(aes(
      x = target_end_date,
      y = as.numeric(value),
      color = quantile
    )) +
    geom_line() + geom_point() + labs(y = "Daily hospital admissions")
  ggsave(plot_path4, p4)
  
  plot_path5 <- file.path(plot_dir, "death-forecasts.png")
  death_fcst <-
    create_forecast_df(
      means = dets$sim_means[3, case_inds],
      vars = dets$sim_cov[3, 3, case_inds],
      target_type = "inc deaths",
      location = forecast_loc,
      fdt = fdt
    )
  
  p5 <- death_fcst %>%
    ggplot(aes(
      x = target_end_date,
      y = as.numeric(value),
      color = quantile
    )) +
    geom_line() + geom_point() + labs(y = "Weekly deaths")
  
  ggsave(filename = plot_path5, plot = p5)
  
  plot_path6 <- file.path(plot_dir, "Rt-time-series.png")
  png(
    plot_path6,
    width = 7.5,
    height = 7,
    units = "in",
    res = 90
  )
  make_rt_plot(fit, cov, no_hosps, wfixed)
  dev.off()
  
  plot_path7 <- file.path(plot_dir, "params-tab-1.png")  
  png(plot_path7,
      width = 10,
      height = 2,
      units = "in",
      res = 90)
  
  if(no_hosps){
    t1inds <- seq_len(8)
  } else{
    t1inds <- seq_len(10)
  }
  
  est_tab <-
    rbind(init = winit,
          MLE = fit$par,
          sd = sqrt(diag(solve(h)))) %>% signif(3)
  gridExtra::grid.table(est_tab[, t1inds])
  dev.off()
  
  plot_path8 <- file.path(plot_dir, "params-tab-2.png")
  png(plot_path8,
      width = 15,
      height = 2,
      units = "in",
      res = 90)
  gridExtra::grid.table(est_tab[, -t1inds])
  dev.off()
}

forecast_date <- Sys.getenv("fdt", unset = "2021-03-29")
forecast_loc <- Sys.getenv("loc", unset = "36")

fit_dir <-
  file.path("fits",
            paste0(forecast_date,
                   "-fips",
                   forecast_loc))

load(file.path(fit_dir, "fit.RData"))

ti <- max(x$target_end_date) + lubridate::ddays(1)
fdtw <- lubridate::wday(forecast_date)
if (fdtw > 2) {
  extra <-
    7 - fdtw ## extra days to ensure we sim up to the end of the 4 ahead target
} else {
  extra <- 0
}
tf <- lubridate::ymd(forecast_date) + lubridate::ddays(28 + extra)
target_end_dates <- seq(from = ti, to = tf, by = "day")
target_end_times <- lubridate::decimal_date(target_end_dates)
target_wday <- lubridate::wday(target_end_dates)
fets <- tibble(target_end_times, target_wday, target_end_dates)

write_forecasts(
  fit = fit1,
  cov = x,
  z = z,
  winit = winit,
  wfixed = wfixed,
  hess = h1,
  p_hsd = p_hsd,
  β_0sd = β_0sd, 
  τ_csd = τ_csd,
  fets = fets,
  fdt = forecast_date,
  forecast_loc = forecast_loc
)
