#!/usr/bin/env R

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

make_params_tab <- function(w, h, cov, z, f, dir, Rtmod){
  ests <- name_params(w, z, cov, f, trans = FALSE)[[1]]
  sevec <- sqrt(diag(solve(h)))
  sd <- name_params(sevec, z, cov, f, trans = FALSE)[[1]]
  ests$ϕ <- coef(Rtmod) %>% unname()
  sd$ϕ <- vcov(Rtmod) %>% sqrt() %>% unname()
  
  pnames <-
    c(
      "doseeffect",
      "residentialeffect",
      "L0",
      "p_d",
      "γ_d12",
      "γ_d34",
      "γ_z17",
      "τ_h",
      "p_hweekend",
      "τ_d",
      "ϕ"
    )
  
  pnames_pretty <-
    c(
      "$\\beta_{\\textrm{dose}}$",
      "$\\beta_{\\textrm{res}}$",
      "$L_0$",
      "$\\operatorname{logit} p_d$",
      "$\\log \\gamma_{d,1} = \\log \\gamma_{d,2}$",
      "$\\log \\gamma_{d,3} = \\log \\gamma_{d,4}$",
      "$\\log \\gamma_{z,1} = \\log \\gamma_{z,7}$",
      "$\\log \\tau_h$",
      "$\\operatorname{logit} p_{h,\\textrm{weekend}}$",
      "$\\log \\tau_d$",
      "$\\phi$"
    )
  if (!"doses_scaled" %in% names(cov)){
    linds <- !grepl("doseeffect", pnames)
    pnames <- pnames[linds]
    pnames_pretty <- pnames_pretty[linds]
  }
  if (all(is.na(z$hospitalizations))){
    linds <- !grepl("p_hweekend", pnames)
    pnames <- pnames[linds]
    pnames_pretty <- pnames_pretty[linds]
  }
  
  tibble(
    parameter = pnames_pretty,
    estimate = unlist(ests[pnames]),
    "standard error" = unlist(sd[pnames]),
  ) %>% knitr::kable(digits = 2, format = "latex", escape = FALSE) %>%
    cat(file = file.path(dir, "param_table.tex"))
  
  q <- qnorm(1 - 0.05 / 2)
  start <- x %>% group_by(p_hmap) %>% slice(1) %>% pull(target_end_date)
  dph <- tibble(x = start, y = ests$p_h, 
                lower = ests$p_h - sd$p_h * q,
                upper = ests$p_h + sd$p_h * q)
  p <- ggplot(dph, aes(x, y)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) + 
    labs(x = "Start of 4 week period", y = expression(paste("logit ", p[h])))
  
  ggsave(filename = file.path(dir, "p_h-plot.png"), plot = p)

  start2 <- x %>% group_by(τ_cmap) %>% slice(1) %>% pull(target_end_date)
  dtc <- tibble(x = start2, y = ests$τ_c, 
                lower = ests$τ_c - sd$τ_c * q,
                upper = ests$τ_c + sd$τ_c * q)
  ptc <- ggplot(dtc, aes(x, y)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) + 
    labs(x = "Start of 4 week period", y = expression(paste("log ", tau[c])))
  ggsave(filename = file.path(dir, "tauc-plot.png"), plot = ptc)
  
  start3 <- x %>% group_by(β_0map) %>% slice(1) %>% pull(target_end_date)
  db0 <- tibble(x = start3, y = ests$β_0, 
                lower = ests$β_0 - sd$β_0 * q,
                upper = ests$β_0 + sd$β_0 * q)
  pb0 <- ggplot(db0, aes(x, y)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper)) + 
    labs(x = "Start of 1 week period", y = expression(paste("log ", beta[0])))
  ggsave(filename = file.path(dir, "b0-plot.png"), plot = pb0)
  
  
}

calc_rt <- function(ft, z, cov, no_hosps, susceptibles, wfixed) {
  nβ_0 <- tail(cov$β_0map, 1)
  np <- length(ft$par)
  inds <- seq(np - nβ_0 + 1, np)
  intercept <- ft$par[inds][cov$β_0map]
  
  p <- name_params(ft$par, z, cov, as.list(wfixed))[[1]]
  
  if ("doses_scaled" %in% names(cov)){
    X <- cbind(cov$residential, cov$doses_scaled)
    effects <- unlist(p[c("residentialeffect", "doseeffect")])
  } else {
    X <- cbind(cov$residential)
    effects <- p[["residentialeffect"]]
  }
  
  (exp(intercept + X %*% effects) / wfixed["γ"] * susceptibles / wfixed["N"]) %>% as.numeric()
}

make_rt_plot <- function(ft, z, cov, no_hosps, susceptibles, wfixed) {
  par(mfrow = c(1, 1))
  nβ_0 <- tail(cov$β_0map, 1)
  np <- length(ft$par)
  inds <- seq(np - nβ_0 + 1, np)
  intercept <- ft$par[inds][cov$β_0map]
  
  p <- name_params(ft$par, z, cov, as.list(wfixed))[[1]]
  
  if ("doses_scaled" %in% names(cov)){
    X <- cbind(cov$residential, cov$doses_scaled)
    effects <- unlist(p[c("residentialeffect", "doseeffect")])
  } else {
    X <- cbind(cov$residential)
    effects <- p[["residentialeffect"]]
  }

  num_all <- exp(intercept + X %*% effects)
  factor <- 1 / wfixed["γ"] * susceptibles / wfixed["N"]
  plot(
    cov$time,
    num_all * factor,
    type = 'l',
    xlab = "Time",
    ylab = expression(paste("Effective reproduction number ", R[e]))
  )
  abline(h = 1000 / wfixed["γ"])
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
          num_nomob * factor,
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
          num_nomob * factor,
          col = "orange")
  }
  lines(cov$time,
        exp(intercept) * factor,
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

make_dashboard_input <- function(dets, cov, z, forecast_loc, fit, wfixed, cov_sim){
  abb <- covidcast::fips_to_abbr(paste0(forecast_loc, "000"))
  snm <- covidcast::fips_to_name(paste0(forecast_loc, "000"))
  pop <- covidHubUtils::hub_locations %>%
    filter(abbreviation == abb &
             geo_type == "state") %>% pull("population")
  
  if (ncol(z) == 2){
    z$hospitalizations <- NA
  }
  ret <- list()
  ret[[1]] <- bind_cols(
    location = snm,
    period = "Past",
    date = cov$target_end_date,
    populationsize = pop,
    state_abr = abb,
    z
  ) %>%
    mutate(cumulative_cases = cumsum(cases),
           cumulative_hosps = cumsum(hospitalizations),
           cumulative_deaths = cumsum(deaths)) %>%
    rename(
      actual_daily_cases = cases,
      actual_daily_hosps = hospitalizations,
      actual_daily_deaths = deaths,
      actual_cumulative_cases = cumulative_cases,
      actual_cumulative_hosps = cumulative_hosps,
      actual_cumulative_deaths = cumulative_deaths
    ) %>%
    pivot_longer(
      actual_daily_cases:actual_cumulative_deaths,
      names_to = "variable",
      values_to = "mean_value"
    ) %>%
    mutate(median_value = mean_value)
  
  cov_all <- rbind(cov[, names(cov_sim)], cov_sim)
  
  gendf <- function(mean, median = mean, sd, vname) {
    df <- bind_cols(
      location = snm,
      period = "Past",
      date = cov_all$target_end_date,
      populationsize = pop,
      state_abr = abb,
      variable = vname,
      mean_value = mean,
      median_value = median,
      lower_80 = qnorm(p = 0.1,
                       mean = mean,
                       sd = sd),
      upper_80 = qnorm(p = 0.9,
                       mean = mean,
                       sd = sd)
    ) %>% mutate(lower_80 = ifelse(lower_80 < 0, 0, lower_80))
  }
  
  ret[[2]] <- gendf(mean = c(dets$xhat_kk["Y",], dets$x_sim[2,]),
                    sd = sqrt(c(dets$P_kk[2, 2,], dets$P_sim[2,2,])),
                    vname = "daily_all_infections")
  
  ret[[3]] <- gendf(
    mean = wfixed["N"] - c(dets$xhat_kk["X",], dets$x_sim[1,]),
    median = NA,
    sd = sqrt(c(dets$P_kk[1, 1, ], dets$P_sim[1, 1, ]  )),
    vname = "cumulative_all_infections"
  )
  
  ret[[4]] <- gendf(mean = c(dets$xhat_kkmo["Z_r",], dets$sim_means[1, ]),
                    sd = sqrt(c(dets$S[1, 1,], dets$sim_cov[1,1,])),
                    vname = "daily_cases")
  
  ret[[5]] <- gendf(
    mean = cumsum(c(dets$xhat_kkmo["Z_r",], dets$sim_means[1, ])),
    median = NA,
    sd = NA,
    vname = "cumulative_cases"
  )
  
  ret[[6]] <- gendf(mean = c(dets$xhat_kkmo["A",], dets$sim_means[2,]),
                    sd = sqrt(c(dets$S[2, 2,], dets$sim_cov[2,2,])),
                    vname = "daily_hosps")
  
  ret[[7]] <- gendf(
    mean = cumsum( c(dets$xhat_kkmo["A",], dets$sim_means[2,] )   )  ,
    median = NA,
    sd = NA,
    vname = "cumulative_hosps"
  )
  
  ret[[8]] <- gendf(mean = c(dets$xhat_kkmo["D_r",], dets$sim_means[3,]),
                    sd = sqrt(c(dets$S[3, 3,], dets$sim_cov[3,3,])),
                    vname = "daily_deaths")
  
  ret[[9]] <- gendf(
    mean = cumsum(c(dets$xhat_kkmo["D_r",], dets$sim_means[3,])),
    median = NA,
    sd = NA,
    vname = "cumulative_deaths"
  )

  
  no_hosps <- all(is.na(dets$ytilde_k[2, ]))
  rt <- calc_rt(fit, z, cov_all, no_hosps, susceptibles = dets$xhat_kk[1,], wfixed)
  
  ret[[10]] <- gendf(
    mean = rt,
    median = NA,
    sd = NA,
    vname = "combined_trend"
  )
  
  if ("doses_scaled" %in% names(cov)){
    ret[[11]] <- gendf(
      mean = fit$par[5] * cov_all$doses_scaled,
      median = NA,
      sd = NA,
      vname = "vaccine_effect"
    )
  } else {
    ret[[11]] <- gendf(
      mean = 0,
      median = NA,
      sd = NA,
      vname = "vaccine_effect"
    )
  }
  
  retall <-
    bind_rows(ret) %>% select(
      location,
      period,
      date,
      variable,
      lower_80,
      mean_value,
      median_value,
      upper_80,
      populationsize,
      state_abr
    )
  retall
}

write_scenarios <- function(fit,
                            cov,
                            z,
                            wfixed,
                            cov_sim,
                            p_hsd,
                            β_0sd,
                            τ_csd,
                            forecast_loc,
                            forecast_date) {

  get_dets <- function(cs){
    calc_kf_nll_r(
      w = fit$par,
      cov = cov,
      z = z,
      p_hsd = p_hsd,
      β_0sd =  β_0sd,
      τ_csd =  τ_csd,
      wfixed = wfixed,
      just_nll = FALSE,
      cov_sim = cs
    )
  }
  dets <- purrr::map(cov_sim, get_dets)
  tmpf <- function(xx, yy){
    make_dashboard_input(
      dets = xx,
      cov = cov,
      z = z,
      forecast_loc = forecast_loc,
      fit = fit,
      wfixed = wfixed,
      cov_sim = yy
    )
  }
  usd <- map2(dets, cov_sim, tmpf) %>% bind_rows(.id = "scenario")
  ## remove duplicated real data time series
  usd2 <- usd %>% 
    filter(!(scenario %in% c("linear_increase_sd", "return_normal") & 
               str_detect(variable, "^actual|^vaccine")))
  test <- str_detect(usd2$variable, "^actual|^vaccine")
  usd2$scenario[test] <- NA
  fcst_dir <-
    file.path("forecasts",
              paste0(forecast_date,
                     "-fips",
                     forecast_loc))
  
  fname <- paste0("dashboard.rds")
  
  path <- file.path(fcst_dir, fname)
  if (!dir.exists(fcst_dir))
    dir.create(fcst_dir, recursive = TRUE)
  write_rds(list(us_dat = usd2), path = path)
}

write_forecasts <-
  function(fit, cov, z, winit, wfixed, hess, cov_sim, p_hsd, β_0sd, τ_csd, fdt, 
           forecast_loc, Rtmod) {
    
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
        fet_zero_cases_deaths = "weekly",
        cov_sim = cov_sim
      )

    case_inds <- which(cov_sim$wday == 7)
    
    if (lubridate::wday(fdt) > 2) {
      case_inds <-
        case_inds[-1] # drop first epiweek, which has been observed too much 
    }
    
    week_start_inds <- case_inds - 6
    
    sum_seq <- function(x, y, varno){
      inds <- seq(from = x, to = y)
      sum(dets$sim_R[varno, varno, inds])
    }
    case_obs_vars <- purrr::map2_dbl(week_start_inds, case_inds, sum_seq, varno = 1)
    case_proc_vars <- dets$sim_P[1,1, case_inds]
    case_Sigma <- case_obs_vars + case_proc_vars
    
    case_fcst <-
      create_forecast_df(
        means = dets$sim_means[1, case_inds],
        vars = case_Sigma,
        location = forecast_loc,
        fdt = fdt
      )
    stopifnot(setequal(
      cov_sim$target_end_dates[case_inds],
      case_fcst$target_end_date %>% unique()
    ))
    
    hosp_inds <- cov_sim$target_end_dates %in%
      (lubridate::ymd(forecast_date) + 1:28)
    hosp_fcst <-
      create_forecast_df(
        means = dets$sim_means[2, hosp_inds],
        vars = dets$sim_Sigma[2, 2, hosp_inds],
        target_type = "hospitalizations",
        location = forecast_loc,
        fdt = fdt
      )
    stopifnot(setequal(
      cov_sim$target_end_dates[hosp_inds],
      hosp_fcst$target_end_date %>% unique()
    ))
    
    death_inds <- case_inds
    week_start_inds2 <- death_inds - 6
    
    death_obs_vars <- purrr::map2_dbl(week_start_inds2, death_inds, sum_seq, varno = 3)
    death_proc_vars <- dets$sim_P[3,3, death_inds]
    death_Sigma <- death_obs_vars + death_proc_vars

    death_fcst <-
      create_forecast_df(
        means = dets$sim_means[3, death_inds],
        vars = death_Sigma,
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
      wfixed = wfixed,
      cov_sim = cov_sim,
      z,
      Rtmod = Rtmod
    )
}

make_fit_plots <- function(dets, cov, winit, h, p_hsd, β_0sd, τ_csd, fdt, 
                           forecast_loc, fit, wfixed, cov_sim, z, Rtmod) {
  plot_dir <-
    file.path("plots",
              paste0(fdt,
                     "-fips",
                     forecast_loc),
              paste0("lambda",
                     sprintf("%06.2f", 1 / β_0sd)))
  
  if (!dir.exists(plot_dir))
    dir.create(plot_dir, recursive = TRUE)
    
  make_params_tab(fit$par, h, cov, z, as.list(wfixed), plot_dir, Rtmod)
  
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

    pred_hosps <- dets$xhat_kkmo["A", ]
    se_hosps <- sqrt(dets$S[2, 2, ])
    plot(x$time,
         pred_hosps,
         type = "l",
         xlab = "Time",
         ylab = "Hospitalizations")
    lines(x$time, se_hosps * 2 + pred_hosps, col = "grey")
    points(x$time, z$hospitalizations)
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
  
  lambda <- 1 / β_0sd
  fcst_dir <-
    file.path("forecasts",
              paste0(fdt,
                     "-fips",
                     forecast_loc))
  
  fcst_name <- paste0("lambda",
                      sprintf("%06.2f", lambda),
                      "-status-quo",
                      "-CEID-InfectionKalman.csv")
  
  fcst <- read_csv(
    file.path(fcst_dir, fcst_name),
    col_types = cols(
      forecast_date = col_date(format = ""),
      target = col_character(),
      target_end_date = col_date(format = ""),
      location = col_character(),
      type = col_character(),
      quantile = col_character(),
      value = col_double()
    )
  )
  
  p3 <- fcst %>% 
    filter(str_detect(target, "wk ahead inc case$")) %>%
    ggplot(aes(
      x = target_end_date,
      y = as.numeric(value),
      color = quantile
    )) +
    geom_line() + geom_point() + labs(y = "Weekly cases")
  ggsave(plot_path3, p3)
  
  plot_path4 <- file.path(plot_dir, "hosp-forecasts.png")

  p4 <- fcst %>%
    filter(str_detect(target, "day ahead inc hosp$")) %>%
    ggplot(aes(
      x = target_end_date,
      y = as.numeric(value),
      color = quantile
    )) +
    geom_line() + geom_point() + labs(y = "Daily hospital admissions")
  ggsave(plot_path4, p4)
  
  plot_path5 <- file.path(plot_dir, "death-forecasts.png")

  p5 <- fcst %>% 
    filter(str_detect(target, "wk ahead inc death$")) %>%
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
  make_rt_plot(fit, z, cov, no_hosps, susceptibles = dets$xhat_kk[1, ], wfixed)
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
  
  plot_path9 <- file.path(plot_dir, "mobility.png")
  pmob <- cov %>% ggplot(aes(x = target_end_date, y = residential_pcb)) + 
    geom_point(color = "grey") + 
    geom_point(aes(y = residential * 100)) + 
    theme_minimal() + 
    labs(x = "Date", y = "Percent increase in time spent\nin residential areas")
  ggsave(plot_path9, pmob,  dpi = 600, width = 5.2)
  
  if ("doses" %in% names(cov)){
    pvac_path <- file.path(plot_dir, "doses.png")
    pvac <- cov %>% 
      filter(target_end_date >= "2020-12-01") %>% 
      ggplot(aes(x = target_end_date, y = doses)) + 
      geom_line() + 
      theme_minimal() + 
      labs(x = "Date", y = "Doses administered")
    ggsave(pvac_path, pvac, dpi = 600, width = 5.2)
  }
  
  prhot <- cov %>% filter(target_end_date < "2020-07-20") %>% 
    ggplot(aes(x = target_end_date, y = ρ)) + geom_line() +
    xlab("Date") +
    ylab(expression(paste(rho[t]," = Pr(Removal is reported)"))) +
    theme_minimal()
  
  prhot_path <- file.path(plot_dir, "rhot-estimate.png")
  ggsave(prhot_path, prhot, dpi = 600, width = 5.2)
}

forecast_date <- Sys.getenv("fdt", unset = "2020-06-29")
forecast_loc <- Sys.getenv("loc", unset = "06")

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
sim_weekly_days <- 28 + extra
tf <- lubridate::ymd(forecast_date) + lubridate::ddays(sim_weekly_days)
target_end_dates <- seq(from = ti, to = tf, by = "day")
target_end_times <- lubridate::decimal_date(target_end_dates)
target_wday <- lubridate::wday(target_end_dates)

if ("doses_scaled" %in% names(x)){
  doses_slope <- mean(diff(tail(x$doses_scaled, n = 7)))
  weekly_sim_doses <- tail(x$doses_scaled, n = 1) + seq_len(length(target_end_dates)) * doses_slope
  
  cov_weekly_fcst <- tibble(
    time = target_end_times,
    wday = target_wday,
    target_end_dates,
    doses_scaled = weekly_sim_doses,
    residential = tail(x$residential, n = 1),
    β_0map = tail(x$β_0map, 1),
    p_hmap = tail(x$p_hmap, 1),
    ρ = tail(x$ρ, 1)
  )
} else {
  cov_weekly_fcst <- tibble(
    time = target_end_times,
    wday = target_wday,
    target_end_dates,
    residential = tail(x$residential, n = 1),
    β_0map = tail(x$β_0map, 1),
    p_hmap = tail(x$p_hmap, 1),
    ρ = tail(x$ρ, 1)
  )
}

# AR-1 modeling of R_t

p <- name_params(fit1$par, z, x, as.list(wfixed))[[1]]
if ("doses_scaled" %in% names(x)) {
  β   <-
    exp(p$β_0[x$β_0map] + p$residentialeffect * x$residential + p$doseeffect * x$doses_scaled)
} else {
  β   <-
    exp(p$β_0[x$β_0map] + p$residentialeffect * x$residential)
}

include <- x$target_end_date > "2020-05-01"
ar1ts <- β[include] - wfixed["γ"]
mod <- arima(ar1ts, order = c(1, 0, 0), include.mean = FALSE)
cov_weekly_fcst$β <- predict(mod, n.ahead = nrow(cov_weekly_fcst))$pred + wfixed["γ"]

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
  cov_sim = cov_weekly_fcst,
  fdt = forecast_date,
  forecast_loc = forecast_loc,
  Rtmod = mod
)

q("no")
sim_days <- 42
tf2 <- ti + lubridate::ddays(sim_days - 1)
target_end_dates2 <- seq(from = ti, to = tf2, by = "day")
target_end_times2 <- lubridate::decimal_date(target_end_dates2)
target_wday2 <- lubridate::wday(target_end_dates2)

doses_slope <- mean(diff(tail(x$doses_scaled, n = 7)))
sim_doses <- tail(x$doses_scaled, n = 1) + seq_len(sim_days) * doses_slope

cov_status_quo <-
  tibble(
    time = target_end_times2,
    wday = target_wday2,
    target_end_date = target_end_dates2,
    doses_scaled = sim_doses,
    residential = tail(x$residential, n = 1),
    β_0map = tail(x$β_0map, 1)
  )

res_increase_slope <- max((30 / 100 - tail(x$residential, n = 1) ) / sim_days, 0)
res_increase_residential <- tail(x$residential, n = 1) + seq_len(sim_days) * res_increase_slope

cov_linear_increase_sd <-
  tibble(
    time = target_end_times2,
    wday = target_wday2,
    target_end_date = target_end_dates2,
    doses_scaled = sim_doses,
    residential = res_increase_residential,
    β_0map = tail(x$β_0map, 1)
  )

return_normal_slope <- min((0 - tail(x$residential, n = 1) ) / sim_days, 0)
return_normal_residential <- tail(x$residential, n = 1) + seq_len(sim_days) * return_normal_slope

cov_return_normal <-
  tibble(
    time = target_end_times2,
    wday = target_wday2,
    target_end_date = target_end_dates2,
    doses_scaled = sim_doses,
    residential = return_normal_residential,
    β_0map = tail(x$β_0map, 1)
  )

write_scenarios(
  fit = fit1,
  cov = x,
  z = z,
  wfixed = wfixed,
  cov_sim = list(status_quo = cov_status_quo, 
                 linear_increase_sd = cov_linear_increase_sd,
                 return_normal = cov_return_normal),
  p_hsd = p_hsd,
  β_0sd = β_0sd,
  τ_csd = τ_csd,
  forecast_loc = forecast_loc,
  forecast_date = forecast_date
)
