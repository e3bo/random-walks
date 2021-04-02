#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

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
  n <- nrow(means)
  stopifnot(nrow(vars) == n)
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
    mutate(q1 = map2_dbl(
      quantile,
      h,
      ~ KScorrect::qmixnorm(
        p = .x,
        mean = means[.y,],
        sd = sqrt(vars[.y,]),
        pro = rep(1, ncol(means)),
        expand = 0.0
      )
    ),
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

write_forecasts <- function(fits, fet, agrid, betagrid, fdt = forecast_date) {
  n1 <- length(agrid)
  n2 <- length(betagrid)
  for (ind1 in 1:n1) {
    for (ind2 in 1:n2) {
      wfit <- fits[[ind1]][[ind2]]$par
      a <- agrid[ind1]
      betasd <- betagrid[ind2]
      dets <-
        kf_nll_details(
          wfit,
          x,
          y,
          param_map,
          betasd = betasd,
          a = a,
          fet = fet
        )
      make_fit_plots(dets, x = x, a = a, betasd = betasd)
      case_inds <- which(fet$target_wday == 7)
      if (lubridate::wday(fdt) > 2){
        case_inds <- case_inds[-1] # drop first epiweek, which has been observed too much to be forecasted
      }
      case_fcst <- create_forecast_df(means = dets$sim_means[1, case_inds,],
                                      vars = dets$sim_cov[1, 1, case_inds,],
                                      location = forecast_loc)
      stopifnot(setequal(
        fet$target_end_dates[case_inds],
        case_fcst$target_end_date %>% unique()
      ))
      hosp_inds <- fet$target_end_dates %in% 
        (lubridate::ymd(forecast_date) + 1:28)
      hosp_fcst <- create_forecast_df(means = dets$sim_means[2, hosp_inds,],
                                      vars = dets$sim_cov[2, 2, hosp_inds,],
                                      target_type = "hospitalizations",
                                      location = forecast_loc)
      stopifnot(setequal(
        fet$target_end_dates[hosp_inds],
        hosp_fcst$target_end_date %>% unique()
      ))
      
      death_inds <- case_inds
      death_fcst <- create_forecast_df(means = dets$sim_means[3, death_inds,],
                                       vars = dets$sim_cov[3, 3, death_inds,],
                                       target_type = "inc deaths",
                                       location = forecast_loc)
      fcst <- bind_rows(case_fcst, hosp_fcst, death_fcst)
      
      lambda <- 1 / betasd
      fcst_dir <-
        file.path(
          "forecasts",
          paste0(
            forecast_date,
            "-fips",
            forecast_loc))
      
      fcst_name <- paste0("lambda",
                          sprintf("%06.2f", lambda),
                          "-a",
                          sprintf("%02.2f", a),
                          "-CEID-InfectionKalman.csv"
      )
      fcst_path <- file.path(fcst_dir, fcst_name)
      if (!dir.exists(fcst_dir))
        dir.create(fcst_dir, recursive = TRUE)
      write_csv(x = fcst, path = fcst_path)
    }
  }
}

kf_nll_details <- function(w, x, y, betasd, a, pm, fet) {
  p <- pm(x, w)
  nll <- kfnll(
    logbeta = p$bpars,
    logE0 = p$logE0,
    logH0 = p$logH0,
    logtauc = p$logtauc,
    logtauh = p$logtauh,
    logtaud = p$logtaud,
    logchp = p$logchp,
    loghfp = p$loghfp,
    loggammahd = p$loggammahd,
    logdoseeffect = p$logdoseeffect,
    logprophomeeffect = p$logprophomeeffect,
    eta = p$eta,
    gamma = p$gamma,
    N = p$N,
    z = y,
    t0 = p$t0,
    times = p$times,
    doses = p$doses,
    prophome = p$prophome,
    fet = fet,
    just_nll = FALSE,
    fet_zero_cases_deaths = "weekly",
    nsim = 20,
    betasd = betasd,
    a = a
  )
  nll
}



kfnll <-
  function(logbeta,
           logE0,
           logH0,
           logtauc,
           logtauh,
           logtaud,
           logchp,
           loghfp,
           loggammahd,
           logdoseeffect,
           logprophomeeffect,
           eta,
           gamma,
           N,
           z,
           t0,
           times,
           doses,
           prophome,
           Phat0 = diag(c(1, 1, 1, 0, 0, 1, 1, 0)),
           fets = NULL,
           fet_zero_cases_deaths = "daily",
           nsim,
           a = .98,
           betasd = 1,
           maxzscore = Inf,
           just_nll = TRUE) {
    E0 <- exp(logE0)
    I0 <- E0 * eta / gamma
    H0 <- exp(logH0)
    gamma_d <- gamma_h <- exp(loggammahd)
    D0 <- H0 * gamma_h / gamma_d
    xhat0 <- c(N - E0 - I0 - H0 - D0, E0, I0, 0, 0, H0, D0, 0)
    names(xhat0) <- c("S", "E", "I", "C", "Hnew", "H", "D", "Drep")
    doseeffect <- exp(logdoseeffect)
    prophomeeffect <- exp(logprophomeeffect)

    if (ncol(z) == 1 && "cases" %in% names(z)){
      z$hospitalizations <- NA
      z$deaths <- NA
    }
    
    z <- data.matrix(z[, c("cases", "hospitalizations", "deaths")])
    is_z_na <- is.na(z)
    T <- nrow(z)
    dobs <- ncol(z)
    dstate <- length(xhat0)
    stopifnot(T > 0)
    
    ytilde_kk <- ytilde_k <- array(NA_real_, dim = c(dobs, T))
    S <- array(NA_real_, dim = c(dobs, dobs, T))
    K <- array(NA_real_, dim = c(dstate, dobs, T))
    rdiagadj <- array(1, dim = c(dobs, T))
    
    xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(dstate, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(dstate, dstate, T))
    
    H <- function(time, t0 = 2020.164){
      day <- (time - t0) * 365.25
      rbind(c(0, 0, 0, detect_frac(day), 1, 0, 0, 0), 
            c(0, 0, 0,    0, 1, 0, 0, 0),
            c(0, 0, 0,    0, 0, 0, 0, 1))
    }
    R <- diag(exp(c(logtauc, logtauh, logtaud)))
    
    
    for (i in seq(1, T)) {
      if (i == 1) {
        xhat_init <- xhat0
        PNinit <- Phat0
      } else {
        xhat_init <- xhat_kk[, i - 1]
        PNinit <- P_kk[, , i - 1]
      }
      
      xhat_init["C"] <- 0
      xhat_init["Hnew"] <- 0
      xhat_init["Drep"] <- 0
      
      PNinit[, 4] <- PNinit[4,] <- 0
      PNinit[, 5] <- PNinit[5,] <- 0
      PNinit[, 8] <- PNinit[8,] <- 0
      
      #if(i == 360) browser()
      XP <- iterate_f_and_P(
        xhat_init,
        PN = PNinit,
        eta = eta,
        gamma = gamma,
        gamma_d = gamma_d,
        gamma_h = gamma_h,
        chp = exp(logchp),
        hfp = exp(loghfp),
        N = N,
        beta_t = exp(logbeta[i] -doseeffect * doses[i] - prophomeeffect * prophome[i]),
      )
      xhat_kkmo[, i] <- XP$xhat
      P_kkmo[, , i] <- XP$PN

      for (j in 1:dstate){
        if (P_kkmo[j,j,i] < 0){
          P_kkmo[j,,i] <- 0
          P_kkmo[,j,i] <- 0
        }
      }
      
      ytilde_k[, i] <- matrix(z[i, ], ncol = 1) - 
        H(x$time[i]) %*% xhat_kkmo[, i, drop = FALSE]     
      S[, , i] <- H(x$time[i]) %*% P_kkmo[, , i] %*% t(H(x$time[i])) + R
      
      for (j in 1:dobs){
        if (is.na(z[i,j])){
          zscore <- 0
          rdiagadj[j,i] <- rdiagadj[j,i] + 0
        } else {
          sd <- sqrt(S[j,j,i])
          zscore <- ytilde_k[j,i] / sd 
          if (abs(zscore) > maxzscore){
            adjzscore <- maxzscore / (1 + abs(zscore) - maxzscore)
            newsd <- abs(ytilde_k[j,i]) / adjzscore
            rdiagadj[j,i] <- rdiagadj[j,i] + (newsd) ^ 2 - sd ^ 2
          } else {
            rdiagadj[j,i] <- rdiagadj[j,i] + 0
          }
        }
        S[j,j,i] <- S[j,j,i] + rdiagadj[j,i]
      }
      K[, , i] <- P_kkmo[, , i] %*% t(H(x$time[i])) %*% solve(S[, , i])
      desel <- is_z_na[i, ]
      K[, desel, i] <- 0
      
      xhat_kk[, i] <-
        xhat_kkmo[, i, drop = FALSE] +
        K[, !desel, i] %*% ytilde_k[!desel, i, drop = FALSE]
      
      xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
      P_kk[, , i] <-
        (diag(dstate) - K[, , i] %*% H(x$time[i])) %*% P_kkmo[, , i]
      ytilde_kk[, i] <- matrix(z[i, ], ncol = 1) - 
        H(x$time[i]) %*% xhat_kk[, i, drop = FALSE]
    }
    
    rwlik <- 0
    for (i in 1:(T - 1)){
      step <- (logbeta[i + 1] - log(gamma)) - a * (logbeta[i] - log(gamma))
      rwlik <- rwlik + dnorm(step, mean = 0, sd = betasd, log = TRUE)
    }

    nll <- 0
    for (i in seq(1, T)){
      sel <- !is_z_na[i, ]
      nll <- nll + 
        t(ytilde_k[sel, i]) %*% solve(S[sel, sel, i]) %*% ytilde_k[sel, i] + 
        log(det(S[,,i][sel, sel, drop = FALSE])) + dobs * log(2 * pi)
    }
    nll <- 0.5 * nll - rwlik
    
    if (!just_nll) {
      if (!is.null(fets)) {
        nsimdays <- nrow(fets)
        sim_means <- array(NA_real_, dim = c(dobs, nsimdays, nsim))
        sim_cov <- array(NA_real_, dim = c(dobs, dobs, nsimdays, nsim))
        for (j in seq_len(nsim)) {
          logbeta_fet <- numeric(nsimdays)
          logbeta_fet[1] <-
            log(gamma) + a * (logbeta[T] - log(gamma)) +  
            rnorm(1, mean = 0, sd = betasd)
          if (length(logbeta_fet) > 1) {
            for (jj in seq(2, length(logbeta_fet))) {
              logbeta_fet[jj] <-
                log(gamma) + a * (logbeta_fet[jj - 1] - log(gamma)) + 
                rnorm(1, mean = 0, sd = betasd)
            }
          }
          xhat_init <- xhat_kk[, T]
          PNinit <- P_kk[, , T]
          if (fet_zero_cases_deaths == "daily" ||
              fets$target_wday[1] == 1) {
            xhat_init["C"] <- 0
            PNinit[, 4] <- PNinit[4,] <- 0
            xhat_init["Drep"] <- 0
            PNinit[, 8] <- PNinit[8,] <- 0 
          }
          xhat_init["Hnew"] <- 0
          PNinit[, 5] <- PNinit[5,] <- 0
          
          XP <- iterate_f_and_P(
            xhat_init,
            PN = PNinit,
            eta = eta,
            gamma = gamma,
            gamma_d = gamma_d,
            gamma_h = gamma_h,
            chp = exp(logchp),
            hfp = exp(loghfp),
            N = N,
            beta_t = exp(logbeta_fet[1]) * exp(-doseeffect * doses[T])
          )
          sim_means[, 1, j] <- H(fets$target_end_times[1]) %*% XP$xhat
          sim_cov[, , 1, j] <- H(fets$target_end_times[1]) %*% XP$PN %*% t(H(fets$target_end_times[1])) + R
          for (i in seq_along(fets$target_end_times[-1])) {
            xhat_init <- XP$xhat
            PNinit <- XP$PN
            if (fet_zero_cases_deaths == "daily" ||
                fets$target_wday[i + 1] == 1) {
              xhat_init["C"] <- 0
              PNinit[, 4] <- PNinit[4,] <- 0
              xhat_init["Drep"] <- 0
              PNinit[, 8] <- PNinit[8,] <- 0
            }
            xhat_init["Hnew"] <- 0
            PNinit[, 5] <- PNinit[5,] <- 0
            
            XP <-
              iterate_f_and_P(
                xhat_init,
                PN = PNinit,
                eta = eta,
                gamma = gamma,
                gamma_d = gamma_d,
                gamma_h = gamma_h,
                chp = exp(logchp),
                hfp = exp(loghfp),
                N = N,
                beta_t = exp(logbeta_fet[i + 1]) * exp(-doseeffect * doses[T])
              )
            sim_means[, i + 1, j] <- H(fets$target_end_times[i + 1]) %*% XP$xhat
            sim_cov[, , i + 1, j] <-
              H(fets$target_end_times[i + 1]) %*% XP$PN %*% t(H(fets$target_end_times[i + 1])) + R
          }
        }
      } else {
        sim_means <- sim_cov <- NULL
      }
      
      list(
        nll = nll,
        xhat_kkmo = xhat_kkmo,
        xhat_kk = xhat_kk,
        P_kkmo = P_kkmo,
        P_kk = P_kk,
        ytilde_k = ytilde_k,
        S = S,
        sim_means = sim_means,
        sim_cov = sim_cov,
        logbeta = logbeta,
        doseeffect = doseeffect,
        gamma = gamma,
        rdiagadj = rdiagadj
      )
    } else {
      nll
    }
  }

make_fit_plots <- function(dets, x, a, betasd) {
  
  plot_dir <-
    file.path(
      "plots",
      paste0(
        forecast_date,
        "-fips",
        forecast_loc),
      paste0(
        "lambda",
        sprintf("%06.2f", 1 / betasd),
        "-a",
        sprintf("%02.2f", a)))
  
  if (!dir.exists(plot_dir))
    dir.create(plot_dir, recursive = TRUE)
  
  plot_path1 <- file.path(plot_dir, "qqplots.png")
  
  png(
    plot_path1,
    width = 7.5,
    height = 10,
    units = "in",
    res = 90
  )
  
  par(mfrow = c(3, 1))
  qqnorm(dets$ytilde_k[1, ] / sqrt(dets$S[1, 1, ]), sub = "Cases")
  abline(0, 1)
  qqnorm(dets$ytilde_k[2, ] / sqrt(dets$S[2, 2, ]), sub = "Hospitalizations")
  abline(0, 1)
  qqnorm(dets$ytilde_k[3, ] / sqrt(dets$S[3, 3, ]), sub = "Deaths")
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
  
  
  rho_t <- detect_frac(seq_len(nrow(y)))
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
  dev.off()
  
  plot_path3 <- file.path(plot_dir, "case-forecasts.png")
  case_inds <- which(fet$target_wday == 7)
  case_fcst <-
    create_forecast_df(
      means = dets$sim_means[1, case_inds, ],
      vars = dets$sim_cov[1, 1, case_inds, ],
      location = forecast_loc
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
    fet$target_end_dates %in% (lubridate::ymd(forecast_date) + 1:28)
  stopifnot(sum(hosp_inds) == 28)
  
  hosp_fcst <-
    create_forecast_df(
      means = dets$sim_means[2, hosp_inds, ],
      vars = dets$sim_cov[2, 2, hosp_inds, ],
      target_type = "hospitalizations",
      location = forecast_loc
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
      means = dets$sim_means[3, case_inds, ],
      vars = dets$sim_cov[3, 3, case_inds, ],
      target_type = "inc deaths",
      location = forecast_loc
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
  
  df <- tibble(
    time = x$time,
    beta_t = exp(dets$logbeta),
    doses = x$doses,
    gamma = dets$gamma,
    de = dets$doseeffect
  ) %>%
    mutate(
      Rtnovacc = beta_t / gamma,
      u = exp(-doses * de),
      Rt = beta_t * u / gamma,
      vacc_reduction = Rtnovacc - Rt
    ) %>%
    select(time, vacc_reduction, Rt) %>%
    pivot_longer(-time) %>%
    mutate(Estimate = factor(
      name,
      levels = c("vacc_reduction", "Rt"),
      labels = c("Reduction in reproduction number due to vaccination", "Value")
    ))
  
  p6 <- ggplot(df, aes(x = time, y = value, fill = Estimate)) +
    geom_area() +
    labs(x = "Time", y = expression(paste("Instantaneous reproduction number, ", R[t]))) +
    theme_light() +
    theme(legend.position = "top")
  ggsave(filename = plot_path6, plot = p6)
}

forecast_date <- Sys.getenv("fdt", unset = "2021-03-20")
forecast_loc <- Sys.getenv("loc", unset = "36")

fit_dir <-
  file.path(
    "fits",
    paste0(
      forecast_date,
      "-fips",
      forecast_loc))

load(file.path(fit_dir, "fit.RData"))
write_forecasts(fits, fet, agrid, betasdgrid)
