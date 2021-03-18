#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

JuliaCall::julia_setup("/opt/julia-1.5.3/bin")
JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")

## main functions

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

iterate_f_and_P <-
  function(xhat,
           PN,
           eta,
           gamma,
           gamma_d,
           gamma_h,
           chp,
           hfp,
           N,
           beta_t,
           time.steps) {
    P <- PN / N
    dt <- diff(time.steps)
    vf <-  with(as.list(c(xhat, beta_t)), {
      F <- rbind(
        c(-beta_t * I / N,    0,    -beta_t * S / N, 0, 0,             0,        0, 0),
        c( beta_t * I / N, -eta,     beta_t * S / N, 0, 0,             0,        0, 0),
        c(              0,  eta,             -gamma, 0, 0,             0,        0, 0),
        c(              0,    0,  (1 - chp) * gamma, 0, 0,             0,        0, 0),
        c(              0,    0,        chp * gamma, 0, 0,             0,        0, 0),
        c(              0,    0,        chp * gamma, 0, 0,      -gamma_h,        0, 0),
        c(              0,    0,                  0, 0, 0, hfp * gamma_h, -gamma_d, 0),
        c(              0,    0,                  0, 0, 0,             0,  gamma_d, 0))
      
      f <-
        c(
          0,
          beta_t * (S / N) * (I / N),
          eta * E / N,
          (1 - chp) * gamma * I / N,
          chp * gamma * I / N,
          (1 - hfp) * gamma_h * H / N,
          hfp * gamma_h * H / N,
          gamma_d * D / N
        )
      
      Q <-
        rbind(
          c(f[1] + f[2],       -f[2],                  0,     0,     0,             0,           0,    0),
          c(      -f[2], f[2] + f[3],              -f[3],     0,     0,             0,           0,    0),
          c(          0,       -f[3], f[3] + f[4] + f[5], -f[4], -f[5],         -f[5],           0,    0),
          c(          0,           0,              -f[4],  f[4],     0,             0,           0,    0),
          c(          0,           0,              -f[5],     0,  f[5],             0,           0,    0),
          c(          0,           0,              -f[5],     0,     0,f[7]+f[5]+f[6],       -f[7],    0),
          c(          0,           0,                  0,     0,     0,         -f[7], f[8] + f[7], -f[8]),
          c(          0,           0,                  0,     0,     0,             0,       -f[8],  f[8])
        )
      
      dS <- -beta_t * S * I / N
      dE <- beta_t * S * I / N - eta * E
      dI <- eta * E -  gamma * I
      dC <- gamma * (1 - chp) * I
      dHnew <- gamma * chp * I
      dH <- gamma * chp * I - gamma_h * H
      dD <- gamma_h * H * hfp - gamma_d * D
      dDrep <- gamma_d * D
      
      dP <-  F %*% P + P %*% t(F) + Q
      
      list(vf = c(dS, dE, dI, dC, dHnew, dH, dD, dDrep),
           dP = dP)
    })
    xhat_new <- xhat + vf$vf * dt
    xhat_new[xhat_new < 0] <- 0
    P_new <- P + vf$dP * dt
    names(xhat_new) <- c("S", "E", "I", "C", "Hnew", "H", "D", "Drep")
    PN_new <- P_new * N
    list(xhat = xhat_new, PN = PN_new)
  }

kfnll <-
  function(bpars,
           logE0,
           logH0,
           rho1,
           logtauc,
           logtauh,
           logtaud,
           logchp,
           loghfp,
           eta,
           gamma,
           gamma_h,
           gamma_d,
           N,
           z,
           t0,
           times,
           Phat0 = diag(c(1, 1, 1, 0, 0, 1, 1, 0)),
           fets = NULL,
           fet_zero_cases_deaths = "daily",
           nsim,
           a = .98,
           betasd = 1,
           maxzscore = Inf,
           just_nll = TRUE,
           logmaxRt = 1.6) {
    E0 = exp(logE0)
    I0 = E0 * eta / gamma
    H0 = exp(logH0)
    D0 = H0 * gamma_h / gamma_d
    xhat0 = c(N - E0 - I0 - H0 - D0, E0, I0, 0, 0, H0, D0, 0)
    names(xhat0) <- c("S", "E", "I", "C", "Hnew", "H", "D", "Drep")
    
    z <- data.matrix(z)
    is_z_na <- is.na(z)
    T <- nrow(z)
    dobs <- ncol(z)
    dstate <- length(xhat0)
    stopifnot(T > 0)
    
    logbeta <- array(NA_real_, dim = c(T))
    ytilde_kk <- ytilde_k <- array(NA_real_, dim = c(dobs, T))
    S <- array(NA_real_, dim = c(dobs, dobs, T))
    K <- array(NA_real_, dim = c(dstate, dobs, T))
    rdiagadj <- array(NA_real_, dim = c(dobs, T))
    
    xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(dstate, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(dstate, dstate, T))
    
    H <- rbind(c(0, 0, 0, rho1, 1, 0, 0, 0), 
               c(0, 0, 0,    0, 1, 0, 0, 0),
               c(0, 0, 0,    0, 0, 0, 0, 1))
    R <- diag(exp(c(logtauc, logtauh, logtaud)))
    logbeta[T] <- bpars[T]
    for (i in seq(T - 1, 1)){
      logbeta[i] <- min((logbeta[i + 1] - log (gamma) - bpars[i]) / a + log(gamma), maxlogRt + log(gamma))
    }

    for (i in seq(1, T)) {
      if (i == 1) {
        xhat_init <- xhat0
        PNinit <- Phat0
        time.steps <- c(t0, times[1])
      } else {
        xhat_init <- xhat_kk[, i - 1]
        PNinit <- P_kk[, , i - 1]
        time.steps <- c(times[i - 1], times[i])
      }

      xhat_init["C"] <- 0
      xhat_init["Hnew"] <- 0
      xhat_init["Drep"] <- 0
      
      PNinit[, 4] <- PNinit[4,] <- 0
      PNinit[, 5] <- PNinit[5,] <- 0
      PNinit[, 8] <- PNinit[8,] <- 0

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
        beta_t = exp(logbeta[i]),
        time.steps = time.steps
      )
      xhat_kkmo[, i] <- XP$xhat
      P_kkmo[, , i] <- XP$PN
      
      ytilde_k[, i] <- matrix(z[i, ], ncol = 1) - 
        H %*% xhat_kkmo[, i, drop = FALSE]     
      S[, , i] <- H %*% P_kkmo[, , i] %*% t(H) + R
      
      for (j in 1:dobs){
        if (is.na(z[i,j])){
          zscore <- 0
          rdiagadj[j,i] <- 0
        } else {
          sd <- sqrt(S[j,j,i])
          zscore <- ytilde_k[j,i] / sd 
          if (abs(zscore) > maxzscore){
            adjzscore <- maxzscore / (1 + abs(zscore) - maxzscore)
            newsd <- abs(ytilde_k[j,i]) / adjzscore
            rdiagadj[j,i] <- (newsd) ^ 2 - sd ^ 2
          } else {
            rdiagadj[j,i] <- 0
          }
        }
        S[j,j,i] <- S[j,j,i] + rdiagadj[j,i]
      }
      K[, , i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, , i])
      desel <- is_z_na[i, ]
      K[, desel, i] <- 0

      xhat_kk[, i] <-
        xhat_kkmo[, i, drop = FALSE] +
        K[, !desel, i] %*% ytilde_k[!desel, i, drop = FALSE]
      
      xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
      P_kk[, , i] <-
        (diag(dstate) - K[, , i] %*% H) %*% P_kkmo[, , i]
      ytilde_kk[, i] <- matrix(z[i, ], ncol = 1) - 
        H %*% xhat_kk[, i, drop = FALSE]
    }
    
    nll <- 0
    for (i in seq(1, T)){
      sel <- !is_z_na[i, ]
      nll <- nll + 
        t(ytilde_k[sel, i]) %*% solve(S[sel, sel, i]) %*% ytilde_k[sel, i] + 
        log(det(S[,,i][sel, sel, drop = FALSE])) + dobs * log(2 * pi)
    }
    nll <- 0.5 * nll - sum(dnorm(bpars[-T], sd = betasd, log = TRUE))

    if (!just_nll) {
      if (!is.null(fets)) {
        nsimdays <- nrow(fets)
        sim_means <- array(NA_real_, dim = c(dobs, nsimdays, nsim))
        sim_cov <- array(NA_real_, dim = c(dobs, dobs, nsimdays, nsim))
        for (j in seq_len(nsim)) {
          logbeta_fet <- numeric(nsimdays)
          logbeta_fet[1] <-
            log(gamma) + a * (bpars[T] - log(gamma)) +  
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
            beta_t = exp(logbeta_fet[1]),
            time.steps = c(times[T], fets$target_end_times[1])
          )
          sim_means[, 1, j] <- H %*% XP$xhat
          sim_cov[, , 1, j] <- H %*% XP$PN %*% t(H) + R
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
                beta_t = exp(logbeta_fet[i + 1]),
                time.steps = c(fets$target_end_times[i], 
                               fets$target_end_times[i + 1])
              )
            sim_means[, i + 1, j] <- H %*% XP$xhat
            sim_cov[, , i + 1, j] <-
              H %*% XP$PN %*% t(H) + R
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
        rdiagadj = rdiagadj
      )
    } else {
      nll
    }
  }

moving_average <- function(x, n = 7) {
  stats::filter(x, rep(1 / n, n), sides = 1)
}

calc_kf_nll <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  pvar <- c(p$logE0, p$logH0, p$logtauc, p$logtauh, p$logtaud, p$logchp, p$loghfp, p$bpars)
  JuliaCall::julia_assign("pvar", pvar)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("N", p$N)
  JuliaCall::julia_assign("ρ", p$rho1)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("γd", p$gamma_d)
  JuliaCall::julia_assign("γh", p$gamma_h)
  JuliaCall::julia_assign("z", data.matrix(y))
  JuliaCall::julia_assign("a", a)
  JuliaCall::julia_assign("betasd", betasd)
  nll <- JuliaCall::julia_eval(paste0(
    "InfectionKalman.obj",
    "(pvar, z; ρ = ρ, N = N, η = η, γ = γ, γd = γd, γh = γh, a = a, betasd = betasd)"
  ))
  nll
}

calc_kf_grad <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  pvar <- c(p$logE0, p$logH0, p$logtauc, p$logtauh, p$logtaud, p$logchp, p$loghfp, p$bpars)
  JuliaCall::julia_assign("pvar", pvar)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("N", p$N)
  JuliaCall::julia_assign("ρ", p$rho1)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("γd", p$gamma_d)
  JuliaCall::julia_assign("γh", p$gamma_h)
  JuliaCall::julia_assign("z", data.matrix(y))
  JuliaCall::julia_assign("a", a)
  JuliaCall::julia_assign("betasd", betasd)
  g <- JuliaCall::julia_eval(paste0(
    "InfectionKalman.grad",
    "(pvar, z; ρ = ρ, N = N, η = η, γ = γ, γd = γd, γh = γh, a = a, betasd = betasd)"
  ))
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

kf_nll_details <- function(w, x, y, betasd, a, pm, fet) {
  p <- pm(x, w)
  nll <- kfnll(
    bpars = p$bpars,
    logE0 = p$logE0,
    logH0 = p$logH0,
    rho1 = p$rho1,
    logtauc = p$logtauc,
    logtauh = p$logtauh,
    logtaud = p$logtaud,
    logchp = p$logchp,
    loghfp = p$loghfp,
    eta = p$eta,
    gamma = p$gamma,
    gamma_d = p$gamma_d,
    gamma_h = p$gamma_h,
    N = p$N,
    z = y,
    t0 = p$t0,
    times = p$times,
    fet = fet,
    just_nll = FALSE,
    fet_zero_cases_deaths = "weekly",
    nsim = 20,
    betasd = betasd,
    a = a
  )
  nll
}

write_forecasts <- function(fits, fet, agrid, betagrid) {
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
      case_inds <- which(fet$target_wday == 7)
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

wind <- obs_data %>% filter(target_end_date >= "2020-03-01")
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

wind <- tail(obs_data, n = wsize)
y <- wind %>% select(cases, hospitalizations, deaths)
x <- matrix(wind$time, ncol = 1)

tau_cases_init <- var(y$cases, na.rm = TRUE)
tau_hosp_init <- var(y$hospitalizations, na.rm = TRUE)
tau_deaths_init <- var(y$deaths, na.rm = TRUE)

E0init <- ((mean(y$cases) / wfixed["rho1"]) * (365.25 / wfixed["eta"])) %>% unname()
H0init <- (mean(y$hospitalizations, na.rm = TRUE) * 365.25 / wfixed["gamma_h"]) %>% unname()

chp_init <- sum(y$hospitalizations, na.rm = TRUE) / sum(y$cases, na.rm = TRUE)
hfp_init <- sum(y$deaths, na.rm = TRUE) / sum(y$hospitalizations, na.rm = TRUE)

binit <- c(rep(0, wsize - 1), log(wfixed["gamma"]))
names(binit) <- paste0("b", seq_len(wsize))
winit <- c(
  logE0 = log(E0init),
  logH0 = log(H0init),
  logtauc = log(tau_cases_init),
  logtauh = log(tau_hosp_init),
  logtaud = log(tau_deaths_init),
  logchp = log(chp_init),
  loghfp = log(hfp_init),
  binit
)

param_map <- function(x, w, fixed = wfixed){
  ret <- list()
  ret$bpars <- w[seq(8, length(w))]
  ret$logE0 <- w[1]
  ret$logH0 <- w[2]
  ret$logtauc <- w[3]
  ret$logtauh <- w[4]
  ret$logtaud <- w[5]
  ret$logchp <- w[6]
  ret$loghfp <- w[7]
  ret$times <- x[, 1]
  ret$rho1 <- fixed["rho1"]
  ret$eta <- fixed["eta"]
  ret$gamma <- fixed["gamma"]
  ret$gamma_h <- fixed["gamma_h"]
  ret$gamma_d <- fixed["gamma_d"]
  ret$N <- fixed["N"]
  ret$t0 <- fixed["t0"]
  ret
}

tictoc::tic("optimization")

fit_over_betagrid <- function(a, betagrid) {
  fits <- list()
  fits[[1]] <- lbfgs::lbfgs(
    calc_kf_nll,
    calc_kf_grad,
    x = x,
    betasd = betagrid[1],
    epsilon = 1e-3,
    a = a,
    y = y,
    pm = param_map,
    winit,
    invisible = 1
  )
  for (i in seq(2, length(betagrid))) {
    fits[[i]] <- lbfgs::lbfgs(
      calc_kf_nll,
      calc_kf_grad,
      x = x,
      betasd = betagrid[i],
      epsilon = 1e-3,
      a = a,
      y = y,
      pm = param_map,
      fits[[i - 1]]$par,
      invisible = 1
    )
  }
  return(fits)
}

betagrid <- seq(0.001, 0.1, len = 10)
agrid <- c(0.94, 0.95)

fits <- map(agrid, fit_over_betagrid, betagrid = betagrid)
tictoc::toc()

ti <- max(wind$target_end_date) + lubridate::ddays(1)
tf <- lubridate::ymd(forecast_date) + lubridate::ddays(28)
target_end_dates <- seq(from = ti, to = tf, by = "day")

target_end_times <- lubridate::decimal_date(target_end_dates)
target_wday <- lubridate::wday(target_end_dates)

fet <- tibble(target_end_times, target_wday, target_end_dates)
write_forecasts(fits, fet, agrid, betagrid)

q("no")
# View diagnostics

fitind1 <- 1
fitind2 <- 2
dets <- kf_nll_details(winit, x = x, y = y, param_map, 
                       betasd = betagrid[fitind2], a = agrid[fitind1],
                       fet)

nll <- calc_kf_nll(winit, x = x, y = y, param_map, 
                       betasd = betagrid[fitind2], a = agrid[fitind1])


fitind1 <- 1
fitind2 <- 2
dets <- kf_nll_details(fits[[fitind1]][[fitind2]]$par, x = x, y = y, param_map, 
                       betasd = betagrid[fitind2], a = agrid[fitind1],
                       fet)
par(mfrow = c(1,1))
qqnorm(dets$ytilde_k[1,] / sqrt(dets$S[1,1,]))
abline(0, 1)
qqnorm(dets$ytilde_k[2,] / sqrt(dets$S[2,2,]))
abline(0, 1)
qqnorm(dets$ytilde_k[3,] / sqrt(dets$S[3,3,]))
abline(0,1 )


plot(x[,1], y$cases, xlab = "Time", ylab = "Cases")
pred_cases <- dets$xhat_kkmo["C",] * wfixed["rho1"] + dets$xhat_kkmo["Hnew",]
se_cases <- sqrt(dets$S[1,1,])
lines(x[,1], se_cases * 2 + pred_cases, col = "grey")
lines(x[,1], pred_cases)
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
