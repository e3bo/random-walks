#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

JuliaCall::julia_setup()
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
  } else{
    stop("not implemented")
    probs <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
    if (target_type == "hosps") {
      maxn <- 131
    } else{
      maxn <- 20
    }
  }
  n <- nrow(means)
  stopifnot(nrow(vars) == n)
  stopifnot(n <= maxn)
  wk_ahead <- seq_len(n)
  targets <- paste0(wk_ahead, targ_string)
  dte <- lubridate::ymd(fdt)
  fdt_wd <- lubridate::wday(dte) # 1 = Sunday, 7 = Sat.
  fdt_sun <- lubridate::ymd(dte) - (fdt_wd - 1)
  take_back_step <- fdt_wd <= 2
  if (take_back_step) {
    week0_sun <- fdt_sun - 7
  } else {
    week0_sun <- fdt_sun
  }
  t1 <- expand_grid(quantile = probs, h = wk_ahead) %>%
    mutate(target = paste0(h, targ_string),
           target_end_date = week0_sun + 6 + h * 7) %>%
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


iterate_f_and_P <- function(xhat, PN, eta, gamma, N, beta_t, time.steps){
  P <- PN / N
  dt <- diff(time.steps)
  vf <-  with(as.list(c(xhat, beta_t)), {
    F <- rbind(
      c(-beta_t * I / N,    0,-beta_t * S / N, 0),
      c(beta_t * I / N,-eta,  beta_t * S / N, 0),
      c(0,  eta,-gamma, 0),
      c(0,    0,             gamma, 0)
    )
    
    f <-
      c(0, beta_t * (S / N) * (I / N), eta * E / N, gamma * I / N)
    Q <- rbind(c(f[1] + f[2],-f[2],           0,     0),
               c(-f[2], f[2] + f[3],-f[3],     0),
               c(0,-f[3], f[3] + f[4],-f[4]),
               c(0,           0,-f[4],  f[4]))
    
    dS <- (-beta_t * S * I / N)
    dE <- (beta_t * S * I / N - eta * E)
    dI <- (eta * E -  gamma * I)
    dC <- (gamma * I)
    dP <-  F %*% P + P %*% t(F) + Q
    
    list(vf = c(dS, dE, dI, dC),
         dP = dP)
  })
  xhat_new <- xhat + vf$vf * dt
  xhat_new[xhat_new < 0] <- 0
  P_new <- P + vf$dP * dt
  names(xhat_new) <- c("S", "E", "I", "C")
  PN_new <- P_new * N
  list(xhat = xhat_new, PN = PN_new)
}

kfnll <-
  function(bpars,
           logE0,
           rho1,
           logtau,
           eta,
           gamma,
           N,
           z,
           t0,
           times,
           Phat0 = diag(c(1, 1, 1, 0)),
           fets = NULL,
           fet_zero_cases = "daily",
           nsim,
           a = .98,
           betasd = 1,
           params_Sigma = matrix(0, 2, 2),
           just_nll = TRUE) {
    E0 = exp(logE0)
    I0 = E0 * eta / gamma
    xhat0 = c(N - E0 - I0, E0, I0, 0)
    names(xhat0) <- c("S", "E", "I", "C")
    
    T <- length(z)
    stopifnot(T > 0)
    
    logbeta <-
      ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
    
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))
    
    H <- matrix(c(0, 0, 0, rho1), ncol = 4)
    
    logbeta[T] <- bpars[T]
    for (i in seq(T - 1, 1)){
      logbeta[i] <- (logbeta[i + 1] - log (gamma) - bpars[i]) / a + log(gamma)
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
      R <- exp(logtau)
      
      xhat_init["C"] <- 0
      PNinit[, 4] <- PNinit[4,] <- 0
      
      XP <- iterate_f_and_P(
        xhat_init,
        PN = PNinit,
        eta = eta,
        gamma = gamma,
        N = N,
        beta_t = exp(logbeta[i]),
        time.steps = time.steps
      )
      xhat_kkmo[, i] <- XP$xhat
      P_kkmo[, , i] <- XP$PN
      
      S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
      K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
      ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
      xhat_kk[, i] <-
        xhat_kkmo[, i, drop = FALSE] +
        K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
      xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
      P_kk[, , i] <-
        (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
      ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
    }
    
    nll <-
      0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi)) -
      sum(dnorm(bpars[-T], sd = betasd, log = TRUE))
    if (!just_nll) {
      if (!is.null(fets)) {
        params_mu <- c(logtau, bpars[T])
        sim_means <- sim_cov <- matrix(NA, nrow = (nrow(fets)), ncol = nsim)
        for (j in seq_len(nsim)) {
          logbeta_fet <- numeric(nrow(fets))
          params_samp <- MASS::mvrnorm(n = 1, mu = params_mu, Sigma = params_Sigma)
          logtau_samp <- params_mu[1] # some estimates of variance are unrealistic
          logbetaT_samp <- params_samp[2]
          logbeta_fet[1] <-
            log(gamma) + a * (logbetaT_samp - log(gamma)) +  
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
          if (fet_zero_cases == "daily" ||
              fets$target_wday[1] == 1) {
            xhat_init["C"] <- 0
            PNinit[, 4] <- PNinit[4,] <- 0
          }
          XP <- iterate_f_and_P(
            xhat_init,
            PN = PNinit,
            eta = eta,
            gamma = gamma,
            N = N,
            beta_t = exp(logbeta_fet[1]),
            time.steps = c(times[T], fets$target_end_times[1])
          )
          sim_means[1, j] <- H %*% XP$xhat
          sim_cov[1, j] <- H %*% XP$PN %*% t(H) + exp(logtau)
          for (i in seq_along(fets$target_end_times[-1])) {
            xhat_init <- XP$xhat
            PNinit <- XP$PN
            if (fet_zero_cases == "daily" ||
                fets$target_wday[i + 1] == 1) {
              xhat_init["C"] <- 0
              PNinit[, 4] <- PNinit[4,] <- 0
            }
            XP <-
              iterate_f_and_P(
                xhat_init,
                PN = PNinit,
                eta = eta,
                gamma = gamma,
                N = N,
                beta_t = exp(logbeta_fet[i + 1]),
                time.steps = c(fets$target_end_times[i], 
                               fets$target_end_times[i + 1])
              )
            sim_means[i + 1, j] <- H %*% XP$xhat
            sim_cov[i + 1, j] <-
              H %*% XP$PN %*% t(H) + exp(logtau_samp)
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
        logbeta = logbeta
      )
    } else {
      nll
    }
  }

moving_average <- function(x, n = 7) {
  stats::filter(x, rep(1 / n, n), sides = 1)
}

calc_kf_nll <- function(w, x, y, betasd, pm) {
  p <- pm(x, w)
  pvar <- c(p$logE0, p$logtau, p$bpars)
  JuliaCall::julia_assign("pvar", pvar)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("N", p$N)
  JuliaCall::julia_assign("ρ", p$rho1)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("z", map(y, list))
  JuliaCall::julia_assign("a", p$a)
  JuliaCall::julia_assign("betasd", betasd)
  nll <- JuliaCall::julia_eval(paste0(
    "InfectionKalman.obj",
    "(pvar, z; ρ = ρ, N = N, η = η, γ = γ, a = a, betasd = betasd)"
  ))
  nll
}

calc_kf_grad <- function(w, x, y, betasd, pm) {
  p <- pm(x, w)
  pvar <- c(p$logE0, p$logtau, p$bpars)
  JuliaCall::julia_assign("pvar", pvar)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("N", p$N)
  JuliaCall::julia_assign("ρ", p$rho1)
  JuliaCall::julia_assign("η", p$eta)
  JuliaCall::julia_assign("γ", p$gamma)
  JuliaCall::julia_assign("z", map(y, list))
  JuliaCall::julia_assign("a", p$a)
  JuliaCall::julia_assign("betasd", betasd)
  g <- JuliaCall::julia_eval(paste0(
    "InfectionKalman.grad",
    "(pvar, z; ρ = ρ, N = N, η = η, γ = γ, a = a, betasd = betasd)"
  ))
  g
}

calc_kf_hess <- function(w, x, y, betasd, pm) {
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
  JuliaCall::julia_assign("z", map(y, list))
  JuliaCall::julia_assign("a", p$a)
  JuliaCall::julia_assign("betasd", betasd)
  g <- JuliaCall::julia_eval(paste0(
    "InfectionKalman.hess",
    "(logE0, logtau, bpars, z; ρ = ρ, N = N, η = η, γ = γ, a = a, betasd = betasd)"
  ))
  g
}

kf_nll_details <- function(w, x, y, betasd, pm, fet, params_Sigma) {
  p <- pm(x, w)
  nll <- kfnll(
    bpars = p$bpars,
    logE0 = p$logE0,
    rho1 = p$rho1,
    logtau = p$logtau,
    eta = p$eta,
    gamma = p$gamma,
    N = p$N,
    z = y,
    t0 = p$t0,
    times = p$times,
    fet = fet,
    just_nll = FALSE,
    fet_zero_cases = "weekly",
    nsim = 20,
    betasd = betasd,
    params_Sigma = params_Sigma,
    a = p$a
  )
  nll
}

write_forecasts <- function(fits, fet, betagrid) {
  nlam <- length(betagrid)
  for (penind in 1:nlam) {
    wfit <- fits[[penind]]$par
    betasd <- betagrid[penind]
    hess <- calc_kf_hess(wfit, x, y, betasd = betasd, param_map)
    params_Sigma <- solve(hess)
    dets <-
      kf_nll_details(wfit, x, y, param_map, 
                     betasd = betasd, fet = fet, params_Sigma = params_Sigma)
    inds <- which(fet$target_wday == 7)
    fcst <- create_forecast_df(
      means = dets$sim_means[inds, ],
      vars = dets$sim_cov[inds, ],
      location = forecast_loc
    )
    
    stopifnot(setequal(fet$target_end_dates[inds], 
                       fcst$target_end_date %>% unique()))
    lambda <- 1 / betagrid[penind]
    fcst_path <-
      file.path(
        "forecasts",
        paste0(
          forecast_date,
          "-fips",
          forecast_loc,
          "-lambda",
          sprintf("%06.2f", lambda),
          "-CEID-InfectionKalman.csv"
        )
      )
    if (!dir.exists("forecasts"))
      dir.create("forecasts")
    write_csv(x = fcst, path = fcst_path)
  }
}
## main script

forecast_date <- Sys.getenv("fdt", unset = "2020-10-12")
forecast_loc <- Sys.getenv("loc", unset = "36")
hopdir <- file.path("hopkins", forecast_date)
tictoc::tic("data loading")
tdat <- load_hopkins(hopdir, weekly = FALSE) 
tictoc::toc()

ltdat <- tdat %>% filter(location == forecast_loc) %>% 
  filter(target_type == "day ahead inc case")

ltdat2 <- ltdat %>% mutate(time = lubridate::decimal_date(target_end_date))
case_data <- ltdat2 %>% ungroup() %>% 
  mutate(wday = lubridate::wday(target_end_date)) %>%
  select(target_end_date, time, wday, value) %>% 
  rename(reports = value)

target_end_dates <- max(case_data$target_end_date) + 
  lubridate::ddays(1) * seq(1, 28)
target_end_times <- lubridate::decimal_date(target_end_dates)
target_wday <- lubridate::wday(target_end_dates)

case_data$smooth <- moving_average(case_data$reports)

wsize <- 60
gamma <- 365.25/9
tau_init <- case_data$reports %>% tail(n = wsize) %>% var()

pvar_df <- tribble(
  ~par, ~init, ~lower, ~upper,
  "logE0", log(1e4), log(10), log(1e5),
  "logtau", log(tau_init), log(tau_init * 1e-8), log(tau_init * 10)
) %>% 
  bind_rows(tibble(par = paste0("b", seq_len(wsize - 1)),
                   init = 0,
                   lower = -10,
                   upper = 10)) %>%
  bind_rows(tibble(par = paste0("b", wsize),
                   init = log(gamma), 
                   lower = log(0.1 * gamma),
                   upper = log(4 * gamma)))

winit <- pvar_df$init
names(winit) <- pvar_df$par

N <- covidHubUtils::hub_locations %>% filter(fips == forecast_loc) %>% 
  pull(population)

wfixed <- c(
  N = N,
  rho1 = 0.4,
  gamma = 365.25 / 9, 
  eta = 365.25 / 4,
  t0 = rev(case_data$time)[wsize + 1],
  a = 0.95,
  betasd = 1
)

param_map <- function(x, w, fixed = wfixed){
  ret <- list()
  ret$bpars <- w[seq(3, length(w))]
  ret$logE0 <- w[1]
  ret$rho1 <- fixed["rho1"]
  ret$logtau <- w[2]
  ret$times <- x[, 1]
  ret$eta <- fixed["eta"]
  ret$gamma <- fixed["gamma"]
  ret$N <- fixed["N"]
  ret$t0 <- fixed["t0"]
  ret$a <- fixed["a"]
  ret
}

wind <- tail(case_data, n = wsize)
y <- wind$smooth
x <- matrix(wind$time, ncol = 1)

tictoc::tic("optimization")
fits <- list()
betagrid <- seq(0.001, 0.1, len = 10)
fits[[1]] <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x,
  betasd = betagrid[1],
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
    y = y,
    pm = param_map,
    fits[[i - 1]]$par,
    invisible = 1
  )
}
tictoc::toc()

fet <- tibble(target_end_times, target_wday, target_end_dates)
write_forecasts(fits, fet, betagrid)

q("no")
# View diagnostics

fitind <- 1
dets <- kf_nll_details(fits[[fitind]]$par, x = x, y = y, param_map, 
                       betasd = betagrid[fitind], fet)
par(mfrow = c(1,1))
qqnorm(dets$ytilde_k / sqrt(dets$S))
abline(0, 1)

tgrid <- tail(case_data$time, n = wsize)
plot(tgrid, tail(case_data$smooth, n = wsize), xlab = "Time", ylab = "Cases")
lines(tgrid, dets$xhat_kkmo["C",] * wfixed["rho1"])

inds <- which(fet$target_wday == 7)

fcst <- create_forecast_df(means = dets$sim_means[inds,],
                           vars = dets$sim_cov[inds,],
                           location = forecast_loc)

fcst %>% 
  ggplot(aes(x = target_end_date, y = as.numeric(value), color = quantile)) + 
  geom_line() + geom_point()