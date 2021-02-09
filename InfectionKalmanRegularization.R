#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")
source("regularization.R")

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
  bind_rows(t1, t2) %>% select(forecast_date,
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
           xhat0,
           rho1, 
           tau,
           eta,
           gamma,
           N,
           z,
           t0,
           times, 
           Phat0 = diag(c(1, 1, 1, 0)),
           fets = NULL,
           fet_zero_cases = "daily",
           lambda,
           nsim,
           just_nll = TRUE) {
    
    T <- length(z)
    stopifnot(T > 0)
    
    ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))
    
    H <-matrix(c(0, 0, 0, rho1), ncol = 4)
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
      R <- tau
      
      xhat_init["C"] <- 0
      PNinit[, 4] <- PNinit[4, ] <- 0
      
      XP <- iterate_f_and_P(
        xhat_init,
        PN = PNinit,
        eta = eta,
        gamma = gamma,
        N = N,
        beta_t = bpars[i],
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
      0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi)) 
    if (!just_nll) {
      if (!is.null(fets)) {
        sim_means <-
          sim_cov <- matrix(NA, nrow = (nrow(fets)), ncol = nsim)
        for (j in seq_len(nsim)) {
          bpars_fet <- numeric(nrow(fets))
          bpars_fet[1] <- bpars[T] * exp(sample(c(-1, 1), 1) * rexp(n = 1, rate = lambda))
          if (length(bpars_fet) > 1){
            for (jj in seq(2, length(bpars_fet))){
              bpars_fet[jj] <- bpars_fet[jj - 1] * exp(sample(c(-1, 1), 1) * rexp(n = 1, rate = lambda))
            }
          }
          xhat_init <- xhat_kk[, T]
          PNinit <- P_kk[, , T]
          if (fet_zero_cases == "daily" || fets$target_wday[1] == 1){
            xhat_init["C"] <- 0
            PNinit[, 4] <- PNinit[4, ] <- 0
          }
          XP <- iterate_f_and_P(
            xhat_init,
            PN = PNinit,
            eta = eta,
            gamma = gamma,
            N = N,
            beta_t = bpars_fet[1],
            time.steps = c(times[T], fets$target_end_times[1])
          )
          sim_means[1, j] <- H %*% XP$xhat
          sim_cov[1, j] <- H %*% XP$PN %*% t(H) + tau
          for (i in seq_along(fets$target_end_times[-1])) {
            xhat_init <- XP$xhat
            PNinit <- XP$PN
            if (fet_zero_cases == "daily" || fets$target_wday[i + 1] == 1){
              xhat_init["C"] <- 0
              PNinit[, 4] <- PNinit[4, ] <- 0
            }
            XP <-
              iterate_f_and_P(
                xhat_init,
                PN = PNinit,
                eta = eta,
                gamma = gamma,
                N = N,
                beta_t = bpars_fet[i + 1],
                time.steps = c(fets$target_end_times[i], fets$target_end_times[i + 1])
              )
            sim_means[i + 1, j] <- H %*% XP$xhat
            sim_cov[i + 1, j] <-
              H %*% XP$PN %*% t(H) + tau
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
        sim_cov = sim_cov
      )
    } else {
      nll
    }
  }

moving_average <- function(x, n = 7) {
  stats::filter(x, rep(1 / n, n), sides = 1)
}

calc_kf_nll <- function(w, x, y, pm) {
  p <- pm(x, w)
  nll <- kfnll(
    bpars = p$bpars,
    xhat0 = p$xhat0,
    rho1 = p$rho1,
    tau = p$tau,
    eta = p$eta,
    gamma = p$gamma,
    N = p$N,
    z = y,
    t0 = p$t0,
    times = p$times,
    just_nll = TRUE
  )
  nll
}

kf_nll_details <- function(w, x, y, pm, lambda, fet) {
  p <- pm(x, w)
  nll <- kfnll(
    bpars = p$bpars,
    xhat0 = p$xhat0,
    rho1 = p$rho1,
    tau = p$tau,
    eta = p$eta,
    gamma = p$gamma,
    N = p$N,
    z = y,
    t0 = p$t0,
    times = p$times,
    fet = fet,
    lambda = lambda,
    just_nll = FALSE,
    fet_zero_cases = "weekly",
    nsim = 10
  )
  nll
}

write_forecasts <- function(gpnet, fet) {
  nlam <- length(gpnet$lambda)
  for (penind in 1:nlam) {
    wfit <- c(rpath$a0[, penind], rpath$beta[, penind])
    names(wfit)[4:length(wfit)] <- paste0("b", 2:wsize)
    dets <-
      kf_nll_details(wfit, x, y, param_map, rpath$lambda[penind], fet)
    inds <- which(fet$target_wday == 7)
    fcst <- create_forecast_df(
      means = dets$sim_means[inds, ],
      vars = dets$sim_cov[inds, ],
      location = forecast_loc
    )
    
    stopifnot(setequal(fet$target_end_dates[inds], 
                       fcst$target_end_date %>% unique()))
    lambda <- rpath$lambda[penind]
    fcst_path <-
      file.path(
        "forecasts",
        paste0(
          forecast_date,
          "-fips",
          forecast_loc,
          "-lambda",
          lambda,
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

target_end_dates <- max(case_data$target_end_date) + lubridate::ddays(1) * seq(1, 28)
target_end_times <- lubridate::decimal_date(target_end_dates)
target_wday <- lubridate::wday(target_end_dates)

case_data$smooth <- moving_average(case_data$reports)

wsize <- 60
gamma <- 365.25/9
tau_init <- case_data$reports %>% tail(n = wsize) %>% var()

pvar_df <- tribble(
  ~par, ~init, ~lower, ~upper,
  "logI_0", log(1e4), log(10), log(1e5),
  "logtau", log(tau_init), log(tau_init * 1e-8), log(tau_init * 10)
) %>% 
  bind_rows(tibble(par = "b1",
                   init = log(gamma), 
                   lower = log(0.1 * gamma),
                   upper = log(4 * gamma))) %>%
  bind_rows(tibble(par = paste0("b", 1 + seq_len(wsize - 1)),
                   init = 0,
                   lower = -10,
                   upper = 10))

winit <- pvar_df$init
names(winit) <- pvar_df$par

wfixed <- c(
  N = 20e6,
  rho1 = 0.4,
  gamma = 365.25 / 9, 
  eta = 365.25 / 4,
  t0 = rev(case_data$time)[wsize + 1]
)

param_map <- function(x, w, fixed = wfixed){
  ret <- list()
  is_spline_par <- grepl("^b[0-9]+$", names(w))
  tmp <- exp(w[is_spline_par]) # tranform from log scale
  bpars_diff <- tmp[order(as.integer(str_remove(names(tmp), "^b")))]
  stopifnot(length(bpars_diff) == nrow(x))
  ret$bpars <- cumprod(bpars_diff)
  
  I_0 <- exp(w["logI_0"])
  E_0 <- I_0 * fixed["gamma"] / fixed["eta"]
  S_0 <- fixed["N"] - E_0 - I_0
  
  ret$xhat0 = c(S_0, E_0, I_0, 0)
  names(ret$xhat0) <- c("S", "E", "I", "C")
  
  ret$rho1 <- fixed["rho1"]
  ret$tau <- exp(w["logtau"])
  ret$times <- x[, 1]
  ret$eta <- fixed["eta"]
  ret$gamma <- fixed["gamma"]
  ret$N <- fixed["N"]
  ret$t0 <- fixed["t0"]
  ret
}

wind <- tail(case_data, n = wsize)
y <- wind$smooth
x <- matrix(wind$time, ncol = 1)

#calc_kf_nll(winit, x, y, param_map)
pen_factor <- rep(1, length(winit))
pen_factor[1:3] <- 0

max_lambda <- 2.2
log_dec <- 1.5
len_lambda <- 10
lambda <-
  10 ^ (seq(max_lambda, max_lambda - log_dec, length.out = len_lambda)) %>% 
  round(2)

tictoc::tic("optimization")
rpath <-
  get_gpnet(
    x = x,
    y = y,
    calc_convex_nll = calc_kf_nll,
    param_map = param_map,
    lambda = lambda[1:8],
    penalty.factor = pen_factor,
    winit = winit,
    make_log = TRUE
  )
tictoc::toc()

fet <- tibble(target_end_times, target_wday, target_end_dates)

write_forecasts(rpath, fet)

q("no")
# View diagnostics

penind <- 10
wfit <- c(rpath$a0[,penind], rpath$beta[,penind])
names(wfit)[4:length(wfit)] <- paste0("b", 2:wsize)


fet <- tibble(target_end_times, target_wday, target_end_dates)


dets <- kf_nll_details(wfit, x, y, param_map, rpath$lambda[penind], fet)
par(mfrow = c(1,1))
qqnorm(dets$ytilde_k / sqrt(dets$S))
abline(0, 1)

par(mfrow = c(4, 1))
tgrid <- tail(case_data$time, n = wsize)

plot(tgrid, tail(case_data$smooth, n = wsize), xlab = "Time", ylab = "Cases")
lines(tgrid, dets$xhat_kkmo["C",] * wfixed["rho1"])

plot(tgrid, dets$S, log = "y", xlab = "Time", ylab = "Variance in smoother")
plot(tgrid, dets$ytilde_k, xlab = "Time", ylab = "Residual in process 1-ahead prediction")
plot(tgrid, dets$ytilde_k / sqrt(dets$S), xlab = "Time", ylab = "Standardized residual")



inds <- which(fet$target_wday == 7)

fcst <- create_forecast_df(means = dets$sim_means[inds,],
                           vars = dets$sim_cov[inds,],
                           location = forecast_loc)


fcst %>% ggplot(aes(x = target_end_date, y = as.numeric(value), color = quantile)) + geom_line() + geom_point()

fcst_path <- file.path("forecasts", paste0(forecast_date, "-", forecast_loc, "-CEID-InfectionKalman.csv"))
if(!dir.exists("forecasts")) dir.create("forecasts")
write_csv(x = fcst, path = fcst_path)

