#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")
source("regularization.R")


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

moving_average <- function(x, n = 7) {
  stats::filter(x, rep(1 / n, n), sides = 1)
}

case_data$smooth <- moving_average(case_data$reports)

wsize <- 60
gamma <- 365.25/9
tau_init <- case_data$reports %>% tail(n = wsize) %>% var()

pvar_df <- tribble(
  ~par, ~init, ~lower, ~upper,
  "I_0", log(1e4), log(10), log(1e5),
  "tau", log(tau_init), log(tau_init * 1e-8), log(tau_init * 10)
) %>% 
  bind_rows(tibble(par = "b1",
                   init = log(gamma), 
                   lower = log(0.1 * gamma),
                   upper = log(4 * gamma))) %>%
  bind_rows(tibble(par = paste0("b", 1 + seq_len(wsize - 1)),
                   init = 0,
                   lower = -10,
                   upper = 10))

the_t0 <- rev(case_data$time)[wsize + 1]

iterate_f_and_P <- function(xhat, PN, pvec, beta_t, time.steps){
  P <- PN / pvec["N"]
  dt <- diff(time.steps)
  vf <-  with(as.list(c(pvec, xhat, beta_t)), {
    eta <- pvec["eta"]
    gamma <- pvec["gamma"]
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
    dI <- (iota + eta * E -  gamma * I)
    dC <- (gamma * I)
    dP <-  F %*% P + P %*% t(F) + Q
    
    list(vf = c(dS, dE, dI, dC),
         dP = dP)
  })
  xhat_new <- xhat + vf$vf * dt
  xhat_new[xhat_new < 0] <- 0
  P_new <- P + vf$dP * dt
  names(xhat_new) <- c("S", "E", "I", "C")
  PN_new <- P_new * pvec["N"]
  list(xhat = xhat_new, PN = PN_new)
}

kfnll <-
  function(pvar,
           pfixed,
           cdata,
           t0,
           Phat0 = diag(c(1, 1, 1, 0)),
           just_nll = TRUE,
           nsim = 10,
           fets = NULL,
           Rzzero = 1e6,
           fet_zero_cases = "daily") {
    p <- c(exp(pvar), pfixed)
    
    is_spline_par <- grepl("^b[0-9]+$", names(p))
    tmp <- p[is_spline_par]
    bpars_diff <- tmp[order(as.integer(str_remove(names(tmp), "^b")))]
    bpars_diff[-1] <- log(bpars_diff[-1]) ## backtransform step sizes so that they can be negative
    stopifnot(length(bpars_diff) == nrow(cdata))
    bpars <- cumsum(bpars_diff)
    
    is_wday_par <- grepl("rho[1-7]", names(p))
    tmp <- p[is_wday_par]
    wpars <- tmp[order(names(tmp))]
    stopifnot(length(wpars) == 1)
    
    xhat0 = c(p[c("S_0", "E_0", "I_0")], 0)
    names(xhat0) <- c("S", "E", "I", "C")
    
    T <- nrow(cdata)
    stopifnot(T > 0)
    z <- cdata$smooth
    times <- cdata$time
    wday <- cdata$wday
    
    ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))
    
    H <-matrix(c(0, 0, 0, p["rho1"]), ncol = 4)
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
      R <- p["tau"]
      if (z[i] < 1){
        R <- Rzzero * p["tau"]
      }
      if (TRUE){
        xhat_init["C"] <- 0
        PNinit[, 4] <- PNinit[4, ] <- 0
      }
      XP <- iterate_f_and_P(
        xhat_init,
        PN = PNinit,
        pvec = p,
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
          bpars_fet[1] <- p["gamma"] + (bpars[T] - p["gamma"] )* p["a"] + rnorm(n = 1, sd = p["betasd"])
          rand <- runif(1)
          if (rand > 0.66){
            bpars_fet[1] <- bpars_fet[1] * 1.5
          } else if (rand < 0.33){
            bpars_fet[1] <- bpars_fet[1] * 0.5
          }
          if (length(bpars_fet) > 1){
            for (jj in seq(2, length(bpars_fet))){
              bpars_fet[jj] <- p["gamma"] + (bpars_fet[jj - 1]  - p["gamma"]) * p["a"] + rnorm(n = 1, sd = p["betasd"])
            }
          }
          bpars_fet[bpars_fet < 0] <- 0
          xhat_init <- xhat_kk[, T]
          PNinit <- P_kk[, , T]
          if (fet_zero_cases == "daily" || fets$target_wday[1] == 1){
            xhat_init["C"] <- 0
            PNinit[, 4] <- PNinit[4, ] <- 0
          }
          XP <- iterate_f_and_P(
            xhat_init,
            PN = PNinit,
            pvec = p,
            beta_t = bpars_fet[1],
            time.steps = c(times[T], fets$target_end_times[1])
          )
          sim_means[1, j] <- H %*% XP$xhat
          sim_cov[1, j] <- H %*% XP$PN %*% t(H) + p["tau"]
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
                pvec = p,
                beta_t = bpars_fet[i + 1],
                time.steps = c(fets$target_end_times[i], fets$target_end_times[i + 1])
              )
            sim_means[i + 1, j] <- H %*% XP$xhat
            sim_cov[i + 1, j] <-
              H %*% XP$PN %*% t(H) + p["tau"]
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

pvar <- pvar_df$init
names(pvar) <- pvar_df$par

pfixed <- c(
  N = 20e6,
  betasd = 2.,
  a = .97,
  rho1 = 0.4,
  gamma = 365.25 / 9, 
  eta = 365.25 / 4,
  iota = 0
)

pfixed["E_0"] <- pvar["I_0"] * pfixed["gamma"] / pfixed["eta"]
pfixed["S_0"] <- pfixed["N"] - pfixed["E_0"] - pvar["I_0"]

if(FALSE){
  tictoc::tic("optimization")
  ans <- optim(
    par = pvar,
    fn = kfnll,
    pfixed = pfixed,
    cdata = tail(case_data, n = wsize),
    t0 = the_t0,
    method = "L-BFGS-B",
    lower = pvar_df$lower,
    upper = pvar_df$upper,
    control = list(trace = 1, maxit = 300)
  )
  tictoc::toc()
}
