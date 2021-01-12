#!/usr/bin/env R

tictoc::tic()

library(tidyverse)
source("covidhub-common.R")

forecast_date <- Sys.getenv("fdt", unset = "2020-10-12")
forecast_loc <- "36"
hopdir <- file.path("hopkins", forecast_date)
tdat <- load_hopkins(hopdir) 

nyc <- tdat %>% filter(location == forecast_loc) %>% 
  filter(target_type == "wk ahead inc case")

nyc2 <- nyc %>% mutate(time = lubridate::decimal_date(target_end_date))
case_data <- nyc2 %>% select(time, value) %>% rename(reports = value)

target_end_dates <- max(nyc2$target_end_date) + lubridate::dweeks(1:4)
target_end_times <- lubridate::decimal_date(target_end_dates)

PsystemSEIR <- function(pvec,
                        beta_t,
                        init.vars,
                        time.steps = c(1995, 1995 + 1 / 52)) {
  PModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      eta <- 365 / 4
      gamma <- 365 / 9
      S <- exp(lS)
      E <- exp(lE)
      I <- exp(lI)
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
      
      P <- rbind(
        c(Pss, Pse, Psi, Psc),
        c(Pse, Pee, Pei, Pec),
        c(Psi, Pei, Pii, Pic),
        c(Psc, Pec, Pic, Pcc)
      )
      
      dlS <- (-beta_t * S * I / N) / S
      dlE <- (beta_t * S * I / N - eta * E) / E
      dlI <- (iota + eta * E -  gamma * I) / I
      dC <- (gamma * I)
      dP <-  F %*% P + P %*% t(F) + Q
      
      list(
        c(
          dlS = dlS,
          dlE = dlE,
          dlI = dlI,
          dC = dC,
          dPss = dP[1, 1],
          dPse = dP[1, 2],
          dPsi = dP[1, 3],
          dPsc = dP[1, 4],
          dPee = dP[2, 2],
          dPei = dP[2, 3],
          dPec = dP[2, 4],
          dPii = dP[3, 3],
          dPic = dP[3, 4],
          dPcc = dP[4, 4]
        )
      )
    })
  }
  deSolve::lsoda(init.vars, time.steps, PModel, c(pvec, beta_t = beta_t))
}

gamma <- 365/9
pvec <- params <- c(
  N = 20e6,
  beta_mu = 50,
  beta_sd = 1,
  b1 = 1 * gamma,
  b2 = 1 * gamma,
  b3 = 1 * gamma,
  b4 = 1 * gamma,
  b5 = 1 * gamma,
  b6 = 1.1 * gamma,
  iota = 2,
  rho = 0.5,
  S_0 = 20e6, 
  E_0 = 26000,
  I_0 = 13425,
  tau = 0.001
)

init.vars <- c(lS=log(3e-2 * 20e6), lE=log(1.6e-4 * 20e6), 
               lI=log(1.6e-4 * 20e6), C = 0, Pss = 1, Pse = 0, Psi = 0, Psc = 0, 
               Pee = 1, Pei = 0, Pec = 0, Pii = 1, Pic = 0, Pcc = 1)

#out <- PsystemSEIR(pvec = pvec, beta_t = 365 / 9, init.vars = init.vars, time.steps = case_data$time)

iterate_f_and_P <- function(xhat, PN, pvec, beta_t, time.steps){
  P <- PN / pvec["N"]
  xhat_trans <- c(log(xhat[c("S", "E", "I")]), xhat["C"])
  if(!all(is.finite(xhat_trans))) {
    browser()
  }
  names(xhat_trans)[1:3] <- c("lS", "lE", "lI")
  init.vars <- c(xhat_trans, Pss = P[1,1], Pse = P[1,2], 
                 Psi = P[1,3], Psc = P[1,4], Pee = P[2,2], Pei = P[2,3], Pec = P[2,4], 
                 Pii = P[3,3], Pic = P[3,4], Pcc = P[4,4])
  ret <- PsystemSEIR(pvec = pvec, init.vars = init.vars, beta_t = beta_t, time.steps)[2, ]
  xhat_new <- c(exp(ret[c("lS", "lE", "lI")]), ret["C"])
  names(xhat_new)[1:3] <- c("S", "E", "I")
  P_new  <- with(as.list(ret),        
                 rbind(c(Pss, Pse, Psi, Psc),
                       c(Pse, Pee, Pei, Pec),
                       c(Psi, Pei, Pii, Pic),
                       c(Psc, Pec, Pic, Pcc)))
  PN_new <- P_new * pvec["N"]
  list(xhat = xhat_new, PN = PN_new)
}


iterate_f_and_P(c(S=20e6, E=26e3, I=13e3, C=0), PN = diag(nrow=4), pvec = pvec, 
                beta_t = 44, time.steps = c(0, 1 / 52))

kfnll <-
  function(cdata,
           pvec,
           logit_I0,
           logit_E0,
           logit_b1,
           logit_b2,
           logit_b3,
           logit_b4,
           logit_b5,
           logit_b6,
           logit_beta_sd,
           logit_tau,
           logit_iota,
           t0,
           xhat0 = structure(c(20e6, 99.2, 99.2, 0),
                             .Dim = c(4L, 1L),
                             .Dimnames = list(c("S", "E", "I", "C"), NULL)),
           Phat0 = diag(c(1, 1, 1, 0)),
           just_nll = TRUE,
           nsim = 10,
           fets = NULL) {
    pvec["b1"] <- scaled_expit(logit_b1, a_bpar, b_bpar)
    pvec["b2"] <- scaled_expit(logit_b2, a_bpar, b_bpar)
    pvec["b3"] <- scaled_expit(logit_b3, a_bpar, b_bpar)
    pvec["b4"] <- scaled_expit(logit_b4, a_bpar, b_bpar)
    pvec["b5"] <- scaled_expit(logit_b5, a_bpar, b_bpar)
    pvec["b6"] <- scaled_expit(logit_b6, a_bpar, b_bpar)
    pvec["beta_sd"] <- scaled_expit(logit_beta_sd, a_beta_sd, b_beta_sd)
    pvec["rho"] <- 0.4 # scaled_expit(logit_rho, a_rho, b_rho)
    pvec["iota"] <- scaled_expit(logit_iota, a_iota, b_iota)
    xhat0["S", 1] <- 20e6
    xhat0["I", 1] <- scaled_expit(logit_I0, a_I0, b_I0)
    xhat0["E", 1] <- scaled_expit(logit_E0, a_E0, b_E0)

    is_spline_par <- grepl("^b[0-9]+$", names(pvec))
    bpars <- pvec[is_spline_par]
    stopifnot(length(bpars) == nrow(cdata) + 1)
    #obs_per_b <- ceiling(length(cdata$reports) / 6)
    #bpars <- c(pvec["b1"], pvec["b2"], pvec["b3"], pvec["b4"], pvec["b5"], pvec["b6"])
    #y <- rep(bpars, each = obs_per_b)
    #beta_fun <- approxfun(x = cdata$time, y[seq_along(cdata$reports)], rule = 2)
    #x <- c(t0, cdata$time)
    #bpars <- rep(pvec["b1"], length(x))
    #beta_fun <- approxfun(x = x, y = bpars, rule = 2)
    
    #print(c("                                     ", pvec["beta_mu"], xhat0["S", 1] / b_S0, pvec["b2"]))
    
    # Initialize
    z_1 <-  cdata$reports[1]
    H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
    R <- max(5, z_1 * pvec["tau"])
    #R <- max(1, z[1] * (1 - pvec["rho"]))
    #R <- max(tau2, z_1 + z_1 ^ 2 /  pvec["tau"])
    
    
    # Predict
    XP_1_0 <- iterate_f_and_P(
      xhat0[, 1],
      PN = Phat0,
      pvec = pvec,
      beta_t = bpars[1],
      time.steps = c(t0, cdata$time[1])
    )
    xhat_1_0 <- XP_1_0$xhat
    P_1_0 <- XP_1_0$PN
    # Update
    
    K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R)
    ytilde_1 <- z_1 - H %*% xhat_1_0
    xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
    xhat_1_1[xhat_1_1 < 0] <- 1e-4
    P_1_1 <- (diag(4) - K_1 %*% H) %*% P_1_0
    
    ## Now calculate for each step in simulation
    
    T <- nrow(cdata)
    z <- cdata$reports
    
    ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))
    
    K[, 1] <- K_1
    xhat_kkmo[, 1] <- xhat_1_0
    xhat_kk[, 1] <- xhat_1_1
    P_kk[, , 1] <- P_1_1
    P_kkmo[, , 1] <- P_1_0
    S[, 1] <- H %*% P_kkmo[, , 1] %*% t(H) + R
    ytilde_kk[, 1] <- z[1] - H %*% xhat_kk[, 1]
    ytilde_k[, 1] <- ytilde_1
    
    if (T > 1) {
      for (i in seq(2, T)) {
        xhat_init <- xhat_kk[, i - 1]
        xhat_init["C"] <- 0
        PNinit <- P_kk[, , i - 1]
        PNinit[, 4] <- PNinit[4,] <- 0
        XP <- iterate_f_and_P(
          xhat_init,
          PN = PNinit,
          pvec = pvec,
          beta_t = bpars[i + 1],
          time.steps = cdata$time[c(i - 1, i)]
        )
        xhat_kkmo[, i] <- XP$xhat
        P_kkmo[, , i] <- XP$PN
        R <- z[i - 1] * pvec["tau"]
        #R <- max(5, z[i - 1] * pvec["tau"])
        #R <- max(1, z[i - 1] * pvec["rho"]))
        #R <- max(tau2, xhat_kkmo["C", i] * pvec["rho"] + (xhat_kkmo["C", i] * pvec["rho"]) ^ 2 / pvec["tau"])
        S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
        K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
        ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
        xhat_kk[, i] <-
          xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
        xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
        P_kk[, , i] <-
          (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
        ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
      }
    }
    
    rwlik <-
      sum(dnorm(diff(bpars), sd = pvec["beta_sd"], log = TRUE))
    nll <-
      0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi)) - rwlik
    if (!just_nll) {
      sim_means <- sim_cov <- matrix(NA, nrow = (length(fets)), ncol = nsim)
      if (!is.null(fets)) {
        for(j in seq_len(nsim)){
          bpars_fet <- bpars[T + 1] + cumsum(rnorm(n = length(fets), mean = 0, sd = pvec["beta_sd"]))
          xhat_init <- xhat_kk[, T]
          xhat_init["C"] <- 0
          PNinit <- P_kk[, , T]
          PNinit[, 4] <- PNinit[4,] <- 0
          XP <- iterate_f_and_P(
            xhat_init,
           PN = PNinit,
            pvec = pvec,
            beta_t = bpars_fet[1],
            time.steps = c(cdata$time[T], fets[1])
          )
          sim_means[1, j] <- H %*% XP$xhat
          sim_cov[1, j] <- H %*% XP$PN %*% t(H)
          for (i in seq_along(fets[-1])) {
            xhat_init <- XP$xhat
            xhat_init["C"] <- 0
            PNinit <- XP$PN
            PNinit[, 4] <- PNinit[4,] <- 0
            XP <-
              iterate_f_and_P(
                xhat_init,
                PN = PNinit,
                pvec = pvec,
                beta_t = bpars_fet[i + 1],
                time.steps = c(fets[i], fets[i + 1])
              )
            sim_means[i + 1, j] <- H %*% XP$xhat
            sim_cov[i + 1, j] <- H %*% XP$PN %*% t(H)
          }
        }
        pred_means <- rowMeans(sim_means)
        pred_cov <- rowMeans(sim_cov)
      } else {
        pred_means <- pred_cov <- NULL
      }
      list(
        nll = nll,
        xhat_kkmo = xhat_kkmo,
        xhat_kk = xhat_kk,
        P_kkmo = P_kkmo,
        P_kk = P_kk,
        ytilde_k = ytilde_k,
        S = S,
        pred_means = pred_means,
        pred_cov = pred_cov
      )
    } else {
      nll
    }
  }

library(bbmle)

scaled_expit <- function(y, a, b){
  (b - a) * exp(y)/ (1 + exp(y)) + a
} 

scaled_logit <- function(x, a, b){
  log((x - a) / (b - x))
}
a_beta_mu <- 5
b_beta_mu <- 1000

a_I0 <- 0
b_I0 <- 1e5
a_E0 <- 0
b_E0 <- 1e5

a_bpar <- 0
b_bpar <- 10 * gamma

a_iota <- 0
b_iota <- 300

a_tau <- 0.0001
b_tau <- 1

a_tau2 <- 0
b_tau2 <- 20

a_beta_sd <- 0
b_beta_sd <- 10

Phat0 <- diag(c(1e4, 1e2, 1e2, 0))


## Find reasonable initial values for transmission terms


#plot(case_data$reports)
Rt <- case_data$reports[-1] / case_data$reports[-nrow(case_data)]

kfret_sample <-  
              kfnll(cdata = tail(case_data, n = 3), pvec = pvec, 
                    logit_E0 = scaled_logit(pvec["E_0"], a_E0, b_E0),
                    logit_I0 = scaled_logit(pvec["I_0"], a_I0, b_I0),
                    logit_b1 = scaled_logit(pvec["b1"], a_bpar, b_bpar),
                    logit_b2 = scaled_logit(pvec["b2"], a_bpar, b_bpar), 
                    logit_b3 = scaled_logit(pvec["b3"], a_bpar, b_bpar),
                    logit_b4 = scaled_logit(pvec["b4"], a_bpar, b_bpar),
                    logit_b5 = scaled_logit(pvec["b5"], a_bpar, b_bpar),
                    logit_b6 = scaled_logit(pvec["b6"], a_bpar, b_bpar),
                    logit_iota = scaled_logit(pvec["iota"], a_iota, b_iota),
                    just_nll = FALSE,
                    fet = target_end_times)

the_n <- 5
the_t0 <- rev(case_data$time)[the_n + 1]
system.time(
  m0 <- mle2(
    minuslogl = kfnll,
    start = list(
      logit_I0 = scaled_logit(14000, a_I0, b_I0),
      logit_E0 = scaled_logit(16000, a_E0, b_E0),
      logit_b1 = scaled_logit(40, a_bpar, b_bpar),
      logit_b2 = scaled_logit(40, a_bpar, b_bpar),
      logit_b3 = scaled_logit(40, a_bpar, b_bpar),
      logit_b4 = scaled_logit(40, a_bpar, b_bpar),
      logit_b5 = scaled_logit(40, a_bpar, b_bpar),
      logit_b6 = scaled_logit(40, a_bpar, b_bpar),
      logit_beta_sd = scaled_logit(0.01, a_beta_sd, b_beta_sd),
      logit_tau = scaled_logit(0.1, a_tau, b_tau),
      logit_iota = scaled_logit(0.6, a_iota, b_iota)
    ),
    method = "Nelder-Mead",
    skip.hessian = TRUE,
    control = list(
      reltol = 1e-4,
      trace = 1,
      maxit = 2000
    ),
    data = list(
      cdata = tail(case_data, n = the_n),
      pvec = pvec,
      Phat0 = diag(c(1, 1, 1, 0)),
      t0 = the_t0
    )
  )
)

kfret <- with(
  as.list(coef(m0)),
  kfnll(
    cdata = tail(case_data, n = the_n),
    pvec = pvec,
    logit_E0 = logit_E0,
    logit_I0 = logit_I0,
    logit_b1 = logit_b1,
    logit_b2 = logit_b2,
    logit_b3 = logit_b3,
    logit_b4 = logit_b4,
    logit_b5 = logit_b5,
    logit_b6 = logit_b6,
    logit_beta_sd = logit_beta_sd,
    logit_iota = logit_iota,
    t0 = the_t0,
    just_nll = FALSE,
    fet = target_end_times,
    Phat0 = diag(c(1, 1, 1, 0))
  )
)

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
  n <- length(means)
  stopifnot(length(vars) == n)
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
    mutate(q1 = qnorm(
      p = quantile,
      mean = means[h],
      sd = sqrt(vars[h])
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

fcst <- create_forecast_df(means = kfret$pred_means,
                           vars = kfret$pred_cov,
                           location = forecast_loc)

fcst_path <- file.path("forecasts", paste0(forecast_date, "-CEID-SIR_EKF.csv"))
if(!dir.exists("forecasts")) dir.create("forecasts")
write_csv(x = fcst, path = fcst_path)

## Produce metrics
time <- tictoc::toc()
walltime <- list(wall = time$toc - time$tic)
if(!dir.exists("metrics")) dir.create("metrics")
mpath <- file.path("metrics", paste0(forecast_date, 
                                     "-forecast-calc-time.json"))
jsonlite::write_json(walltime, mpath, auto_unbox = TRUE)


q("no")

scaled_expit(coef(m0)["logit_I0"], a_I0, b_I0)
scaled_expit(coef(m0)["logit_E0"], a_E0, b_E0)
#(rho_hat <- scaled_expit(coef(m0)["logit_rho"], a_rho, b_rho))
scaled_expit(coef(m0)["logit_iota"], a_iota, b_iota)
scaled_expit(coef(m0)["logit_tau"], a_tau, b_tau)
scaled_expit(coef(m0)["logit_beta_sd"], a_beta_sd, b_beta_sd)
#scaled_expit(coef(m0)["logit_tau2"], a_tau2, b_tau2)

is_spline_par <- grepl("^logit_b[1-6]$", names(coef(m0)))
bhat <- scaled_expit(coef(m0)[is_spline_par], a_bpar, b_bpar)
R0hat <- bhat / gamma



par(mfrow = c(1, 1))
test <- case_data$time > 1990
qqnorm(kfret$ytilde_k[test]/ kfret$S[test]) # evalutate departure from normality
abline(0, 1)

rho_hat <- 0.4
test <- case_data$time >= 2020.18
par(mfrow = c(4, 1))
plot(case_data$time[test], kfret$xhat_kkmo["C",] * rho_hat, ylim = c(0, 1e5))
points(case_data$time[test], kfret$xhat_kk["C",] * rho_hat, col = 2, pch = 2)
lines(case_data$time[test], case_data$reports[test])
plot(case_data$time[test], kfret$S, log = "y")
plot(case_data$time[test], kfret$ytilde_k)
plot(case_data$time[test], kfret$ytilde_k / kfret$S)

par(mfrow = c(2, 1))
plot(case_data$time[test], kfret$xhat_kkmo["I",])
points(case_data$time[test], kfret$xhat_kk["I",], col = 2, pch = 2)

plot(case_data$time[test], kfret$xhat_kkmo["S",])
points(case_data$time[test], kfret$xhat_kk["S",], col = 2, pch = 2)


tgrid <- c(case_data$time, target_end_times)


matplot(tgrid, ximat)

is_spline_par <- grepl("^logit_b[1-6]$", names(coef(m0)))
bhat <- scaled_expit(coef(m0)[is_spline_par], a_bpar, b_bpar)
R0hat <- bhat / gamma

plot(tgrid, R0grid, log = "y")
## R0 stays constant at last value in forecast, seems like a good starting point

