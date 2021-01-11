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

cov_data0 <- tibble(time = case_data$time)

bspline_basis <- pomp::periodic.bspline.basis(
  cov_data0$time,
  nbasis = 6,
  degree = 3,
  period = 1,
  names = "xi%d"
) %>%
  as_tibble()

cov_data <- bind_cols(cov_data0, bspline_basis)

PsystemSEIR <- function(pvec, covf,
                        init.vars,
                        time.steps = c(1995, 1995 + 1 / 52)){
  
  PModel <- function(t, x, parms) {
    with(as.list(c(parms, x)), {
      xi1t <- covf$xi1(t)
      xi2t <- covf$xi2(t)
      xi3t <- covf$xi3(t)
      xi4t <- covf$xi4(t)
      xi5t <- covf$xi5(t)
      xi6t <- covf$xi6(t)
      beta_t <- beta_mu * (1 + exp(xi1t * b1 + xi2t * b2 + xi3t * b3 + xi4t * b4 + xi5t * b5 + xi6t * b6))
      eta <- 365 / 4  
      gamma <- 365 / 9
      S <- exp(lS)
      E <- exp(lE)
      I <- exp(lI)
      F <- rbind(c(-beta_t * I / N,    0, -beta_t * S / N, 0),
                 c( beta_t * I / N, -eta,  beta_t * S / N, 0),
                 c(                0,  eta,            -gamma, 0),
                 c(                0,    0,             gamma, 0))
      
      
      f <- c(0, beta_t * (S / N) * (I / N), eta * E / N, gamma * I / N)
      Q <- rbind(c(f[1] + f[2],       -f[2],           0,     0),
                 c(      -f[2], f[2] + f[3],       -f[3],     0),
                 c(          0,       -f[3], f[3] + f[4], -f[4]),
                 c(          0,           0,       -f[4],  f[4]))
      
      P <- rbind(c(Pss, Pse, Psi, Psc),
                 c(Pse, Pee, Pei, Pec),
                 c(Psi, Pei, Pii, Pic),
                 c(Psc, Pec, Pic, Pcc))
      
      dlS <- (- beta_t * S * I / N) / S
      dlE <- (beta_t * S * I / N - eta * E) / E
      dlI <- (iota + eta * E -  gamma * I) / I
      dC <- (gamma * I)
      dP <-  F %*% P + P %*% t(F) + Q
      
      list(c(dlS = dlS, dlE = dlE, dlI = dlI, dC = dC, dPss = dP[1,1], dPse = dP[1,2], 
             dPsi = dP[1,3], dPsc = dP[1,4], dPee = dP[2,2], dPei = dP[2,3], 
             dPec = dP[2,4], dPii = dP[3,3], dPic = dP[3,4], dPcc = dP[4,4]))
    })
  }
  deSolve::lsoda(init.vars, time.steps, PModel, pvec)
}

genfun <- function(y) {
  approxfun(cov_data$time, y, rule = 2)
}
covf <- apply(cov_data[, -1], 2, genfun)

pvec <-   params <- c(
  N = 20e6,
  beta_mu = 50,
  beta_sd = 0.001,
  b1 = 3,
  b2 = 3,
  b3 = 6,
  b4 = 5,
  b5 = 2,
  b6 = 0,
  iota = 2,
  rho = 0.5,
  S_0 = 20e6, 
  E_0 = 10,
  I_0 = 10,
  tau = 0.001
)

init.vars <- c(lS=log(3e-2 * 20e6), lE=log(1.6e-4 * 20e6), 
               lI=log(1.6e-4 * 20e6), C = 0, Pss = 1, Pse = 0, Psi = 0, Psc = 0, 
               Pee = 1, Pei = 0, Pec = 0, Pii = 1, Pic = 0, Pcc = 1)

#out <- PsystemSEIR(pvec = pvec, covf = covf, init.vars = init.vars, time.steps = case_data$time)

iterate_f_and_P <- function(xhat, PN, pvec, covf, time.steps){
  P <- PN / pvec["N"]
  xhat_trans <- c(log(xhat[c("S", "E", "I")]), xhat["C"])
  if(!all(is.finite(xhat_trans))) {
    browser()
  }
  names(xhat_trans)[1:3] <- c("lS", "lE", "lI")
  init.vars <- c(xhat_trans, Pss = P[1,1], Pse = P[1,2], 
                 Psi = P[1,3], Psc = P[1,4], Pee = P[2,2], Pei = P[2,3], Pec = P[2,4], 
                 Pii = P[3,3], Pic = P[3,4], Pcc = P[4,4])
  ret <- PsystemSEIR(pvec = pvec, init.vars = init.vars, covf, time.steps)[2, ]
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

xhat <-  c(S=(3e-2 * 20e6), E=(1.6e-4 * 20e6), I=(1.6e-4 * 20e6), C = 0)
P <- with(as.list(init.vars[-c(1:4)]),        
          rbind(c(Pss, Pse, Psi, Psc),
                c(Pse, Pee, Pei, Pec),
                c(Psi, Pei, Pii, Pic),
                c(Psc, Pec, Pic, Pcc)))

if (FALSE){ # useful only for debugging
iterate_f_and_P(xhat = xhat, PN = P, pvec = pvec, covf = covf, 
                time.steps = c(2020.1, 2020.1 + 1 / 52))

# Initialize
xhat0 <- matrix(xhat, ncol = 1)
rownames(xhat0) <- names(xhat)
Phat0 <- P
z_1 <-  case_data$reports[1]
H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
R <- z_1 * (1 - pvec["rho"])

# Predict
XP_1_0 <- iterate_f_and_P(xhat0[, 1], PN = Phat0, pvec = pvec, covf = covf,
                          time.steps = case_data$time[c(1, 2)])
xhat_1_0 <- XP_1_0$xhat
P_1_0 <- XP_1_0$PN
# Update

K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R[1])
ytilde_1 <- z_1 - H %*% xhat_1_0
xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
P_1_1 <- (diag(4) - K_1 %*% H) %*% P_1_0

## Now calculate for each step in simulation


T <- nrow(case_data)
z <- case_data$reports

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

for (i in seq(2, T)){
  xhat_init <- xhat_kk[, i - 1]
  xhat_init["C"] <- 0
  PNinit <- P_kk[,,i - 1]
  PNinit[, 4] <- PNinit[4, ] <- 0
  XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                        time.steps = case_data$time[c(i - 1, i)])
  xhat_kkmo[, i] <- XP$xhat
  P_kkmo[, , i] <- XP$PN
  #R <- xhat_kkmo["C", i] * pvec["rho"] * (1 - pvec["rho"])
  R <- max(5, z[i - 1] * (1 - pvec["rho"]))
  S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
  K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
  ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
  xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
  xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
  P_kk[, , i] <- (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
  ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
}

log_lik <- function(Sigma, resids){
  -0.5 * sum(resids ^ 2 / Sigma + log(Sigma) + log(2 * pi))
}

log_lik(S, ytilde_k)

}

kfnll <-
  function(cdata,
           pvec,
           logit_beta_mu,
           logit_I0,
           logit_E0,
           logit_b1,
           logit_b2,
           logit_b3,
           logit_b4, 
           logit_b5,
           logit_b6,
           logit_tau,
           logit_iota,
           xhat0 = structure(c(20e6, 99.2, 99.2, 0), .Dim = c(4L, 1L), 
                             .Dimnames = list(c("S", "E", "I", "C"), NULL)),
           Phat0 = diag(c(1, 1, 1, 0)),
           just_nll = TRUE,
           fets = NULL) {
    
    pvec["beta_mu"] <- scaled_expit(logit_beta_mu, a_beta_mu, b_beta_mu)
    pvec["b1"] <- scaled_expit(logit_b1, a_bpar, b_bpar)
    pvec["b2"] <- scaled_expit(logit_b2, a_bpar, b_bpar)
    pvec["b3"] <- scaled_expit(logit_b3, a_bpar, b_bpar)
    pvec["b4"] <- scaled_expit(logit_b4, a_bpar, b_bpar)
    pvec["b5"] <- scaled_expit(logit_b5, a_bpar, b_bpar)
    pvec["b6"] <- scaled_expit(logit_b6, a_bpar, b_bpar)
    pvec["rho"] <- 0.4 # scaled_expit(logit_rho, a_rho, b_rho)
    pvec["iota"] <- scaled_expit(logit_iota, a_iota, b_iota)
    xhat0["S", 1] <- 20e6
    xhat0["I", 1] <- scaled_expit(logit_I0, a_I0, b_I0)
    xhat0["E", 1] <- scaled_expit(logit_E0, a_E0, b_E0)

    
    #print(c("                                     ", pvec["beta_mu"], xhat0["S", 1] / b_S0, pvec["b2"]))
    
    # Initialize
    z_1 <-  cdata$reports[1]
    H <- matrix(c(0, 0, 0, pvec["rho"]), ncol = 4)
    R <- max(5, z_1 * pvec["tau"])
    #R <- max(1, z[1] * (1 - pvec["rho"]))
    #R <- max(tau2, z_1 + z_1 ^ 2 /  pvec["tau"]) 
    
    
    # Predict
    XP_1_0 <- iterate_f_and_P(xhat0[, 1], PN = Phat0, pvec = pvec, covf = covf,
                              time.steps = cdata$time[c(2, 3)])
    xhat_1_0 <- XP_1_0$xhat
    P_1_0 <- XP_1_0$PN
    # Update
    
    K_1 <- P_1_0 %*% t(H) %*% solve(H %*% P_1_0 %*% t(H) + R)
    ytilde_1 <- z_1 - H %*% xhat_1_0
    xhat_1_1 <- xhat_1_0 + K_1 %*% ytilde_1
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
    
    for (i in seq(2, T)){
      xhat_init <- xhat_kk[, i - 1]
      xhat_init["C"] <- 0
      PNinit <- P_kk[,,i - 1]
      PNinit[, 4] <- PNinit[4, ] <- 0
      XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                            time.steps = cdata$time[c(i - 1, i)])
      xhat_kkmo[, i] <- XP$xhat
      P_kkmo[, , i] <- XP$PN
      R <- max(5, z[i - 1] * pvec["tau"])
      #R <- max(1, z[i - 1] * pvec["rho"]))
      #R <- max(tau2, xhat_kkmo["C", i] * pvec["rho"] + (xhat_kkmo["C", i] * pvec["rho"]) ^ 2 / pvec["tau"])
      S[, i] <- H %*% P_kkmo[, , i] %*% t(H) + R
      K[, i] <- P_kkmo[, , i] %*% t(H) %*% solve(S[, i])
      ytilde_k[, i] <- z[i] - H %*% xhat_kkmo[, i, drop = FALSE]
      xhat_kk[, i] <- xhat_kkmo[, i, drop = FALSE] + K[, i, drop = FALSE] %*% ytilde_k[, i, drop = FALSE]
      xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
      P_kk[, , i] <- (diag(4) - K[, i, drop = FALSE] %*% H) %*% P_kkmo[, , i]
      ytilde_kk[i] <- z[i] - H %*% xhat_kk[, i, drop = FALSE]
    }
    
    nll <- 0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi))
    if (!just_nll){
      if (!is.null(fets)){
        xhat_init <- xhat_kk[, T]
        xhat_init["C"] <- 0
        PNinit <- P_kk[,, T]
        PNinit[, 4] <- PNinit[4, ] <- 0
        XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                             time.steps = c(cdata$time[T], fets[1]))
        pred_means <- pred_cov <- numeric(length(fets))
        pred_means[1] <- H %*% XP$xhat 
        pred_cov[1] <- H %*% XP$PN %*% t(H)
        for(i in seq_along(fets[-1])){
          xhat_init <- XP$xhat
          xhat_init["C"] <- 0
          PNinit <- XP$PN
          PNinit[, 4] <- PNinit[4, ] <- 0
          XP <- iterate_f_and_P(xhat_init, PN = PNinit, pvec = pvec, covf = covf,
                                time.steps = c(fets[i], fets[i + 1]))
          pred_means[i + 1] <- H %*% XP$xhat
          pred_cov[i + 1] <- H %*% XP$PN %*% t(H)
        }
      } else {
        pred_means <- pred_cov <- NULL
      }
      list(nll = nll, xhat_kkmo = xhat_kkmo, xhat_kk = xhat_kk, 
           P_kkmo = P_kkmo, P_kk = P_kk, 
           ytilde_k = ytilde_k, S = S, pred_means = pred_means, 
           pred_cov = pred_cov)
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
b_I0 <- 100
a_E0 <- 0
b_E0 <- 200

a_bpar <- -11
b_bpar <- 11

a_iota <- 0
b_iota <- 300

a_tau <- 0.0001
b_tau <- 2

a_tau2 <- 0
b_tau2 <- 20

pvec2 <- pvec
pvec2["rho"] <- 0.1
pvec2["b1"] <- 0.3
pvec2["b2"] <- 0.3
pvec2["b3"] <- 0.6
pvec2["b4"] <- 0.5
pvec2["b5"] <- 0.2
pvec2["b6"] <- 0

Phat0 <- diag(c(1e4, 1e2, 1e2, 0))


## Find reasonable initial values for transmission terms


#plot(case_data$reports)
Rt <- case_data$reports[-1] / case_data$reports[-nrow(case_data)]

lhs <- log(Rt / 0.4 - 1)
rhs <- cov_data[-1, -1]

rdata <- cbind(lhs, rhs)

bm <- lm(lhs~0 + xi1 + xi2 + xi3 + xi4 + xi5 + xi6, data = rdata %>% 
           filter(is.finite(lhs)))

system.time(m0 <- mle2(minuslogl = kfnll, 
                       start = list(logit_beta_mu = scaled_logit(22, a_beta_mu, b_beta_mu), 
                                    logit_I0 = scaled_logit(85, a_I0, b_I0),
                                    logit_E0 = scaled_logit(197, a_E0, b_E0),
                                    logit_b1 = scaled_logit(10, a_bpar, b_bpar),
                                    logit_b2 = scaled_logit(4.6, a_bpar, b_bpar),
                                    logit_b3 = scaled_logit(-1.7, a_bpar, b_bpar),
                                    logit_b4 = scaled_logit(0.99, a_bpar, b_bpar),
                                    logit_b5 = scaled_logit(0.062, a_bpar, b_bpar),
                                    logit_b6 = scaled_logit(0.87, a_bpar, b_bpar),
                                    logit_tau = scaled_logit(0.005, a_tau, b_tau),
                                    logit_iota = scaled_logit(0.6, a_iota, b_iota)),
                       method = "Nelder-Mead",
                       skip.hessian = TRUE,
                       control = list(reltol = 1e-4, trace = 1, maxit = 1000),
                       data = list(cdata = case_data[-c(1,2,3,4,5),], pvec = pvec2, Phat0 = Phat0)))

kfret <- with(as.list(coef(m0)), 
              kfnll(cdata = case_data[-seq(1,5),], pvec = pvec2, 
                    logit_beta_mu = logit_beta_mu, 
                    logit_E0 = logit_E0,
                    logit_I0 = logit_I0,
                    logit_b1 = logit_b1,
                    logit_b2 = logit_b2, 
                    logit_b3 = logit_b3,
                    logit_b4 = logit_b4,
                    logit_b5 = logit_b5,
                    logit_b6 = logit_b6,
                    logit_iota = logit_iota,
                    just_nll = FALSE,
                    fet = target_end_times))

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
(rho_hat <- scaled_expit(coef(m0)["logit_rho"], a_rho, b_rho))
scaled_expit(coef(m0)["logit_iota"], a_iota, b_iota)
scaled_expit(coef(m0)["logit_tau"], a_tau, b_tau)
scaled_expit(coef(m0)["logit_tau2"], a_tau2, b_tau2)

par(mfrow = c(1, 1))
test <- case_data$time > 1990
qqnorm(kfret$ytilde_k[test]/ kfret$S[test]) # evalutate departure from normality
abline(0, 1)

rho_hat <- 0.4
test <- case_data$time >= 2020.18
par(mfrow = c(4, 1))
plot(case_data$time[test], kfret$xhat_kkmo["C",] * rho_hat)
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


tgrid <- c(case_data$time, forecast_times)
ximat <- cbind(covf$xi1(tgrid),
               covf$xi2(tgrid),
               covf$xi3(tgrid),
               covf$xi4(tgrid),
               covf$xi5(tgrid),
               covf$xi6(tgrid))

matplot(tgrid, ximat)

is_spline_par <- grepl("^logit_b[1-6]$", names(coef(m0)))
bhat <- scaled_expit(coef(m0)[is_spline_par], a_bpar, b_bpar)
seasgrid <- 1 + exp(ximat %*% bhat)
beta_mu_hat <- scaled_expit(coef(m0)["logit_beta_mu"], a_beta_mu, b_beta_mu)
R0grid <- beta_mu_hat * seasgrid / (365 / 5)
plot(tgrid, R0grid, log = "y")
## R0 stays constant at last value in forecast, seems like a good starting point

