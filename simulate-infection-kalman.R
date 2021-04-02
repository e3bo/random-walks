#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

kf_nll_suff_stats <- function(w, x, pm, Phat0 = diag(c(1, 1, 1, 0, 0, 1, 1, 0))){
  p <- pm(x, w)
  E0 <- exp(p$logE0)
  I0 <- E0 * p$eta / p$gamma
  H0 <- exp(p$logH0)
  gamma_d <- gamma_h <- exp(p$loggammahd)
  D0 <- H0 * gamma_h / gamma_d
  xhat0 <- c(N - E0 - I0 - H0 - D0, E0, I0, 0, 0, H0, D0, 0)
  names(xhat0) <- c("S", "E", "I", "C", "Hnew", "H", "D", "Drep")
  doseeffect <- exp(p$logdoseeffect)
  prophomeeffect <- exp(p$logprophomeeffect)
  logbeta <- p$bpars
  
  T <- nrow(x)
  dobs <- 3
  dstate <- length(xhat0)
  stopifnot(T > 0)
  
  S <- array(NA_real_, dim = c(dobs, dobs, T))
  ybar <- array(NA_real_, dim = c(dobs, T))
  rdiagadj <- array(1, dim = c(dobs, T))
  
  
  xhat_kkmo <- array(NA_real_, dim = c(dstate, T))
  rownames(xhat_kkmo) <- names(xhat0)
  P_kkmo <- array(NA_real_, dim = c(dstate, dstate, T))
  
  H <- function(time, t0 = 2020.164){
    day <- (time - t0) * 365.25
    rbind(c(0, 0, 0, detect_frac(day), 1, 0, 0, 0), 
          c(0, 0, 0,    0, 1, 0, 0, 0),
          c(0, 0, 0,    0, 0, 0, 0, 1))
  }
  R <- diag(exp(c(p$logtauc, p$logtauh, p$logtaud)))
  
  for (i in seq(1, T)) {
    if (i == 1) {
      xhat_init <- xhat0
      PNinit <- Phat0
    } else {
      xhat_init <- xhat_kkmo[, i - 1]
      PNinit <- P_kkmo[, , i - 1]
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
      eta = p$eta,
      gamma = p$gamma,
      gamma_d = gamma_d,
      gamma_h = gamma_h,
      chp = exp(p$logchp),
      hfp = exp(p$loghfp),
      N = p$N,
      beta_t = exp(logbeta[i] -doseeffect * x$doses[i] - prophomeeffect * x$prophome[i]),
    )
    xhat_kkmo[, i] <- XP$xhat
    P_kkmo[, , i] <- XP$PN
    
    for (j in 1:dstate){
      if (P_kkmo[j,j,i] < 0){
        P_kkmo[j,,i] <- 0
        P_kkmo[,j,i] <- 0
      }
    }
    ybar[, i ] <- H(x$time[i]) %*% xhat_kkmo[, i]
    S[, , i] <- H(x$time[i]) %*% P_kkmo[, , i] %*% t(H(x$time[i])) + R + diag(rdiagadj[, i])
  }
  ret <- list(X = xhat_kkmo, P = P_kkmo, ybar = ybar, S = S)
  ret
}

stats <- kf_nll_suff_stats(w = winit, x = x, pm = param_map)
stats$ybar

