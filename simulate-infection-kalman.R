#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
source("covidhub-common.R")

kf_nll_suff_stats <- function(w, x, pm, Phat0 = diag(c(1, 1, 1, 0, 0, 1, 1, 0))){
  p <- pm(x, w)
  E0 <- exp(p$logE0)
  I0 <- E0 * p$eta / p$gamma
  H0 <- exp(p$logH0)
  gamma_d <- gamma_h <- exp(p$loggammahd)
  D0 <- H0 * gamma_h / gamma_d * exp(p$loghfp)
  xhat0 <- c(p$N - E0 - I0 - H0 - D0, E0, I0, 0, 0, H0, D0, 0)
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

wsim <- c(logE0 = 8, logH0 = log(116), logtauc = 4, 
  logtauh = 3, logtaud =2, logchp = -2.54689972836648, 
  loghfp = -1.64459615121798, loggammahd = 5.20743487007086, logdoseeffect.N = -16.7835406339017, 
  logprophomeeffect = -20, b1 = 4.51428768951092, b2 = 4.51428768951092, 
  b3 = 4.51428768951092, b4 = 4.51428768951092, b5 = 4.51428768951092, 
  b6 = 4.51428768951092, b7 = 4.51428768951092, b8 = 4.51428768951092, 
  b9 = 4.51428768951092, b10 = 4.51428768951092, b11 = 4.51428768951092, 
  b12 = 4.51428768951092, b13 = 4.51428768951092, b14 = 4.51428768951092, 
  b15 = 4.51428768951092, b16 = 4.51428768951092, b17 = 4.51428768951092, 
  b18 = 4.51428768951092, b19 = 4.51428768951092, b20 = 4.51428768951092, 
  b21 = 4.51428768951092, b22 = 4.51428768951092, b23 = 4.51428768951092, 
  b24 = 4.51428768951092, b25 = 4.51428768951092, b26 = 4.51428768951092, 
  b27 = 4.51428768951092, b28 = 4.51428768951092, b29 = 4.51428768951092, 
  b30 = 4.51428768951092, b31 = 4.51428768951092, b32 = 4.51428768951092, 
  b33 = 4.51428768951092, b34 = 4.51428768951092, b35 = 4.51428768951092, 
  b36 = 4.51428768951092, b37 = 4.51428768951092, b38 = 4.51428768951092, 
  b39 = 4.51428768951092, b40 = 4.51428768951092, b41 = 4.51428768951092, 
  b42 = 4.51428768951092, b43 = 4.51428768951092, b44 = 4.51428768951092, 
  b45 = 4.51428768951092, b46 = 4.51428768951092, b47 = 4.51428768951092, 
  b48 = 4.51428768951092, b49 = 4.51428768951092, b50 = 4.51428768951092, 
  b51 = 4.51428768951092, b52 = 4.51428768951092, b53 = 4.51428768951092, 
  b54 = 4.51428768951092, b55 = 4.51428768951092, b56 = 4.51428768951092, 
  b57 = 4.51428768951092, b58 = 4.51428768951092, b59 = 4.51428768951092, 
  b60 = 4.51428768951092, b61 = 4.51428768951092, b62 = 4.51428768951092, 
  b63 = 4.51428768951092)

x <-
  structure(
    list(
      time = c(
        2021.06575342466,
        2021.06849315068,
        2021.07123287671,
        2021.07397260274,
        2021.07671232877,
        2021.07945205479,
        2021.08219178082,
        2021.08493150685,
        2021.08767123288,
        2021.0904109589,
        2021.09315068493,
        2021.09589041096,
        2021.09863013699,
        2021.10136986301,
        2021.10410958904,
        2021.10684931507,
        2021.1095890411,
        2021.11232876712,
        2021.11506849315,
        2021.11780821918,
        2021.12054794521,
        2021.12328767123,
        2021.12602739726,
        2021.12876712329,
        2021.13150684931,
        2021.13424657534,
        2021.13698630137,
        2021.1397260274,
        2021.14246575342,
        2021.14520547945,
        2021.14794520548,
        2021.15068493151,
        2021.15342465753,
        2021.15616438356,
        2021.15890410959,
        2021.16164383562,
        2021.16438356164,
        2021.16712328767,
        2021.1698630137,
        2021.17260273973,
        2021.17534246575,
        2021.17808219178,
        2021.18082191781,
        2021.18356164384,
        2021.18630136986,
        2021.18904109589,
        2021.19178082192,
        2021.19452054795,
        2021.19726027397,
        2021.2,
        2021.20273972603,
        2021.20547945205,
        2021.20821917808,
        2021.21095890411,
        2021.21369863014,
        2021.21643835616,
        2021.21917808219,
        2021.22191780822,
        2021.22465753425,
        2021.22739726027,
        2021.2301369863,
        2021.23287671233,
        2021.23561643836
      ),
      doses = c(
        -0.71,
        0.973982485380633,
        1.86416977197159,
        -2.22658922511262,
        0.889918452522017,
        0.470906570981503,
        -0.257419590116702,
        0.212307132840038,-0.83105881614593,
        -0.108248101367014,
        -0.949788953456508,
        -1.14504374421476,
        1.31491920545564,
        -1.24216899752115,
        -0.981977494846581,
        -1.2053168906173,
        1.44418964293532,
        1.94277898101826,
        -0.353139777832902,
        0.367082692667085,-1.09926029921498,
        1.79958799618293,
        -0.3574216005876,
        0.356776246672285,-1.12704638051469,
        1.38566834746188,
        -0.163738651220553,
        1.58310634495189,
        0.765189406153406,
        -0.290133548867813,
        -0.582679670538113,
        0.72278941862414,-0.583195688572441,
        -1.10686916332658,
        -1.85899284818982,
        0.699407033593285,-0.0545276522155156,
        -1.11665774767763,
        -0.160557695397889,
        0.456771756486734,-0.472282466937697,
        0.352091820781306,
        0.119165392829192,
        -0.829498252598463,
        1.01149671738213,
        -0.36521007303057,
        0.76935506180555,
        2.2405303818552,
        0.0552331675316245,
        -0.256496833768724,
        -0.677368965545138,
        -0.782734798555358,
        0.0604586906367204,
        -0.947417250795686,
        0.720873215154343,
        0.153626435834447,
        1.04437422220652,
        -0.94325049230429,
        -1.95947587238177,
        0.313183163096498,-0.803523778148557,
        -0.398439110190028,
        0.653765025886135
      ),
      prophome = c(
        0.38817,
        0.3901391,
        0.3902349,
        0.389917,
        0.3794764,
        0.3907247,
        0.3920108,
        0.3901786,
        0.3890768,
        0.3873175,
        0.3932854,
        0.3705576,
        0.3575584,
        0.3549622,
        0.3557152,
        0.3548081,
        0.3533463,
        0.3411294,
        0.3410463,
        0.3395478,
        0.3355191,
        0.342282,
        0.3446422,
        0.3447609,
        0.3471406,
        0.3460999,
        0.3435818,
        0.3438205,
        0.3340976,
        0.3260415,
        0.3263853,
        0.3288008,
        0.325047,
        0.3243972,
        0.3217037,
        0.3210908,
        0.3233132,
        0.3195206,
        0.3149208,
        0.316503,
        0.316354,
        0.3190463,
        0.3157146,
        0.3149865,
        0.3163087,
        0.3171952,
        0.3158711,
        0.3139757,
        0.3112235,
        0.3135903,
        0.3120052,
        0.3078487,
        0.2992052,
        0.2928327,
        0.2883606,
        0.2886316,
        0.2842454,
        0.2811453,
        0.2811453,
        0.2811453,
        0.2811453,
        0.2811453,
        0.2811453
      )
    ),
    row.names = c(NA,-63L),
    class = c("tbl_df",
              "tbl", "data.frame")
  )

wfixed <- c(
  N = 1e7,
  rho1 = 0.4,
  gamma = 365.25 / 4,
  gamma_h = 365.25 / 10,
  gamma_d = 365.25 / 10,
  eta = 365.25 / 4
)
stats <- kf_nll_suff_stats(w = wsim, x = x, pm = param_map)

sampts <- function(stats){
  y <- array(NA_real_, dim = dim(t(stats$ybar)), 
             dimnames = list(NULL, c("cases", "hospitalizations", "deaths"))) %>% 
    as_tibble()
  for (i in 1:nrow(y)){
      y[i, ] <- mvtnorm::rmvnorm(mean = stats$ybar[, i], sigma = stats$S[,,i], n = 1) 
  }
  y
}

ysim <- replicate(10, sampts(stats), simplify = FALSE)

JuliaCall::julia_setup("/opt/julia-1.5.3/bin")
JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
JuliaCall::julia_eval("using DataFrames")

#winit <- initialize_estimates(y = ysim, wfixed = wfixed)

tmpf <- function(ys){
  winit <- initialize_estimates(y = ys, wfixed = wfixed)
fit <- lbfgs::lbfgs(
  calc_kf_nll,
  calc_kf_grad,
  x = x,
  betasd = 0.0001,
  epsilon = 1e-4,
  max_iterations = 1e2,
  a = 0.9,
  y = ys,
  pm = param_map,
  winit,
  invisible = 0
)
fit
}

fits <- lapply(ysim, tmpf)

system.time(fit <- tmpf(ysim[[1]]))

dets <- kf_nll_details(w=fit$par, x=x, y=ysim[[1]], betasd = .0001, a = 0.9, pm = param_map, fet = NULL)

par(mfrow = c(3, 1))
qqnorm(dets$ytilde_k[1, ] / sqrt(dets$S[1, 1, ]), sub = "Cases")
abline(0, 1)
qqnorm(dets$ytilde_k[2, ] / sqrt(dets$S[2, 2, ]), sub = "Hospitalizations")
abline(0, 1)
qqnorm(dets$ytilde_k[3, ] / sqrt(dets$S[3, 3, ]), sub = "Deaths")
abline(0, 1)

rho_t <- detect_frac(365.25 * (x$time - 2020.164))
plot(x$time, ysim[[1]]$cases, xlab = "Time", ylab = "Cases")
pred_cases <- dets$xhat_kkmo["C", ] * rho_t + dets$xhat_kkmo["Hnew", ]
est_cases <- dets$xhat_kkmo["C", ] + dets$xhat_kkmo["Hnew", ]
se_cases <- sqrt(dets$S[1, 1, ])
lines(x$time, se_cases * 2 + pred_cases, col = "grey")
lines(x$time, pred_cases)
lines(x$time, est_cases, lty = 2)
lines(x$time,-se_cases * 2 + pred_cases, col = "grey")

plot(x$time, ysim[[1]]$hospitalizations, xlab = "Time",
     ylab = "Hospitalizations", ylim = c(0, 100))
pred_hosps <- dets$xhat_kkmo["Hnew", ]
se_hosps <- sqrt(dets$S[2, 2, ])
lines(x$time, se_hosps * 2 + pred_hosps, col = "grey")
lines(x$time, pred_hosps)
lines(x$time,-se_hosps * 2 + pred_hosps, col = "grey")

plot(x$time, ysim[[1]]$deaths, xlab = "Time",
     ylab = "Deaths", ylim = c(0, 30))
pred_deaths <- dets$xhat_kkmo["Drep", ]
se_deaths <- sqrt(dets$S[3, 3, ])
lines(x$time, se_deaths * 2 + pred_deaths, col = "grey")
lines(x$time, pred_deaths)
lines(x$time,-se_deaths * 2 + pred_deaths, col = "grey")

