#!/usr/bin/env R

tictoc::tic()

suppressPackageStartupMessages(library(tidyverse))

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

forecast_date <- Sys.getenv("fdt", unset = "2020-10-12")
forecast_loc <- Sys.getenv("loc", unset = "36")
data_fname <- paste0("data/", forecast_date, "--", forecast_loc, ".csv")
case_data <- read_csv(data_fname)

target_end_dates <- max(case_data$target_end_date) + lubridate::ddays(1) * seq(1, 28)
target_end_times <- lubridate::decimal_date(target_end_dates)
target_wday <- lubridate::wday(target_end_dates)

res_fname <- paste0("minimizer/", forecast_date, "--", forecast_loc, ".csv")
pvar_df <- read_csv(res_fname)

wsize <- nrow(pvar_df) - 2L
the_t0 <- rev(case_data$time)[wsize + 1]

iterate_f_and_P <- function(xhat, PN, pvec, beta_t, time.steps, vif = 10){
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
        c(0, beta_t * (S / N) * (I / N) * vif, eta * E / N, gamma * I / N)
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
    p <- c(pvar, pfixed)
    
    is_spline_par <- grepl("^b[0-9]+$", names(p))
    tmp <- p[is_spline_par]
    bpars <- tmp[order(as.integer(str_remove(names(tmp), "^b")))]
    stopifnot(length(bpars) == nrow(cdata))
    
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
    diff_ar <- (bpars[-1] - p["gamma"]) - p["a"] * (bpars[-length(bpars)] - p["gamma"])
    stat_sd = sqrt(p["betasd"]^2 / (1 - p["a"]^2))
    
    rwlik <-
      dnorm(bpars[1], mean = p["gamma"], sd = stat_sd, log = TRUE) + sum(dnorm(diff_ar, sd = p["betasd"], log = TRUE))
    nll <-
      0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi)) - rwlik
    if (!just_nll) {
      if (!is.null(fets)) {
        sim_means <-
          sim_cov <- matrix(NA, nrow = (nrow(fets)), ncol = nsim)
        for (j in seq_len(nsim)) {
          bpars_fet <- numeric(nrow(fets))
          bpars_fet[1] <- p["gamma"] + (bpars[T] - p["gamma"] )* p["a"] + rnorm(n = 1, sd = p["betasd"])
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

pvar <- pvar_df$minimizer
names(pvar) <- pvar_df$par

pfixed <- unlist(read_csv(paste0("initial-pars/", forecast_loc, ".csv"))[1,])

pfixed <- c(
  pfixed, 
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

fet <- tibble(target_end_times, target_wday, target_end_dates)

kfret <- kfnll(pvar = pvar,
      pfixed = pfixed,
      cdata = tail(case_data, n = wsize),
      t0 = the_t0,
      just_nll = FALSE,
      fet = fet, 
      fet_zero_cases = "weekly")

inds <- which(fet$target_wday == 7)

fcst <- create_forecast_df(means = kfret$sim_means[inds,],
                           vars = kfret$sim_cov[inds,],
                           location = forecast_loc)

stopifnot(setequal(fet$target_end_dates[inds], fcst$target_end_date %>% unique()))

fcst_path <- file.path("forecasts", paste0(forecast_date, "-", forecast_loc, "-CEID-InfectionKalman.csv"))
if(!dir.exists("forecasts")) dir.create("forecasts")
write_csv(x = fcst, path = fcst_path)

q("no")

fcst %>% ggplot(aes(x = target_end_date, y = as.numeric(value), color = quantile)) + geom_line() + geom_point()


kfret2 <- kfnll(pvar = pvar,
               pfixed = pfixed,
               cdata = tail(case_data, n = wsize),
               t0 = the_t0,
               just_nll = FALSE,
               fet = NULL, 
               fet_zero_cases = "daily")

par(mfrow = c(1,1))
qqnorm(kfret2$ytilde_k / sqrt(kfret2$S))
abline(0, 1)

par(mfrow = c(4, 1))
tgrid <- tail(case_data$time, n = wsize)

plot(tgrid, tail(case_data$smooth, n = wsize), xlab = "Time", ylab = "Cases")
lines(tgrid, kfret$xhat_kkmo["C",] * pfixed["rho1"])

plot(tgrid, kfret$S, log = "y", xlab = "Time", ylab = "Variance in smoother")
plot(tgrid, kfret$ytilde_k, xlab = "Time", ylab = "Residual in process 1-ahead prediction")
plot(tgrid, kfret$ytilde_k / sqrt(kfret$S), xlab = "Time", ylab = "Standardized residual")

par(mfrow = c(4, 1))

plot(tgrid, pvar[-c(1,2)] / pfixed["gamma"], xlab = "Time", ylab = expression(R[t]))

plot(tgrid, kfret$xhat_kkmo["E",], xlab = "Time", ylab = "Predicted exposed")
points(tgrid, kfret$xhat_kk["E", ], xlab = "Time", col = 2)

plot(tgrid, kfret$xhat_kkmo["I",], xlab = "Time", ylab = "Predicted infected")
points(tgrid, kfret$xhat_kk["I", ], xlab = "Time", col = 2)

plot(tgrid, kfret$xhat_kkmo["S",], xlab = "Time", ylab = "Predicted susceptible")
points(tgrid, kfret$xhat_kk["S", ], xlab = "Time", col = 2)


