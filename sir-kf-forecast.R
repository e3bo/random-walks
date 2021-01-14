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

forecast_date <- Sys.getenv("fdt", unset = "2020-10-12")
forecast_loc <- "36"
hopdir <- file.path("hopkins", forecast_date)
tictoc::tic("data loading")
tdat <- load_hopkins(hopdir, weekly = FALSE) 
tictoc::toc()

nys <- tdat %>% filter(location == forecast_loc) %>% 
  filter(target_type == "day ahead inc case")

nys2 <- nys %>% mutate(time = lubridate::decimal_date(target_end_date))
case_data <- nys2 %>% ungroup() %>% 
  mutate(wday = lubridate::wday(target_end_date)) %>%
  select(time, wday, value) %>% 
  rename(reports = value)

target_end_dates <- max(nys2$target_end_date) + lubridate::dweeks(1:4)
target_end_times <- lubridate::decimal_date(target_end_dates)

iterate_f_and_P <- function(xhat, PN, pvec, beta_t, time.steps){
  P <- PN / pvec["N"]
  dt <- diff(time.steps)
  vf <-  with(as.list(c(pvec, xhat, beta_t)), {
      eta <- 365 / 4
      gamma <- 365 / 9
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
           fets = NULL) {
    p <- c(pvar, pfixed)
    
    is_spline_par <- grepl("^b[0-9]+$", names(p))
    tmp <- p[is_spline_par]
    bpars <- tmp[order(names(tmp))]
    stopifnot(length(bpars) == nrow(cdata) + 1)
    
    is_wday_par <- grepl("rho[1-7]", names(p))
    tmp <- p[is_wday_par]
    wpars <- tmp[order(names(tmp))]
    stopifnot(length(wpars) == 7)
    
    xhat0 = c(p[c("S_0", "E_0", "I_0")], 0)
    names(xhat0) <- c("S", "E", "I", "C")
    
    T <- nrow(cdata)
    stopifnot(T > 0)
    z <- cdata$reports
    times <- cdata$time
    wday <- cdata$wday
    
    ytilde_kk <- ytilde_k <- S <- array(NA_real_, dim = c(1, T))
    K <- xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(4, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(4, 4, T))
    Hfun <- function(w, par = p){
      name <- paste0("rho", w)
      matrix(c(0, 0, 0, par[name]), ncol = 4)
    }
    for (i in seq(1, T)) {
      if (i == 1) {
        xhat_init <- xhat0
        PNinit <- Phat0
        time.steps <- c(t0, times[1])
        R <- max(5, z[1] * p["tau"])
      } else {
        xhat_init <- xhat_kk[, i - 1]
        PNinit <- P_kk[, , i - 1]
        PNinit[, 4] <- PNinit[4, ] <- 0
        time.steps <- c(times[i - 1], times[i])
        R <- max(5, z[i - 1] * p["tau"])
      }
      xhat_init["C"] <- 0
      H <- Hfun(wday[i])
      XP <- iterate_f_and_P(
        xhat_init,
        PN = PNinit,
        pvec = p,
        beta_t = bpars[i + 1],
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
    
    rwlik <-
      sum(dnorm(diff(bpars), sd = p["beta_sd"], log = TRUE))
    nll <-
      0.5 * sum(ytilde_k ^ 2 / S + log(S) + log(2 * pi)) - rwlik
    if (!just_nll) {
      sim_means <-
        sim_cov <- matrix(NA, nrow = (length(fets)), ncol = nsim)
      if (!is.null(fets)) {
        for (j in seq_len(nsim)) {
          bpars_fet <-
            bpars[T + 1] + cumsum(rnorm(
              n = length(fets),
              mean = 0,
              sd = p["beta_sd"]
            ))
          xhat_init <- xhat_kk[, T]
          xhat_init["C"] <- 0
          PNinit <- P_kk[, , T]
          PNinit[, 4] <- PNinit[4, ] <- 0
          XP <- iterate_f_and_P(
            xhat_init,
            PN = PNinit,
            pvec = p,
            beta_t = bpars_fet[1],
            time.steps = c(times[T], fets[1])
          )
          sim_means[1, j] <- H %*% XP$xhat
          sim_cov[1, j] <- H %*% XP$PN %*% t(H)
          for (i in seq_along(fets[-1])) {
            xhat_init <- XP$xhat
            xhat_init["C"] <- 0
            PNinit <- XP$PN
            PNinit[, 4] <- PNinit[4, ] <- 0
            XP <-
              iterate_f_and_P(
                xhat_init,
                PN = PNinit,
                pvec = p,
                beta_t = bpars_fet[i + 1],
                time.steps = c(fets[i], fets[i + 1])
              )
            sim_means[i + 1, j] <- H %*% XP$xhat
            sim_cov[i + 1, j] <-
              H %*% XP$PN %*% t(H) + sim_means[i + 1, j] * p["tau"]
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

the_n <- 80
the_t0 <- rev(case_data$time)[the_n + 1]
gamma <- 365/9

pvar_df <- tribble(
  ~par, ~init, ~lower, ~upper,
  "beta_sd",  1, 0.01, 10,
  "E_0", 1e4, 10, 1e5,
  "I_0", 1e4, 10, 1e5,
  "tau", 1e-2, 1e-4, 1e2
) %>% 
  bind_rows(tibble(par = paste0("rho", seq(2, 7)),
            init = 0.4,
            lower = 0,
            upper = 1)) %>%
  bind_rows(tibble(par = paste0("b", seq_len(the_n + 1)),
                init = gamma,
                lower = 0.1 * gamma,
                upper = 4 * gamma))

pvar <- pvar_df$init
names(pvar) <- pvar_df$par

pfixed <- c(
  N = 20e6,
  S_0 = 19e6,
  rho1 = 0.4,
  iota = 2
)

tictoc::tic("optimization")
ans <- optim(
  par = pvar,
  fn = kfnll,
  pfixed = pfixed,
  cdata = tail(case_data, n = the_n),
  t0 = the_t0,
  method = "L-BFGS-B",
  lower = pvar_df$lower,
  upper = pvar_df$upper,
  control = list(trace = 1, maxit = 200)
)
tictoc::toc()

kfret <- kfnll(pvar = ans$par,
      pfixed = pfixed,
      cdata = tail(case_data, n = the_n),
      t0 = the_t0,
      just_nll = FALSE,
      fet = target_end_times)

fcst <- create_forecast_df(means = kfret$pred_means,
                           vars = kfret$pred_cov,
                           location = forecast_loc)

fcst_path <- file.path("forecasts", paste0(forecast_date, "-CEID-SIR_KF.csv"))
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

is_spline_par <- grepl("^b[0-9]+$", names(ans$par))
bhat <- ans$par[is_spline_par]
R0hat <- bhat / gamma

qqnorm(kfret$ytilde_k / sqrt(kfret$S))
abline(0, 1)

par(mfrow = c(4, 1))
tgrid <- tail(case_data$time, n = the_n)
plot(tgrid, tail(case_data$reports, n = the_n), xlab = "Time", ylab = "Cases")
lines(tgrid, kfret$xhat_kkmo["C",] * pfixed["rho"])
lines(tgrid, kfret$xhat_kk["C",] * pfixed["rho"], col = 2)

plot(tgrid, kfret$S, log = "y", xlab = "Time", ylab = "Variance in smoother")
plot(tgrid, kfret$ytilde_k, xlab = "Time", ylab = "Residual in process 1-ahead prediction")
plot(tgrid, kfret$ytilde_k / sqrt(kfret$S), xlab = "Time", ylab = "Standardized residual")

par(mfrow = c(3, 1))

plot(tgrid, kfret$xhat_kkmo["E",], xlab = "Time", ylab = "Predicted exposed")
points(tgrid, kfret$xhat_kk["E", ], xlab = "Time", col = 2)

plot(tgrid, kfret$xhat_kkmo["I",], xlab = "Time", ylab = "Predicted infected")
points(tgrid, kfret$xhat_kk["I", ], xlab = "Time", col = 2)

plot(tgrid, kfret$xhat_kkmo["S",], xlab = "Time", ylab = "Predicted susceptible")
points(tgrid, kfret$xhat_kk["S", ], xlab = "Time", col = 2)
