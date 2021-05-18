state_abb_fips <-
  readr::read_csv(
    file = "state,state_code,state_name
AK,02,Alaska
AL,01,Alabama
AR,05,Arkansas
AS,60,American Samoa
AZ,04,Arizona
CA,06,California
CO,08,Colorado
CT,09,Connecticut
DC,11,District of Columbia
DE,10,Delaware
FL,12,Florida
GA,13,Georgia
GU,66,Guam
HI,15,Hawaii
IA,19,Iowa
ID,16,Idaho
IL,17,Illinois
IN,18,Indiana
KS,20,Kansas
KY,21,Kentucky
LA,22,Louisiana
MA,25,Massachusetts
MD,24,Maryland
ME,23,Maine
MI,26,Michigan
MN,27,Minnesota
MO,29,Missouri
MP,69,Northern Mariana Islands
MS,28,Mississippi
MT,30,Montana
NC,37,North Carolina
ND,38,North Dakota
NE,31,Nebraska
NH,33,New Hampshire
NJ,34,New Jersey
NM,35,New Mexico
NV,32,Nevada
NY,36,New York
OH,39,Ohio
OK,40,Oklahoma
OR,41,Oregon
PA,42,Pennsylvania
PR,72,Puerto Rico
RI,44,Rhode Island
SC,45,South Carolina
SD,46,South Dakota
TN,47,Tennessee
TX,48,Texas
UM,74,U.S. Minor Outlying Islands
UT,49,Utah
VA,51,Virginia
VI,78,Virgin Islands
VT,50,Vermont
WA,53,Washington
WI,55,Wisconsin
WV,54,West Virginia
WY,56,Wyoming"
  )
covidhub_locations <- c("District of Columbia", state.name, "All")

process_hopkins <- function(path, weekly = TRUE){
  is_deaths <- str_detect(path, "deaths_US.csv$")
  if (is_deaths){
    targ_suffix <- " death"
    if (weekly){
      targ_pattern <- "^wk ahead inc|^wk ahead cum"
    } else {
      targ_pattern <- "^day ahead inc|^day ahead cum"
    }
  } else {
    targ_suffix <- " case"
    if (weekly){
      targ_pattern <- "^wk ahead inc"
    } else {
      targ_pattern <- "^day ahead inc"
    }
  }

  ## construct a col spec to read in only the state and date columns
  colnms <- readLines(path, n = 1) %>% str_split(",") %>% "[["(1)
  is_date_col <- colnms %>% str_detect("^\\d{1,2}/\\d{1,2}/\\d{2}$")
  date_cols <- colnms[is_date_col]
  colspec <- sapply(date_cols, function(x)
    "i")
  col_types <-
    do.call(cols_only, c(list(FIPS = col_double(), 
                              Province_State = col_character()), 
                         as.list(colspec)))
  
  hpd_raw <-
    read_csv(path, col_types = col_types) %>%
    pivot_longer(-c(Province_State, FIPS), names_to = "date_string", 
                 values_to = "day ahead cum") %>%
    mutate(target_end_date = lubridate::mdy(date_string)) %>%
    mutate(location = sprintf("%05d", FIPS)) %>%
    mutate(location = str_remove(location, "^000")) %>%
    select(-FIPS)
  
  hpd_county <- hpd_raw %>% filter(location !="   NA") %>%
    filter(!str_detect(location, "^900|^999|^800|^888"))
  
  hpd_state <- hpd_raw %>%  group_by(Province_State, target_end_date) %>%
    mutate(has_id = str_detect(location, "^[0-9][0-9]")) %>%
    summarise(`day ahead cum` = sum(`day ahead cum`),
              location = substring(location[has_id][1], 1, 2), 
              n = n(), .groups = "drop") %>%
    filter(n > 1) %>%
    select(-n)
  
  hpd_us <- hpd_raw %>%  group_by(target_end_date) %>%
    summarise(`day ahead cum` = sum(`day ahead cum`),
              location = "US", .groups = "drop")
  
  hpd <- bind_rows(hpd_county, hpd_state, hpd_us) %>% 
    arrange(location, target_end_date) %>%
    group_by(location) %>% 
    mutate(`day ahead inc`= c(NA, diff(`day ahead cum`)))    
  
  if(weekly){
  hpd2 <- 
    hpd %>% mutate(week = lubridate::epiweek(target_end_date)) %>% 
    group_by(location, week) %>% 
    arrange(location, week, target_end_date) %>% 
    summarise(`wk ahead inc` = sum(`day ahead inc`),
              `wk ahead cum` = tail(`day ahead cum`, n = 1),
              target_end_date = tail(target_end_date, n = 1), 
              n = n(), .groups = "drop") %>% 
    filter(n == 7) %>%
    select(-n, -week)
  } else {
    hpd2 <- hpd %>% select(location, `day ahead inc`, `day ahead cum`, 
                           target_end_date)
  }
  
  hpd3 <- 
    hpd2 %>% pivot_longer(-c(target_end_date, location), 
                          names_to = "target_type", values_to = "value") %>%
    mutate(target_type = str_c(target_type, targ_suffix)) %>%
    filter(str_detect(target_type, targ_pattern)) %>%
    filter(!str_detect(location, "72[0-9]{3}")) %>% #remove PR counties
    filter(location != "11001") # remove duplicated location of DC as county
  
  hpd3 
}

load_google_mobility <- function(issue_date, rootdir = ".") {
  if (issue_date > "2021-12-31") {
    stop("reading of reports issued after 2021 is not implemented")
  }
  path <- file.path(
    rootdir, 
    "google-mobility-reports",
    issue_date,
    c(
      "2020_US_Region_Mobility_Report.csv",
      "2021_US_Region_Mobility_Report.csv"
    )
  )
  rcsv <- function(p) {
    read_csv(
      p,
      col_types = cols_only(
        country_region = col_character(),
        sub_region_1 = col_character(),
        sub_region_2 = col_character(),
        iso_3166_2_code = col_character(),
        date = col_date(format = "%Y-%m-%d"),
        residential_percent_change_from_baseline = col_double()
      )
    )
  }
  purrr::map(path, rcsv) %>% bind_rows() %>% filter(is.na(sub_region_2)) %>%
    select(-sub_region_2)
}

load_hopkins <- function(loc, weekly = TRUE) {
  # load and clean data
  hpp <- dir(loc, full.names = TRUE)
  hpd <- str_subset(hpp, "time_series_covid19_deaths_US.csv$")
  hpc <- str_subset(hpp, "time_series_covid19_confirmed_US.csv$")
  
  ddat <- process_hopkins(hpd, weekly = weekly)
  cdat <- process_hopkins(hpc, weekly = weekly)
  bind_rows(ddat, cdat)
}

read_forecast <- function(file) {
  read_csv(
    file,
    col_types =
      cols(
        forecast_date = col_date(format = "%Y-%m-%d"),
        target = col_character(),
        target_end_date = col_date(format = "%Y-%m-%d"),
        location = col_character(),
        type = col_character(),
        quantile = col_double(),
        value = col_double()
      )
  )
}

quant <- function(x, p){
  quantile(x, prob = p, names = FALSE, type = 8, na.rm = TRUE)
}

vardf <- function(var, samp){
  cname <- switch(var,
                  "inc hosp" = "hosps",
                  "inc death" = "deaths",
                  "inc case" = "cases",
                  "cum death" = "cum_deaths")
  if (var != "inc case") {
    prob <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  } else {
    prob <- c(0.025, 0.100, 0.250, 0.500, 0.750, 0.900, 0.975)
  }
  t1 <- tibble(
    target = var,
    type  = "quantile",
    quantile = prob,
    value = quant(samp[[cname]], prob)
  )
  t2 <-
    tibble(
      target = var,
      type = "point",
      quantile = NA,
      value = quant(samp[[cname]], 0.5)
      
    )
  bind_rows(t1, t2)
}

samp_to_df <-
  function(sampdf,
           vars = c("inc hosp", "inc death", "cum death", "inc case")) {
  purrr::map_dfr(vars, vardf, samp = sampdf)
}

# Take simulation trajectories and output a data frame in the format described
# here: https://github.com/reichlab/covid19-forecast-hub/blob/6a7e5624ef540a55902770b7c17609d19e1f593a/data-processed/README.md
paths_to_forecast <- function(out, loc = "13", wks_ahead = 1:6, hop, fdt) {
  if(any(wks_ahead > 20)){
    stop("Max weeks ahead accepted is 20", .call = FALSE)
  }

  out2 <- 
    out %>% group_by(Rep) %>% arrange(Date) %>% 
    mutate(cum_deaths = cumsum(deaths),
           day_diff = c(NA, diff(Date)))
  
  prior_deaths <- hop %>% 
    filter(target_end_date < fdt & location == loc & 
             target_type == "wk ahead cum death") %>%
    arrange(target_end_date) %>%
    pull("value") %>%
    tail(n = 1)
  out2$cum_deaths <- out2$cum_deaths + prior_deaths
  
  take_back_step <- lubridate::wday(fdt, label = TRUE) %in% c("Sun", "Mon")
  fdt_yr <- lubridate::year(fdt)
  fdt_wk <- lubridate::epiweek(fdt)
  fdt_sun <- MMWRweek::MMWRweek2Date(fdt_yr, fdt_wk)
  if (take_back_step){
    week0_sun <- fdt_sun - 7
  } else {
    week0_sun <- fdt_sun
  }
  forecast_epiweeks <- (week0_sun + wks_ahead * 7) %>% lubridate::epiweek()
  if(any(na.omit(out2$day_diff) > 7)){
    stop("Non-continuous series of weeks, unable to compute cumulative forecasts", 
         .call = FALSE)
  }

  weekly <- out2 %>%
    as_tibble() %>%
    mutate(epiweek = lubridate::epiweek(Date)) %>%
    filter(epiweek %in% forecast_epiweeks) %>%
    rename("target_end_date" = Date) %>%
    nest(data = c(Rep, cases, deaths, cum_deaths)) %>%
    mutate(pred_df = purrr::map(data, samp_to_df, 
                                vars = c("inc death", "cum death", "inc case"))) %>%
    select(-data) %>%
    unnest(pred_df) %>%
    mutate(target = paste((target_end_date - (week0_sun + 6)) / lubridate::ddays(7), 
                          "wk ahead", target)) %>%
    add_column(location = loc) %>%
    mutate(quantile = round(quantile, digits = 3),
           value = round(value, digits = 3)) %>%
    add_column(forecast_date = fdt) %>%
    select(forecast_date, target, target_end_date, location, type, quantile, 
           value)

  weekly %>% 
    filter(!is.na(value)) %>%
    filter((nchar(location)) <= 2 | str_detect(target, "inc case$")) ## only case forecasts accepted for counties
}

julia_assign2 <- function(w, cov, z, wfixed, p_hsd, β_0sd,  τ_csd){
  JuliaCall::julia_assign("w", w)
  JuliaCall::julia_assign("η", wfixed["η"])
  JuliaCall::julia_assign("N", wfixed["N"])
  JuliaCall::julia_assign("γ", wfixed["γ"])
  JuliaCall::julia_assign("cov", cov)
  JuliaCall::julia_assign("z", z)
  JuliaCall::julia_assign("p_h", wfixed["p_h"])
  JuliaCall::julia_assign("τ_h", wfixed["τ_h"])
  JuliaCall::julia_assign("p_hsd", p_hsd)
  JuliaCall::julia_assign("β_0sd", β_0sd)
  JuliaCall::julia_assign("τ_csd", τ_csd)  
  JuliaCall::julia_assign("γ_d", wfixed["γ_d"])
  JuliaCall::julia_assign("γ_h", wfixed["γ_h"])
}

julia_assign3 <- function(w, cov, z, wfixed, p_hsd, β_0sd, τ_csd){
  JuliaCall::julia_assign("w", w)
  JuliaCall::julia_assign("η", wfixed["η"])
  JuliaCall::julia_assign("N", wfixed["N"])
  JuliaCall::julia_assign("γ", wfixed["γ"])
  JuliaCall::julia_assign("cov", cov)
  JuliaCall::julia_assign("z", z)
  JuliaCall::julia_assign("p_hsd", p_hsd)
  JuliaCall::julia_assign("β_0sd", β_0sd)
  JuliaCall::julia_assign("τ_csd", τ_csd)
  JuliaCall::julia_assign("γ_d", wfixed["γ_d"])
  JuliaCall::julia_assign("γ_h", wfixed["γ_h"])
  
}

ncol2tuple <- "(w, cov, z; N = N, η = η, γ = γ, p_h = p_h, γ_h = γ_h, γ_d = γ_d,  p_hsd = p_hsd, β_0sd = β_0sd, τ_csd = τ_csd, just_nll = true)"
ncol3tuple <- "(w, cov, z; N = N, η = η, γ = γ, γ_h = γ_h, γ_d = γ_d, p_hsd = p_hsd, β_0sd = β_0sd, τ_csd = τ_csd, just_nll = true)"

calc_kf_nll <- function(w, cov, z,  p_hsd,  β_0sd,   τ_csd, wfixed) {
  if (ncol(z) == 3) {
    julia_assign3(w, cov, z, wfixed,  p_hsd,  β_0sd,   τ_csd)
    nll <-
      JuliaCall::julia_eval(paste0("InfectionKalman.obj", ncol3tuple))
  } else if (ncol(z) == 2 && !"hospitalizations" %in% names(z)) {
    julia_assign2(w, cov, z, wfixed,  p_hsd,  β_0sd,   τ_csd)
    nll <- JuliaCall::julia_eval(paste0("InfectionKalman.obj", ncol2tuple))
  }
  nll
}

calc_kf_grad <- function(w, cov, z, p_hsd,  β_0sd,   τ_csd, wfixed) {
  if (ncol(z) == 3) {
    julia_assign3(w, cov, z, wfixed, p_hsd,  β_0sd,   τ_csd)
    g <- JuliaCall::julia_eval(paste0("InfectionKalman.grad", ncol3tuple))
  } else if (ncol(z) == 2 && !"hospitalizations" %in% names(z)) {
    julia_assign2(w, cov, z, wfixed, p_hsd,  β_0sd,   τ_csd)
    g <- JuliaCall::julia_eval(paste0("InfectionKalman.grad", ncol2tuple))
  }
  g
}

calc_kf_hess <- function(w, cov, z, p_hsd,  β_0sd,   τ_csd, wfixed) {
  if (ncol(z) == 3) {
    julia_assign3(w, cov, z, wfixed, p_hsd,  β_0sd,   τ_csd)
    g <- JuliaCall::julia_eval(paste0("InfectionKalman.hess", ncol3tuple))
  } else if (ncol(z) == 2 && !"hospitalizations" %in% names(z)) {
    julia_assign2(w, cov, z, wfixed, p_hsd,  β_0sd,   τ_csd)
    g <- JuliaCall::julia_eval(paste0("InfectionKalman.hess", ncol2tuple))
  }
  g
}

initialize_estimates <- function(x, y, wfixed, dt = 0.00273224) {
  τ_c_init <- max(var(y$cases, na.rm = TRUE) / mean(y$cases, na.rm = TRUE), 1)
  nτ_c <- tail(x$τ_cmap, n = 1)
  np_h <- tail(x$p_hmap, n = 1)
  
  τ_d_init <- max(var(y$deaths, na.rm = TRUE), 1) #/ mean(y$deaths)
  wsize <- nrow(y)

  est_true_cases_hosps <- y$cases / x$ρ
  Yest <- est_true_cases_hosps / (wfixed["γ"] * dt)
  
  Y0init <- mean(head(Yest, n = 7), na.rm = TRUE)
  L0init <- (Y0init * wfixed["γ"] / wfixed["η"])  %>% unname()
  
  Hest <- Yest  * wfixed["γ"] * wfixed["p_h"] / wfixed["γ_h"]
  
  p_d_init <-
    max(sum(y$deaths, na.rm = TRUE) / sum(Hest, na.rm = TRUE),
        0.01)
  
  l <- 7 
  Rt <- lead(y$cases / x$ρ, n = l) / (y$cases / x$ρ)
  Rt[Rt > 4] <- 4
  Rt[Rt < 0] <- 0.1
  m <- 10
  Rtfilt <- signal::sgolayfilt(na.omit(Rt), p = 2, n = 2 * m  + 1) # edge values seem biased

  Rtfilt[1:m] <- NA
  lf <- length(Rtfilt)
  Rtfilt[(lf - m + 1):lf] <- NA
  Rtfilt[Rtfilt < 0] <- 0.1
  
  H0init <- wfixed["p_h"] * wfixed["γ"] * Y0init / wfixed["γ_h"]
  D0 <- H0init * wfixed["γ_h"] / wfixed["γ_d"] * p_d_init
  X0init <- wfixed["N"] - L0init - Y0init - H0init - D0
  Xdecrement <- cumsum(lead(y$cases / x$ρ, n = l))
  Xt <- X0init - Xdecrement
  β <- Rtfilt * wfixed["γ"] * wfixed["N"] / na.omit(Xt)
  β[1:m] <- β[m + 1]
  β[(lf - m + 1):lf] <- β[lf - m]
  β2 <- c(β, rep(β[lf - m], l))
  
  df <- data.frame(logβ = log(β2), residential = x$residential)
  mod <- lm(logβ ~ residential, data = df)
  print(summary(mod))
  residentialeffect_init <- coef(mod)["residential"] %>% unname()
  intercept_init <- coef(mod)["(Intercept)"] + residuals(mod)
  intercept_split <- split(intercept_init, x$β_0map)
  β_0 <- sapply(intercept_split, mean)
  names(β_0) <- paste0("β_0_", seq_along(β_0))

  winit <- c(
    logL0 = log(L0init + 1),
    logτ_d = log(τ_d_init),
    residentialeffect = residentialeffect_init,
    logitp_d = qlogis(p_d_init),
    logγ_d12 = log(365.25 / 10),
    logγ_d34 = log(365.25 / 10),
    logγ_z17 = log(365.25 / 1),
    logτ_c = rep(log(τ_c_init), times = nτ_c),
    β_0
  )
  winit
}

calc_kf_nll_r <-
  function(w,
           cov,
           z,
           p_hsd,
           β_0sd,
           τ_csd,
           wfixed,
           fets = NULL,
           fet_zero_cases_deaths = "daily",
           just_nll = TRUE) {
    diffeqr::diffeq_setup("/opt/julia-1.5.3/bin")
    JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
    
    η  <- unname(wfixed["η"])
    γ  <- unname(wfixed["γ"])
    N <- unname(wfixed["N"])
    γ_h <- unname(wfixed["γ_h"])
    γ_d <- unname(wfixed["γ_d"])
    γ_z <- unname(wfixed["γ_z"])
    nτ_c <- tail(cov$τ_cmap, 1)
    np_h <- tail(cov$p_hmap, 1)
    if (ncol(z) == 3) {
      τ_h <- exp(w[2])
      if ("doses_scaled" %in% names(cov)){
        doseeffect <- w[5]
        w2 <- c(w[1], w[3:4], w[6:length(w)])
      } else{
        w2 <- c(w[1], w[3:length(w)])
      }
      p_h <- plogis(w2[(8 + nτ_c):(8 + nτ_c + np_h - 1)])
    } else if(ncol(z) == 2 && !"hospitalizations" %in% names(z)) {
      τ_h <- wfixed["τ_h"]
      p_h <- wfixed["p_h"]
      np_h <- 0
      z$hospitalizations <- NA
      w2 <- w
    }
    L0 <- exp(w2[1])
    τ_d <- exp(w2[2])
    residentialeffect <- w2[3]
    p_d <- plogis(w2[4])
    γ_d12 <- exp(w2[5])
    γ_d34 <- exp(w2[6])
    γ_z17 <- exp(w2[7])
    τ_c <- exp(w2[seq(8, 8 + nτ_c - 1)])
    β_0 <- w2[seq(8 +  nτ_c + np_h, length(w2))]
    
    zloc <-
      data.matrix(z[, c("cases", "hospitalizations", "deaths")])
    zmiss <- is.na(zloc)
    
    Y0 <- L0 * η / γ
    H0 <- p_h[1] * γ * Y0 / γ_h
    D0 <- H0 *  γ_h /  γ_d * p_d

    x0 <- c(
      max(N - L0 - Y0 - H0 - D0, N * 0.1),
      #X
      min(Y0, N),
      #Y
      min(L0, N),
      #L
      min(cov$ρ[1] * Y0 *  γ  * (1 - p_h[cov$p_hmap[1]]) /  γ_z, N),
      #Z
      0,
      #Z_r
      min(H0, N),
      #H
      0,
      #A
      min(D0, N),
      #D
      0
    ) #D_r
    
    names(x0) <- c("X", "Y", "L", "Z", "Z_r", "H", "A", "D", "D_r")
    
    p0 <- diag(c(1, 1, 1, 1, 0, 1, 0, 1, 0))
    
    dstate = length(x0)
    dobs = ncol(zloc)
    
    nobs <- nrow(zloc)
    stopifnot(nobs > 0)
    
    ytilde_kk <- ytilde_k <- array(NA_real_, dim = c(dobs, nobs))
    S <- array(NA_real_, dim = c(dobs, dobs, nobs))
    K <- array(NA_real_, dim = c(dstate, dobs, nobs))
    rdiagadj <- array(1, dim = c(dobs, nobs))
    
    xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(dstate, nobs))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(x0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(dstate, dstate, nobs))
    
    hmat <-
      rbind(c(0, 0, 0, 0, 1, 0, 0, 0, 0),
            c(0, 0, 0, 0, 0, 0, 1, 0, 0),
            c(0, 0, 0, 0, 0, 0, 0, 0, 1))
    JuliaCall::julia_eval("vectorfield = InfectionKalman.genvectorfield()")
    zerovars <- c(5, 7, 9)
    
    for (i in seq(1, nobs)) {
      if (i == 1) {
        xhat_init <- x0
        PNinit <- p0
      } else {
        xhat_init <- xhat_kk[, i - 1]
        PNinit <- P_kk[, , i - 1]
      }
      for (zv in zerovars) {
        xhat_init[zv] <- 0
        PNinit[, zv] <- PNinit[zv,] <- 0
      }
      
      u0 <- cbind(xhat_init, PNinit)
      if ("doses_scaled" %in% names(cov)) {
        β   <-
          exp(β_0[cov$β_0map[i]] + residentialeffect * cov$residential[i] + doseeffect * cov$doses_scaled[i])
      } else {
        β   <-
          exp(β_0[cov$β_0map[i]] + residentialeffect * cov$residential[i])
      }
      if(β > 1000){
        β <- 1000
      }
      par <- c(β, N,  η,  γ,  γ_d,  γ_z,  γ_h, p_h[cov$p_hmap[i]], p_d, cov$ρ[i])
      if (cov$wday[i] == 1) {
        par[5] <-  γ_d12
        par[6] <- γ_z17
      } else if (cov$wday[i] == 2){
        par[5] <-  γ_d12
      } else if (cov$wday[i] %in% c(3, 4)) {
        par[5] <-  γ_d34
      } else if (cov$wday[i] == 7){
        par[6] <- γ_z17
      }
      JuliaCall::julia_assign("u0", u0)
      JuliaCall::julia_assign("tspan", c(0, 0.00273224))
      JuliaCall::julia_assign("par", par)
      JuliaCall::julia_eval("prob = ODEProblem(vectorfield,u0,tspan,par)")
      XP <-
        JuliaCall::julia_eval("solve(prob, Tsit5(), saveat = tspan[2]).u[2]")
      xhat_kkmo[, i] <- XP[, 1]
      P_kkmo[, , i] <- XP[,-1]
      
      for (j in 1:dstate) {
        if (P_kkmo[j, j, i] < 0) {
          P_kkmo[j, , i] <- 0
          P_kkmo[, j, i] <- 0
        }
      }
      r <- diag(c(τ_c[cov$τ_cmap[i]] * xhat_init[3],  τ_h,   τ_d))
      ytilde_k[, i] <-
        matrix(zloc[i,], ncol = 1) - hmat %*% xhat_kkmo[, i, drop = FALSE]
      S[, , i] <- hmat %*% P_kkmo[, , i] %*% t(hmat) + r
      
      for (j in 1:dobs) {
        S[j, j, i] <- S[j, j, i] + rdiagadj[j, i]
      }
      K[, , i] <-
        P_kkmo[, , i] %*% t(hmat) %*% solve(S[, , i])
      desel <- zmiss[i,]
      K[, desel, i] <- 0
      
      xhat_kk[, i] <-
        xhat_kkmo[, i, drop = FALSE] +
        K[,!desel, i] %*% ytilde_k[!desel, i, drop = FALSE]
      
      xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
      P_kk[, , i] <-
        (diag(dstate) - K[, , i] %*% hmat) %*% P_kkmo[, , i]
      
      for (j in 1:dstate) {
        if (P_kk[j, j, i] < 0) {
          P_kk[j, , i] <- 0
          P_kk[, j, i] <- 0
        }
      }
      ytilde_kk[, i] <- matrix(zloc[i,], ncol = 1) -
        hmat %*% xhat_kk[, i, drop = FALSE]
    }
    
    rwlik <- 0
    for (i in seq_len(length(β_0) - 1)) {
      step <-  β_0[i + 1] -  β_0[i]
      rwlik <-
        rwlik + dnorm(step,
                      mean = 0,
                      sd =  β_0sd,
                      log = TRUE)
    }
    
    for (i in seq_len(length(τ_c) - 1)) {
      τ_cstep <- log(τ_c[i + 1]) - log(τ_c[i])
      rwlik <-
        rwlik + dnorm(τ_cstep,
                        mean = 0,
                        sd =  τ_csd,
                        log = TRUE)
    }
    
    for (i in seq_len(length(p_h) - 1)) {
      p_hstep <- qlogis(p_h[i + 1]) - qlogis(p_h[i])
      rwlik <-
        rwlik + dnorm(p_hstep,
                        mean = 0,
                        sd =  p_hsd,
                        log = TRUE)
    }
    
    nll <- 0
    for (i in seq(1, nobs)) {
      sel <- !zmiss[i,]
      nll <- nll +
        t(ytilde_k[sel, i]) %*% solve(S[sel, sel, i]) %*% ytilde_k[sel, i] +
        log(det(S[, , i][sel, sel, drop = FALSE])) + dobs * log(2 * pi)
    }
    nll <- 0.5 * nll - rwlik
    
    if (!just_nll) {
      if (!is.null(fets)) {
        nsimdays <- nrow(fets)
        sim_means <- array(NA_real_, dim = c(dobs, nsimdays))
        sim_cov <- array(NA_real_, dim = c(dobs, dobs, nsimdays))
        
        x_sim <- array(NA_real_, dim = c(dstate, nsimdays))
        rownames(x_sim) <- names(x0)
        P_sim <- array(NA_real_, dim = c(dstate, dstate, nsimdays))
        
        
        for (j in seq_len(nsimdays)) {
          if (j == 1) {
            xhat_init <- xhat_kk[, nobs]
            PNinit <- P_kk[, , nobs]
          } else {
            xhat_init <- XP[, 1]
            PNinit <- XP[,-1]
          }
          
          if (fet_zero_cases_deaths == "daily" ||
              fets$target_wday[j] == 1) {
            xhat_init[5] <- 0
            PNinit[, 5] <- PNinit[5,] <- 0
            xhat_init[9] <- 0
            PNinit[, 9] <- PNinit[9,] <- 0
          }
          xhat_init[7] <- 0
          PNinit[, 7] <- PNinit[7,] <- 0
          
          u0 <- cbind(xhat_init, PNinit)
          
          if ("doses_scaled" %in% names(cov)) {
            β   <-
              exp(β_0[cov$β_0map[i]] + residentialeffect * cov$residential[i] + doseeffect * cov$doses_scaled[i])
          } else {
            β   <-
              exp(β_0[cov$β_0map[i]] + residentialeffect * cov$residential[i])
          }
          if(β > 1000){
            β <- 1000
          }
          par <- c(β, N,  η,  γ,  γ_d,  γ_z,  γ_h, p_h, p_d, cov$ρ[i])
          if (cov$wday[i] %in% c(1, 2)) {
            par[5] <-  γ_d12
          } else if (cov$wday[i] %in% c(3, 4)) {
            par[5] <-  γ_d34
          }
          JuliaCall::julia_assign("u0", u0)
          JuliaCall::julia_assign("tspan", c(0, 0.00273224))
          JuliaCall::julia_assign("par", par)
          JuliaCall::julia_eval("prob = ODEProblem(vectorfield,u0,tspan,par)")
          XP <-
            JuliaCall::julia_eval("solve(prob, Tsit5(), saveat = tspan[2]).u[2]")
          
          r <-
            diag(c(τ_c[cov$τ_cmap[nobs]] * xhat_init[3],  τ_h,   τ_d))
          
          sim_means[, j] <- hmat %*% XP[, 1]
          sim_cov[, , j] <- hmat %*% XP[, -1] %*% t(hmat) + r
          
          x_sim[, j] <- XP[, 1]
          P_sim[, , j] <- XP[,-1]
          
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
        x_sim = x_sim,
        P_sim = P_sim,
        rdiagadj = rdiagadj
      )
    } else {
      nll
    }
  }
