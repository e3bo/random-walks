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

load_hopkins <- function(loc, weekly = TRUE) {
  # load and clean data
  hpp <- dir(loc, full.names = TRUE)
  hpd <- str_subset(hpp, "time_series_covid19_deaths_US.csv$")
  hpc <- str_subset(hpp, "time_series_covid19_confirmed_US.csv$")
  
  ddat <- process_hopkins(hpd, weekly = weekly)
  cdat <- process_hopkins(hpc, weekly = weekly)
  bind_rows(ddat, cdat)
}

load_covidtracking <- function(loc) {
  sn <- state_abb_fips$state_name
  names(sn) <- state_abb_fips$state
  
  tfile <-
    dir(loc, pattern = "^covid19us-.*\\.csv$", full.names = TRUE) %>%
    tail(1)
  
  tdf <- read_csv(
    tfile,
    col_types = cols_only(
      date = col_date(),
      state = col_character(),
      hospitalized_increase = col_integer()
    )
  ) %>% 
    mutate(Province_State = sn[state]) %>%
    select(-state)
  
  fips <- state_abb_fips$state_code
  names(fips) <- state_abb_fips$state_name 
  
  tdf2 <-
    tdf %>%
    add_column(target_type = "day ahead inc hosp") %>%
    rename(target_end_date = date,
            value = hospitalized_increase) %>%
    mutate(location = fips[Province_State])
  
  tdf3 <- 
    tdf2 %>% group_by(target_end_date) %>%
    summarise(value = sum(value)) %>%
    mutate(Province_State = "All", 
           location = "US") %>%
    add_column(target_type = "day ahead inc hosp")
  
  tdf4 <- bind_rows(tdf2, tdf3)
  tdf4 %>% filter(Province_State %in% covidhub_locations)
}

pull_data <- function(compid, dir){
  tstamp <- format(Sys.time(), "%Y-%m-%d--%H-%M-%S")
  if (compid == "state"){
    idt <- cdcfluview::ilinet(region = c("state"))
    stem <- "state-ILInet"
  } else if (compid == "national") {
    idt_nat <- cdcfluview::ilinet(region = c("national")) %>% 
      mutate(region = as.character(region))
    idt_reg <- cdcfluview::ilinet(region = c("hhs")) %>% 
      mutate(region = as.character(region))
    idt <- bind_rows(idt_nat, idt_reg)
    stem <- "national-regional-ILInet"
  } else if (compid == "hosp") {
    idt <- covid19us::get_states_daily(state = "GA")
    stem <- "covid19us"
  }
  file <- paste0(stem, "-", tstamp, ".csv")
  path <- file.path(dir, file)
  if(!dir.exists(dir)) dir.create(dir)
  write_csv(idt, path) %>% tail()
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

param_map <- function(x, w, fixed = wfixed){
  ret <- list()
  ret$logE0 <- w[1]
  ret$logH0 <- w[2]
  ret$logtauc <- w[3]
  ret$logtauh <- w[4]
  ret$logtaud <- w[5]
  ret$logchp <- w[6]
  ret$loghfp <- w[7]
  ret$loggammahd <- log(fixed["gamma_h"])
  ret$noiseeffect <- w[8]
  ret$prophomeeffect <- w[9]
  ret$bvec <- w[seq(10, length(w))]
  ret$times <- x$time
  ret$dosesiqr <- x$dosesiqr
  ret$prophomeiqr <- x$prophomeiqr
  ret$betanoise <- x$betanoise
  ret$eta <- fixed["eta"]
  ret$gamma <- fixed["gamma"]
  ret$N <- fixed["N"]
  ret$t0 <- fixed["t0"]
  ret
}

detect_frac <- function(t,
                        max_detect_par = 0.4,
                        detect_inc_rate = 1.1,
                        half_detect = 30,
                        base_detect_frac = 0.1) {
  max_detect_par * (t ^ detect_inc_rate)  / ((half_detect ^ detect_inc_rate) + (t ^ detect_inc_rate)) + base_detect_frac
}

calc_kf_nll <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  if (ncol(y) == 3) {
    pvar <-
      c(
        p$logE0,
        p$logH0,
        p$logtauc,
        p$logtauh,
        p$logtaud,
        p$logchp,
        p$loghfp,
        p$noiseeffect,
        p$prophomeeffect,
        p$bvec
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    nll <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.obj",
      "(pvar, cov, z; N = N, η = η, γ = γ, γh = γh, γd = γd, a = a, betasd = betasd, just_nll = true)"
    ))
  } else if(ncol(y) == 1 && "cases" %in% names(y)) {
    pvar <-
      c(
        p$logE0,
        p$logtauc,
        p$logdoseeffect,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    JuliaCall::julia_assign("h0", exp(p$logH0))
    JuliaCall::julia_assign("τh", exp(p$logtauh))
    JuliaCall::julia_assign("τd", exp(p$logtaud))
    JuliaCall::julia_assign("chp", exp(p$logchp))
    JuliaCall::julia_assign("hfp", exp(p$loghfp))
    nll <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.obj",
      "(pvar, cov, z; N = N, η = η, γ = γ, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τh, chp = chp, hfp = hfp, a = a, betasd = betasd, just_nll = true)"
    ))
  }
  
  nll
}

calc_kf_grad <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  if (ncol(y) == 3) {
    pvar <-
      c(
        p$logE0,
        p$logH0,
        p$logtauc,
        p$logtauh,
        p$logtaud,
        p$logchp,
        p$loghfp,
        p$noiseeffect,
        p$prophomeeffect,
        p$bvec
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    g <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.grad",
      "(pvar, cov, z; N = N, η = η, γ = γ, γd = γd, γh = γh, a = a, betasd = betasd, just_nll = true)"
    ))
  } else if(ncol(y) == 1 && "cases" %in% names(y)) {
    pvar <-
      c(
        p$logE0,
        p$logtauc,
        p$logdoseeffect,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    JuliaCall::julia_assign("h0", exp(p$logH0))
    JuliaCall::julia_assign("τh", exp(p$logtauh))
    JuliaCall::julia_assign("τd", exp(p$logtaud))
    JuliaCall::julia_assign("chp", exp(p$logchp))
    JuliaCall::julia_assign("hfp", exp(p$loghfp))
    g <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.grad",
      "(pvar, cov, z; N = N, η = η, γ = γ, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τd, chp = chp, hfp = hfp, a = a, betasd = betasd, just_nll = true)"
    ))
  }
  g
}

calc_kf_hess <- function(w, x, y, betasd, a, pm) {
  p <- pm(x, w)
  if (ncol(y) == 3) {
    pvar <-
      c(
        p$logE0,
        p$logH0,
        p$logtauc,
        p$logtauh,
        p$logtaud,
        p$logchp,
        p$loghfp,
        p$noiseeffect,
        p$prophomeeffect,
        p$bvec
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    g <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.hess",
      "(pvar, cov, z; N = N, η = η, γ = γ, γd = γd, γh = γh, a = a, betasd = betasd, just_nll = true)"
    ))
  } else if(ncol(y) == 1 && "cases" %in% names(y)) {
    pvar <-
      c(
        p$logE0,
        p$logtauc,
        p$logdoseeffect,
        p$bpars
      )
    JuliaCall::julia_assign("pvar", pvar)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("N", p$N)
    JuliaCall::julia_assign("η", p$eta)
    JuliaCall::julia_assign("γ", p$gamma)
    JuliaCall::julia_assign("z", y)
    JuliaCall::julia_assign("cov", x)
    JuliaCall::julia_assign("a", a)
    JuliaCall::julia_assign("betasd", betasd)
    JuliaCall::julia_assign("γd", exp(p$loggammahd))
    JuliaCall::julia_assign("γh", exp(p$loggammahd))
    JuliaCall::julia_assign("h0", exp(p$logH0))
    JuliaCall::julia_assign("τh", exp(p$logtauh))
    JuliaCall::julia_assign("τd", exp(p$logtaud))
    JuliaCall::julia_assign("chp", exp(p$logchp))
    JuliaCall::julia_assign("hfp", exp(p$loghfp))
    g <- JuliaCall::julia_eval(paste0(
      "InfectionKalman.hess",
      "(pvar, cov, z; N = N, η = η, γ = γ, γd = γd, γh = γh, h0 = h0, τh = τh, τd = τd, chp = chp, hfp = hfp, a = a, betasd = betasd, just_nll = true)"
    ))
  }
  g
}

initialize_estimates <- function(x, y, wfixed, t0 = 2020.164, dt = 0.00273224) {
  tau_cases_init <- max(var(y$cases, na.rm = TRUE), 1)
  tau_hosp_init <- max(var(y$hospitalizations, na.rm = TRUE), 1)
  tau_deaths_init <- max(var(y$deaths, na.rm = TRUE), 1)
  wsize <- nrow(y)
  stopifnot(!wsize %% 7)
  rhot <- detect_frac((x$time - t0) * 365.25)
  
  hfill <- zoo::na.fill(y$hospitalizations, "extend")
  est_true_cases_hosps <- (y$cases - hfill) / rhot + hfill
  Iest <- est_true_cases_hosps / (wfixed["gamma"] * dt)
  
  I0init <- mean(head(Iest, n = 7), na.rm = TRUE)
  E0init <- I0init * wfixed["gamma"] / wfixed["eta"] 
  
  hosp_obs <- !is.na(y$hospitalizations)
  chp_init <-
    sum(y$hospitalizations, na.rm = TRUE) / sum(est_true_cases_hosps)
  
  Hest <- Iest  * wfixed["gamma"] * chp_init / wfixed["gamma_h"]
  
  H0init <- mean(head(Hest, n = 7), na.rm = TRUE)
  hfp_init <-
    max(sum(y$deaths[hosp_obs], na.rm = TRUE) / sum(y$hospitalizations, na.rm = TRUE),
        0.01)
  
  l <- 7 
  Rt <- lead(y$cases / rhot + hfill, n = l) / (y$cases / rhot + hfill)
  m <- 10
  Rtfilt <- signal::sgolayfilt(na.omit(Rt), p = 2, n = 2 * m  + 1) # edge effects seem biased
  Rtfilt[1:m] <- NA
  lf <- length(Rtfilt)
  Rtfilt[(lf - m + 1):lf] <- NA
  
  D0 <- H0init * wfixed["gamma_h"] / wfixed["gamma_d"] * hfp_init
  S0init <- wfixed["N"] - E0init - I0init - H0init - D0
  Sdecrement <- cumsum(lead(y$cases / rhot + hfill, n = l))
  St <- S0init - Sdecrement
  betat <- Rtfilt * wfixed["gamma"] * wfixed["N"] / na.omit(St)
  betat[1:m] <- betat[m + 1]
  betat[(lf - m + 1):lf] <- betat[lf - m]
  betat2 <- c(betat, rep(betat[lf - m], l))
  
  df <- data.frame(logbeta = log(betat2), dosesiqr = x$dosesiqr, prophomeiqr = x$prophomeiqr)
  mod <- lm(logbeta~ prophomeiqr, data = df)
  print(summary(mod))
  #doseeffect_init <- coef(mod)["dosesiqr"] %>% unname()
  noiseeffect_init <- sd(residuals(mod))
  prophomeeffect_init <- coef(mod)["prophomeiqr"] %>% unname()
  binit <- coef(mod)["(Intercept)"] + residuals(mod)
  binit_weekly <- matrix(binit, nrow = 7) %>% colMeans()
  names(binit_weekly) <- paste0("b", seq_along(binit_weekly))

  winit <- c(
    logE0 = log(E0init),
    logH0 = log(H0init),
    logtauc = log(tau_cases_init),
    logtauh = log(tau_hosp_init),
    logtaud = log(tau_deaths_init),
    logchp = log(chp_init),
    loghfp = log(hfp_init),
    noiseeffect = noiseeffect_init,
    prophomeeffect = prophomeeffect_init,
    binit_weekly
  )
  winit
}

kf_nll_details <- function(w, x, y, betasd, a, pm, fet) {
  p <- pm(x, w)
  nll <- kfnll(
    bvec = p$bvec,
    logE0 = p$logE0,
    logH0 = p$logH0,
    logtauc = p$logtauc,
    logtauh = p$logtauh,
    logtaud = p$logtaud,
    logchp = p$logchp,
    loghfp = p$loghfp,
    loggammahd = p$loggammahd,
    noiseeffect = p$noiseeffect,
    prophomeeffect = p$prophomeeffect,
    eta = p$eta,
    gamma = p$gamma,
    N = p$N,
    z = y,
    t0 = p$t0,
    cov = x,
    fet = fet,
    just_nll = FALSE,
    fet_zero_cases_deaths = "weekly",
    nsim = 20,
    betasd = betasd,
    a = a
  )
  nll
}

kfnll <-
  function(bvec,
           logE0,
           logH0,
           logtauc,
           logtauh,
           logtaud,
           logchp,
           loghfp,
           loggammahd,
           noiseeffect,
           prophomeeffect,
           eta,
           gamma,
           N,
           z,
           t0,
           cov, 
           Phat0 = diag(c(1, 1, 1, 0, 0, 1, 1, 0)),
           fets = NULL,
           fet_zero_cases_deaths = "daily",
           nsim,
           a = .98,
           betasd = 1,
           maxzscore = Inf,
           just_nll = TRUE) {
    diffeqr::diffeq_setup()
    JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
    
    E0 <- exp(logE0)
    I0 <- E0 * eta / gamma
    H0 <- exp(logH0)
    gamma_d <- gamma_h <- exp(loggammahd)
    D0 <- H0 * gamma_h / gamma_d * exp(loghfp)
    xhat0 <- c(N - E0 - I0 - H0 - D0, E0, I0, 0, 0, H0, D0, 0)
    names(xhat0) <- c("S", "E", "I", "C", "Hnew", "H", "D", "Drep")

    if (ncol(z) == 1 && "cases" %in% names(z)){
      z$hospitalizations <- NA
      z$deaths <- NA
    }
    
    z <- data.matrix(z[, c("cases", "hospitalizations", "deaths")])
    is_z_na <- is.na(z)
    T <- nrow(z)
    dobs <- ncol(z)
    dstate <- length(xhat0)
    stopifnot(T > 0)
    
    ytilde_kk <- ytilde_k <- array(NA_real_, dim = c(dobs, T))
    S <- array(NA_real_, dim = c(dobs, dobs, T))
    K <- array(NA_real_, dim = c(dstate, dobs, T))
    rdiagadj <- array(1, dim = c(dobs, T))
    
    xhat_kk <- xhat_kkmo <- array(NA_real_, dim = c(dstate, T))
    rownames(xhat_kk) <- rownames(xhat_kkmo) <- names(xhat0)
    P_kk <- P_kkmo <- array(NA_real_, dim = c(dstate, dstate, T))
    
    H <- function(time, t0 = 2020.164){
      day <- (time - t0) * 365.25
      rbind(c(0, 0, 0, detect_frac(day), 1, 0, 0, 0), 
            c(0, 0, 0,    0, 1, 0, 0, 0),
            c(0, 0, 0,    0, 0, 0, 0, 1))
    }
    R <- diag(exp(c(logtauc, logtauh, logtaud)))
    
    for (i in seq(1, T)) {
      if (i == 1) {
        xhat_init <- xhat0
        PNinit <- Phat0
      } else {
        xhat_init <- xhat_kk[, i - 1]
        PNinit <- P_kk[, , i - 1]
      }
      xhat_init["C"] <- 0
      xhat_init["Hnew"] <- 0
      xhat_init["Drep"] <- 0
      
      PNinit[, 4] <- PNinit[4,] <- 0
      PNinit[, 5] <- PNinit[5,] <- 0
      PNinit[, 8] <- PNinit[8,] <- 0
      
      u0 <- cbind(xhat_init, PNinit)
      beta_t <- min(exp(bvec[(i - 1) %/% 7 + 1] + cov$betanoise[i] * noiseeffect + prophomeeffect * cov$prophomeiqr[i]), 4 * gamma)
      par <- c(beta_t, N,  0, eta, gamma, gamma_d, gamma_h, exp(logchp), exp(loghfp))
      
      JuliaCall::julia_assign("u0", u0)
      JuliaCall::julia_assign("tspan", c(0, 0.00273224))
      JuliaCall::julia_assign("par", par)
      JuliaCall::julia_eval("prob = ODEProblem(InfectionKalman.vectorfield,u0,tspan,par)")
      XP <- JuliaCall::julia_eval("solve(prob, Tsit5(), saveat = tspan[2]).u[2]")
      xhat_kkmo[, i] <- XP[, 1]
      P_kkmo[, , i] <- XP[, -1]
      
      for (j in 1:dstate){
        if (P_kkmo[j,j,i] < 0){
          P_kkmo[j,,i] <- 0
          P_kkmo[,j,i] <- 0
        }
      }
      
      ytilde_k[, i] <- matrix(z[i, ], ncol = 1) - 
        H(cov$time[i]) %*% xhat_kkmo[, i, drop = FALSE]     
      S[, , i] <- H(cov$time[i]) %*% P_kkmo[, , i] %*% t(H(cov$time[i])) + R
      
      for (j in 1:dobs){
        if (is.na(z[i,j])){
          zscore <- 0
          rdiagadj[j,i] <- rdiagadj[j,i] + 0
        } else {
          sd <- sqrt(S[j,j,i])
          zscore <- ytilde_k[j,i] / sd 
          if (abs(zscore) > maxzscore){
            adjzscore <- maxzscore / (1 + abs(zscore) - maxzscore)
            newsd <- abs(ytilde_k[j,i]) / adjzscore
            rdiagadj[j,i] <- rdiagadj[j,i] + (newsd) ^ 2 - sd ^ 2
          } else {
            rdiagadj[j,i] <- rdiagadj[j,i] + 0
          }
        }
        S[j,j,i] <- S[j,j,i] + rdiagadj[j,i]
      }
      K[, , i] <- P_kkmo[, , i] %*% t(H(cov$time[i])) %*% solve(S[, , i])
      desel <- is_z_na[i, ]
      K[, desel, i] <- 0
      
      xhat_kk[, i] <-
        xhat_kkmo[, i, drop = FALSE] +
        K[, !desel, i] %*% ytilde_k[!desel, i, drop = FALSE]
      
      xhat_kk[xhat_kk[, i] < 0, i] <- 1e-4
      P_kk[, , i] <-
        (diag(dstate) - K[, , i] %*% H(cov$time[i])) %*% P_kkmo[, , i]
      
      for (j in 1:dstate){
        if (P_kk[j,j,i] < 0){
          P_kk[j,,i] <- 0
          P_kk[,j,i] <- 0
        }
      }
      ytilde_kk[, i] <- matrix(z[i, ], ncol = 1) - 
        H(cov$time[i]) %*% xhat_kk[, i, drop = FALSE]
    }
    
    rwlik <- 0
    for (i in 1:((T %/% 7) - 1)){
      step <- bvec[i + 1] - bvec[i]
      rwlik <- rwlik + dnorm(step, mean = 0, sd = betasd, log = TRUE)
    }
    nll <- 0
    for (i in seq(1, T)){
      sel <- !is_z_na[i, ]
      nll <- nll + 
        t(ytilde_k[sel, i]) %*% solve(S[sel, sel, i]) %*% ytilde_k[sel, i] + 
        log(det(S[,,i][sel, sel, drop = FALSE])) + dobs * log(2 * pi)
    }
    nll <- 0.5 * nll - rwlik
    
    if (!just_nll) {
      if (!is.null(fets)) {
        nsimdays <- nrow(fets)
        sim_means <- array(NA_real_, dim = c(dobs, nsimdays, nsim))
        sim_cov <- array(NA_real_, dim = c(dobs, dobs, nsimdays, nsim))
        for (j in seq_len(nsim)) {
          logbeta_fet <- numeric(nsimdays)
          logbeta_fet[1] <-
            log(gamma) + a * (logbeta[T] - log(gamma)) +  
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
          if (fet_zero_cases_deaths == "daily" ||
              fets$target_wday[1] == 1) {
            xhat_init["C"] <- 0
            PNinit[, 4] <- PNinit[4,] <- 0
            xhat_init["Drep"] <- 0
            PNinit[, 8] <- PNinit[8,] <- 0 
          }
          xhat_init["Hnew"] <- 0
          PNinit[, 5] <- PNinit[5,] <- 0
          
          XP <- iterate_f_and_P(
            xhat_init,
            PN = PNinit,
            eta = eta,
            gamma = gamma,
            gamma_d = gamma_d,
            gamma_h = gamma_h,
            chp = exp(logchp),
            hfp = exp(loghfp),
            N = N,
            beta_t = exp(logbeta_fet[1]) * exp(-doseeffect * doses[T])
          )
          sim_means[, 1, j] <- H(fets$target_end_times[1]) %*% XP$xhat
          sim_cov[, , 1, j] <- H(fets$target_end_times[1]) %*% XP$PN %*% t(H(fets$target_end_times[1])) + R
          for (i in seq_along(fets$target_end_times[-1])) {
            xhat_init <- XP$xhat
            PNinit <- XP$PN
            if (fet_zero_cases_deaths == "daily" ||
                fets$target_wday[i + 1] == 1) {
              xhat_init["C"] <- 0
              PNinit[, 4] <- PNinit[4,] <- 0
              xhat_init["Drep"] <- 0
              PNinit[, 8] <- PNinit[8,] <- 0
            }
            xhat_init["Hnew"] <- 0
            PNinit[, 5] <- PNinit[5,] <- 0
            
            XP <-
              iterate_f_and_P(
                xhat_init,
                PN = PNinit,
                eta = eta,
                gamma = gamma,
                gamma_d = gamma_d,
                gamma_h = gamma_h,
                chp = exp(logchp),
                hfp = exp(loghfp),
                N = N,
                beta_t = exp(logbeta_fet[i + 1]) * exp(-doseeffect * doses[T])
              )
            sim_means[, i + 1, j] <- H(fets$target_end_times[i + 1]) %*% XP$xhat
            sim_cov[, , i + 1, j] <-
              H(fets$target_end_times[i + 1]) %*% XP$PN %*% t(H(fets$target_end_times[i + 1])) + R
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
        bvec = bvec,
        gamma = gamma,
        rdiagadj = rdiagadj
      )
    } else {
      nll
    }
  }
