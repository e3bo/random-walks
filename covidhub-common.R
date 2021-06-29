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
    hpd %>% mutate(week = lubridate::epiweek(target_end_date),
                   year = lubridate::epiyear(target_end_date)) %>% 
    group_by(location, week, year) %>% 
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

read_us_mob_report <- function(path, locpat = "US-CA"){
  readr::read_csv(
    path,
    col_types = readr::cols_only(
      country_region = col_character(),
      sub_region_1 = col_character(),
      sub_region_2 = col_character(),
      iso_3166_2_code = col_character(),
      date = col_date(format = "%Y-%m-%d"),
      residential_percent_change_from_baseline = col_double()
    )
  ) %>% filter(is.na(sub_region_2)) %>%
    select(-sub_region_2) %>%
    filter(str_detect(iso_3166_2_code, locpat))
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

initialize_estimates <- function(x, y, wfixed, dt = 1 /365.25) {
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

name_params <- function(w, z, cov, f, trans = TRUE){
  p <- list()
  p$nτ_c <- tail(cov$τ_cmap, 1)
  p$np_h <- tail(cov$p_hmap, 1)
  
  if (trans){
    tfun <- plogis
    tfun2 <- exp
  } else{
    tfun <- tfun2 <- identity
  }
  
  if (ncol(z) == 3) {
    p$τ_h <- tfun2(w[2])
    p$p_hweekend <- tfun(w[3])
    if ("doses_scaled" %in% names(cov)){
      p$doseeffect <- w[6]
      w2 <- c(w[1], w[4:5], w[7:length(w)])
    } else{
      w2 <- c(w[1], w[4:length(w)])
    }

    p$p_h <- tfun(w2[(8 + p$nτ_c):(8 + p$nτ_c + p$np_h - 1)])
  } else if(ncol(z) == 2 && !"hospitalizations" %in% names(z)) {
    p$τ_h <- f$τ_h
    p$p_h <- f$p_h
    p$np_h <- 0
    z$hospitalizations <- NA
    w2 <- w
  }

  p$L0 <- tfun2(w2[1])
  p$τ_d <- tfun2(w2[2])
  p$residentialeffect <- w2[3]
  p$p_d <- tfun(w2[4])
  p$γ_d12 <- tfun2(w2[5])
  p$γ_d34 <- tfun2(w2[6])
  p$γ_z17 <- tfun2(w2[7])
  p$τ_c <- tfun2(w2[seq(8, 8 + p$nτ_c - 1)])
  p$β_0 <- w2[seq(8 +  p$nτ_c + p$np_h, length(w2))]
  zloc <- data.matrix(z[, c("cases", "hospitalizations", "deaths")])
  list(p, zloc)
}

advance_proc <- function(p, f, cov, xhat_init, PNinit, t){
  
  u0 <- cbind(xhat_init, PNinit)
  if ("doses_scaled" %in% names(cov)) {
    β   <-
      exp(p$β_0[cov$β_0map[t]] + p$residentialeffect * cov$residential[t] + p$doseeffect * cov$doses_scaled[t])
  } else {
    β   <-
      exp(p$β_0[cov$β_0map[t]] + p$residentialeffect * cov$residential[t])
  }
  if(β > 1000){
    β <- 1000
  }
  par <- c(β, f$N,  f$η,  f$γ,  f$γ_d,  f$γ_z, f$γ_h, p$p_h[cov$p_hmap[t]], p$p_d, cov$ρ[t])
  if (cov$wday[t] == 1) {
    par[5] <- p$γ_d12
    par[6] <- p$γ_z17
    par[8] <- par[8] * p$p_hweekend
  } else if (cov$wday[t] == 2){
    par[5] <- p$γ_d12
  } else if (cov$wday[t] %in% c(3, 4)) {
    par[5] <- p$γ_d34
  } else if (cov$wday[t] == 7){
    par[6] <- p$γ_z17
    par[8] <- par[8] * p$p_hweekend
  }
  JuliaCall::julia_assign("u0", u0)
  JuliaCall::julia_assign("tspan", c(0, 1 / 365.25))
  JuliaCall::julia_assign("par", par)
  JuliaCall::julia_eval("prob = ODEProblem(vectorfield,u0,tspan,par)")
  XP <-
    JuliaCall::julia_eval("solve(prob, Tsit5(), saveat = tspan[2]).u[2]")
  list(xhat_kkmo = XP[, 1], P_kkmo = XP[,-1])
}

calc_kf_nll_r <-
  function(w,
           cov,
           z,
           p_hsd,
           β_0sd,
           τ_csd,
           wfixed,
           cov_sim = NULL,
           fet_zero_cases_deaths = "daily",
           just_nll = TRUE) {
    diffeqr::diffeq_setup("/opt/julia-1.5.3/bin")
    JuliaCall::julia_eval("include(\"InfectionKalman.jl\")")
    
    f <- as.list(wfixed)
    npout <- name_params(w, z, cov, f)
    p <- npout[[1]]
    zloc <- npout[[2]]
    zmiss <- is.na(zloc)
    
    Y0 <- p$L0 * f$η / f$γ
    H0 <- p$p_h[1] * f$γ * Y0 / f$γ_h
    D0 <- H0 *  f$γ_h /  f$γ_d * p$p_d
    Z0 <- (p$p_h[cov$p_hmap[1]] + cov$ρ[1] * (1 - p$p_h[cov$p_hmap[1]])) * Y0 *  f$γ / f$γ_z

    x0 <- c(
      max(f$N - p$L0 - Y0 - H0 - D0, f$N * 0.1),
      #X
      min(Y0, f$N),
      #Y
      min(p$L0, f$N),
      #L
      min(Z0, f$N),
      #Z
      0,
      #Z_r
      min(H0, f$N),
      #H
      0,
      #A
      min(D0, f$N),
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
      
      XPL <-
        advance_proc(
          p = p,
          f = f,
          cov = cov,
          xhat_init = xhat_init,
          PNinit = PNinit,
          t = i
        )
      xhat_kkmo[, i] <- XPL$xhat_kkmo
      P_kkmo[,, i] <- XPL$P_kkmo
      
      for (j in 1:dstate) {
        if (P_kkmo[j, j, i] < 0) {
          P_kkmo[j, , i] <- 0
          P_kkmo[, j, i] <- 0
        }
      }
      r <- diag(c(p$τ_c[cov$τ_cmap[i]] * xhat_init[2],  p$τ_h,  p$τ_d))
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
    for (i in seq_len(length(p$β_0) - 1)) {
      step <-  p$β_0[i + 1] -  p$β_0[i]
      rwlik <-
        rwlik + dnorm(step,
                      mean = 0,
                      sd =  β_0sd,
                      log = TRUE)
    }
    
    for (i in seq_len(length(p$τ_c) - 1)) {
      τ_cstep <- log(p$τ_c[i + 1]) - log(p$τ_c[i])
      rwlik <-
        rwlik + dnorm(τ_cstep,
                        mean = 0,
                        sd =  τ_csd,
                        log = TRUE)
    }
    
    for (i in seq_len(length(p$p_h) - 1)) {
      p_hstep <- qlogis(p$p_h[i + 1]) - qlogis(p$p_h[i])
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
        log(det(S[, , i][sel, sel, drop = FALSE])) + sum(sel) * log(2 * pi)
    }
    nll <- 0.5 * nll - rwlik
    if (!just_nll) {
      if (!is.null(cov_sim)) {
        nsimdays <- nrow(cov_sim)
        sim_means <- array(NA_real_, dim = c(dobs, nsimdays))
        sim_Sigma <- sim_P <- sim_R <- array(NA_real_, dim = c(dobs, dobs, nsimdays))
        
        x_sim <- array(NA_real_, dim = c(dstate, nsimdays))
        rownames(x_sim) <- names(x0)
        P_sim <- array(NA_real_, dim = c(dstate, dstate, nsimdays))
        
        if (fet_zero_cases_deaths != "daily" && cov_sim$wday[1] > 2){
          warning("first weekly forecast is for incomplete week")
        }
        
        for (j in seq_len(nsimdays)) {
          if (j == 1) {
            xhat_init <- xhat_kk[, nobs]
            PNinit <- P_kk[, , nobs]
          } else {
            xhat_init <- XPL$xhat_kkmo
            PNinit <- XPL$P_kkmo
          }
          
          if (fet_zero_cases_deaths == "daily" ||
              cov_sim$wday[j] == 1) {
            xhat_init[5] <- 0
            PNinit[, 5] <- PNinit[5,] <- 0
            xhat_init[9] <- 0
            PNinit[, 9] <- PNinit[9,] <- 0
          }
          xhat_init[7] <- 0
          PNinit[, 7] <- PNinit[7,] <- 0
          
          XPL <-
            advance_proc(
              p = p,
              f = f,
              cov = cov_sim,
              xhat_init = xhat_init,
              PNinit = PNinit,
              t = j
            )

          r <-
            diag(c(p$τ_c[cov$τ_cmap[nobs]] * xhat_init[2],  p$τ_h,   p$τ_d))
          
          x_sim[, j] <- XPL$xhat_kkmo
          P_sim[, , j] <- XPL$P_kkmo
          
          sim_means[, j] <- hmat %*% x_sim[, j]
          sim_P[, ,j] <- hmat %*% P_sim[,, j] %*% t(hmat)
          sim_R[, ,j] <- r
          sim_Sigma[, , j] <- hmat %*% P_sim[,, j] %*% t(hmat) + r
        }
      } else {
        sim_means <- sim_Sigma <- sim_R <- sim_P <- x_sim <- P_sim <- NULL
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
        sim_P = sim_P,
        sim_R = sim_R,
        sim_Sigma = sim_Sigma,
        x_sim = x_sim,
        P_sim = P_sim,
        rdiagadj = rdiagadj
      )
    } else {
      nll
    }
  }

load_health_data <-
  function(forecast_date, forecast_loc, cov_thresh = 0.9) {
    healthd <-
      file.path("healthdata", forecast_date, forecast_loc, "epidata.csv")
    tdat2 <- read_csv(
      healthd,
      col_types = cols_only(
        date = col_date("%Y%m%d"),
        previous_day_admission_adult_covid_confirmed = col_integer(),
        previous_day_admission_pediatric_covid_confirmed = col_integer(),
        previous_day_admission_adult_covid_confirmed_coverage = col_integer()
      )
    )
    most_recent_coverage <- tdat2 %>% arrange(date) %>%
      pull(previous_day_admission_adult_covid_confirmed_coverage) %>%
      tail(n = 1)
    min <- most_recent_coverage * cov_thresh
    max <- most_recent_coverage / cov_thresh
    tdat3 <- tdat2 %>%
      filter(
        previous_day_admission_adult_covid_confirmed_coverage >
          min
      ) %>%
      filter(
        previous_day_admission_adult_covid_confirmed_coverage <
          max 
      ) %>%
      mutate(
        hospitalizations = previous_day_admission_adult_covid_confirmed +
          previous_day_admission_pediatric_covid_confirmed,
        target_end_date = date - lubridate::ddays(1)
      ) %>%
      select(target_end_date, hospitalizations)
    tdat3
  }
