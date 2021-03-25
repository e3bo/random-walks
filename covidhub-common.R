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
  ret$loggammahd <- w[8]
  ret$logdoseeffect <- w[9]
  ret$bpars <- w[seq(10, length(w))]
  ret$times <- x[, 1]
  ret$doses <- x[, 2]
  ret$eta <- fixed["eta"]
  ret$gamma <- fixed["gamma"]
  ret$N <- fixed["N"]
  ret$t0 <- fixed["t0"]
  ret
}