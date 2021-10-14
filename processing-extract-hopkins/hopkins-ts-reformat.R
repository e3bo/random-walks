
suppressPackageStartupMessages(library(tidyverse))

process_hopkins <- function(path, weekly = FALSE){
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

input_dir <- "/opt/ml/processing/input" 
hpd <- file.path(input_dir, "deaths", "time_series_covid19_deaths_US.csv")
hpc <- file.path(input_dir, "cases", "time_series_covid19_confirmed_US.csv")

print("Reformating deaths time series into standard format")
ddat <- process_hopkins(hpd)
print("Reformating cases time series into standard format")
cdat <- process_hopkins(hpc)
tdat <- bind_rows(ddat, cdat) %>% ungroup()

output_path <- file.path("/opt/ml/processing/output", "time-series-covid19-US.parquet")
print(paste("Saving all time series to ", output_path))
arrow::write_parquet(tdat, output_path)
