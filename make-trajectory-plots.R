#!/usr/bin/env Rscript

source("covidhub-common.R")
library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

## load forecasts
lambda <- 20
par2name <- function(lambda){
  paste0("lambda", sprintf("%06.2f", lambda), 
         "-status-quo",
         "-CEID-InfectionKalman")
}
dirnames <- purrr::map(lambda, par2name)

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}

fdat2 <- map_dfr(dirnames, load_from_dir) %>% 
  bind_rows(readRDS("other-model-forecasts.rds")) %>%
  filter(model != "CEID-Walk") %>%
  filter(model != "COVIDhub-baseline")

## split forecasts by location
fdat3 <- fdat2 %>% filter(target_variable == "inc case")
splt_cases <- split(fdat3, fdat3$location)

fdat4 <- fdat2 %>% filter(target_variable == "inc hosp") %>%
  filter(forecast_date >= "2020-12-07") %>%
  filter(target_end_date <= "2021-05-01")
splt_hosp <- split(fdat4, fdat4$location)

fdat5 <- fdat2 %>% filter(target_variable == "inc death")
splt_death <- split(fdat5, fdat5$location)

## load JHU data
data_date <- Sys.getenv("ddt", unset = "2021-05-24")
hopdir <- file.path("hopkins", data_date)
tdat <- load_hopkins(hopdir, weekly = TRUE)
tdat2 <- tdat %>% rename(true_value=value) %>% 
  filter(target_type == "wk ahead inc case") %>% 
  select(location, target_end_date, true_value)

tdat3 <- tdat %>% rename(true_value=value) %>% 
  filter(target_type == "wk ahead inc death") %>% 
  select(location, target_end_date, true_value)


## load hhs data
healthdirs <-
  dir(file.path("healthdata", data_date), full.names = TRUE)
names(healthdirs) <- basename(healthdirs)
load_healthd <- function(dirpath) {
  path <- file.path(dirpath, "epidata.csv")
  tdat <- read_csv(
    path,
    col_types = cols_only(
      date = col_date("%Y%m%d"),
      previous_day_admission_adult_covid_confirmed = col_integer(),
      previous_day_admission_pediatric_covid_confirmed = col_integer()
    )
  ) %>%
    mutate(
      hospitalizations = previous_day_admission_adult_covid_confirmed +
        previous_day_admission_pediatric_covid_confirmed,
      target_end_date = date - lubridate::ddays(1)
    )
  tdat %>% select(target_end_date, hospitalizations)
}
healthdata <-
  map(healthdirs, load_healthd) %>% bind_rows(.id = "location") %>%
  rename(true_value = hospitalizations)

## make the plots
plot_forecast_grid <- function(locdata, tv = "inc case"){
  xname <- "Year-Epiweek"
  db <- "2 weeks"
  labf <- function(x) {
    wk <- sprintf("%02d", lubridate::epiweek(x))
    yr <- lubridate::epiyear(x)
    paste(yr, wk, sep = '-')
  }
  if (tv == "inc case"){
    yname <- "Incident cases"
    truth <- tdat2
  } else if (tv == "inc death") {
    yname <- "Incident deaths"
    truth <- tdat3
  } else {
    xname <- "Date"
    yname <- "Incident hospitalizations"
    db <- "7 days"
    labf <- waiver()
    truth <- healthdata
  }
  locdata %>% left_join(truth, by = c("location", "target_end_date")) %>%
  filter(quantile %in% c(0.025, 0.5, 0.975)) %>%
  ggplot(aes(x = target_end_date, y = value, group = interaction(forecast_date, model), color = model)) + 
  geom_line(aes(y = true_value), color = "grey", size = 2, alpha = 0.5) + 
  geom_point(aes(y = true_value), color = "grey", size = 3, alpha = 0.5) +
  geom_line() + 
  geom_point() + 
  scale_x_date(
    name = xname,
    date_breaks = db,
    labels = labf,
    expand = expansion()
  ) +
  labs(y = yname) + 
  facet_grid(quantile~., scales = "free_y") + 
  ggthemes::scale_colour_colorblind(name = "Model") +
  theme_minimal() + 
  theme(legend.position = "top") + 
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

plots_cases <- map(splt_cases, plot_forecast_grid)
plots_hosp <- map(splt_hosp, plot_forecast_grid, tv = "inc hosp")
plots_deaths <- map(splt_death, plot_forecast_grid, tv = "inc death")

plot_writer <- function(p, loc, pref){
  plot_name <- paste0("trajectories-all/", pref, loc, ".png")
  ggsave(plot_name, p, width = 4, height = 10)
}
if (!dir.exists("trajectories-all")){
  dir.create("trajectories-all")
}
map2(plots_cases, names(plots_cases), plot_writer, pref = "cases")
map2(plots_hosp, names(plots_hosp), plot_writer, pref = "hosp")
map2(plots_deaths, names(plots_deaths), plot_writer, pref = "death")
