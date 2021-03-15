#!/usr/bin/env Rscript

source("covidhub-common.R")
library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

lambda <- 1 / seq(0.001, 0.1, length.out = 10)
agrid <- c(0.94, 0.95)
par2name <- function(lambda, a){
  paste0("lambda", sprintf("%06.2f", lambda), 
         "-a", sprintf("%02.2f", a),
         "-CEID-InfectionKalman")
}
dirnames <- outer(lambda, agrid, par2name)

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat2 <- map_dfr(dirnames, load_from_dir)

fdat2$lambda <- fdat2$model %>% str_extract("^lambda\\d+.\\d{2}-") %>% 
  str_remove("^lambda") %>% str_remove("-$") %>% as.numeric()
fdat2$a <- fdat2$model %>% str_extract("-a\\d{1}.\\d{2}-") %>% 
  str_remove("^-a") %>% str_remove("-$") %>% as.numeric()

fdat3 <- fdat2 %>% filter(target_variable == "inc case")
splt <- split(fdat3, fdat3$location)
data_date <- Sys.getenv("ddt", unset = "2021-02-01")
hopdir <- file.path("hopkins", data_date)
tdat <- load_hopkins(hopdir, weekly = TRUE)
tdat2 <- tdat %>% rename(true_value=value) %>% 
  filter(target_type == "wk ahead inc case") %>% 
  select(location, target_end_date, true_value)


healthdirs <- dir(file.path("healthdata", data_date), full.names = TRUE)
names(healthdirs) <- basename(healthdirs)

load_healthd <- function(dirpath){
  path <- file.path(dirpath, "epidata.csv")
  tdat <- read_csv(path,
                  col_types = cols_only(date = col_date("%Y%m%d"),
                                        previous_day_admission_adult_covid_confirmed = col_integer(),
                                        previous_day_admission_pediatric_covid_confirmed = col_integer())) %>%
    mutate(
      hospitalizations = previous_day_admission_adult_covid_confirmed +
        previous_day_admission_pediatric_covid_confirmed,
      target_end_date = date - lubridate::ddays(1)
    )
  tdat %>% select(target_end_date, hospitalizations)
}

healthdata <- map(healthdirs, load_healthd) %>% bind_rows(.id = "loc")


plot_forecast_grid <- function(locdata){
  locdata %>% left_join(tdat2, by = c("location", "target_end_date")) %>%
  filter(type == "point") %>%
  ggplot(aes(x = target_end_date, y = value, group = forecast_date)) + 
  geom_line() + 
  geom_point() + 
  geom_line(aes(y = true_value), color = "red") + 
  facet_grid(a~lambda) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_date(
    name = "Year-Epiweek",
    date_breaks = "1 weeks",
    labels = function(x) {
      wk <- strftime(x, format = "%U") %>% as.integer() + 1
      yr <- strftime(x, format = "%y")
      paste(yr, wk, sep = '-')
      }
    ,
    expand = expansion()
  ) +
  labs(y = "Incident cases")
}

plots <- map(splt, plot_forecast_grid)

plot_writer <- function(p, loc){
  plot_name <- paste0("trajectories-all/", loc, ".png")
  ggsave(plot_name, p, height = 7.5, width = 24)
}
if (!dir.exists("trajectories-all")){
  dir.create("trajectories-all")
}
map2(plots, names(plots), plot_writer)
