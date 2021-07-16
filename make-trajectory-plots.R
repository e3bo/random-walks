#!/usr/bin/env Rscript

source("covidhub-common.R")
library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

## load forecasts
lambda <- 20
loc <- "06"
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
  filter(model != "COVIDhub-baseline") %>%
  mutate(model = ifelse(model == "lambda020.00-status-quo-CEID-InfectionKalman", 
                        "GISST", model))

GISSTstart <- fdat2 %>% filter(model == "GISST") %>% pull(forecast_date) %>% 
  min()
date_ranges <-
  fdat2 %>% filter(model == "COVIDhub-ensemble") %>%
  group_by(target_variable) %>%
  summarise(start = max(min(forecast_date), GISSTstart),
            stop = max(forecast_date), .groups = "drop")

fdat22 <- fdat2 %>% left_join(date_ranges, by = "target_variable")
## split forecasts by location

fdat3 <- fdat22 %>% filter(target_variable == "inc case") %>%
  filter(forecast_date >= start & forecast_date <= stop)
splt_cases <- fdat3 %>% filter(location == loc)

fdat4 <- fdat22 %>% filter(target_variable == "inc hosp") %>%
  filter(forecast_date >= start & forecast_date <= stop)

splt_hosp <- fdat4 %>% filter(location == loc)

fdat5 <- fdat22 %>% filter(target_variable == "inc death") %>% 
  filter(forecast_date >= start & forecast_date <= stop) 
splt_death <- fdat5 %>% filter(location == loc)

## load JHU data
data_date <- Sys.getenv("ddt", unset = "2021-06-21")
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
  db <- "4 weeks"
  labf <- function(x) {
    wk <- sprintf("%02d", lubridate::epiweek(x))
    yr <- lubridate::epiyear(x)
    paste(yr, wk, sep = '-')
  }
  if (tv == "inc case"){
    yname <- "Cases"
    truth <- tdat2
  } else if (tv == "inc death") {
    yname <- "Deaths"
    truth <- tdat3
  } else {
    xname <- "Day"
    yname <- "Hospital admissions"
    db <- "14 days"
    labf <- waiver()
    truth <- healthdata
  }
  locdata %>% left_join(truth, by = c("location", "target_end_date")) %>%
  filter(quantile %in% c(0.025, 0.5, 0.975)) %>%
  ggplot(aes(x = target_end_date, y = value, 
             group = interaction(forecast_date, model), color = model)) + 
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

plots_cases <- plot_forecast_grid(splt_cases)
plots_hosp <- plot_forecast_grid(splt_hosp, tv = "inc hosp")
plots_deaths <- plot_forecast_grid(splt_death, tv = "inc death")

if (!dir.exists("trajectories-all")){
  dir.create("trajectories-all")
}

prow <- cowplot::plot_grid(
  plots_cases + theme(legend.position="none"),
  plots_hosp + theme(legend.position="none"),
  plots_deaths + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  nrow = 1
)
prow

leg <- cowplot::get_legend(
  plots_cases + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "top")
)

pall <- cowplot::plot_grid(leg, prow, ncol = 1, rel_heights = c(0.1, 1))
ggsave(file.path("trajectories-all", "traj.png"), plot = pall, width = 7.3, 
       dpi = 600)