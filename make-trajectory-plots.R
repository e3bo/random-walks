#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")
lambda <- c(158.49, 107.98, 73.56, 50.12, 34.15, 23.26, 15.85, 10.8)
dirnames <- paste0("lambda", sprintf("%06.2f", lambda), "-CEID-InfectionKalman")

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat1 <- map_dfr(dirnames, load_from_dir)
fdat2 <- bind_rows(fdat0, fdat1)

plotter <- function(dat) {
  p <- plot_forecast(
    dat,
    target_variable = "inc case",
    truth_source = "JHU",
    intervals = c(.5, .95),
    facet = location ~ model,
    facet_scales = "free_y",
    fill_by_model = TRUE,
    plot = FALSE
  )
  
  p2 <- p +
    scale_x_date(
      name = NULL,
      date_breaks = "1 months",
      date_labels = "%b",
      limits = lubridate::ymd(c("2020-07-01", "2020-12-31")),
      expand = expansion()
    ) +
    theme(
      axis.ticks.length.x = unit(0.5, "cm"),
      axis.text.x = element_text(vjust = 7, hjust = -0.2)
    )
  p2
}


pall <- plotter(fdat2)

ggsave("trajectories-all.png", pall, height = 18, width = 24)

fdat3 <- fdat2 %>% mutate(ewk = lubridate::epiweek(forecast_date),
                          wkcyc = ewk %% 4)

splt <- split(fdat3, fdat3$wkcyc)
plots_staggered <- map(splt, plotter)

pnames <- paste0("trajectories-", names(plots_staggered), ".png")
map2(pnames, plots_staggered, ggsave, height = 18, width = 24)
