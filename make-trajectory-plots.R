#!/usr/bin/env Rscript


library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

fdat0 <- readRDS("other-model-forecasts.rds")
fdat1 <-
  dir("CEID-InfectionKalman", full.names = TRUE) %>% load_forecast_files_repo()
fdat2 <- bind_rows(fdat0, fdat1)

plotter <- function(dat) {
  p <- plot_forecast(
    dat,
    target_variable = "inc case",
    truth_source = "JHU",
    intervals = c(.5, .95),
    facet = location ~ model,
    fill_by_model = TRUE,
    plot = FALSE
  )
  
  p2 <- p +
    scale_x_date(
      name = NULL,
      date_breaks = "1 months",
      date_labels = "%b",
      limits = lubridate::ymd(c("2020-09-01", "2020-12-31")),
      expand = expansion()
    ) +
    theme(
      axis.ticks.length.x = unit(0.5, "cm"),
      axis.text.x = element_text(vjust = 7, hjust = -0.2)
    )
  p2
}


pall <- plotter(fdat2)

ggsave("trajectories-all.png", pall, width = 8)

fdat3 <- fdat2 %>% mutate(ewk = lubridate::epiweek(forecast_date),
                          wkcyc = ewk %% 4)

splt <- split(fdat3, fdat3$wkcyc)
plots_staggered <- map(splt, plotter)

pnames <- paste0("trajectories-", names(plots_staggered), ".png")
map2(pnames, plots_staggered, ggsave, width = 8)
