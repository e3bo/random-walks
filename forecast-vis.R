#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

source("covidhub-common.R")
plotd <- "visuals"
fdt <- Sys.getenv("fdt")
ddt <- Sys.getenv("ddt")
targ_dir <- file.path("hopkins", ddt)
output_dir <- "forecasts"
datf <- load_hopkins(targ_dir) %>% rename(true_value = "value")
forecast_path <- file.path("forecasts", paste0(fdt, "-CEID-Walk.csv"))
fcsts <- read_forecast(forecast_path) %>% 
  mutate(target_type = str_remove(target, "^\\d+ "))

plot_fcsts <- function(targ, datf, fcsts) {
  if (targ == "inc death") {
    tt <- "wk ahead inc death"
    ylab <- "Weekly incident deaths"
    maxc <- 2
    nrow <- 10
  } else if (targ == "cum death"){
    tt <- "wk ahead cum death"
    ylab <- "Weekly cumulative deaths"
    maxc <- 2
    nrow <- 10
  } else if (targ == "inc case"){
    tt <- "wk ahead inc case"
    ylab <- "Weekly incident cases"
    maxc <- 5
    nrow <- 10 * 50
  }
  pdata <-
    datf %>%
    filter(target_type == tt) %>%
    filter(nchar(location) <= maxc)
  fdata <-
    fcsts %>% filter(target_type == tt) %>%
    filter(nchar(location) <= maxc)
  g2 <-
    ggplot(fdata, aes(
      x = target_end_date,
      y = value,
      color = as.factor(quantile)
    )) +
    geom_line() +
    geom_line(data = pdata, aes(y = true_value), color = "black") +
    labs(x = "Date", y = ylab) +
    guides(color = guide_legend(title = "Probability for\nquantile")) +
    facet_wrap(vars(location), scales = "free_y", nrow = nrow)
  g2
}

dir.create(plotd)
inc_death_g <- plot_fcsts("inc death", datf = datf, fcsts = fcsts)
date_spec <- paste0("fdt", fdt, "-ddt", ddt, "-")

file.path(plotd, paste0(date_spec, "inc-death-forecasts.png")) %>% 
  ggsave(plot = inc_death_g, width = 9, height = 15)

cum_death_g <- plot_fcsts("cum death", datf = datf, fcsts = fcsts)
file.path(plotd, paste0(date_spec, "cum-death-forecasts.png")) %>% 
  ggsave(plot = cum_death_g, width = 9, height = 15)

inc_case_g <- plot_fcsts("inc case", datf = datf, fcsts = fcsts)
file.path(plotd, paste0(date_spec, "inc-case-forecasts.pdf")) %>% 
  ggsave(plot = inc_case_g, width = 9, height = 900, limitsize = FALSE)
