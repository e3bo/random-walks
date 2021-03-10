#!/usr/bin/env Rscript

library(covidHubUtils)
suppressPackageStartupMessages(library(tidyverse))

lambda <- 1 / seq(0.001, 0.1, length.out = 10)
agrid <- c(0.94, 0.95)
par2name <- function(lambda, a){
  paste0("lambda", sprintf("%06.2f", lambda), 
         "-a", sprintf("%02.2f", a),
         "-CEID-InfectionKalmanEmp")
}
dirnames <- outer(lambda, agrid, par2name)

load_from_dir <- function(dname){
  dir(dname, full.names = TRUE) %>% load_forecast_files_repo() 
}
fdat2 <- map_dfr(dirnames, load_from_dir)

plotter <- function(dat) {
  p <- plot_forecast(
    dat,
    target_variable = "inc case",
    truth_source = "JHU",
    intervals = c(.5, .95),
    facet = a ~ lambda,
    facet_scales = "free_y",
    fill_by_model = TRUE,
    plot = FALSE
  )
  
  p2 <- p +
    scale_x_date(
      name = NULL,
      date_breaks = "1 months",
      date_labels = "%b",
      limits = lubridate::ymd(c("2020-07-01", "2021-03-01")),
      expand = expansion()
    ) +
    theme(
      axis.ticks.length.x = unit(0.5, "cm"),
      axis.text.x = element_text(vjust = 7, hjust = -0.2)
    )
  p2
}

make_page_plot <- function(page, pobj) {
  pobj + ggforce::facet_wrap_paginate(location~model, nrow = 2, ncol = 12, 
                                      page = page, scales = "free_y")
}

save_page_plot <- function(pp, page_num, dname){
  pstring <- sprintf("%02d", page_num)
  fname <- paste0(dname, "/page", pstring, ".png")
  ggsave(fname, pp, height = 7.5, width = 24) 
}

save_multi_page_plot <- function(gg, dirname){
 page_plots <- lapply(1:19, make_page_plot, pobj = gg)
 if(!dir.exists(dirname)){
   dir.create(dirname)
 }
 mapply(save_page_plot, page_plots, seq_along(page_plots), dname = dirname)
}

fdat2$lambda <- fdat2$model %>% str_extract("^lambda\\d+.\\d{2}-") %>% 
  str_remove("^lambda") %>% str_remove("-$") %>% as.numeric()
fdat2$a <- fdat2$model %>% str_extract("-a\\d{1}.\\d{2}-") %>% 
  str_remove("^-a") %>% str_remove("-$") %>% as.numeric()

splt <- split(fdat2, fdat2$location)
plots <- map(splt, plotter)




%>% group_by(location) %>% nest() %>% 
  mutate(plot = purrr::map(data, plotter))

pall <- plotter(fdat2)
save_multi_page_plot(pall, "trajectories-all")

#fdat3 <- fdat2 %>% mutate(ewk = lubridate::epiweek(forecast_date),
#                          wkcyc = ewk %% 4)
#splt <- split(fdat3, fdat3$wkcyc)
#plots_staggered <- map(splt, plotter)

#pnames <- paste0("trajectories-", names(plots_staggered))
#map2(plots_staggered, pnames, save_multi_page_plot)
