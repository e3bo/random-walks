#!/usr/bin/env Rscript

library(tidyverse)
source("covidhub-common.R")

fdts <- seq.Date(as.Date("2020-06-29"), as.Date("2021-04-26"), by = "7 days")
fitdirs <- paste0(fdts, "-fips06")
paths <- file.path("fits", fitdirs, "fit.RData")

loadpar <- function(p) {
  load(p)
  f <- as.list(wfixed)
  ests <- name_params(fit1$par, z, x, f, trans = FALSE)[[1]]
  sevec <- sqrt(diag(solve(h1)))
  sd <- name_params(sevec, z, x, f, trans = FALSE)[[1]]
  list(ests, sd)
}
pars <- map(paths, loadpar)

resests <- map(pars, 1) %>% map_dbl("residentialeffect")
ressds <- map(pars, 2) %>% map_dbl("residentialeffect")

q <- qnorm(1 - 0.05 / 2)

dph <- tibble(x = lubridate::ymd(fdts), y = resests, 
              lower = resests - ressds * q,
              upper = resests + ressds * q)

p <- ggplot(dph, aes(x, y)) + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) + 
  labs(x = "Forecast date", y = expression(beta[res])) + 
  theme_minimal()

ggsave("res-effect-by-forecast-date.png", p, width=5.2, height=4, dpi=600)
