#!/usr/bin/env R

source("covidhub-common.R")
ddt <- Sys.getenv("ddt")
targ_dir <- file.path("hopkins", ddt)
output_dir <- "forecasts"
datf <- load_hopkins(targ_dir) %>% rename(true_value = "value")
forecast_path <- file.path("forecasts", paste0(ddt, "-CEID-Walk.csv"))
fcsts <- read_forecast(forecast_path) %>% mutate(target_type = str_remove(target, "^\\d+ "))

g1 <- ili %>% filter(year == 2020) %>% ggplot(aes(
  x = date,
  y = cili_encount,
)) + geom_line()+
  facet_grid(location ~ ., scales = "free_y") +
  labs(y = "Covid+ILI encounters") +
  theme(strip.text.y = element_text(angle = 0)) + 
  geom_line(aes(y = value), data = fcst, col = "blue") + 
  geom_ribbon(aes(y = NULL, ymin = llim, ymax = ulim), data = fcst, alpha = 0.3, fill = "blue")


pdata <- 
  datf %>% 
  filter(target_type == "wk ahead inc death") %>%
  filter(nchar(location) == 2)

fdata <- 
  fcsts %>% filter(target_type == "wk ahead inc death") %>%
  filter(nchar(location) == 2)
  
g2 <-
  ggplot(fdata, aes(x = target_end_date, y = value, color = as.factor(quantile))) + 
  geom_line() + 
  facet_wrap(vars(location), scales = "free_y", ncol = 10) +
  geom_line(data=pdata, aes(y = true_value), color = "black") +
  labs(x = "Date", y = "Weekly incident deaths") + 
  guides(color = guide_legend(title = "Probability for\nquantile"))
g2
