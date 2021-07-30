#!/usr/bin/env Rscript

library(tidyverse)

system("dvc metrics show > metrics.csv")

mets <- read.table("metrics.csv", header = TRUE, na.strings = "-") %>%
  mutate(date = str_extract(Path, "\\d{4}-\\d{2}-\\d{2}"),
         loc = str_extract(Path, "fips\\d{2}"),
         loc2 = str_remove(loc, "^fips"),
         wday = lubridate::wday(date, label = TRUE)) 

m2 <- mets %>% 
  filter(wday == "Mon" & loc2 == "06") %>% 
  filter(date >= "2020-06-29" & date <= "2021-04-26") %>%
  filter(fit1hessposdef == "True")

m3 <- m2 %>%
  select(date, fit1mase_cases, fit1mase_cases_w, 
         fit1mase_hosps, fit1mase_hosps_w,
         fit1mase_death, fit1mase_death_w) %>%
  pivot_longer(-date) %>%
  mutate(Weekly = str_detect(name, "_w$"),
         var = str_extract(name, "_[a-z]{5}"),
         var2 = str_remove(var, "^_"),
         var3 = case_when(var2 == "cases" ~ "Cases", 
                          var2 == "death" ~ "Deaths", 
                          var2 == "hosps" ~ "Hospital admissions"),
         var4 = factor(var3, 
                       levels = c("Cases", "Hospital admissions", "Deaths")))

p <- ggplot(filter(m3, !is.na(value)), 
       aes(x = lubridate::ymd(date), y = value, shape = Weekly)) + 
  geom_point() + 
  ylim(c(0, 1)) + 
  facet_wrap(~var4, ncol = 1) + 
  labs(x = "Forecast date", y = "MASE") +
  theme(legend.position="top") +
  theme_minimal()

ggsave("california-mase.png", p, width = 5.2, dpi = 600, units = "in")

m4 <- m2 %>% select(date, fit1walltime) %>% 
  mutate(date = lubridate::ymd(date),
         hours = fit1walltime / 60 / 60)

p2 <- m4 %>% ggplot(aes(x = date, y = hours)) + 
  geom_point() +  
  labs(x = "Forecast date", y = "Fitting time (hours)") +
  theme_minimal()

ggsave("california-fitting-time.png", p2, dpi = 600, width = 5.2, units = "in",
       height = 4)
