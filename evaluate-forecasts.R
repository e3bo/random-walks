library(dplyr)
library(evalcast)
library(magrittr)
library(ggplot2)
library(purrr)

# collect prediction cards
cpred <- list()
cpred$infkalm <- load_local_covidhub("CEID-InfectionKalman")

fdts <- map(cpred[[1]], attr, "forecast_date") %>% 
  map_chr(as.character) %>% unique()

fdtsrw <- c(paste0("2020-12-", c("07", "14", "21")),
            paste0("2020-11-", c("30", "23", "16", "09", "01")),
            paste0("2020-10-", c("25", "18", "11", "04")),
            paste0("2020-09-", c("27", "20", "13")))
cpred$rwf <- get_covidhub_predictions("CEID-Walk",
                                      forecast_dates = fdtsrw)
cpred$cef <- get_covidhub_predictions("COVIDhub-ensemble",
                                      forecast_dates = fdts)

cpred2 <- map(cpred, filter_predictions, geo_type = "state")

## trajectory_plots

ikselc <- cpred$infkalm %>% 
  filter_predictions(response_signal = "confirmed_incidence_num", 
                     ahead = c(1, 2))
rwselc <- cpred$rwf %>% 
  filter_predictions(response_signal = "confirmed_incidence_num", 
                     geo_type = "state", ahead = c(1, 2))

ceselc <- cpred$cef %>%
  filter_predictions(response_signal = "confirmed_incidence_num", 
                     geo_type = "state", ahead = c(1, 2))

evalcast:::intersect_locations(c(ikselc, rwselc)) %>% 
  plot_trajectory(first_day = "2020-09-06", last_day = "2020-12-31")

evalcast:::intersect_locations(c(ikselc, ceselc)) %>% 
  plot_trajectory(first_day = "2020-09-06", last_day = "2020-12-31")

# case forecast evaluations
pull_cases <- function(x, h = c(1, 2)){
  filter_predictions(x, 
                     response_signal = "confirmed_incidence_num", 
                     ahead = h)
}

cpred3 <- cpred2[-5] %>% map(pull_cases)

case_evals <- cpred3 %>% map(evaluate_predictions, backfill_buffer = 4)

ci <- list()
ci$h1 <- evalcast:::intersect_locations(map(case_evals, 1))
ci$h2 <- evalcast:::intersect_locations(map(case_evals, 2))

add_ci <- function(p){
  # confidence interval for median calculated by `boxplot.stats`
  # idea from koshke at https://stackoverflow.com/a/8135865
  
  f <- function(x) {
    ans <- boxplot.stats(x)
    data.frame(ymin = ans$conf[1], ymax = ans$conf[2], y = ans$stat[3])
  }
  p + stat_summary(fun.data = f, geom = "crossbar", color = NA, 
                   fill = "skyblue", size = 5, alpha = 0.5, width = 0.75) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

plot_measure(ci$h1, "ae") %>% add_ci()
plot_measure(ci$h1, "wis") %>% add_ci()

plot_measure(ci$h2, "ae") %>% add_ci()
plot_measure(ci$h2, "wis") %>% add_ci()

my_dotplot <- function(sc, err_name = "wis"){
  nm <- attr(sc, "name_of_forecaster")
  rsp <- attr(sc, "response")$signal
  sc2 <- sc %>% mutate(location = covidcast::fips_to_abbr(location))
  ordered_levels <- sc2 %>% group_by(.data$location) %>% 
    summarize(avg_err = mean(.data[[err_name]])) %>% 
    arrange(.data$avg_err)
  sc2$location <- factor(sc2$location, levels = ordered_levels$location)

  sc2 %>% ggplot(aes(x = wis, y = location)) + geom_point(alpha = 0.5) + 
    scale_x_continuous(trans = "log1p") + 
    ggtitle(nm, subtitle = rsp) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}

my_dotplot(ci$h1[[1]], "wis")
my_dotplot(ci$h1[[3]], "wis")
  
plot_width(ci$h1, levels = 0.9)
plot_width(ci$h2, levels = 0.9)

map(ci$h1, plot_calibration)

## regression modeling

fit_mod <- function(sc){
  df <- dplyr::bind_rows(sc, .id = "forecaster")
  m <- lme4::lmer(log(wis + 1) ~ forecaster + (1|end), data = df)
  ci <- confint(m)
  list(model = m, conf = ci)  
}

(ci_mods <- map(ci, fit_mod))