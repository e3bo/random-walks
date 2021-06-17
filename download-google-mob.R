#!/usr/bin/env Rscript

library(magrittr)

resp <- httr::GET("http://timetravel.mementoweb.org/timemap/json/http://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv")

parsed <-
  httr::content(resp, as = "parsed", simplifyVector = TRUE)

if (httr::http_error(resp)) {
  stop(
    sprintf(
      "memento API request failed [%s]\nresult : %d\n%s\n<%s>",
      httr::status_code(resp),
      parsed$result,
      parsed$message,
      "http://mementoweb.org/guide/howto/"
    ),
    call. = FALSE
  )
}

is_archiv <-
  purrr::map_lgl(parsed$mementos[[1]]$uri,
                 ~ stringr::str_detect(.x, "^https://web.archive.org"))

times_avail <- parsed$mementos[[1]]$datetime[is_archiv]

download_past_issue <- function(dname, csvuri) {
  destdir <- file.path("google-mobility-reports-wayback", dname)
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  curdir <- setwd(destdir)
  on.exit(setwd(curdir))
  download.file(url = csvuri, destfile = "Global_Mobility_Report.csv")
  system("head -n1 Global_Mobility_Report.csv > US-states-mobility.csv")
  system("grep \"US-[A-Z][A-Z]\" Global_Mobility_Report.csv >> US-states-mobility.csv")
  unlink("Global_Mobility_Report.csv")
  cat(csvuri, file = "uri")
}

timestamps <- seq.Date(lubridate::ymd("2020-06-29"),
                       lubridate::ymd("2021-04-26"),
                       by = 7) %>% lubridate::as_datetime()

find_previous_ind <- function(x){
  match(TRUE, x < lubridate::as_datetime(parsed$mementos[[1]]$datetime[is_archiv])) - 1
}

inds <- purrr::map_dbl(timestamps, find_previous_ind)
uris <- parsed$mementos[[1]]$uri[is_archiv][inds]
dirnames <- timestamps %>% as.character()

purrr::map2(dirnames, uris, download_past_issue)
