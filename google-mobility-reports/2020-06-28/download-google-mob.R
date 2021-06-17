#!/usr/bin/env Rscript

library(magrittr)

download_past_issue <- function(tstamp = "20200629") {
  base_url <- "http://timetravel.mementoweb.org/api/json"
  moburl <- "http://www.google.com/covid19/mobility"
  
  url <- paste(base_url, tstamp, moburl, sep = "/")
  
  resp <- httr::GET(url)
  if (httr::http_type(resp) != "application/json") {
    stop("API did not return json", call. = FALSE)
  }
  
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
  indexuri <- parsed$mementos$prev$uri
  csvuri <- stringr::str_replace(
    indexuri,
    "https://(www.)?google.com/covid19/mobility/?",
    "https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv"
  )
  
  destdir <- file.path("google-mobility-reports-wayback", tstamp)
  dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  curdir <- setwd(destdir)
  on.exit(setwd(curdir))
  download.file(url = csvuri, destfile = "Global_Mobility_Report.csv")
  system("head -n1 Global_Mobility_Report.csv > US-states-mobility.csv")
  system("grep \"US-[A-Z][A-Z]\" Global_Mobility_Report.csv >> US-states-mobility.csv")
  unlink("Global_Mobility_Report.csv")
  jsonlite::write_json(parsed, "response.json")
}

timestamps <- seq.Date(lubridate::ymd("2020-06-29"),
                       lubridate::ymd("2021-04-26"),
                       by = 7) %>%
  as.character() %>% stringr::str_remove_all("-")

purrr::map(timestamps[1:3], download_past_issue)
