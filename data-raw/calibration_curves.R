
library(tidyverse)
# http://www.radiocarbon.org/IntCal13.htm

files <- tibble(
  urls = c(
    "http://www.radiocarbon.org/IntCal13%20files/intcal13.14c",
    "http://www.radiocarbon.org/IntCal13%20files/marine13.14c",
    "http://www.radiocarbon.org/IntCal13%20files/shcal13.14c",
    "http://www.radiocarbon.org/IntCal09%20files/intcal09.14c",
    "http://www.radiocarbon.org/IntCal09%20files/marine09.14c",
    "http://www.radiocarbon.org/IntCal04%20files/intcal04.14c",
    "http://www.radiocarbon.org/IntCal04%20files/marine04.14c",
    "http://www.radiocarbon.org/IntCal04%20files/shcal04.14c"
  )
) %>%
  extract(urls, "file_name", "([^/]*?)$", remove = FALSE) %>%
  mutate(file_loc = file.path("data-raw", file_name)) %>%
  mutate(object_name = str_replace(file_name, ".14c$", ""))

#CAL BP, 14C age,Error,Delta 14C,Sigma
#      , YR BP  ,YR BP,per mil  ,per mil

for(item in transpose(files)) {
  curl::curl_download(item$url, item$file_loc)
  df <- read_csv(
    item$file_loc,
    comment = "#",
    col_names = c("cal_bp", "age_14C", "error", "delta_14C", "sigma"),
    col_types = cols(.default = col_double())
  )
  class(df) <- c("age_calibration_curve", class(df))
  attr(df, "spec") <- NULL
  attr(df, "calibration") <- list(
    cal_age = "cal_bp",
    measured_age = "age_14C",
    measured_age_error = "error"
  )
  assign(item$object_name, df)
  do.call(devtools::use_data, list(rlang::sym(item$object_name), overwrite = TRUE))
}
