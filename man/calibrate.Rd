% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibrate.R
\name{calibrate}
\alias{calibrate}
\title{Calibrate a Carbon-14 Date}
\usage{
calibrate(.data = NULL, measured_age, measured_age_error, df = Inf,
  curve = intcal13, name = NULL, cal_age_type = c("Calibrated BP",
  "Year AD"))
}
\arguments{
\item{.data}{An optional data frame from which age_14C, sd_14C, curve, and name could be sourced}

\item{measured_age, measured_age_error}{The age and error of the date. For most calibration
curves, this is in Radiocarbon years BP, but could also be a calendar age/error with
curve = "identity".}

\item{df}{Degrees of freedom for the t-distribution of the measured age.
Use Inf for the normal distribution.}

\item{curve}{A calibration curve (e.g. \link{intcal13}), or a vector of calibration curves.
Can be a character vector (e.g., "intcal13") or a calibration curve created using
\link{calibration_curve}. Use NULL for previously calibrated dates.}

\item{name}{Sample names to pass to the output.}

\item{cal_age_type}{The type of the calibrated age (must be
"Calibrated BP" or "Year AD"). See \link{in_year_ad} or \link{in_cal_bp} to switch
an existing calibration curve from one unit to the other.}
}
\value{
A data.frame
}
\description{
Calibrate a Carbon-14 Date
}
