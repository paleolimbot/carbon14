% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration_curve.R
\name{null_calibration_curve}
\alias{null_calibration_curve}
\title{Create an identity calibration curve}
\usage{
null_calibration_curve(age_type = c("Calibrated BP", "Year AD"),
  range = c(-70000, 70000))
}
\arguments{
\item{age_type}{The input age type (Calibrated BP or Year AD)}

\item{range}{Range of ages that are valid}
}
\value{
A calibration curve
}
\description{
Create an identity calibration curve
}
