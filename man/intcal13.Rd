% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{intcal13}
\alias{intcal13}
\alias{marine13}
\alias{shcal13}
\alias{intcal09}
\alias{marine09}
\alias{intcal04}
\alias{marine04}
\alias{shcal04}
\title{Calibration curves}
\format{An object of class \code{age_calibration_curve} (inherits from \code{tbl_df}, \code{tbl}, \code{data.frame}) with 5141 rows and 5 columns.}
\source{
\url{http://www.radiocarbon.org/IntCal04.htm}
\url{http://www.radiocarbon.org/IntCal09.htm}
\url{http://www.radiocarbon.org/IntCal13.htm}
}
\usage{
intcal13

marine13

shcal13

intcal09

marine09

intcal04

marine04

shcal04
}
\description{
Data frame objects for the IntCal, Marine, and SHCal series
}
\details{
All contain the columns \code{cal_bp} (number of years before 1950),
\code{age_14C} (radiocarbon years BP), \code{error} (estimated error in radiocarbon years BP),
\code{delta14C}, and \code{sigma}.
}
\references{
Reimer PJ, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk Ramsey
C, Buck CE, Cheng H, Edwards RL, Friedrich M, Grootes PM, Guilderson TP,
Haflidason H, Hajdas I, Hatté C, Heaton TJ, Hoffmann DL, Hogg AG, Hughen KA,
Kaiser KF, Kromer B, Manning SW, Niu M, Reimer RW, Richards DA, Scott EM,
Southon JR, Staff RA, Turney CSM, van der Plicht J. 2013. IntCal13 and
Marine13 radiocarbon age calibration curves 0–50,000 years cal BP.
Radiocarbon 55(4):1869–1887.

Reimer PJ, Baillie MGL, Bard E, Bayliss A, Beck JW, Blackwell PG, Bronk
Ramsey C, Buck CE, Burr GS, Edwards RL, Friedrich M, Grootes PM, Guilderson
TP, Hajdas I, Heaton TJ, Hogg AG, Hughen KA, Kaiser KF, Kromer B, McCormac
FG, Manning SW, Reimer RW, Richards DA, Southon JR, Talamo S, Turney CSM, van
der Plicht J, Weyhenmeyer CE. 2009. IntCal09 and Marine09 radiocarbon age
calibration curves, 0–50,000 years cal BP. Radiocarbon 51(4):1111–50.

Reimer, Paula J., Mike G. L. Baillie, Edouard Bard, Alex Bayliss, J. Warren
Beck, Chanda J. H. Bertrand, Paul G. Blackwell, et al. "IntCal04 Terrestrial
Radiocarbon Age Calibration, 0-26 Cal Kyr BP." Radiocarbon 46, no. 3
(November 20, 2004): 1029–58.
}
\keyword{datasets}
