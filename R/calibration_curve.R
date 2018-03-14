
#' Create a calibration curve
#'
#' @param .data An optional data frame
#' @param cal_age The age in calibrated years BP (before 1950)
#' @param measured_age The measured age (usually 14C years  measurements)
#' @param measured_age_error Error on the measured age, plus or minus
#' @param name a name for the calibration curve
#' @param cal_age_type The type of the calibrated age (must be
#'   "Calibrated BP" or "Year AD"). See \link{in_year_ad} or \link{in_cal_bp} to switch
#'   an existing calibration curve from one unit to the other.
#' @param measured_age_type The type of measured age (usually 'Radiocarbon Years BP')
#'
#' @return A tibble classed as an age_calibration_curve
#' @export
#'
calibration_curve <- function(.data = NULL, cal_age, measured_age,
                              measured_age_error = 0, name = NULL,
                              cal_age_type = c("Calibrated BP", "Year AD"),
                              measured_age_type = NULL) {
  data <- data_eval(
    .data, cal_age = !!rlang::enquo(cal_age), measured_age = !!rlang::enquo(measured_age),
    measured_age_error = !!rlang::enquo(measured_age_error)
  )

  attr(data, "calibration") <- list(
    cal_age = "cal_age",
    measured_age = "measured_age",
    measured_age_error = "measured_age_error"
  )
  if(is.null(name)) {
    name <- "<unnamed>"
  } else {
    if(!is.character(name) || (length(name) != 1)) stop("name must be a character vector of length 1")
  }
  if(!is.null(measured_age_type)) {
    if(!is.character(measured_age_type) || (length(measured_age_type) != 1)) {
      stop("measured_age_type must be a character vector of length 1")
    }
  }

  attr(data, "curve_name") <- name %||% "<unnamed>"
  attr(data, "measured_age_type")  <- measured_age_type
  attr(data, "cal_age_type") <- match.arg(cal_age_type)

  class(data) <- c("age_calibration_curve", class(data))

  data
}

#' Create an identity calibration curve
#'
#' @param age_type The input age type (Calibrated BP or Year AD)
#' @param range Range of ages that are valid
#'
#' @return A calibration curve
#' @export
#'
null_calibration_curve <- function(
  age_type = c("Calibrated BP", "Year AD"),
  range = c(-70000, 70000)
) {
  age_type <- match.arg(age_type)
  curve <- calibration_curve(
    measured_age = range,
    cal_age = range,
    measured_age_type = age_type,
    cal_age_type = age_type,
    name = "identity"
  )
  # classify the curve such that it can be handled specially
  class(curve) <- c("null_calibration_curve", class(curve))
  curve
}

#' Convert an age calibration curve to a different calibrated age type
#'
#' @param x An existing age calibration curve (\link{calibration_curve}) or numeric vector
#'
#' @return A new age calibration curve
#' @export
#'
in_cal_bp <- function(x) {
  UseMethod("in_cal_bp")
}

#' @rdname in_cal_bp
#' @export
in_year_ad <- function(x) {
  UseMethod("in_year_ad")
}

#' @rdname in_cal_bp
#' @export
in_year_ad.default <- function(x) {
  if(is.null(x)) return(NULL)
  1950 - x
}

#' @rdname in_cal_bp
#' @export
in_cal_bp.default <- function(x) {
  if(is.null(x)) return(NULL)
  1950 - x
}

#' @rdname in_cal_bp
#' @export
in_cal_bp.age_calibration_curve <- function(x) {
  type <- attr(x, "cal_age_type")
  if(is.null(type)) stop("NULL calibrated age type")

  if(type == "Calibrated BP") {
    x
  } else if(type == "Year AD") {
    cal_col <- attr(x, "calibration")$cal_age
    x[[cal_col]] <- 1950 - x[[cal_col]]
    attr(x, "cal_age_type") <- "Calibrated BP"
    x
  } else {
    stop("Unknown calibrated age type: ", type)
  }
}

#' @rdname in_cal_bp
#' @export
in_year_ad.age_calibration_curve <- function(x) {
  type <- attr(x, "cal_age_type")
  if(is.null(type)) stop("NULL calibrated age type")

  if(type == "Year AD") {
    x
  } else if(type == "Calibrated BP") {
    cal_col <- attr(x, "calibration")$cal_age
    x[[cal_col]] <- 1950 - x[[cal_col]]
    attr(x, "cal_age_type") <- "Year AD"
    x
  } else {
    stop("Unknown calibrated age type: ", type)
  }
}

#' Print, plot an age-depth calibration object
#'
#' @param x An age-depth calibration
#' @param type,xlab,ylab,... Passed to \link[graphics]{plot}
#' @param measured_age_limits Plot high and low values of the measured age
#' @param title A title for the plot
#'
#' @return The input, invisibly
#' @export
#'
print.age_calibration_curve <- function(x, ...) {
  name <- attr(x, "curve_name")
  if(!is.null(name)) {
    cat(sprintf("<age_calibration_curve: '%s'>\n", name))
  } else {
    cat("<age_calibration_curve: <unnamed>>\n")
  }
  mapping <- unlist(attr(x, "calibration"))
  cat(paste(names(mapping), mapping, sep = " = ", collapse = ", "))
  cat("\n")
  cat(attr(x, "measured_age_type") %||% "<unknown>", "=>", attr(x, "cal_age_type"))
  cat("\n\n")
  print(tibble::as_tibble(unclass(x)), ...)
  invisible(x)
}

#' @rdname print.age_calibration_curve
#' @export
plot.age_calibration_curve <- function(x, type = "l",
                                       xlab = attr(x, "cal_age_type"),
                                       ylab = attr(x, "measured_age_type") %||% "Measured Age",
                                       title = attr(x, "curve_name"),
                                       ...,
                                       measured_age_limits = TRUE) {
  cols <- attr(x, "calibration")
  meas <- x[[cols$measured_age]]
  graphics::plot(x[[cols$cal_age]], meas, type = type,
                 xlab = xlab, ylab = ylab, ...)
  if(measured_age_limits && any(x[[cols$measured_age_error]] != 0)) {
    meas_max <- meas + x[[cols$measured_age_error]]
    meas_min <- meas - x[[cols$measured_age_error]]

    graphics::lines(x[[cols$cal_age]], meas_max, lty = 3)
    graphics::lines(x[[cols$cal_age]], meas_min, lty = 3)
  }

  if(!is.null(title)) {
    graphics::title(title)
  }
}

#' @export
as.character.age_calibration_curve <- function(x, ...) {
  xlab <- attr(x, "cal_age_type")
  ylab <- attr(x, "measured_age_type") %||% "Measured Age"
  sprintf("<age_calibration_curve: '%s'>", attr(x, "curve_name"))
}

#' @export
format.age_calibration_curve <- function(x, ...) {
  as.character(x, ...)
}

#' @export
as.character.age_calibration_curve_list <- function(x, ...) {
  purrr::map_chr(x, as.character, ...)
}

#' @export
format.age_calibration_curve_list <- function(x, ...) {
  purrr::map_chr(x, format, ...)
}

#' @export
print.age_calibration_curve_list <- function(x, ...) {
  cat("<age_calibration_curve_list>\n")
  print(format(x), quote = FALSE, ...)
  invisible(x)
}

#' @export
`[.age_calibration_curve_list` <- function(x, i, ...) {
  result <- NextMethod()
  class(result) <- "age_calibration_curve_list"
  result
}
