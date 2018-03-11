
#' Create a calibration curve
#'
#' @param .data An optional data frame
#' @param cal_age The age in calibrated years BP (before 1950)
#' @param measured_age The measured age (usually 14C years  measurements)
#' @param measured_age_error,measured_age_min,measured_age_max Error on the measured age, plus or minus
#'
#' @return A tibble classed as an age_calibration_curve
#' @export
#'
calibration_curve <- function(.data = NULL, cal_age, measured_age,
                              measured_age_error = NULL,
                              measured_age_min = NA_real_, measured_age_max = NA_real_) {
  data <- data_eval(
    .data, cal_age = rlang::enquo(cal_age), measured_age = rlang::enquo(cal_age),
    measured_age_min = rlang::enquo(measured_age_min),
    measured_age_max = rlang::enquo(measured_age_max),
    measured_age_error = rlang::enquo(measured_age_error)
  )

  if("measured_age_error" %in% colnames(data)) {
    attr(data, "calibration") <- list(
      cal_age = "cal_age",
      measured_age = "measured_age",
      measured_age_min = "measured_age_min",
      measured_age_max = "measured_age_max",
      measured_age_error = "measured_age_error"
    )
  } else {
    attr(data, "calibration") <- list(
      cal_age = "cal_age",
      measured_age = "measured_age",
      measured_age_min = "measured_age_min",
      measured_age_max = "measured_age_max"
    )
  }

  class(data) <- c("age_calibration_curve", class(data))

  data
}

#' Print, plot an age-depth calibration object
#'
#' @param x An age-depth calibration
#' @param type,xlab,ylab,... Passed to \link[graphics]{plot}
#' @param measured_age_limits Plot high and low values of the measured age
#'
#' @return The input, invisibly
#' @export
#'
print.age_calibration_curve <- function(x, ...) {
  cat("<age_calibration_curve>\n")
  mapping <- unlist(attr(x, "calibration"))
  cat(paste(names(mapping), mapping, sep = " = ", collapse = ", "))
  cat("\n\n")
  print(tibble::as_tibble(unclass(x)), ...)
  invisible(x)
}

#' @rdname print.age_calibration_curve
#' @export
plot.age_calibration_curve <- function(x, type = "l", xlab = "Age (cal BP)",
                                       ylab = "Measured Age", ..., measured_age_limits = TRUE) {
  cols <- attr(x, "calibration")
  meas <- x[[cols$measured_age]]
  graphics::plot(x[[cols$cal_age]], meas, type = type,
                 xlab = xlab, ylab = ylab, ...)
  if(measured_age_limits) {
    if(!is.null(cols$measured_age_error)) {
      meas_max <- meas + x[[cols$measured_age_error]]
      meas_min <- meas - x[[cols$measured_age_error]]
    } else {
      meas_max <- x[[cols$measured_age_max]]
      meas_min <- x[[cols$measured_age_min]]
    }

    graphics::lines(x[[cols$cal_age]], meas_max, lty = 3)
    graphics::lines(x[[cols$cal_age]], meas_min, lty = 3)
  }
}
