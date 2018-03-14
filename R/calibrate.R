
#' Calibrate a Carbon-14 Date
#'
#' @param .data An optional data frame from which age_14C, sd_14C, curve, and name could be sourced
#' @param measured_age,measured_age_error The age and error of the date. For most calibration
#'   curves, this is in Radiocarbon years BP, but could also be a calendar age/error with
#'   curve = "identity".
#' @param df Degrees of freedom for the t-distribution of the measured age.
#'   Use Inf for the normal distribution.
#' @param curve A calibration curve (e.g. \link{intcal13}), or a vector of calibration curves.
#'   Can be a character vector (e.g., "intcal13") or a calibration curve created using
#'   \link{calibration_curve}. Use NULL for previously calibrated dates.
#' @param name Sample names to pass to the output.
#' @param cal_age_type The type of the calibrated age (must be
#'   "Calibrated BP" or "Year AD"). See \link{in_year_ad} or \link{in_cal_bp} to switch
#'   an existing calibration curve from one unit to the other.
#'
#' @return A data.frame
#' @export
#'
calibrate <- function(.data = NULL, measured_age, measured_age_error, df = Inf,
                      curve = intcal13, name = NULL,
                      cal_age_type = c("Calibrated BP", "Year AD")) {

  curve_quo <- enquo(curve)
  curve_eval <- try(rlang::eval_tidy(curve_quo, rlang::as_data_mask(.data)), silent = TRUE)
  cal_age_type <- match.arg(cal_age_type)

  if(is.null(curve_eval)) {
    curve <- list(null_calibration_curve(cal_age_type))
  } else if(is.data.frame(curve_eval)) {
    if(!inherits(curve_eval, "age_calibration_curve")) {
      stop("curve must be a calibration_curve()")
    }
    curve <- list(curve_eval)
  } else if(is.vector(curve_eval)) {
    curve <- lapply(curve_eval, resolve_curve, null_age_type = cal_age_type, env = parent.frame())
  } else {
    stop("curve must be a data frame, a vector, or NULL")
  }

  if(cal_age_type == "Calibrated BP") {
    curve <- lapply(curve, in_cal_bp)
  } else {
    curve <- lapply(curve, in_year_ad)
  }

  data <- data_eval(
    .data,
    name = !!enquo(name),
    measured_age = !!enquo(measured_age),
    measured_age_error = !!enquo(measured_age_error),
    df = !!enquo(df),
    curve = !!curve
  )

  data$cal_age <- as_cdist(purrr::pmap(data, calibrate_single))
  if("name" %in% colnames(data)) {
    data$cal_age <- purrr::set_names(data$cal_age, data$name)
  }
  data$curve_name <- purrr::map_chr(data$curve, attr, "curve_name")
  data$measured_age_type <- purrr::map_chr(data$curve, attr, "measured_age_type")
  data$cal_age_type <- purrr::map_chr(data$curve, attr, "cal_age_type")
  class(data$curve) <- "age_calibration_curve_list"

  class(data) <- c("calibrate_result", class(data))
  data
}

#' @export
#' @importFrom dplyr slice
slice.calibrate_result <- function(.data, ...) {
  result <- NextMethod()
  class(result) <- c("calibrate_result", class(result))
  result
}

#' @export
#' @importFrom dplyr filter
filter.calibrate_result <- function(.data, ...) {
  result <- NextMethod()
  class(result) <- c("calibrate_result", class(result))
  result
}

#' @export
dplyr::filter

#' @export
`[.calibrate_result` <- function(x, i, j, ...) {
  result <- NextMethod()
  if(inherits(result, "data.frame")) {
    class(result) <- c("calibrate_result", class(result))
  }
  result
}


#' Plot a calibration result
#'
#' @param x An age calibration
#' @param max_plot Maximum number of items to plot
#' @param n_col Number of columns
#' @param ... Passed to \link[graphics]{plot}
#'
#' @export
#'
plot.calibrate_result <- function(x, ..., max_plot = 9, n_col = 3) {
  if(nrow(x) > max_plot) {
    message("Plotting first ", max_plot, " dates. Use max_plot = Inf to plot all dates.")
  }
  x <- utils::head(x, max_plot)
  n <- nrow(x)
  rows <- (n + n_col - 1) %/% n_col
  if(rows == 1) {
    n_col <- n
  }
  withr::with_par(list(mfrow = c(rows, n_col)), {
    purrr::map(purrr::transpose(x), plot_single_result, ..., env = parent.frame())
  })
  invisible(NULL)
}

plot_single_result <- function(x, ..., xlim = NULL, ylim = NULL, eps = 1e-4, env = parent.frame()) {
  curve <- x$curve
  cols <- attr(curve, "calibration")

  # get xlim
  range_cal <- range(x$cal_age, eps = eps)

  # get ylim
  cal_age <- curve[[cols$cal_age]]
  age_filter <- (cal_age >= range_cal[1]) & (cal_age <= range_cal[2])
  meas <- curve[[cols$measured_age]]

  if(!any(age_filter)) {
    # no points on the calibration curve in range
    # (probably the null calibration curve)
    range_meas <- stats::approx(cal_age, meas, xout = range_cal)$y
  } else {
    meas <- curve[[cols$measured_age]]
    meas_max <- meas + curve[[cols$measured_age_error]]
    meas_min <- meas - curve[[cols$measured_age_error]]
    range_meas <- range(meas_max[age_filter], meas_min[age_filter])
  }

  # get plot title
  if(!is.null(x$name)) {
    title <- sprintf("%s (%s)", x$name, x$curve_name)
  } else {
    title <- x$curve_name
  }

  # plot curve
  graphics::plot(
    curve,
    xlim = xlim %||% range_cal,
    ylim = ylim %||% range_meas,
    title = title,
    ...
  )

  # plot distributions
  plot_range <- graphics::par("usr")

  test_bp <- seq(range_cal[1], range_cal[2], length.out = 1024)
  cal_dens <- density(x$cal_age, test_bp)
  cal_dens <- cal_dens / max(cal_dens)
  graphics::lines(test_bp, cal_dens * (plot_range[4] - plot_range[3]) * 0.25 + plot_range[3],
                  col = "red")

  date <- dist_item_parameterized("t", list(df = x$df, m = x$measured_age, s = x$measured_age_error))
  range_14c <- quantile(date, c(eps, 1 - eps))
  test_14c <- seq(range_14c[1], range_14c[2], length.out = 1024)
  date_dens <- density(date, test_14c)
  date_dens <- date_dens / max(date_dens)
  graphics::lines(date_dens * (plot_range[2] - plot_range[1]) * 0.25 + plot_range[1], test_14c,
                  col = "blue")

  # calibrated range text
  graphics::mtext(
    sprintf("%0.0f %0.0f-%0.0f (95%%)", weighted.mean(x$cal_age),
            quantile(x$cal_age, 0.05), quantile(x$cal_age, 0.95)),
    cex = 0.5,
    side = 3
  )

  # measured range text
  txt <- sprintf("%0.0f %%+-%% %0.0f", x$measured_age, x$measured_age_err)
  graphics::mtext(
    parse(text = txt),
    cex = 0.5,
    side = 4
  )
}


#' Resolve a calibration curve
#'
#' @param curve A calibration curve, a character vector with an object name in this package
#'   or the calling environment, or NULL for the identity \link{null_calibration_curve}.
#' @param null_age_type The input age type for the NULL curve.
#' @param env The environment to look for character names of calibration curves
#'
#' @return A calibration curve
#' @export
#'
resolve_curve <- function(curve, null_age_type = "Calibrated BP", env = parent.frame()) {

  if(is.null(curve)) {
    curve <- null_calibration_curve(null_age_type)
  } else if(is.character(curve)) {
    if(curve == "identity") {
      curve <- null_calibration_curve(null_age_type)
    } else {
      curve <- get0(curve, envir = env, ifnotfound = NULL) %||%
        get0(curve, envir = environment(), ifnotfound = NULL)
    }
    if(is.null(curve)) stop("Could not resolve curve '", curve, "'")
  }

  if(!inherits(curve, "age_calibration_curve")) {
    stop("Curve '", curve, "' is not a calibration_curve()")
  }

  curve
}

calibrate_single <- function(measured_age, measured_age_error, df, curve, name = NULL) {
  dist <- dist_item_parameterized("t", list(m = measured_age, s = measured_age_error, df = df))

  if(!inherits(curve, "age_calibration_curve")) stop("curve is not an age_calibration_curve")

  translate_distribution(age_calibration_curve = curve, dist = dist)
}

#' Translate a distribution using a calibration curve
#'
#' @param age_calibration_curve An age_calibration_curve (\link{calibration_curve}) or NULL
#'   for the identity curve
#' @param dist A \link{cdist} containing source densities
#' @param eps Initial densities smaller than this number will not be considered
#' @param n The number of equally-spaced x values at which density should be calculated
#'
#' @importFrom rlang !!
#' @importFrom rlang enquo
#' @importFrom rlang .data
#' @export
translate_distribution <- function(age_calibration_curve, dist, eps = 1e-8, n = 512) {

  if(inherits(dist, "cdist")) {
    # vectorize
    return(
      new_cdist(purrr::map(
        dist, translate_distribution,
        age_calibration_curve = age_calibration_curve,
        eps = eps,
        n = n
      ))
    )
  }

  # NULL age calibration curve == identity
  if(inherits(age_calibration_curve, "null_calibration_curve")) {
    if(attr(age_calibration_curve, "cal_age_type") ==
       attr(age_calibration_curve, "measured_age_type")) {
      # no transformation
      return(dist)
    } else {
      # 1950 - transformation
      return(1950 - dist)
    }
  }

  # map columns
  cols <- attr(age_calibration_curve, "calibration")
  data <- tibble::tibble(
    x = age_calibration_curve[[cols$cal_age]],
    y = age_calibration_curve[[cols$measured_age]],
    y_sd = age_calibration_curve[[cols$measured_age_error]]
  )

  # calculate density along the measured axis
  data$density_y <- density(dist, data$y)

  # tries to use parameterized values, defaults to the middle 68% of data
  nominal_sd <- dist$dist_info$sd %||%
    dist$dist_info$s %||%
    (diff(quantile(dist, c(0.16, 0.84))) / 2)

  xrange <- scales::expand_range(range(data$x[data$density_y > eps], na.rm = TRUE), mul = 1)

  dist_values <- tibble::tibble(
    x = seq(xrange[1], xrange[2], length.out = n)
  )
  dist_values$y <- stats::approx(data$x, data$y, xout = dist_values$x)$y
  dist_values$y_sd <- stats::approx(data$x, data$y_sd, xout = dist_values$x)$y
  dist_values$tau <- dist_values$y_sd + nominal_sd ^ 2
  dist_values$density <- density(dist, values = dist_values$y) / sqrt(dist_values$tau)
  dist_values$density <- dist_values$density / sum(dist_values$density, na.rm = TRUE)
  dist_values <- dplyr::filter(dist_values, !is.na(.data$density), !is.na(.data$x))

  dist_item_custom(values = dist_values$x, densities = dist_values$density)
}
