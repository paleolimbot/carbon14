
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
      stop("curve must be an age_calibration_curve()")
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
    name = !!enquo(name),
    measured_age = !!enquo(measured_age),
    measured_age_error = !!enquo(measured_age_error),
    df = !!enquo(df),
    curve = !!curve
  )

  data$age_bp_distribution <- as_cdist(purrr::pmap(data, calibrate_single))
  if("name" %in% colnames(data)) {
    data$age_bp_distribution <- purrr::set_names(data$age_bp_distribution, data$name)
  }
  data$curve_name <- purrr::map_chr(data$curve, attr, "curve_name")
  data$measured_age_type <- purrr::map_chr(data$curve, attr, "measured_age_type")
  data$cal_age_type <- purrr::map_chr(data$curve, attr, "cal_age_type")
  data$curve <- NULL

  data
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
    curve <- get0(curve, envir = env, ifnotfound = NULL) %||%
      get0(curve, envir = environment(), ifnotfound = NULL)
    if(is.null(curve)) stop("Could not resolve curve '", curve, "'")
  }

  if(!inherits(curve, "age_calibration_curve")) {
    stop("Curve '", curve, "' is not an age_calibration_curve()")
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
