
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
#'   \link{calibration_curve}.
#' @param name Sample names to pass to the output.
#'
#' @return A data.frame
#' @export
#'
calibrate <- function(.data = NULL, measured_age, measured_age_error, df = Inf,
                      curve = "intcal13", name = NULL) {
  curve_quo <- enquo(curve)
  curve_eval <- try(rlang::eval_tidy(curve_quo, rlang::as_data_mask(.data)), silent = TRUE)
  if(is.data.frame(curve_eval)) {
    if(!inherits(curve_eval, "age_calibration_curve")) stop("curve must be an age_calibration_curve()")
    curve <- list(curve_eval)
  } else if(is.vector(curve_eval)) {
    curve <- lapply(curve_eval, resolve_curve, env = parent.frame())
  } else {
    stop("curve must be a data frame or a vector")
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
  data$curve <- NULL

  data
}

resolve_curve <- function(curve, env = parent.frame()) {

  if(is.character(curve)) {
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
  dist <- cdist_item(m = measured_age, s = measured_age_error, df = df, dist = "t")

  if(!inherits(curve, "age_calibration_curve")) stop("curve is not an age_calibration_curve")

  cols <- attr(curve, "calibration")
  translate_distribution(
    x = curve[[cols$cal_age]],
    y = curve[[cols$measured_age]],
    y_sd = curve[[cols$measured_age_error]],
    dist = dist
  )
}
