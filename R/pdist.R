
#' Create custom probability distributions
#'
#' @param ... Parameters to pass to density and quantile functions
#' @param dist Suffix for distribution (e.g., "norm", "t", etc.)
#' @param density_function Function like \link[stats]{dnorm}. Overrides dist argument.
#' @param quantile_function Function like \link[stats]{qnorm}. Overrides dist argument.
#' @param dist_info Meta information about the distribution
#'
#' @importFrom purrr %||%
#' @export
cdist_item <- function(..., dist = "norm", density_function = NULL,
                       quantile_function = NULL, dist_info = list()) {
  params <- tibble::lst(...)
  stats_ns <- getNamespace("stats")
  if(!is.null(dist)) {
    density_function_def <- stats_ns[[paste0("d", dist)]]
    quantile_function_def <- stats_ns[[paste0("q", dist)]]
    dist_info <- c(dist_info, list(dist = dist), params)
  } else {
    density_function_def <- NULL
    quantile_function_def <- NULL
  }

  item <- new_cdist_item(list(
    density_function = density_function %||% density_function_def,
    quantile_function = quantile_function %||% quantile_function_def,
    dist_info = dist_info,
    params = params
  ))

  validate_cdist_item(item)

  item
}

#' Create custom probability distributions from discrete densities
#'
#' @param .data A data frame from which to source values and densities
#' @param values,densities Values and relative probability densities
#'
#' @importFrom rlang !!
#' @importFrom rlang enquo
#' @export
custom_cdist_item <- function(.data = NULL, values, densities) {
  data <- data_eval(.data, values = !!enquo(values), densities = !!enquo(densities))
  cdist_item(
    density_function = function(values) {
      stats::approx(data$values, data$densities, xout = values)$y
    },
    quantile_function = function(q) {
      cdf <- cumsum(data$densities)
      stats::approx(cdf / max(cdf), values, xout = q)$y
    }
  )
}


#' Summarise, print, characterify distributions
#'
#' @param object,x A \link{cdist_item} object.
#' @param digits The number of digits to display
#' @param alpha The confidence level to which print when coercing to character
#' @param ... Passed to parent methods
#'
#' @export
#'
summary.cdist_item <- function(object, ...) {
  quantile_vals <- c(0.01, 0.05, 0.16, 0.5, 0.84, 0.95, 0.99)
  quants <- purrr::set_names(
    stats::quantile(object, quantile_vals),
    paste0("quantile_", quantile_vals * 100)
  )
  tibble::as.tibble(as.list(quants))
}

#' @rdname summary.cdist_item
#' @export
as.character.cdist_item <- function(x, digits = 3, alpha = 0.05, ...) {
  conf_level <- 1 - alpha
  quants <- format(quantile(x, c(alpha, 0.5, conf_level)), digits = digits, ...)
  sprintf("%s (%s%% CI: %s-%s)", quants[2], conf_level * 100, quants[1], quants[3])
}

#' @rdname summary.cdist_item
#' @export
print.cdist_item <- function(x, ...) {
  info <- x$dist_info
  if(length(info) > 0) {
    desc <- paste(names(info), purrr::map_chr(info, format), sep = " = ", collapse = ", ")
    desc <- paste0(" ", desc)
  } else {
    desc <- ""
  }
  chr <- format.cdist_item(x, ...)
  cat(sprintf("<continuous distribution%s>\n", desc))
  print(chr, quote = FALSE)
  invisible(x)
}

#' @rdname summary.cdist_item
#' @export
format.cdist_item <- function(x, ...) {
  as.character(x, ...)
}

new_cdist_item <- function(x) {
  if(!is.list(x)) stop("class cpdist must be a list")
  structure(x, class = "cdist_item")
}

validate_cdist_item <- function(x) {
  if(!inherits(x, "cdist_item")) stop("x is not a of class cdist_item")
  if(!is.function(x$density_function)) stop("x$density_function is not a function")
  if(!is.function(x$quantile_function)) stop("x$quantile_function is not a function")
  if(!is.list(x$params)) stop("x$params is not a list")

  for(func in x[c("density_function", "quantile_function")]) {
    if(!all(names(x$params) %in% names(formals(func)))) stop("function missing arguments for params")
  }

  invisible(x)
}

#' @importFrom stats quantile
#' @export
quantile.cdist_item <- function(x, probs = seq(0, 1, 0.25), names = FALSE, ...) {
  q <- do.call(x$quantile_function, c(list(probs), x$params))
  if(names) {
    purrr::set_names(q, format(probs))
  } else {
    q
  }
}

#' @importFrom stats density
#' @export
density.cdist_item <- function(x, values, ...) {
  do.call(x$density_function, c(list(values), x$params))
}

#' Translate a distribution from y to x along a discrete function
#'
#' @param .data A data frame from which to source x, y, and y_sd
#' @param x A vector of destination values
#' @param y A vector of source values (x->y must be unique)
#' @param y_sd Optional standard deviation of source value curve
#' @param dist A \link{cdist_item} containing source densities
#' @param eps Initial densities smaller than this number will not be considered
#' @param n The number of equally-spaced x values at which density should be calculated
#'
#' @importFrom rlang !!
#' @importFrom rlang enquo
#' @export
translate_distribution <- function(.data = NULL, x, y, y_sd = 0, dist, eps = 1e-8, n = 512) {
  data <- data_eval(.data, x = !!enquo(x), y = !!enquo(y), y_sd = !!enquo(y_sd))
  data$density_y <- density(dist, data$y)

  nominal_sd <- diff(quantile(dist, c(0.16, 0.84))) / 2
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

  custom_cdist_item(values = dist_values$x, densities = dist_values$density)
}
