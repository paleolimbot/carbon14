
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
cdist_item <- function(..., dist = NULL, density_function = NULL,
                       quantile_function = NULL, dist_info = list()) {
  params <- tibble::lst(...)
  stats_ns <- getNamespace("stats")
  if(!is.null(dist)) {
    density_function_def <- stats_ns[[paste0("d", dist)]]
    quantile_function_def <- stats_ns[[paste0("q", dist)]]
    dist_info <- c(dist_info, list(dist = dist), params)
    dist_info <- dist_info[unique(names(dist_info))]
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
cdist_item_from_densities <- function(.data = NULL, values, densities) {
  data <- data_eval(.data, values = !!enquo(values), densities = !!enquo(densities))
  cdist_item(
    density_function = function(values) {
      density_from_data(data, values)
    },
    quantile_function = function(q) {
      quantile_from_data(data, q)
    },
    dist_info = list(dist = "custom", data = data)
  )
}

density_from_data <- function(data, values) {
  dens <- stats::approx(data$values, data$densities, xout = values)$y
  dplyr::if_else(is.na(dens), 0, dens)
}

quantile_from_data <- function(data, q) {
  cdf <- cumsum(data$densities)
  cdf <- cdf / max(cdf)
  vals <- stats::approx(cdf, data$values, xout = q)$y
  dplyr::if_else(q < min(cdf), min(data$values), vals)
}

#' Summarise, print, characterify distributions
#'
#' @param object,x A \link{cdist_item} object.
#' @param quantiles Quantiles to include in the summary object
#' @param digits The number of digits to display
#' @param alpha The confidence level to which print when coercing to character
#' @param ... Passed to parent methods
#'
#' @export
#'
summary.cdist_item <- function(object, quantiles = c(0.05, 0.25, 0.5, 0.75, 0.95), ...) {
  quants <- purrr::set_names(
    stats::quantile(object, quantiles),
    paste0("quantile_", quantiles * 100)
  )

  info <- object$dist_info
  if(length(info) > 0) {
    is_atomic <- purrr::map_lgl(info, is.atomic)
    len_1 <- purrr::map_lgl(info, function(item) length(item) == 1)
    info <- info[is_atomic & len_1]
  }

  dplyr::bind_cols(
    tibble::tibble(weighted_mean = weighted.mean(object)),
    tibble::as_tibble(as.list(quants)),
    tibble::as_tibble(info)
  )
}

#' @rdname summary.cdist_item
#' @importFrom tibble as_tibble
#' @export
as_tibble.cdist_item <- function(x, ...) {
  summary.cdist_item(x, ...)
}

#' @rdname summary.cdist_item
#' @export
as.data.frame.cdist_item <- function(x, ...) {
  as.data.frame(summary.cdist_item(x, ...))
}

#' @rdname summary.cdist_item
#' @export
as.character.cdist_item <- function(x, digits = 3, alpha = 0.05, ...) {
  conf_level <- 1 - alpha
  quants <- format(quantile(x, c(alpha, 0.5, conf_level)), digits = digits, trim = TRUE, ...)
  sprintf("%s (%s%% CI: %s-%s)", quants[2], conf_level * 100, quants[1], quants[3])
}

#' @rdname summary.cdist_item
#' @export
print.cdist_item <- function(x, ...) {
  info <- x$dist_info
  if(length(info) > 0) {
    is_atomic <- purrr::map_lgl(info, is.atomic)
    len_1 <- purrr::map_lgl(info, function(item) length(item) == 1)
    info <- info[is_atomic & len_1]
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
  format(as.character(x), quote = FALSE, ...)
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
  if(!is.list(x$dist_info)) stop("x$dist_info is not a list")

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

#' @importFrom stats density
#' @export
range.cdist_item <- function(..., na.rm = FALSE, eps = 0) {
  items <- list(...)
  low <- purrr::map_dbl(items, quantile, eps)
  high <- purrr::map_dbl(items, quantile, 1 - eps)
  range(c(low, high))
}

#' @importFrom stats weighted.mean
#' @export
weighted.mean.cdist_item <- function(x, w = NULL, eps = 1e-8, n = 512, ...) {
  if(!is.null(w)) stop("Custom weights not supported")
  rge <- range(x, eps = eps)
  vals <- seq(rge[1], rge[2], length.out = n)
  dens <- density(x, vals)
  weighted.mean(vals, w = dens / sum(dens))
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
  if(inherits(dist, "cdist")) {
    # vectorize
    return(
      new_cdist(purrr::map(
        dist, translate_distribution,
        .data = NULL,
        x = data$x,
        y = data$y,
        y_sd = data$y_sd,
        eps = eps,
        n = n
      ))
    )
  }

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

  cdist_item_from_densities(values = dist_values$x, densities = dist_values$density)
}

#' Construct parameterized distribution vectors
#'
#' @param .data An optional data frame from which parameters should be sourced
#' @param mean,sd Parameters of the normal distribution
#' @param ncp,df Parameters of the t distribution
#' @param names Names for output objects
#'
#' @export
cdist_norm <- function(.data = NULL, mean, sd, names = NULL) {
  data <- data_eval(.data, mean = !!enquo(mean), sd = !!enquo(sd), names = !!enquo(names))
  cd <- new_cdist(purrr::pmap(data[c("mean", "sd")], cdist_item, dist = "norm"))
  if("names" %in% colnames(data)) {
    purrr::set_names(cd, data$names)
  } else {
    cd
  }
}

#' @rdname cdist_norm
#' @export
cdist_t <- function(.data = NULL, df, ncp, names = NULL) {
  data <- data_eval(.data, df = !!enquo(df), ncp = !!enquo(ncp), names = !!enquo(names))
  cd <- new_cdist(purrr::pmap(data[c("df", "ncp")], cdist_item, dist = "t"))
  if("names" %in% colnames(data)) {
    purrr::set_names(cd, data$names)
  } else {
    cd
  }
}

#' Construct a continuous distribution vector
#'
#' @param ... Items created with \link{cdist_item} or \link{cdist_item_from_densities}.
#'
#' @return A list-ish object of type cdist.
#' @export
#'
cdist <- function(...) {
  as_cdist.list(list(...))
}

#' Coerce to a continuous distribution vector
#'
#' @param x An object
#' @param ... Passed to other methods
#'
#' @return A continuous distribution vector
#' @export
#'
as_cdist <- function(x, ...) {
  UseMethod("as_cdist")
}

#' @rdname as_cdist
#' @export
as_cdist.cdist <- function(x, ...) {
  x
}

#' @rdname as_cdist
#' @export
as_cdist.cdist_item <- function(x, ...) {
  as_cdist.list(list(x))
}

#' @rdname as_cdist
#' @export
as_cdist.list <- function(x, ...) {
  cd <- new_cdist(x)
  validate_cdist(cd)
  cd
}

#' Create a continuous distribution vector
#'
#' @param x A list
#'
#' @return A continuous distribution vector
#' @export
#'
new_cdist <- function(x) {
  if(!is.list(x)) stop("x must be a list")
  structure(x, class = "cdist")
}

#' Validate a continuous distribution vector
#'
#' @param x A continuous distribution vector
#'
#' @return The input, invisibly
#' @export
#'
validate_cdist <- function(x) {
  if(!inherits(x, "cdist")) stop("x must be of class cdist")
  item_result <- lapply(x, function(item) try(validate_cdist_item(item), silent = TRUE))
  errors <- purrr::map_lgl(item_result, inherits, "try-error")
  if(any(errors)) {
    stop(
      "items at positions ",
      paste(which(errors), collapse = ", "),
      " are not valid cdist objects"
    )
  }
  invisible(x)
}

#' Summarise, coerce continuous distribution vectors
#'
#' @param object,x A \link{cdist} vector
#' @param alpha level of confience for character representation of objects
#' @param digits Number of digits to display
#' @param quote Whether or not to quote formatted objects
#' @param ... Passed to parent functions
#' @param .id column in which names should be placed, if present
#'
#' @export
#'
summary.cdist <- function(object, ..., .id = ".name") {
  purrr::map_dfr(object, summary, ..., .id = .id)
}

#' @rdname summary.cdist
#' @export
as_tibble.cdist <- function(x, ...) {
  summary.cdist(x, ...)
}

#' @rdname summary.cdist
#' @export
as.data.frame.cdist <- function(x, ...) {
  as.data.frame(summary.cdist(x, ...))
}

#' @rdname summary.cdist
#' @export
as.character.cdist <- function(x, alpha = 0.05, digits = NULL, ...) {
  sum <- summary(x, quantiles = c(alpha, 1 - alpha))
  sum_quant <- dplyr::select(sum, dplyr::starts_with("quantile"))
  sprintf(
    "%s (%s%% CI: %s/%s)",
    format(sum$weighted_mean, trim = TRUE, digits = digits, ...),
    (1 - alpha) * 100,
    format(sum_quant[[1]], trim = TRUE, digits = digits, ...),
    format(sum_quant[[2]], trim = TRUE, digits = digits, ...)
  )
}

#' @rdname summary.cdist
#' @export
format.cdist <- function(x, quote = FALSE, alpha = 0.05, digits = 3, ...) {
  format(as.character(x, alpha = alpha, digits = digits), quote = quote, ...)
}

#' @rdname summary.cdist
#' @export
print.cdist <- function(x, alpha = 0.05, digits = 3, ...) {
  conf <- 1 - alpha
  chr <- as.character(x, alpha = alpha, digits = digits, ...)
  cat("<continuous distribution vector>\n")
  print(chr, quote = FALSE)
  invisible(x)
}

#' @export
`[[.cdist<-` <- function(x, i, value) {
  validate_cdist_item(value)
  new_cdist(NextMethod())
}

#' @export
`[.cdist<-` <- function(x, i, value, ...) {
  if(!is.list())
  lapply(value, validate_cdist_item(value))
  new_cdist(NextMethod())
}

#' @export
`[.cdist` <- function(x, i, ...) {
  new_cdist(NextMethod())
}

#' @export
c.cdist <- function(...) {
  new_cdist(NextMethod())
}
