

#' Construct individual continuous distribution items
#'
#' @param params Parameters of the distribution (e.g. mean, sd)
#' @param dist The distribution (e.g., norm)
#' @param values Values for custom density distributions
#' @param densities Paired density values for values
#'
#' @return A continuous distribution object
#' @export
#'
dist_item_parameterized <- function(dist, params) {
  if(!is.character(dist) || (length(dist) != 1)) stop("dist must be a character vector of length 1")
  if(!is.list(params)) stop("params must be a list")
  structure(
    list(params = params, dist = dist),
    class = c(paste("dist_item", dist, sep = "_"), "dist_item_parameterized", "cdist_item")
  )
}

#' @rdname dist_item_parameterized
#' @export
dist_item_custom <- function(values, densities) {
  data <- tibble::tibble(values = !!enquo(values), densities = !!enquo(densities))
  if(!is.numeric(data$values) || !is.numeric(data$densities)) {
    stop("values and densities must be numeric")
  }
  structure(
    list(data = data, dist = "custom"),
    class = c("dist_item_custom", "cdist_item")
  )
}



#' Quantile functions for continuous distribution items
#'
#' @param x The distribution
#' @param probs The probabilities
#' @param names Whether or not to include names in output
#' @param ... Unused
#'
#' @return A numeric vector (possibly named) of quantiles
#' @export
#' @importFrom purrr %||%
#' @importFrom stats quantile
#'
quantile.dist_item_parameterized <- function(x, probs = seq(0, 1, by = 0.25), names = FALSE, ...) {
  stats_ns <- getNamespace("stats")
  this_ns <- getNamespace("carbon14")
  quantile_function <- this_ns[[paste0("q", x$dist)]] %||% stats_ns[[paste0("q", x$dist)]]
  if(is.null(quantile_function)) {
    stop("could not resolve quantile function q", x$dist, "()")
  }
  q <- do.call(quantile_function, c(list(probs), x$params))
  if(names) {
    purrr::set_names(q, format(probs))
  } else {
    q
  }
}

#' @rdname quantile.dist_item_parameterized
#' @export
quantile.dist_item_custom <- function(x, probs = seq(0, 1, by = 0.25), names = FALSE, ...) {
  cdf <- cumsum(x$data$densities)
  cdf <- cdf / max(cdf)
  vals <- stats::approx(cdf, x$data$values, xout = probs)$y
  values <- dplyr::if_else(probs < min(cdf), min(x$data$values), vals)
  if(names) {
    purrr::set_names(values, format(probs))
  } else {
    values
  }
}

#' Density functions for continuous distribution items
#'
#' @param x A continuous distribution item
#' @param values Values at which densities should be calculated
#' @param ... Unused
#'
#' @return A numeric vector of densities
#' @export
#' @importFrom purrr %||%
#' @importFrom stats density
#'
density.dist_item_parameterized <- function(x, values, ...) {
  stats_ns <- getNamespace("stats")
  this_ns <- getNamespace("carbon14")
  density_function <- this_ns[[paste0("d", x$dist)]] %||% stats_ns[[paste0("d", x$dist)]]
  if(is.null(density_function)) {
    stop("could not resolve density function d", x$dist, "()")
  }
  do.call(density_function, c(list(values), x$params))
}

#' @rdname density.dist_item_parameterized
#' @export
density.dist_item_custom <- function(x, values, ...) {
  dens <- stats::approx(x$data$values, x$data$densities, xout = values)$y
  dplyr::if_else(is.na(dens), 0, dens)
}


#' Alternative t distribution
#'
#' T distribution but with mean and standard deviation (same params as the MASS implementation)
#'
#' @param values Values at which densities should be calculated
#' @param df Degrees of freedom
#' @param m "location" (mean)
#' @param s "scale" (sd)
#'
#' @noRd
#'
dt <- function(values, df, m = 0, s = 1) {
  stats::dt((values - m) / s, df = df)
}

qt <- function(q, df, m = 0, s = 1) {
  stats::qt(q, df) * s + m
}

#' Summarise, print, characterify distributions
#'
#' @param object,x A continuous distribution object.
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

  info <- c(list(dist = object$dist), object$params)
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
  info <- c(list(dist = x$dist), x$params)
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
  quantiles <- try(quantile(x, c(0.05, 0.95)), silent = TRUE)
  if(inherits(quantiles, "try-error")) stop("quantile() is not properly defined for x: ",
                                            as.character(quantiles))
  densities <- try(density(x, quantiles), silent = TRUE)
  if(inherits(densities, "try-error")) stop("density() is not properly defined for x: ",
                                            as.character(densities))

  invisible(x)
}

#' @importFrom stats density
#' @export
range.cdist_item <- function(..., na.rm = FALSE, eps = 0) {
  items <- list(...)
  low <- purrr::map_dbl(items, quantile, eps)
  high <- purrr::map_dbl(items, quantile, 1 - eps)
  range(c(low, high))
}

#' @export
min.cdist_item <- function(..., na.rm = FALSE, eps = 0) {
  range(..., na.rm = na.rm, eps = eps)[1]
}

#' @export
max.cdist_item <- function(..., na.rm = FALSE, eps = 0) {
  range(..., na.rm = na.rm, eps = eps)[2]
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

#' @importFrom stats weighted.mean
#' @export
weighted.mean.dist_item_norm <- function(x, w = NULL, ...) {
  x$params$mean
}

#' @importFrom stats weighted.mean
#' @export
weighted.mean.dist_item_t <- function(x, w = NULL, ...) {
  x$params$m
}

#' @export
mean.cdist_item <- function(x, trim = 0, ...) {
  weighted.mean(x, eps = trim)
}

#' @export
mean.dist_item_norm <- function(x, ...) {
  x$params$mean
}

#' @export
mean.dist_item_t <- function(x, ...) {
  x$params$m
}

#' @export
is.na.cdist_item <- function(x, ...) {
  all(is.na(quantile(c(0.05, 0.5, 0.95))))
}

#' @export
`+.dist_item_norm` <- function(x, y) {
  if(missing(y)) return(x)
  if(inherits(x, "cdist_item")) {
    x$params$mean <- x$params$mean + y
    x
  } else {
    y$params$mean <- y$params$mean + x
    y
  }
}

#' @export
`*.dist_item_norm` <- function(x, y) {
  if(inherits(x, "cdist_item")) {
    dist <- x
    operand <- y
  } else {
    dist <- y
    operand <- x
  }
  dist$params$mean <- dist$params$mean * operand
  dist$params$sd <- dist$params$sd * abs(operand)
  dist
}

#' @export
`+.dist_item_t` <- function(x, y) {
  if(missing(y)) return(x)
  if(inherits(x, "cdist_item")) {
    x$params$m <- x$params$m + y
    x
  } else {
    y$params$m <- y$params$m + x
    y
  }

}

#' @export
`*.dist_item_t` <- function(x, y) {
  if(inherits(x, "cdist_item")) {
    dist <- x
    operand <- y
  } else {
    dist <- y
    operand <- x
  }
  dist$params$m <- dist$params$m * operand
  dist$params$s <- dist$params$s * abs(operand)
  dist
}

#' @export
`+.dist_item_custom` <- function(x, y) {
  if(missing(y)) return(x)
  if(inherits(x, "cdist_item")) {
    x$data$values <- x$data$values + y
    x
  } else {
    y$data$values <- y$data$values + x
    y
  }

}

#' @export
`*.dist_item_custom` <- function(x, y) {
  if(inherits(x, "cdist_item")) {
    dist <- x
    operand <- y
  } else {
    dist <- y
    operand <- x
  }
  dist$data$values <- dist$data$values * operand
  dist
}

#' @export
`-.cdist_item` <- function(x, y) {
  if(missing(y)) return(x * -1)
  x + (y * -1)
}

#' @export
`/.cdist_item` <- function(x, y) {
  if(inherits(x, "cdist_item")) {
    x * (1 / y)
  } else {
    stop("division by a distribution is not defined")
  }
}

#' Construct parameterized distribution vectors
#'
#' @param .data An optional data frame from which parameters should be sourced
#' @param mean,sd Parameters of the normal distribution
#' @param m,s,df Parameters of the t distribution (location of m, scale of s)
#' @param names Names for output objects
#'
#' @export
dist_norm <- function(.data = NULL, mean, sd, names = NULL) {
  data <- data_eval(.data, mean = !!enquo(mean), sd = !!enquo(sd), names = !!enquo(names))
  cd <- new_cdist(purrr::map(purrr::transpose(data[c("mean", "sd")]),
                              dist_item_parameterized, dist = "norm"))
  if("names" %in% colnames(data)) {
    purrr::set_names(cd, data$names)
  } else {
    cd
  }
}

#' @rdname dist_norm
#' @export
dist_t <- function(.data = NULL, m, s, df, names = NULL) {
  data <- data_eval(.data, df = !!enquo(df), m = !!enquo(m), s = !!enquo(s),
                    names = !!enquo(names))
  cd <- new_cdist(purrr::map(purrr::transpose(data[c("m", "s", "df")]),
                              dist_item_parameterized, dist = "t"))
  if("names" %in% colnames(data)) {
    purrr::set_names(cd, data$names)
  } else {
    cd
  }
}

#' Construct a continuous distribution vector
#'
#' @param ... Items created with \link{dist_item_parameterized} or \link{dist_item_custom}.
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

#' @export
quantile.cdist <- function(x, probs = seq(0, 1, by = 0.25), names = FALSE, ...) {
  sum <- dplyr::select(summary(x, quantiles = probs), dplyr::starts_with("quantile"))
  mat <- as.matrix(sum)
  if(names) {
    rownames(mat) <- names(x)
    colnames(mat) <- format(probs)
  } else {
    rownames(mat) <- NULL
    colnames(mat) <- NULL
  }

  mat
}

#' @export
mean.cdist <- function(x, na.rm = FALSE, ...) {
  purrr::map_dbl(x, mean, na.rm = na.rm, ...)
}

#' @export
weighted.mean.cdist <- function(x, w, ...) {
  purrr::map_dbl(x, weighted.mean, ...)
}

#' @export
#' @importFrom stats median
median.cdist <- function(x, ...) {
  quantile(x, 0.5, ...)
}

#' @export
is.na.cdist <- function(x, ...) {
  purrr::map_lgl(x, is.na, ...)
}

#' @export
`+.cdist` <- function(x, y) {
  if(missing(y)) return(x)
  if(inherits(x, "cdist")) {
    new_cdist(lapply(x, "+", y = y))
  } else {
    new_cdist(lapply(y, "+", x = x))
  }
}

#' @export
`-.cdist` <- function(x, y) {
  if(missing(y)) return(x * -1)
  x + (y * -1)
}

#' @export
`*.cdist` <- function(x, y) {
  if(inherits(x, "cdist")) {
    new_cdist(lapply(x, "*", y = y))
  } else {
    new_cdist(lapply(y, "*", x = x))
  }
}

#' @export
`/.cdist` <- function(x, y) {
  if(inherits(x, "cdist")) {
    new_cdist(lapply(x, "/", y = y))
  } else {
    stop("division by a distribution is not defined")
  }
}

#' Boxplot a continuous distribution vector using base graphics
#'
#' @param x A continuous distribution vector
#' @param whisker_quantiles Quantiles for the ends of the whiskers
#' @param box_quantiles Quantiles for the ends of the box
#' @param mid_quantile Quantile for the middle marker
#' @param ... Passed to \link[graphics]{bxp}
#'
#' @export
#' @importFrom graphics boxplot
#'
boxplot.cdist <- function(x, whisker_quantiles = c(0.01, 0.99), box_quantiles = c(0.25, 0.75),
                          mid_quantile = 0.5, ...) {

  stats <- sapply(
    x, quantile,
    c(whisker_quantiles[1], box_quantiles[1], mid_quantile, box_quantiles[2], whisker_quantiles[2])
  )

  boxplot_info <- list(
    stats = stats,
    n = rep(Inf, length(x)),
    conf = matrix(NA_real_, nrow = 2, ncol = length(x)),
    out = numeric(0),
    names = names(x) %||% as.character(seq_len(length(x)))
  )

  graphics::bxp(boxplot_info, ...)
}
