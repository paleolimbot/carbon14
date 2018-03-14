context("test-pdist.R")

test_that("dist_item works for parameterized distributions", {
  cd <- dist_item_parameterized("norm", list(mean = 10, sd = 1))
  expect_equal(quantile(cd, 0.5), 10)
  expect_equal(density(cd, 5:15), dnorm(5:15, mean = 10, sd = 1))
})

test_that("custom specified distributions work", {
  vals <- seq(-5, 5, 0.01)
  densities <- dnorm(vals, mean = 0, sd = 1)
  cd <- dist_item_custom(values = vals, densities = densities)
  expect_true(abs(quantile(cd, 0.5)) <= 0.01)
  expect_true(
    all(
      na.omit(density(cd, vals + 0.005) - dnorm(vals + 0.005, mean = 0, sd = 1)) %>%
        abs() %>%
        dplyr::near(0, tol = 1e-4)
    )
  )
  expect_is(summary(cd), "tbl_df")
})

test_that("distribution math is close enough", {
  truth <- dist_item_parameterized("norm", list(mean = 0, sd = 1))
  approx <- dist_item_custom(
    values = seq(-5, 5, 0.01),
    densities = dnorm(values, mean = 0, sd = 1)
  )

  expect_identical(range(truth), c(-Inf, Inf))
  expect_identical(range(approx), c(-5, 5))
  expect_true(dplyr::near(weighted.mean(truth), weighted.mean(approx), tol = 1e-4))
})

test_that("summary, as_tibble, as.data.frame functions returns a 1-row data frame", {
  cd <- dist_item_parameterized("norm", list(mean = 10, sd = 1))
  expect_is(summary(cd), "tbl_df")
  expect_equal(nrow(summary(cd)), 1)
  expect_is(tibble::as_tibble(cd), "tbl_df")
  expect_is(as.data.frame(cd), "data.frame")
})

test_that("as.character returns a character vector", {
  cd <- dist_item_parameterized("norm", list(mean = 10, sd = 1))
  expect_is(as.character(cd), "character")
  expect_length(as.character(cd), 1)
})

test_that("print prints things", {
  cd <- dist_item_parameterized("norm", list(mean = 10, sd = 1))
  expect_output(expect_is(print(cd), "cdist_item"), "^<continuous distribution")
})

test_that("vectorized normal distributions work according to plan", {
  cdv <- cdist(
    dist_item_parameterized("norm", list(mean = 10, sd = 1)),
    dist_item_parameterized("norm", list(mean = 0, sd = 1))
  )

  expect_is(cdv, "cdist")
  expect_is(head(cdv), "cdist")
  expect_is(tail(cdv), "cdist")
  expect_output(expect_is(print(cdv), "cdist"), "^<continuous distribution vector")

  cdv2 <- cdist(
    dist_item_parameterized("norm", list(mean = 5, sd = 1)),
    dist_item_parameterized("norm", list(mean = 15, sd = 1))
  )

  expect_is(c(cdv, cdv2), "cdist")
  expect_length(c(cdv, cdv2), 4)

  cdv[1] <- list(dist_item_parameterized("norm", list(mean = 22, sd = 1)))
  expect_is(cdv, "cdist")
  expect_equal(cdv[[1]]$params$mean, 22)

  cdv[[1]] <- dist_item_parameterized("norm", list(mean = 12, sd = 1))
  expect_is(cdv, "cdist")
  expect_equal(cdv[[1]]$params$mean, 12)

  expect_is(summary(cdv), "tbl_df")
  expect_is(tibble::as_tibble(cdv), "tbl_df")
  expect_is(as.data.frame(cdv), "data.frame")
})

test_that("vectorized distributions look ok in print output", {
  cdv <- cdist(
    name1 = dist_item_parameterized("norm", list(mean = 10, sd = 1)),
    name2 = dist_item_parameterized("norm", list(mean = 0, sd = 1))
  )

  expect_is(format(cdv), "character")
  # check print output of df and tibble
  df <- data.frame(name = names(cdv), I(cdv))
  tbl <- tibble::tibble(name = names(cdv), cdv)

})

test_that("normal, t shortcuts work as expected", {
  ndist <- dist_norm(mean = c(0, 5, 10), sd = 1)
  expect_length(ndist, 3)
  expect_is(ndist, "cdist")
  expect_true(all(dplyr::near(summary(ndist)$weighted_mean, c(0, 5, 10))))

  ndist_named <- dist_norm(mean = c(0, 5, 10), sd = 1, names = c("zero", "five", "ten"))
  expect_identical(names(ndist_named), c("zero", "five", "ten"))

  tdist <- dist_t(m = c(0, 5, 10), s = 1, df = Inf)
  expect_identical(
    summary(tdist)$weighted_mean,
    summary(ndist)$weighted_mean
  )

  tdist_named <- dist_t(m = c(0, 5, 10), s = 1, df = Inf, names = c("zero", "five", "ten"))
  expect_identical(names(tdist_named), c("zero", "five", "ten"))
})

test_that("boxplotting of cdist vectors works", {
  # visual test
  expect_true(TRUE)

  tdist_named <- dist_t(m = c(0, 5, 10), s = c(1, 2, 3), df = 10,
                        names = c("zero", "five", "ten"))
  boxplot(tdist_named)

  lognorm_sample <- rlnorm(1000)
  lognorm_dens <- density(lognorm_sample)
  lognorm_custom <- dist_item_custom(values = lognorm_dens$x, densities = lognorm_dens$y)
  boxplot(lognorm_sample)
  boxplot(cdist(lognorm_custom))
})

test_that("distribution math checks out", {
  norm <- dist_norm(mean = 10, sd = 1)
  t <- dist_t(m = 10, s = 1, df = 12)
  cust <- cdist(dist_item_custom(
    seq(-10, 30, by = 0.01),
    dnorm(seq(-10, 30, by = 0.01), mean = 10, sd = 1)
  ))

  # normal distribution arithmetic
  expect_identical(mean(norm + 5), 10 + 5); expect_identical((norm + 5)[[1]]$params$sd, 1)
  expect_identical(mean(norm - 5), 10 - 5); expect_identical((norm - 5)[[1]]$params$sd, 1)
  expect_identical(mean(norm * 5), 10 * 5); expect_identical((norm * 5)[[1]]$params$sd, 1 * 5)
  expect_identical(mean(norm / 5), 10 / 5); expect_identical((norm / 5)[[1]]$params$sd, 1 / 5)

  # t distribution arithmetic
  expect_identical(mean(t + 5), 10 + 5); expect_identical((t + 5)[[1]]$params$s, 1)
  expect_identical(mean(t - 5), 10 - 5); expect_identical((t - 5)[[1]]$params$s, 1)
  expect_identical(mean(t * 5), 10 * 5); expect_identical((t * 5)[[1]]$params$s, 1 * 5)
  expect_identical(mean(t / 5), 10 / 5); expect_identical((t / 5)[[1]]$params$s, 1 / 5)

  # custom distribution arithmetic
  tol <- 1e-7
  expect_true(dplyr::near(mean(cust + 5), 10 + 5, tol = tol));
  expect_true(dplyr::near(mean(cust - 5), 10 - 5, tol = tol));
  expect_true(dplyr::near(mean(cust * 5), 10 * 5, tol = tol));
  expect_true(dplyr::near(mean(cust / 5), 10 / 5, tol = tol));

  # normal distribution arithmetic properties
  expect_identical(+norm, norm)
  expect_identical(-norm, norm * -1)
  expect_identical(norm + 5, 5 + norm)
  expect_identical(5 - norm, -1 * (norm - 5))
  expect_identical(norm * 5, 5 * norm)
  expect_error(5 / norm, "division by a distribution is not defined")

  # t distribution arithmetic properties
  expect_identical(+t, t)
  expect_identical(-t, t * -1)
  expect_identical(t + 5, 5 + t)
  expect_identical(5 - t, -1 * (t - 5))
  expect_identical(t * 5, 5 * t)
  expect_error(5 / t, "division by a distribution is not defined")

  # custom distribution arithmetic properties
  expect_identical(+cust, cust)
  expect_identical(-cust, cust * -1)
  expect_identical(cust + 5, 5 + cust)
  expect_identical(5 - cust, -1 * (cust - 5))
  expect_identical(cust * 5, 5 * cust)
  expect_error(5 / cust, "division by a distribution is not defined")
})
