context("test-pdist.R")

test_that("cdist_item works for parameterized distributions", {
  cd <- cdist_item(mean =  10, sd = 1, dist = "norm")
  expect_equal(quantile(cd, 0.5), 10)
  expect_equal(density(cd, 5:15), dnorm(5:15, mean = 10, sd = 1))
})

test_that("invalid args for cdist_item for parameterized distributions are caught", {
  expect_error(cdist_item(not_a_param_of_norm = 4, dist = "norm"),
               "function missing arguments for params")
})

test_that("custom specified distributions work", {
  vals <- seq(-5, 5, 0.01)
  densities <- dnorm(vals, mean = 0, sd = 1)
  cd <- cdist_item_from_densities(values = vals, densities = densities)
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
  truth <- cdist_item(mean = 0, sd = 1, dist = "norm")
  approx <- cdist_item_from_densities(
    values = seq(-5, 5, 0.01),
    densities = dnorm(values, mean = 0, sd = 1)
  )

  expect_identical(range(truth), c(-Inf, Inf))
  expect_identical(range(approx), c(-5, 5))
  expect_true(dplyr::near(weighted.mean(truth), weighted.mean(approx), tol = 1e-4))
})

test_that("summary, as_tibble, as.data.frame functions returns a 1-row data frame", {
  cd <- cdist_item(mean = 10, sd = 1, dist = "norm")
  expect_is(summary(cd), "tbl_df")
  expect_equal(nrow(summary(cd)), 1)
  expect_is(tibble::as_tibble(cd), "tbl_df")
  expect_is(as.data.frame(cd), "data.frame")
})

test_that("as.character returns a character vector", {
  cd <- cdist_item(mean = 10, sd = 1, dist = "norm")
  expect_is(as.character(cd), "character")
  expect_length(as.character(cd), 1)
})

test_that("print prints things", {
  cd <- cdist_item(mean = 10, sd = 1, dist = "norm")
  expect_output(expect_is(print(cd), "cdist_item"), "^<continuous distribution")
})

test_that("translate pdist works according to plan", {
  cd <- cdist_item(mean = 10, sd = 1, dist = "norm")
  dta <- translate_distribution(
    x = seq(-10, 30, length.out = 10),
    y = seq(-10, 30, length.out = 10),
    dist = cd
  )

  test_quantiles <- c(0.05, 0.5, 0.95)
  expect_true(
    all(
      abs(quantile(dta, test_quantiles) -
            qnorm(test_quantiles, mean = 10, sd = 1)) <
        0.5
    )
  )
})

test_that("vectorized normal distributions work according to plan", {
  cdv <- cdist(
    cdist_item(mean = 10, sd = 1, dist = "norm"),
    cdist_item(mean = 0, sd = 1, dist = "norm")
  )

  expect_is(cdv, "cdist")
  expect_is(head(cdv), "cdist")
  expect_is(tail(cdv), "cdist")
  expect_output(expect_is(print(cdv), "cdist"), "^<continuous distribution vector")

  cdv2 <- cdist(
    cdist_item(mean = 5, sd = 1, dist = "norm"),
    cdist_item(mean = 15, sd = 1, dist = "norm")
  )

  expect_is(c(cdv, cdv2), "cdist")
  expect_length(c(cdv, cdv2), 4)

  cdv[1] <- list(cdist_item(mean = 22, sd = 1, dist = "norm"))
  expect_is(cdv, "cdist")
  expect_equal(cdv[[1]]$dist_info$mean, 22)

  cdv[[1]] <- cdist_item(mean = 12, sd = 1, dist = "norm")
  expect_is(cdv, "cdist")
  expect_equal(cdv[[1]]$dist_info$mean, 12)

  expect_is(summary(cdv), "tbl_df")
  expect_is(tibble::as_tibble(cdv), "tbl_df")
  expect_is(as.data.frame(cdv), "data.frame")
})

test_that("vectorized distributions look ok in print output", {
  cdv <- cdist(
    name1 = cdist_item(mean = 10, sd = 1, dist = "norm"),
    name2 = cdist_item(mean = 0, sd = 1, dist = "norm")
  )

  expect_is(format(cdv), "character")
  # check print output of df and tibble
  df <- data.frame(name = names(cdv), I(cdv))
  tbl <- tibble::tibble(name = names(cdv), cdv)

})

test_that("normal, t shortcuts work as expected", {
  ndist <- cdist_norm(mean = c(0, 5, 10), sd = 1)
  expect_length(ndist, 3)
  expect_is(ndist, "cdist")
  expect_true(all(dplyr::near(summary(ndist)$weighted_mean, c(0, 5, 10))))

  ndist_named <- cdist_norm(mean = c(0, 5, 10), sd = 1, names = c("zero", "five", "ten"))
  expect_identical(names(ndist_named), c("zero", "five", "ten"))

  tdist <- cdist_t(m = c(0, 5, 10), s = 1, df = Inf)
  expect_identical(
    summary(tdist)$weighted_mean,
    summary(ndist)$weighted_mean
  )

  tdist_named <- cdist_t(m = c(0, 5, 10), s = 1, df = Inf, names = c("zero", "five", "ten"))
  expect_identical(names(tdist_named), c("zero", "five", "ten"))
})

test_that("boxplotting of cdist vectors works", {
  # visual test
  expect_true(TRUE)

  tdist_named <- cdist_t(m = c(0, 5, 10), s = c(1, 2, 3), df = 10,
                         names = c("zero", "five", "ten"))
  boxplot(tdist_named)

  lognorm_sample <- rlnorm(1000)
  lognorm_dens <- density(lognorm_sample)
  lognorm_custom <- cdist_item_from_densities(values = lognorm_dens$x, densities = lognorm_dens$y)
  boxplot(lognorm_sample)
  boxplot(cdist(lognorm_custom))


})
