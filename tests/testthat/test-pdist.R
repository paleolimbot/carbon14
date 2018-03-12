context("test-pdist.R")

test_that("cdist works for parameterized distributions", {
  cd <- cdist(mean =  10, sd = 1, dist = "norm")
  expect_equal(quantile(cd, 0.5), 10)
  expect_equal(density(cd, 5:15), dnorm(5:15, mean = 10, sd = 1))
})

test_that("invalid args for cdist for parameterized distributions are caught", {
  expect_error(cdist(not_a_param_of_norm = 4, dist = "norm"),
               "function missing arguments for params")
})

test_that("custom specified distributions work", {
  vals <- seq(-5, 5, 0.01)
  densities <- dnorm(vals, mean = 0, sd = 1)
  cd <- custom_cdist(values = vals, densities = densities)
  expect_true(abs(quantile(cd, 0.5)) <= 0.01)
  expect_true(
    all(
      na.omit(density(cd, vals + 0.005) - dnorm(vals + 0.005, mean = 0, sd = 1)) %>%
        abs() %>%
        dplyr::near(0, tol = 1e-4)
    )
  )
})

test_that("summary function returns a 1-row tibble of quantiles", {
  cd <- cdist(mean = 10, sd = 1, dist = "norm")
  expect_is(summary(cd), "tbl_df")
  expect_true(all(grepl("^quantile", colnames(summary(cd)))))
  expect_equal(nrow(summary(cd)), 1)
})

test_that("as.character returns a character vector", {
  cd <- cdist(mean = 10, sd = 1, dist = "norm")
  expect_is(as.character(cd), "character")
  expect_length(as.character(cd), 1)
})

test_that("print prints things", {
  cd <- cdist(mean = 10, sd = 1, dist = "norm")
  expect_output(expect_is(print(cd), "cdist"), "^<continuous distribution")
})

test_that("translate pdist works according to plan", {
  cd <- cdist(mean = 10, sd = 1, dist = "norm")
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

test_that("translate pdist works with a calibration curve", {
  date <- cdist(mean = 340, sd = 30, dist = "norm")
  cal <- translate_distribution(x = intcal13$cal_bp, y = intcal13$age_14C, dist = date)
  bchron_cal <- Bchron::BchronCalibrate(340, 30, "intcal13")
  bchron_dist <- custom_cdist(
    values = bchron_cal$Date1$ageGrid,
    densities = bchron_cal$Date1$densities
  )

  # quantiles are all within 2 yrs of estimate from bchron values
  expect_true(
    all(
      abs(
        quantile(cal, seq(0.05, 0.95, 0.1)) -
          quantile(bchron_dist, seq(0.05, 0.95, 0.1))
      ) < 2
    )
  )

  # visual comparison
  range_cal <- quantile(cal, c(0.00001, 0.99999))

  plot(intcal13, ylim = c(225, 450), xlim = range_cal)
  plot_range <- par('usr')
  test_bp <- seq(range_cal[1], range_cal[2], length.out = 512)
  test_dens <- density(cal, test_bp)
  test_dens <- test_dens / max(test_dens)
  lines(test_bp, test_dens * (plot_range[4] - plot_range[3]) * 0.25 + plot_range[3], col = "red")

  range_14c <- quantile(date, c(0.001, 0.999))
  test_14c <- seq(range_14c[1], range_14c[2], length.out = 512)
  test_dens <- density(date, test_14c)
  test_dens <- test_dens / max(test_dens)
  lines(test_dens * (plot_range[2] - plot_range[1]) * 0.25 + plot_range[1], test_14c,  col = "blue")

  bchron_age <- bchron_cal$Date1$ageGrid
  bchron_dens <- bchron_cal$Date1$densities
  bchron_dens <- bchron_dens / max(bchron_dens)
  lines(bchron_age, bchron_dens * (plot_range[4] - plot_range[3]) * 0.25 + plot_range[3],
        col = "purple")
})
