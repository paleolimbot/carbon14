context("test-calibrate.R")

test_that("calibrate works with curves specified as curves, vectors, character vectors, or NULL", {
  expect_identical(
    summary(calibrate(measured_age = 330, measured_age_error = 30,
                      curve = intcal04)$age_bp_distribution),
    summary(calibrate(measured_age = 330, measured_age_error = 30,
                      curve = "intcal04")$age_bp_distribution)
  )
  expect_identical(
    summary(calibrate(measured_age = 330, measured_age_error = 30,
                      curve = intcal04)$age_bp_distribution),
    summary(calibrate(measured_age = 330, measured_age_error = 30,
                      curve = list(intcal04))$age_bp_distribution)
  )

  # should be close to input values
  cal_null <- calibrate(measured_age = 50, measured_age_error = 10, curve = NULL)
  expect_equal(quantile(cal_null$age_bp_distribution[[1]], c(0.05, 0.95)),
                                qnorm(c(0.05, 0.95), mean = 50, sd = 10))

  # uncomment this when arithmetic on distributions is implemented
  # cal_null_ad <- calibrate(
  #   measured_age = 50,
  #   measured_age_error = 10,
  #   curve = null_calibration_curve(),
  #   cal_age_type = "Year AD"
  # )
  # expect_equal(quantile(cal_null_ad$age_bp_distribution[[1]], c(0.05, 0.95)),
  #              qnorm(c(0.05, 0.95), mean = 1900, sd = 10))

  expect_error(calibrate(measured_age = 330, measured_age_error = 30, curve = environment()),
               "curve must be a data frame, a vector, or NULL")
})




test_that("translate pdist works according to plan", {
  cd <- cdist_item(mean = 10, sd = 1, dist = "norm")
  curve <- calibration_curve(measured_age = seq(-10, 30, length.out = 10),
                             cal_age = seq(-10, 30, length.out = 10))
  dta <- translate_distribution(curve, dist = cd)

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
  date <- cdist_item(m = 340, s = 30, df = 100, dist = "t")
  cal <- translate_distribution(intcal13, dist = date)
  expect_is(cal, "cdist_item")
  bchron_cal <- Bchron::BchronCalibrate(340, 30, "intcal13")
  bchron_dist <- cdist_item_from_densities(
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


test_that("vectorized version of translate_distribution works", {

  date <- cdist_item(mean = 340, sd = 30, dist = "norm")
  date_vec <- cdist(date)
  cal <- translate_distribution(intcal13, dist = date)
  cal_vec <- translate_distribution(intcal13, dist = date_vec)

  expect_is(cal_vec, "cdist")

  expect_identical(
    summary(cal),
    summary(cal_vec[[1]])
  )
})

