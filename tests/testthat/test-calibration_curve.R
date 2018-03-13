context("calibration_curve")

test_that("built-in calibration curves are calibration curves", {
  expect_is(intcal04, "age_calibration_curve")
  expect_is(intcal09, "age_calibration_curve")
  expect_is(intcal13, "age_calibration_curve")
  expect_is(marine04, "age_calibration_curve")
  expect_is(marine09, "age_calibration_curve")
  expect_is(marine13, "age_calibration_curve")
  expect_is(shcal04, "age_calibration_curve")
  expect_is(shcal13, "age_calibration_curve")
})

test_that("custom calibration curves work according to plan", {
  cc1 <- calibration_curve(cal_age = 0:50, measured_age = 1950:1900)
  expect_is(cc1, "age_calibration_curve")
  expect_output(expect_identical(print(cc1), cc1), "^<age_calibration_curve")
  # visual test
  plot(cc1)
})

test_that("printing of named calibration curves works", {
  expect_output(
    print(intcal13),
    "'intcal13'"
  )
})

test_that("naming a calibration curve works", {
  cc2 <- calibration_curve(cal_age = 0:50, measured_age = 1950:1900, name = "falafel")
  expect_output(print(cc2), "'falafel'")
  expect_error(calibration_curve(cal_age = 0:50, measured_age = 1950:1900, name = list("falafel")),
               "name must be a character vector")
})

test_that("switching curves to different cal age works", {
  expect_identical(intcal13, in_cal_bp(intcal13))
  expect_identical(in_cal_bp(in_year_ad(intcal13)), intcal13)
  expect_identical(in_year_ad(in_year_ad(intcal13)), in_year_ad(intcal13))
  expect_equal(min(intcal13$cal_bp), 0)
  expect_equal(max(in_year_ad(intcal13)$cal_bp), 1950)
})

test_that("switching numeric vectors to different cal works", {
  expect_equal(in_year_ad(c(0, 50)), c(1950, 1900))
  expect_equal(in_cal_bp(c(1950, 1900)), c(0, 50))
})

test_that("plotting of a calibration curve works", {
  expect_true(TRUE) # visual tests
  cc1 <- calibration_curve(cal_age = 0:50, measured_age = 1950:1900)
  cc2 <- calibration_curve(cal_age = 0:50, measured_age = 1950:1900, name = "falafel")
  plot(cc1)
  plot(cc1, title = "fish")
  plot(intcal13)
  plot(intcal13, measured_age_limits = FALSE)
})

test_that("translate pdist works with a calibration curve", {
  date <- cdist_item(m = 340, s = 30, df = 100, dist = "t")
  cal <- translate_distribution(x = intcal13$cal_bp, y = intcal13$age_14C, dist = date)
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
  cal <- translate_distribution(x = intcal13$cal_bp, y = intcal13$age_14C, dist = date)
  cal_vec <- translate_distribution(x = intcal13$cal_bp, y = intcal13$age_14C, dist = date_vec)

  expect_is(cal_vec, "cdist")

  expect_identical(
    summary(cal),
    summary(cal_vec[[1]])
  )
})
