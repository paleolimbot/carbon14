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
