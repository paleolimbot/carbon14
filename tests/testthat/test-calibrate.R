context("test-calibrate.R")

test_that("calibrate works with curves specified as curves, vectors, or character vectors", {
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

  expect_error(calibrate(measured_age = 330, measured_age_error = 30, curve = NULL),
               "curve must be a data frame or a vector")
})
