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

test_that("formatting a calibration curve looks ok in tibble/data frame output", {
  expect_identical(format(intcal13), "<age_calibration_curve: 'intcal13'>")
  expect_true(
    all(
      grepl(
        "^<age_calibration_curve",
        format(structure(list(intcal13, intcal09, intcal04), class = "age_calibration_curve_list"))
      )
    )
  )
  expect_true(
    all(
      grepl(
        "^<age_calibration_curve",
        format(
          head(structure(list(intcal13, intcal09, intcal04), class = "age_calibration_curve_list"))
        )
      )
    )
  )

  tibble::tibble(
    i = 1:3,
    curves = structure(list(intcal13, intcal09, intcal04), class = "age_calibration_curve_list")
  )
})

test_that("reading of calibration curves from .14c files works", {
  curve_dir <- system.file("curves", package = "carbon14")
  curves <- list.files(curve_dir, full.names = TRUE)
  for(curve in curves) {
    expect_is(read_14c(curve), "age_calibration_curve")
  }
  expect_identical(read_14c(file.path(curve_dir, "intcal13.14c")), intcal13)

  # check read from url
  expect_is(read_14c("http://www.radiocarbon.org/IntCal13%20files/intcal13.14c"),
            "age_calibration_curve")
})
