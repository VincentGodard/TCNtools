# Reference values ------------------------------------------------------------
# stone2000: US Standard Atmosphere 1976 (NOAA/NASA/USAF), 1013.25 hPa at sea
# level, 898.76 hPa at 1000 m, 701.12 hPa at 3000 m. Stone (2000) uses
# gM/R = 0.03417 K/m, which rounds the US76 exponent: agreement to 0.1 hPa.
# era40: values computed by hand from the grid nodes of
# data-raw/era40/era40_grid.csv and equation 1 of Stone (2000).

test_that("stone2000 matches the US Standard Atmosphere 1976", {
  p <- atm_pressure(c(0, 1000, 3000), model = "stone2000")
  expect_equal(p[1], 1013.25)
  expect_equal(p[2:3], c(898.76, 701.12), tolerance = 0.15 / 900)
})

test_that("stone2000 matches equation 1 of Stone (2000)", {
  # closed form: Ps (1 - xi z / Ts)^(gM / (R xi))
  z <- c(-400, 2500, 5500)
  expected <- 1013.25 * (1 - 0.0065 * z / 288.15)^(0.03417 / 0.0065)
  expect_equal(atm_pressure(z, model = "stone2000"), expected)
})

test_that("stone2000 ignores latitude and longitude", {
  expect_equal(
    atm_pressure(2000, latitude = 60, longitude = 10, model = "stone2000"),
    atm_pressure(2000, model = "stone2000")
  )
})

test_that("era40 at a grid node matches a hand calculation", {
  # node 45 N, 5 E: slp = 1017.36367131 hPa, T1000 = 286.961724101 K,
  # CIRA-86 lapse rate at 45 degrees = 0.006041679289610938 K/m
  expect_equal(atm_pressure(0, 45, 5, model = "era40"), 1017.36367131)
  expect_equal(atm_pressure(1000, 45, 5, model = "era40"), 902.008039201894)
})

test_that("era40 interpolates bilinearly between nodes", {
  # centre of the cell 45-47.5 N, 5-7.5 E: mean of the four corner values
  expect_equal(atm_pressure(0, 46.25, 6.25, model = "era40"),
               1017.3008616075)
})

test_that("era40 wraps around the Greenwich meridian", {
  # cell -2.5-0 N, 357.5-360 E, centre at longitude -1.25 = 358.75
  expected <- 1012.517775455
  expect_equal(atm_pressure(0, -1.25, -1.25, model = "era40"), expected)
  expect_equal(atm_pressure(0, -1.25, 358.75, model = "era40"), expected)
  expect_equal(atm_pressure(0, 30, 0, model = "era40"),
               atm_pressure(0, 30, 360, model = "era40"))
})

test_that("era40 works at the poles", {
  expect_equal(atm_pressure(0, 90, 123, model = "era40"), 1014.33806118)
  expect_equal(atm_pressure(0, -90, 0, model = "era40"), 1005.27893578)
})

test_that("era40 follows the Antarctic relationship of Stone (2000)", {
  # Stone (2000), equation 5: P = 989.1 exp(-z / 7588), RMS error 3 hPa
  # on 67 stations; South Pole, Vostok, Dome C
  z <- c(2835, 3488, 3233)
  p <- atm_pressure(z, c(-90, -78.46, -75.1), c(0, 106.8, 123.3),
                    model = "era40")
  expect_equal(p, 989.1 * exp(-z / 7588), tolerance = 0.015)
  # the standard atmosphere is 35 to 40 hPa too high there
  expect_true(all(atm_pressure(z, model = "stone2000") - p > 30))
})

test_that("pressure decreases with elevation", {
  z <- seq(-400, 6000, by = 100)
  expect_true(all(diff(atm_pressure(z, model = "stone2000")) < 0))
  expect_true(all(diff(atm_pressure(z, -20, 290, model = "era40")) < 0))
})

test_that("arguments are recycled and match a scalar loop", {
  z <- c(0, 1500, 4500)
  lat <- c(-16, 32, 60)
  lon <- c(-68, 78, 10)
  vec <- atm_pressure(z, lat, lon, model = "era40")
  loop <- vapply(1:3, function(i) atm_pressure(z[i], lat[i], lon[i],
                                               model = "era40"), numeric(1))
  expect_equal(vec, loop)
  expect_length(atm_pressure(z, 45, 5, model = "era40"), 3)
  expect_length(atm_pressure(1000, lat, lon, model = "era40"), 3)
})

test_that("invalid input gives an error", {
  expect_error(atm_pressure(1000), "`model` must be given")
  expect_error(atm_pressure(1000, model = "isa"), "`model` must be")
  expect_error(atm_pressure(1000, model = c("era40", "stone2000")),
               "`model` must be")
  expect_error(atm_pressure(1000, model = "era40"), "required")
  expect_error(atm_pressure(1000, latitude = 45, model = "era40"), "required")
  expect_error(atm_pressure("1000", model = "stone2000"), "`elevation`")
  expect_error(atm_pressure(NA_real_, model = "stone2000"), "`elevation`")
  expect_error(atm_pressure(1000, 91, 0, model = "era40"), "`latitude`")
  expect_error(atm_pressure(1000, 45, -190, model = "era40"), "`longitude`")
  expect_error(atm_pressure(1:3, c(10, 20), 0, model = "era40"),
               "length 1 or 3")
  expect_error(atm_pressure(50000, model = "stone2000"), "outside the range")
})
