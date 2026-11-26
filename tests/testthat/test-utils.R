test_that("angle conversions are exact and inverse", {
  expect_equal(deg2rad(c(0, 90, 180, -360)), c(0, pi / 2, pi, -2 * pi))
  expect_equal(rad2deg(c(0, pi / 2, pi)), c(0, 90, 180))
  x <- c(-123.4, 0.5, 89.9)
  expect_equal(rad2deg(deg2rad(x)), x)
})

test_that("angdist matches CRONUScalc angdist.m", {
  # Reference values: CRONUScalc common/production/angdist.m (commit
  # cdbc8f2) run under GNU Octave 11.1.
  lat1 <- c(45.2, -33.9, 64.1, 0, 10, 89.5)
  lon1 <- c(5.7, 18.4, -21.9, 0, -170, 120)
  lat2 <- c(80.1, 80.1, -75, 0, -10, -89.5)
  lon2 <- c(-72.6, -72.6, 130, 180, 175, -60)
  ref <- c(43.649828265412, 123.499549167909, 165.611451346516, 180,
           24.954120252232, 180)
  expect_equal(angdist(lat1, lon1, lat2, lon2), ref, tolerance = 1e-10)
})

test_that("angdist has the expected geometric limits", {
  # Identical points: 0, never NaN from rounding
  lat <- c(0, 12.345, -67.89, 90)
  lon <- c(0, 98.765, -123.4, 0)
  expect_equal(angdist(lat, lon, lat, lon), rep(0, 4), tolerance = 1e-6)
  expect_false(anyNA(angdist(lat, lon, lat, lon)))
  # Pole to any point: 90 - latitude
  expect_equal(angdist(90, 0, c(-30, 0, 45), c(10, -100, 170)),
               c(120, 90, 45))
  # Along the equator: longitude difference, wrapped to [0, 180]
  expect_equal(angdist(0, -170, 0, c(-160, 170, 10)), c(10, 20, 180))
  # Symmetry and recycling
  expect_equal(angdist(10, 20, 30, 40), angdist(30, 40, 10, 20))
  expect_length(angdist(0, 0, c(1, 2, 3), 0), 3)
})

test_that("recycled_length accepts length 1 or n only", {
  expect_equal(recycled_length(a = 1, b = 1:4), 4L)
  expect_equal(recycled_length(a = 1), 1L)
  expect_error(recycled_length(a = 1:2, b = 1:4), "`a`")
})

test_that("interp_bilinear reproduces nodes and plane functions", {
  x <- c(0, 1, 3)
  y <- c(10, 20)
  f <- function(x, y) 2 + 3 * x - 0.5 * y
  z <- outer(x, y, f)
  expect_equal(interp_bilinear(x, y, z, c(0, 3, 1), c(10, 20, 20)),
               f(c(0, 3, 1), c(10, 20, 20)))
  expect_equal(interp_bilinear(x, y, z, c(0.3, 2.2), c(12, 17.5)),
               f(c(0.3, 2.2), c(12, 17.5)))
  expect_error(interp_bilinear(x, y, z, 4, 15), "outside the grid")
})
