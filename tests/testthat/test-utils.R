test_that("angle conversions are exact and inverse", {
  expect_equal(deg2rad(c(0, 90, 180, -360)), c(0, pi / 2, pi, -2 * pi))
  expect_equal(rad2deg(c(0, pi / 2, pi)), c(0, 90, 180))
  x <- c(-123.4, 0.5, 89.9)
  expect_equal(rad2deg(deg2rad(x)), x)
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
