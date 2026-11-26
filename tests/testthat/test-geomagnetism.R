history <- function(time) geomag_history(time, model = "lifton2014",
                                         scaling = "Lm")
rc_lm <- function(...) cutoff_rigidity(..., scaling = "Lm")

# geomag_history() ----------------------------------------------------------

test_that("geomag_history returns the record at its nodes", {
  h <- history(geomag$time)
  expect_named(h, c("time", "moment_ratio", "pole_latitude",
                    "pole_longitude"))
  expect_equal(h$time, geomag$time)
  expect_equal(h$moment_ratio, geomag$moment_ratio, tolerance = 1e-12)
  expect_equal(h$pole_latitude, geomag$pole_latitude, tolerance = 1e-12)
  expect_equal(h$pole_longitude, geomag$pole_longitude, tolerance = 1e-12)
})

test_that("geomag_history starts at the 2010 field", {
  h <- history(0)
  expect_equal(h$moment_ratio, 1)
  expect_equal(c(h$pole_latitude, h$pole_longitude), c(80.01, 287.78),
               tolerance = 1e-12)
})

test_that("geomag_history interpolates between nodes", {
  # Moment ratio: linear
  t0 <- 160
  t1 <- 260
  h <- history(c(t0, t1, (t0 + t1) / 2))
  expect_equal(h$moment_ratio[3], mean(h$moment_ratio[1:2]))
  # Pole: midpoint of the great circle arc
  d01 <- angdist(h$pole_latitude[1], h$pole_longitude[1],
                 h$pole_latitude[2], h$pole_longitude[2])
  d0m <- angdist(h$pole_latitude[1], h$pole_longitude[1],
                 h$pole_latitude[3], h$pole_longitude[3])
  dm1 <- angdist(h$pole_latitude[3], h$pole_longitude[3],
                 h$pole_latitude[2], h$pole_longitude[2])
  expect_equal(c(d0m, dm1), c(d01, d01) / 2, tolerance = 1e-9)
  # Interpolated values stay within the bracketing nodes
  t <- seq(t0, t1, by = 5)
  m <- history(t)$moment_ratio
  expect_true(all(m >= min(h$moment_ratio[1:2]) &
                    m <= max(h$moment_ratio[1:2])))
})

test_that("geomag_history uses an axial dipole from 7060 yr", {
  h <- history(c(7060, 1e4, 1e6, 1e7))
  expect_equal(h$pole_latitude, rep(90, 4))
  expect_equal(h$pole_longitude, rep(0, 4))
  # Between the last pole node (6960 yr) and 7060 yr the pole moves
  # towards the geographic pole
  lat <- history(seq(6960, 7060, by = 20))$pole_latitude
  expect_true(all(diff(lat) > 0))
})

test_that("geomag_history is constant beyond 10 Ma", {
  h <- history(c(1e7, 2e7, 1e9))
  expect_equal(h$moment_ratio, rep(0.7048102, 3), tolerance = 1e-6)
  expect_equal(nrow(h), 3L)
})

test_that("geomag_history validates its arguments", {
  expect_error(history(-1), ">= 0")
  expect_error(history(NA_real_), "finite")
  expect_error(history(numeric(0)), "non-empty")
  expect_error(geomag_history(0, scaling = "Lm"), "`model` must be given")
  expect_error(geomag_history(0, "other", "Lm"), "`model` must be one of")
  expect_error(geomag_history(0, "lifton2014"), "`scaling` must be given")
  expect_error(geomag_history(0, "lifton2014", "St"), "does not depend")
  expect_error(geomag_history(0, "lifton2014", "Sf"), "v1.1")
  expect_error(geomag_history(0, "lifton2014", "Sa"), "v1.1")
  expect_error(geomag_history(0, "lifton2014", "lm"), "must be one of")
})

# cutoff_rigidity() ---------------------------------------------------------

test_that("cutoff_rigidity matches CRONUScalc get_tdsf.m (Lm)", {
  ref <- utils::read.csv(test_path("fixtures", "cronuscalc_lmrc.csv"),
                         comment.char = "#")
  h <- history(ref$time)
  rc <- rc_lm(ref$latitude, ref$longitude, h$moment_ratio,
              h$pole_latitude, h$pole_longitude)
  expect_equal(rc, ref$rc, tolerance = 1e-10)
})

test_that("cutoff_rigidity has the expected dipole limits", {
  # Axial dipole: 14.9 M/M0 cos^4(latitude)
  lat <- c(-90, -60, -30, 0, 30, 60, 90)
  expect_equal(rc_lm(lat, 0, 1), 14.9 * cos(lat * pi / 180)^4)
  expect_equal(rc_lm(0, 123, 1), 14.9)
  expect_equal(rc_lm(c(90, -90), 0, 1), c(0, 0), tolerance = 1e-12)
  # Proportional to the moment ratio, zero without field
  expect_equal(rc_lm(40, 10, c(0, 0.5, 2)), c(0, 0.5, 2) * rc_lm(40, 10, 1))
  # Symmetric about the geomagnetic equator, independent of longitude
  # for an axial dipole
  expect_equal(rc_lm(35, 0:3 * 90, 1), rep(rc_lm(-35, 0, 1), 4))
  # Site at the geomagnetic pole, or at the antipode: zero
  expect_equal(rc_lm(c(80, -80), c(-70, 110), 1, 80, 290), c(0, 0),
               tolerance = 1e-12)
  # Site 90 degrees from the pole: geomagnetic equator
  expect_equal(rc_lm(-10, 290, 0.8, 80, 290), 0.8 * 14.9)
})

test_that("cutoff_rigidity accepts both longitude conventions", {
  expect_equal(rc_lm(45, -82, 1, 80, 288), rc_lm(45, 278, 1, 80, 288))
  expect_equal(rc_lm(45, 5, 1, 80, -72), rc_lm(45, 5, 1, 80, 288))
})

test_that("cutoff_rigidity recycles arguments", {
  h <- history(c(0, 3000, 6000))
  rc <- rc_lm(45, 5, h$moment_ratio, h$pole_latitude, h$pole_longitude)
  expect_length(rc, 3)
  expect_equal(rc[2], rc_lm(45, 5, h$moment_ratio[2], h$pole_latitude[2],
                            h$pole_longitude[2]))
  expect_length(rc_lm(c(0, 10, 20), 0, 1), 3)
  expect_error(rc_lm(c(0, 10), 0, c(1, 1, 1)), "must have length 1 or 3")
})

test_that("cutoff_rigidity validates its arguments", {
  expect_error(cutoff_rigidity(45, 5, 1), "`scaling` must be given")
  expect_error(cutoff_rigidity(45, 5, 1, scaling = "St"), "does not depend")
  expect_error(cutoff_rigidity(45, 5, 1, scaling = "Sa"), "v1.1")
  expect_error(rc_lm(91, 5, 1), "`latitude`")
  expect_error(rc_lm(45, 400, 1), "`longitude`")
  expect_error(rc_lm(45, 5, -0.1), "`moment_ratio`")
  expect_error(rc_lm(45, 5, 1, pole_latitude = 95), "`pole_latitude`")
  expect_error(rc_lm(45, 5, NA_real_), "finite")
})
