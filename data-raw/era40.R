# ERA-40 climatology grid used by atm_pressure(model = "era40").
#
# Source: ERA40.mat of the CRONUScalc code (Marrero et al. 2016,
# marrero2016cosmogenic), https://bitbucket.org/cronusearth/cronus-calc,
# file common/production/ERA40.mat (commit cdbc8f2, 2021-09-14, md5
# bd9f751f13f05bb1edf743bca7078ba8, MATLAB file created 2011-03-31). The
# file dev/data/data_internal/ERA40.mat at tag v0.1.0 is byte-identical.
# It holds the long-term mean sea level pressure (hPa) and mean 1000 hPa
# temperature (K) of the ERA-40 reanalysis (Uppala et al. 2005,
# uppala2005era40) on a 2.5 degree grid. It is read by ERA40atm.m, written
# by G. Balco (2007) for NCEP-NCAR fields and modified by N. Lifton
# (March 2011) to use ERA-40 and a latitude-dependent lapse rate. The same
# grid is used by CREp (Martin et al. 2017, martin2017crep), which credits
# it to the LSD code of Lifton et al. (2014).
# The v0.1.0 build script credited Balco et al. (2008), but that paper
# uses NCEP-NCAR fields, not ERA-40.
#
# Converted to data-raw/era40/era40_grid.csv by
# data-raw/era40/era40_from_mat.py (see that script for details).
#
# Result: a list with
#   latitude   73 values, -90 to 90 by 2.5 (degrees)
#   longitude  145 values, 0 to 360 by 2.5 (degrees east); the 360 column
#              repeats the 0 column so that interpolation wraps around
#   slp        73 x 145 matrix, mean sea level pressure (hPa)
#   t1000      73 x 145 matrix, mean 1000 hPa temperature (K)
# Matrices are indexed [latitude, longitude].

era40 <- local({
  g <- utils::read.csv("data-raw/era40/era40_grid.csv")
  latitude <- sort(unique(g$latitude))
  longitude <- sort(unique(g$longitude))
  stopifnot(
    all.equal(latitude, seq(-90, 90, by = 2.5)),
    all.equal(longitude, seq(0, 357.5, by = 2.5)),
    nrow(g) == length(latitude) * length(longitude),
    !anyNA(g)
  )
  grid <- function(v) {
    m <- matrix(NA_real_, length(latitude), length(longitude))
    m[cbind(match(g$latitude, latitude), match(g$longitude, longitude))] <- v
    cbind(m, m[, 1]) # wrap: 360 = 0
  }
  list(
    latitude = latitude,
    longitude = c(longitude, 360),
    slp = grid(g$slp_hpa),
    t1000 = grid(g$t1000_k)
  )
})
