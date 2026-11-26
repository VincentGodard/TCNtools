# Internal helpers shared across the package.

#' Degrees to radians
#' @param x Angle (degrees).
#' @return Angle (radians).
#' @noRd
deg2rad <- function(x) x * pi / 180

#' Radians to degrees
#' @param x Angle (radians).
#' @return Angle (degrees).
#' @noRd
rad2deg <- function(x) x * 180 / pi

#' Angular distance between points on a sphere
#'
#' Great-circle angle between two points given by latitude and longitude,
#' from the dot product of the two unit vectors, as in `angdist.m` of
#' CRONUScalc (Balco, 2001). The cosine is clamped to \[-1, 1\] so that
#' rounding cannot produce `NaN` for identical or antipodal points (a
#' deliberate departure from `angdist.m`). Used to compute the geomagnetic
#' latitude of a site, `90 - angdist(site, paleopole)`.
#'
#' Arguments are recycled (length 1 or a common length).
#'
#' @param lat1,lon1 Latitude and longitude of the first point(s) (degrees).
#' @param lat2,lon2 Latitude and longitude of the second point(s) (degrees).
#' @return Angular distance (degrees, between 0 and 180).
#' @noRd
angdist <- function(lat1, lon1, lat2, lon2) {
  lat1 <- deg2rad(lat1)
  lon1 <- deg2rad(lon1)
  lat2 <- deg2rad(lat2)
  lon2 <- deg2rad(lon2)
  cos_theta <- cos(lat1) * cos(lat2) *
    (cos(lon1) * cos(lon2) + sin(lon1) * sin(lon2)) +
    sin(lat1) * sin(lat2)
  rad2deg(acos(pmin(pmax(cos_theta, -1), 1)))
}

#' Check that an argument is a non-empty vector of finite numbers
#' @param x Value to check.
#' @param name Argument name used in the error message.
#' @return `x`, invisibly.
#' @noRd
check_finite_numeric <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    stop("`", name, "` must be a non-empty numeric vector of finite values.",
         call. = FALSE)
  }
  invisible(x)
}

#' Common length of recycled arguments
#'
#' Each argument must have length 1 or the length of the longest one.
#'
#' @param ... Named vectors.
#' @return The common length.
#' @noRd
recycled_length <- function(...) {
  lengths <- lengths(list(...))
  n <- max(lengths)
  bad <- lengths != 1L & lengths != n
  if (any(bad)) {
    stop("Arguments ", paste0("`", names(lengths)[bad], "`", collapse = ", "),
         " must have length 1 or ", n, ".", call. = FALSE)
  }
  n
}

#' Bilinear interpolation on a regular grid
#'
#' @param x,y Increasing grid coordinates (rows and columns of `z`).
#' @param z Matrix of values, `z[i, j]` at `(x[i], y[j])`.
#' @param xout,yout Points to interpolate at, within the grid.
#' @return Interpolated values, one per point.
#' @noRd
interp_bilinear <- function(x, y, z, xout, yout) {
  i <- findInterval(xout, x, rightmost.closed = TRUE)
  j <- findInterval(yout, y, rightmost.closed = TRUE)
  if (any(i == 0L | i == length(x) | j == 0L | j == length(y))) {
    stop("Interpolation point outside the grid.", call. = FALSE)
  }
  u <- (xout - x[i]) / (x[i + 1L] - x[i])
  v <- (yout - y[j]) / (y[j + 1L] - y[j])
  (1 - u) * (1 - v) * z[cbind(i, j)] + u * (1 - v) * z[cbind(i + 1L, j)] +
    (1 - u) * v * z[cbind(i, j + 1L)] + u * v * z[cbind(i + 1L, j + 1L)]
}
