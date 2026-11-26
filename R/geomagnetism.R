#' Geomagnetic field history
#'
#' Dipole moment and geomagnetic pole position through time, as used by
#' the Lm scaling scheme (time-dependent Lal 1991 / Stone 2000 scaling) in
#' the CRONUScalc code \insertCite{marrero2016cosmogenic}{TCNtools}, with
#' which the production rates of
#' \insertCite{borchers2016geological;textual}{TCNtools} returned by
#' [tcn_params()] were calibrated.
#'
#' The history is the framework of
#' \insertCite{lifton2014scaling;textual}{TCNtools}, taken from the file
#' `pmag_consts.mat` of CRONUScalc:
#' \itemize{
#'   \item 0 to 60 yr (2010 to 1950 CE): DGRF / IGRF-11
#'     \insertCite{finlay2010international}{TCNtools};
#'   \item 60 yr to 3 ka: CALS3k.3 \insertCite{korte2009geomagnetic}{TCNtools};
#'   \item 3 to 7 ka: CALS7k.2 \insertCite{korte2005continuous}{TCNtools};
#'   \item 7 to 18 ka: GLOPIS-75 \insertCite{laj2004high}{TCNtools};
#'   \item 18 ka to 2 Ma: PADM2M \insertCite{ziegler2011padm2m}{TCNtools};
#'   \item 2 to 10 Ma: constant moment ratio 0.7048, the value used by
#'     CRONUScalc (its origin is not documented, it is not the mean of the
#'     record).
#' }
#' Up to 6960 yr the pole is that of the geocentric dipole fitted to the
#' field models. From 7060 yr the field is taken as a geocentric axial
#' dipole and the pole is the geographic north pole.
#'
#' Between the nodes of the record (10 yr apart before 60 yr, 100 yr to
#' 50 ka, then 1 kyr to 2 Ma), the moment ratio is interpolated linearly
#' and the pole along the great circle (linear interpolation of the unit
#' vectors, then normalisation). Beyond 10 Ma the last values are kept.
#' CRONUScalc only evaluates the history at the nodes, so the two agree
#' there.
#'
#' @param time Numeric vector. Time before 2010 CE (years, >= 0). All
#'   times in the package are counted from 2010, the reference year of the
#'   dipole moment \eqn{M_0}.
#' @param scaling Character scalar. Scaling scheme, `"Lm"` only for now.
#'   No default. `"St"` does not depend on the geomagnetic field, and the
#'   `"Sf"` and `"Sa"` schemes (with their own history of the cutoff
#'   rigidity and solar modulation) arrive in version 1.1.
#'
#' @section Vectorisation:
#' `time` is a vector of any length, the result has one row per element.
#' `scaling` is a scalar.
#'
#' @return Data frame with one row per element of `time` and columns
#'   \describe{
#'     \item{`time`}{time before 2010 (years).}
#'     \item{`moment_ratio`}{dipole moment relative to the 2010 DGRF value
#'       \eqn{M_0 = 7.746 \times 10^{22}}{M0 = 7.746e22} A m², unitless.}
#'     \item{`pole_latitude`}{latitude of the geomagnetic north pole
#'       (degrees).}
#'     \item{`pole_longitude`}{longitude of the geomagnetic north pole
#'       (degrees east, 0 to 360).}
#'   }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [cutoff_rigidity()], which takes these columns as input.
#'
#' @examples
#' # Last 20 kyr
#' h <- geomag_history(seq(0, 20000, by = 2000), scaling = "Lm")
#' h
#'
#' # Dipole moment over the last 2 Myr
#' t <- seq(0, 2e6, by = 1000)
#' plot(t / 1000, geomag_history(t, scaling = "Lm")$moment_ratio,
#'      type = "l", xlab = "Time before 2010 (ka)", ylab = "M / M0")
#' @export
geomag_history <- function(time, scaling) {
  check_geomag_scaling(scaling)
  check_finite_numeric(time, "time")
  if (any(time < 0)) {
    stop("`time` must be >= 0 (years before 2010).", call. = FALSE)
  }

  nodes <- geomag
  moment_ratio <- stats::approx(nodes$time, nodes$moment_ratio, xout = time,
                                rule = 2, ties = "ordered")$y

  # Pole: interpolate the unit vectors, then normalise.
  i <- findInterval(time, nodes$time, rightmost.closed = TRUE)
  i[i == nrow(nodes)] <- nrow(nodes) - 1L
  u <- (time - nodes$time[i]) / (nodes$time[i + 1L] - nodes$time[i])
  u <- pmin(u, 1)
  v0 <- unit_vector(nodes$pole_latitude[i], nodes$pole_longitude[i])
  v1 <- unit_vector(nodes$pole_latitude[i + 1L], nodes$pole_longitude[i + 1L])
  v <- (1 - u) * v0 + u * v1
  v <- v / sqrt(rowSums(v^2))
  pole_latitude <- rad2deg(asin(pmin(pmax(v[, 3], -1), 1)))
  pole_longitude <- rad2deg(atan2(v[, 2], v[, 1])) %% 360
  # Longitude is undefined at the geographic pole: use 0, as in the data.
  at_pole <- abs(pole_latitude) > 90 - 1e-9
  pole_latitude[at_pole] <- sign(pole_latitude[at_pole]) * 90
  pole_longitude[at_pole] <- 0

  data.frame(time = time, moment_ratio = moment_ratio,
             pole_latitude = pole_latitude, pole_longitude = pole_longitude)
}

#' Geomagnetic cutoff rigidity
#'
#' Effective vertical cutoff rigidity of a site for a geocentric dipole
#' field, as in the Lm scaling scheme of CRONUScalc
#' \insertCite{marrero2016cosmogenic}{TCNtools}, from equation 1 of
#' \insertCite{dunai2001influence;textual}{TCNtools}:
#' \deqn{R_c = 14.9 \, \frac{M}{M_0} \cos^4 \lambda_m}{Rc = 14.9 (M/M0) cos^4(lambda_m)}
#' where \eqn{\lambda_m}{lambda_m} is the geomagnetic latitude of the site,
#' \eqn{90 - \theta}{90 - theta}, \eqn{\theta}{theta} being the angular
#' distance between the site and the geomagnetic north pole. With the
#' default pole (geographic north pole, axial dipole) the geomagnetic
#' latitude is the geographic latitude.
#'
#' The function only does the physics: the moment ratio and the pole
#' position at a given time come from [geomag_history()].
#'
#' @param latitude Numeric vector. Site latitude (degrees, -90 to 90).
#' @param longitude Numeric vector. Site longitude (degrees east, -180 to
#'   360).
#' @param moment_ratio Numeric vector. Dipole moment relative to the 2010
#'   value \eqn{M_0}{M0} (unitless, >= 0).
#' @param pole_latitude Numeric vector. Latitude of the geomagnetic north
#'   pole (degrees, -90 to 90). Default 90 (axial dipole).
#' @param pole_longitude Numeric vector. Longitude of the geomagnetic north
#'   pole (degrees east). Default 0.
#' @param scaling Character scalar. Scaling scheme, `"Lm"` only for now.
#'   No default. `"St"` does not use the cutoff rigidity, and the `"Sf"`
#'   and `"Sa"` schemes (which use another relation between moment and
#'   cutoff rigidity) arrive in version 1.1.
#'
#' @section Vectorisation:
#' `latitude`, `longitude`, `moment_ratio`, `pole_latitude` and
#' `pole_longitude` are either of length 1 or of the same length n. The
#' result has length n. This covers one site through time (site of length
#' 1, history of length n) as well as several sites at one time. Other
#' lengths give an error. `scaling` is a scalar.
#'
#' @return Numeric vector of cutoff rigidity (GV).
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [geomag_history()].
#'
#' @examples
#' # Present-day axial dipole, from the equator to the pole
#' cutoff_rigidity(c(0, 30, 60, 90), 0, moment_ratio = 1, scaling = "Lm")
#'
#' # One site through the last 10 kyr
#' h <- geomag_history(seq(0, 10000, by = 1000), scaling = "Lm")
#' cutoff_rigidity(45.2, 5.7, h$moment_ratio, h$pole_latitude,
#'                 h$pole_longitude, scaling = "Lm")
#' @export
cutoff_rigidity <- function(latitude, longitude, moment_ratio,
                            pole_latitude = 90, pole_longitude = 0,
                            scaling) {
  check_geomag_scaling(scaling)
  check_finite_numeric(latitude, "latitude")
  check_finite_numeric(longitude, "longitude")
  check_finite_numeric(moment_ratio, "moment_ratio")
  check_finite_numeric(pole_latitude, "pole_latitude")
  check_finite_numeric(pole_longitude, "pole_longitude")
  if (any(abs(latitude) > 90)) {
    stop("`latitude` must be between -90 and 90 degrees.", call. = FALSE)
  }
  if (any(longitude < -180 | longitude > 360)) {
    stop("`longitude` must be between -180 and 360 degrees.", call. = FALSE)
  }
  if (any(abs(pole_latitude) > 90)) {
    stop("`pole_latitude` must be between -90 and 90 degrees.",
         call. = FALSE)
  }
  if (any(moment_ratio < 0)) {
    stop("`moment_ratio` must be >= 0.", call. = FALSE)
  }
  recycled_length(latitude = latitude, longitude = longitude,
                  moment_ratio = moment_ratio, pole_latitude = pole_latitude,
                  pole_longitude = pole_longitude)

  geomag_latitude <- 90 - angdist(latitude, longitude,
                                  pole_latitude, pole_longitude)
  14.9 * moment_ratio * cos(deg2rad(geomag_latitude))^4
}

#' Check the scaling scheme given to the geomagnetic functions
#' @param scaling Value of the `scaling` argument.
#' @return `scaling`, invisibly.
#' @noRd
check_geomag_scaling <- function(scaling) {
  if (missing(scaling)) {
    stop("`scaling` must be given: \"Lm\".", call. = FALSE)
  }
  check_choice(scaling, scaling_schemes)
  if (scaling %in% c("Sf", "Sa")) {
    stop("Sf/Sa scaling not available before v1.1.", call. = FALSE)
  }
  if (scaling == "St") {
    stop("St scaling does not depend on the geomagnetic field.",
         call. = FALSE)
  }
  invisible(scaling)
}

#' Unit vectors from latitude and longitude
#' @param latitude,longitude Angles (degrees), same length.
#' @return Matrix with columns x, y, z.
#' @noRd
unit_vector <- function(latitude, longitude) {
  lat <- deg2rad(latitude)
  lon <- deg2rad(longitude)
  cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
}
