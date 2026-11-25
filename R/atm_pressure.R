#' Mean atmospheric pressure at a site
#'
#' Converts site elevation into mean atmospheric pressure, the quantity
#' that controls the altitude dependence of cosmogenic production in the
#' scaling schemes.
#'
#' Two models are available.
#'
#' `"stone2000"` is the standard atmosphere of
#' \insertCite{stone2000air;textual}{TCNtools} (equation 1), with fixed sea
#' level pressure (1013.25 hPa), sea level temperature (288.15 K) and lapse
#' rate (0.0065 K/m):
#' \deqn{P = P_s \exp\left(-\frac{g M}{R \xi}
#'   \left[\ln T_s - \ln(T_s - \xi z)\right]\right)}{P = Ps exp(-gM/(R xi) [ln Ts - ln(Ts - xi z)])}
#' with \eqn{gM/R = 0.03417} K/m. It depends on elevation only.
#'
#' `"era40"` uses the same equation, but takes the sea level pressure and
#' the temperature (at 1000 hPa) from the long-term means of the ERA-40
#' reanalysis \insertCite{uppala2005era40}{TCNtools}, interpolated
#' bilinearly at the site on a 2.5 degree grid, and a lapse rate that
#' depends on latitude (sixth-degree polynomial fit to the COSPAR CIRA-86
#' reference atmosphere below 10 km). It reproduces `ERA40atm.m` of the
#' CRONUScalc code \insertCite{marrero2016cosmogenic}{TCNtools}, written by
#' G. Balco for NCEP-NCAR fields \insertCite{balco2008complete}{TCNtools}
#' and adapted to ERA-40 by N. Lifton for the LSD code
#' \insertCite{lifton2014scaling}{TCNtools}, with the same data file. The
#' lapse-rate fit is only published in that code. The same scheme is used by
#' CREp \insertCite{martin2017crep}{TCNtools}. The production rates
#' of \insertCite{borchers2016geological;textual}{TCNtools} returned by
#' [tcn_params()] were calibrated with ERA-40 pressures, so `"era40"` is
#' the consistent choice with these rates. For the calibration samples,
#' \insertCite{marrero2016cosmogenic;textual}{TCNtools} assigned an
#' uncertainty of 5 hPa to these pressures, from the largest differences
#' between the NCEP-NCAR and ERA-40 reanalyses.
#'
#' In Antarctica, ground level pressures are 20 to 40 hPa lower than the
#' standard atmosphere at all altitudes \insertCite{stone2000air}{TCNtools},
#' so `"stone2000"` must not be used there. `"era40"` captures this
#' anomaly and needs no separate Antarctic relationship
#' \insertCite{marrero2016cosmogenic}{TCNtools}: CRONUScalc uses it
#' everywhere. It stays within about 1% of the Antarctic relationship
#' \eqn{P = 989.1 \exp(-z/7588)} of \insertCite{stone2000air;textual}{TCNtools}
#' (equation 5) at the South Pole, Vostok and Dome C.
#'
#' @param elevation Numeric vector. Site elevation above sea level (m).
#'   Negative values (below sea level) are accepted.
#' @param latitude Numeric vector. Site latitude (degrees, -90 to 90).
#'   Required for `model = "era40"`, ignored by `"stone2000"`.
#' @param longitude Numeric vector. Site longitude (degrees east, -180 to
#'   360). Required for `model = "era40"`, ignored by `"stone2000"`.
#' @param model Character scalar. `"era40"` or `"stone2000"`. No default:
#'   the model must be chosen explicitly.
#'
#' @section Vectorisation:
#' `elevation`, `latitude` and `longitude` are either of length 1 (same
#' value for all sites) or of the same length n (one value per site). The
#' result has length n. Other lengths give an error. `model` is a scalar.
#'
#' @return Numeric vector of mean atmospheric pressure (hPa), one value
#'   per site.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Standard atmosphere
#' atm_pressure(c(0, 1000, 3000), model = "stone2000")
#'
#' # ERA-40 for the same elevations at a mid-latitude site
#' atm_pressure(c(0, 1000, 3000), latitude = 45, longitude = 6,
#'              model = "era40")
#'
#' # Several sites at once
#' atm_pressure(elevation = c(4500, 4500), latitude = c(-16, 32),
#'              longitude = c(-68, 78), model = "era40")
#' @export
atm_pressure <- function(elevation, latitude = NULL, longitude = NULL,
                         model) {
  if (missing(model)) {
    stop("`model` must be given: \"era40\" or \"stone2000\".", call. = FALSE)
  }
  if (!is.character(model) || length(model) != 1L ||
      !model %in% c("era40", "stone2000")) {
    stop("`model` must be \"era40\" or \"stone2000\".", call. = FALSE)
  }
  check_finite_numeric(elevation, "elevation")

  if (model == "stone2000") {
    return(pressure_standard(elevation, slp = 1013.25, t0 = 288.15,
                             lapse_rate = 0.0065))
  }

  if (is.null(latitude) || is.null(longitude)) {
    stop("`latitude` and `longitude` are required for model = \"era40\".",
         call. = FALSE)
  }
  check_finite_numeric(latitude, "latitude")
  check_finite_numeric(longitude, "longitude")
  if (any(abs(latitude) > 90)) {
    stop("`latitude` must be between -90 and 90 degrees.", call. = FALSE)
  }
  if (any(longitude < -180 | longitude > 360)) {
    stop("`longitude` must be between -180 and 360 degrees.", call. = FALSE)
  }
  n <- recycled_length(elevation = elevation, latitude = latitude,
                       longitude = longitude)
  elevation <- rep_len(elevation, n)
  latitude <- rep_len(latitude, n)
  longitude <- rep_len(longitude, n) %% 360

  slp <- interp_bilinear(era40$latitude, era40$longitude, era40$slp,
                         latitude, longitude)
  t1000 <- interp_bilinear(era40$latitude, era40$longitude, era40$t1000,
                           latitude, longitude)
  pressure_standard(elevation, slp = slp, t0 = t1000,
                    lapse_rate = lapse_rate_cira86(latitude))
}

#' Standard atmosphere equation (Stone 2000, equation 1)
#'
#' @param elevation Elevation (m).
#' @param slp Sea level pressure (hPa).
#' @param t0 Sea level temperature (K).
#' @param lapse_rate Temperature lapse rate (K/m, positive when temperature
#'   decreases upward).
#' @return Pressure (hPa).
#' @noRd
pressure_standard <- function(elevation, slp, t0, lapse_rate) {
  g_m_over_r <- 0.03417 # K/m
  t_site <- t0 - lapse_rate * elevation
  if (any(t_site <= 0)) {
    stop("`elevation` is outside the range of the atmosphere model.",
         call. = FALSE)
  }
  slp * exp(-g_m_over_r / lapse_rate * (log(t0) - log(t_site)))
}

#' Lapse rate as a function of latitude, fit to CIRA-86 below 10 km
#'
#' Polynomial of N. Lifton used with the ERA-40 grid in the CRONUS-Earth
#' and CREp calculators. Latitude is signed, so the fit is not symmetric
#' about the equator.
#'
#' @param latitude Latitude (degrees).
#' @return Lapse rate (K/m, positive).
#' @noRd
lapse_rate_cira86 <- function(latitude) {
  coef <- c(-6.1517e-03, -3.1831e-06, -1.5014e-07, 1.8097e-09,
            1.1791e-10, -6.5359e-14, -9.5209e-15)
  -drop(outer(latitude, 0:6, `^`) %*% coef)
}
