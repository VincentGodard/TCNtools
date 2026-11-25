# Geomagnetic history of the Lm scaling scheme, used by geomag_history().
#
# Source: common/production/pmag_consts.mat of the CRONUScalc code (Marrero
# et al. 2016, marrero2016cosmogenic), https://bitbucket.org/cronusearth/cronus-calc
# (commit cdbc8f2, md5 db7e74b494b2ee6f40127ce3622ca72a, MATLAB file created
# 2012-09-14). Its geomagnetic fields are identical to
# dev/data/data_internal/PMag_Sep12.mat at tag v0.1.0. The framework is that
# of Lifton et al. (2014, lifton2014scaling):
#   2010 to 1950   DGRF (finlay2010international)
#   0 to 3 ka      CALS3k.3 (korte2009geomagnetic)
#   3 to 7 ka      CALS7k.2 (korte2005continuous)
#   7 to 18 ka     GLOPIS-75 (laj2004high)
#   18 ka to 2 Ma  PADM2M (ziegler2011padm2m)
#   2 to 10 Ma     constant M/M0 = 0.7048, value of CRONUScalc (origin not
#                  documented: it is not the mean of the record)
# Converted to data-raw/geomagnetic/geomag_poles.csv and geomag_moment.csv
# by data-raw/geomagnetic/geomag_from_mat.py (see that script for details).
#
# Result: a data frame `geomag`, one row per node of the history:
#   time            years before 2010 (t = 0 is 2010 CE)
#   moment_ratio    dipole moment relative to the 2010 DGRF value M0
#                   (7.746e22 A m2)
#   pole_latitude   latitude of the geomagnetic north pole (degrees)
#   pole_longitude  longitude of the geomagnetic north pole (degrees east,
#                   0 to 360)
# Up to 6960 yr the pole is that of the geocentric dipole approximating the
# field models (as used by CRONUScalc for Lm). From 7060 yr the field is a
# geocentric axial dipole: the pole is set to the geographic north pole
# (90, 0), which makes the geomagnetic latitude equal to the geographic one.

geomag <- local({
  poles <- utils::read.csv("data-raw/geomagnetic/geomag_poles.csv")
  moment <- utils::read.csv("data-raw/geomagnetic/geomag_moment.csv")
  stopifnot(
    identical(names(poles), c("time", "moment_ratio", "pole_latitude",
                              "pole_longitude")),
    identical(names(moment), c("time", "moment_ratio")),
    !anyNA(poles), !anyNA(moment),
    poles$time[1] == 0,
    all(diff(poles$time) > 0), all(diff(moment$time) > 0),
    max(poles$time) < min(moment$time),
    all(poles$moment_ratio > 0), all(moment$moment_ratio > 0),
    all(abs(poles$pole_latitude) <= 90),
    all(poles$pole_longitude >= 0 & poles$pole_longitude < 360)
  )
  moment$pole_latitude <- 90
  moment$pole_longitude <- 0
  out <- rbind(poles, moment)
  out$time <- as.numeric(out$time)
  rownames(out) <- NULL
  out
})
