#' A set of reference parameters for SLHL production rates and decay constant
#'
#' 10Be, 26Al and 14C spallation production rates from Borchers et al. (2016).
#' 10Be and 26Al muons production rates are from Braucher et al. (2011).
#' 14C muons production rate is from Lupker et al. (2015).
#'
#' Borchers, B., Marrero, S., Balco, G., Caffee, M., Goehring, B., Lifton, N., Nishiizumi, K., Phillips, F., Schaefer, J., & Stone, J. (2016). Geological calibration of spallation production rates in the CRONUS-Earth project. Quaternary Geochronology, 31, 188–198. https://doi.org/10.1016/j.quageo.2015.01.009
#'
#' Braucher, R., Merchel, S., Borgomano, J., & Bourlès, D. L. (2011). Production of cosmogenic radionuclides at great depth: A multi element approach. Earth and Planetary Science Letters, 309(1–2), 1–9. http://dx.doi.org/10.1016/j.epsl.2011.06.036
#'
#' Lupker, M., Hippe, K., Wacker, L., Kober, F., Maden, C., Braucher, R., Bourlès, D., Romani, J. R. V., & Wieler, R. (2015). Depth-dependence of the production rate of in situ14C in quartz from the Leymon High core, Spain. Quaternary Geochronology, 28, 80–87. https://doi.org/10.1016/j.quageo.2015.04.004
#'
#' @format A matrix with nuclides in columns and production rates or decay constant in rows
#' \describe{
#'   \item{Pspal}{SLHL spallation production rate for scaling scheme st (at/g/a)}
#'   \item{Pstop}{SL stopping muons production rate (at/g/a)}
#'   \item{Pfast}{SL fast muons production rate (at/g/a)}
#'   \item{lambda}{Decay constant (1/a)}
#' }
#'
"prm"


#' A set of attenuation lengths for cosmogenic nuclides production pathways
#'
#'
#' @format A vector with 3 elements
#' \describe{
#'   \item{Lspal}{Spallogenic neutrons attenuation length (g/cm2)}
#'   \item{Lstop}{Stopping muons attenuation length (g/cm2)}
#'   \item{Lfast}{Fast muons attenuation length (g/cm2)}
#' }
#'
"Lambda"


#' A dataset of depth profiles (10Be)
#'
#' The general characteristics of each sites and studies are only reported for the first item (label=1)
#'
#' Siame, L., Bellier, O., Braucher, R., Sébrier, M., Cushing, M., Bourlès, D., Hamelin, B., Baroux, E., de Voogd, B., Raisbeck, G., & Yiou, F. (2004). Local erosion rates versus active tectonics: cosmic ray exposure modelling in Provence (south-east France). Earth and Planetary Science Letters, 220(3–4), 345–364. http://dx.doi.org/10.1016/S0012-821X(04)00061-5
#'
#' Hidy, A. J., Gosse, J. C., Pederson, J. L., Mattern, J. P., & Finkel, R. C. (2010). A geologically constrained Monte Carlo approach to modeling exposure ages from profiles of cosmogenic nuclides: An example from Lees Ferry, Arizona. Geochemistry Geophysics Geosystems, 11, Q0AA10. https://doi.org/10.1029/2010GC003084
#'
#' Laloy, E., Beerten, K., Vanacker, V., Christl, M., Rogiers, B., & Wouters, L. (2017). Bayesian inversion of a CRN depth profile to infer Quaternary erosion of the northwestern Campine Plateau (NE Belgium). Earth Surface Dynamics, 5(3), 331–345. https://doi.org/10.5194/esurf-5-331-2017
#'
#' @format A dataframe
#' \describe{
#'   \item{study}{Study key word}
#'   \item{sample}{Sample name}
#'   \item{label}{Sample number}
#'   \item{depth}{Depth below surface (cm)}
#'   \item{C}{10Be concentration (at/g)}
#'   \item{C_e}{Uncertainty on 10Be concentration (at/g)}
#'   \item{latitude}{Site latitude (°) only reported for first sample (label=1)}
#'   \item{altitude}{Site altitude (m) only reported for first sample (label=1)}
#' }
#'
"tcn_depth_profiles"

