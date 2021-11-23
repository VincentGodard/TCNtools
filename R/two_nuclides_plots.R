#' Compute uncertainty ellipses for two-nuclides plot
#'
#' Two nuclides A and B plotted as B/A on the Y-axis and A on the X-axis
#' Return a list with given confidence level ellipses coordinates for each pair of concentrations
#'
#' @param Ca Concentration of first nuclide (usually 10Be) (at)
#' @param Ca_e 1-sigma uncertainty on the concentration of first nuclide  (at)
#' @param Cb Concentration of second nuclide (usually 26Al) (at)
#' @param Cb_e 1-sigma uncertainty on the concentration of second nuclide  (at)
#' @param confidence confidence level (default 0.98)
#'
#' @return
#' @export
#'
#' @examples
tnp_ellipse <- function(Ca, Ca_e, Cb, Cb_e,confidence=0.98) {
  if (!( (length(Ca)==(length(Ca_e))) & (length(Cb)==(length(Cb_e))) & (length(Ca)==(length(Cb))) )) {stop("All inputs should be of same size")}
  n = 10000
  ll = list()
  for (i in 1:length(Ca)) {
    a = rnorm(n, Ca[i], Ca_e[i])
    b = rnorm(n, Cb[i], Cb_e[i])
    ratio = b / a
    cc = cor(ratio, a)
    corr = matrix(c(1, cc, cc, 1), nrow = 2)
    #
    r_e = sqrt((Cb_e[i] / Cb[i])^2 + (Ca_e[i] / Ca[i])^2) * (Cb[i] / Ca[i])
    ll[[i]] <-
      ellipse::ellipse(
        corr,
        level = confidence,
        scale = c(Ca_e[i], r_e),
        centre = c(Ca[i], Cb[i] / Ca[i])
      )
  }
  return(ll)
}





#' Compute steady state denudation and constant exposure curves for two-nuclides plots
#'
#'
#' @param P1 production and decay parameters for the first nuclide (4 elements vector)
#' p[1] -> unscaled spallation production rate (at/g/a)
#' p[2] -> unscaled stopped muons production rate (at/g/a)
#' p[3] -> unscaled fast muons production rate (at/g/a)
#' p[4] -> decay constant (1/a)
#' @param P2 production and decay parameters for the first nuclide (4 elements vector)
#' @param L Attenuation length (3 elements vector)
#' L[1] -> neutrons
#' L[2] -> stopped muons
#' L[3] -> fast muons
#' @param S scaling factors (2 elements vector)
#' S[1] -> scaling factor for spallation
#' S[2] -> scaling factor for muons
#' @param rho density (g/cm3)
#' @param dlim range of denudation rates (m/Ma, optional, default 0.1-1000 m/Ma), two elements vector
#' @param alim range of exposure ages (a, optional, default 100-10e6 a), two elements vector
#' @param n number of evaluation points (optional default 100)
#'
#' @return a list with two dataframes containing the concentrations for the two nuclides as a function of denudation rate or exposure age
#' @export
#'
#' @examples
tnp_curves <- function(P1,P2,Lambda,S,rho,dlim=c(0.1,1000),alim=c(100,10e6),n=100) {
  ero = 10^seq(log10(dlim[1]),log10(dlim[2]),length.out = n) / 1e6*100*rho # log-spaced denudation rate vector  (m/Ma converted to g/cm2/a)
  C1 = solv_conc_eul(0,ero,Inf,0,P1,S,Lambda)
  C2 = solv_conc_eul(0,ero,Inf,0,P2,S,Lambda)
  ss_den = as.data.frame(cbind(ero,C1,C2))
  colnames(ss_den) <- c("denudation_rate","C1","C2")
#
  time = 10^seq(log10(alim[1]),log10(alim[2]),length.out = n) # log-spaced age vector (a)
  C1 = solv_conc_eul(0,0,time,0,P1,S,Lambda)
  C2 = solv_conc_eul(0,0,time,0,P2,S,Lambda)
  cst_exp = as.data.frame(cbind(time,C1,C2))
  colnames(cst_exp) <- c("exposure age","C1","C2")
#
 ll = list(ss_den,cst_exp)
 names(ll) <- c("ss_denudation","cst_exposure")
  return(ll)
}





