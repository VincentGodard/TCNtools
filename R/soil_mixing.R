

#' Average nuclide production rate in soil
#'
#' Compute average nuclide production rates for vertically mixed soil column
#'
#' @param h soil thickness (in g/cm2)
#' @param p production and decay parameters (4 elements vector)
#' p[1] -> unscaled spallation production rate (at/g/a)
#' p[2] -> unscaled stopped muons production rate (at/g/a)
#' p[3] -> unscaled fast muons production rate (at/g/a)
#' p[4] -> decay constant (1/a)
#' @param S scaling factors (2 elements vector)
#' S[1] -> scaling factor for spallation
#' S[2] -> scaling factor for muons
#' @param L attenuation lengths (3 elements vector)
#' L[1] -> neutrons
#' L[2] -> stopped muons
#' L[3] -> fast muons
#'
#' @return
#' @export
#'
#' @examples
depth_averaged_prod <- function(h,p,L,S){
  p = as.numeric(p)
  L = as.numeric(L)
  S = as.numeric(S)
    P = p[1]*S[1]*L[1]/h*(1 - exp(-1*h/L[1])) +
        p[2]*S[2]*L[2]/h*(1 - exp(-1*h/L[2])) +
        p[3]*S[2]*L[3]/h*(1 - exp(-1*h/L[3]))
  P[h<=0]<-NA
  return(P)
}


#' Nuclide concentration in a vertically mixed soil
#'
#' Compute nuclide concentration in a vertically-mixed steady-state soil, using the approach of Foster et al. (2015)
#'
#' @param h soil thickness (cm)
#' @param E erosion rate (m/Ma)
#' @param rhos soil density (g/cm3)
#' @param rhob bedrock density (g/cm3)
#' @param p production and decay parameters (4 elements vector)
#' p[1] -> unscaled spallation production rate (at/g/a)
#' p[2] -> unscaled stopped muons production rate (at/g/a)
#' p[3] -> unscaled fast muons production rate (at/g/a)
#' p[4] -> decay constant (1/a)
#' @param S scaling factors (2 elements vector)
#' S[1] -> scaling factor for spallation
#' S[2] -> scaling factor for muons
#' @param L attenuation lengths (3 elements vector)
#' L[1] -> neutrons
#' L[2] -> stopped muons
#' L[3] -> fast muons
#'
#'
#' @return Concentration in at/g
#'
#'
#' @export
#'
#' @examples
conc_soil_mixing <- function(h,E,rhos,rhob,p,L,S){
  p = as.numeric(p)
  L = as.numeric(L)
  S = as.numeric(S)

  # note that input h is in cm and E is in cm/a
  h = h * rhos
  E = E*100/1e6 # m/Ma to cm/a
  Es = E * rhos #  g/cm2/a : erosion rate soil
  Eb = E * rhob #  g/cm2/a : erosion rate bedrock
  beta = rhob/rhos

  # bedrock concentration
  Cb = solv_conc_eul(h,Eb,Inf,0,p,S,L)

  # average soil nuclide production rate
  Ps = depth_averaged_prod(h,p,L,S)

  # Foster et al. (2015)
  h = h/rhos # back to h in cm
  C = ( (Ps*h/E/beta) + Cb ) / (1 + p[4]*h/E/beta)

  return(C)
}


tnp_soil_mixing <-function(h,E,rhos,rhob,p1,p2,L,S,n=100){

  for (i in 1:length(h))
   res = data.frame(h=rep(NA,n),E=seq(min(E),max(E),length.out = n))
   res$C1 = conc_soil_mixing(h[i],res$E*100/1e6,rhos,rhob,p1,L,S)
   res$C2 = conc_soil_mixing(h[i],res$E*100/1e6,rhos,rhob,p1,L,S)
   if(i==1){res_h=res}else{res_h=rbind(res,res_h)}

}

