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


#' Compute burial age and denudation rate curves array for two-nuclides plots
#'
#'
#'
#' @param A burial age vector (a)
#' @param E erosion rate vector (m/Ma)
#' @param p1 production and decay parameters for the first nuclide (4 elements vector)
#' \itemize{
#'  \item p1[1] unscaled spallation production rate (at/g/a)
#'  \item p1[2] unscaled stopped muons production rate (at/g/a)
#'  \item p1[3] unscaled fast muons production rate (at/g/a)
#'  \item p1[4] decay constant (1/a)
#' }
#' @param p2 production and decay parameters for the second nuclide (4 elements vector, same as p1)
#' @param L Attenuation length (3 elements vector in g/cm2)
#' \itemize{
#' \item L[1]  neutrons
#' \item L[2]  stopped muons
#' \item L[3]  fast muons
#' }
#' @param S scaling factors (2 elements vector)
#' \itemize{
#' \item S[1]  scaling factor for spallation
#' \item S[2]  scaling factor for muons
#' }
#' @param rho material density (g/cm3)
#' @param n number of along-curve evaluation points (optional default 100)
#'
#' @return a list with two dataframes containing the concentrations for the two nuclides as a function of burial age and denudation
#' @export
#'
#' @examples
#' data("prm") # production and decay data
#' p = prm
#' data("Lambda") # attenuation length data
#' L = Lambda
#' altitude = 1000 # site elevation in m
#' latitude = 20 # site latitude in degrees
#' P = atm_pressure(alt=altitude,model="stone2000") # atmospheric pressure at site
#' S = scaling_st(P,latitude) # Stone 2000 scaling parameters
#' rho = 2.65 # bedrock density (g/cm3)
#' N1 = "Be10" # longer half-life
#'N2 = "Al26" # shorter half-life
#'res = tnp_curves(prm[,N1],prm[,N2],Lambda,S,rho)
#'plot(NA,xlim=c(0.25,10),ylim=c(0.5,7),log="x",
#'     xlab=paste(N1,"(x10^6 at/g)"),ylab=paste(N2,"/",N1))
#'lines(res[[1]]$C1/1e6,res[[1]]$C2/res[[1]]$C1,lty=2,lwd=2,col="khaki4") # constant exposure
#'lines(res[[2]]$C1/1e6,res[[2]]$C2/res[[2]]$C1,lty=1,lwd=2,col="khaki4") # steady-state erosion
#' A = c(1,2,3,4)*1e6 # increments in age (a)
#' E = c(0.1,0.2,0.5,1,2,5,10) # increments in denudation rate (m/Ma)
#' res = tnp_burial(A,E,prm[,N1],prm[,N2],L,S,rho,n=100) # compute array
#' lines(res[[1]]$C1/1e6,res[[1]]$C2/res[[1]]$C1,col="grey") # constant age
#' lines(res[[2]]$C1/1e6,res[[2]]$C2/res[[2]]$C1,col="black") # constant denudation
tnp_burial <-function(A,E,p1,p2,L,S,rho,n=100){
  p1 = as.numeric(p1)
  p2 = as.numeric(p2)
  L = as.numeric(L)
  S = as.numeric(S)
  A = as.numeric(A)
  E = as.numeric(E)
  # courbes iso age
  for (i in 1:length(A)){
    res = data.frame(A=rep(A[i],n),E=10^seq(log10(min(E)),log10(max(E)),length.out = n))
    C1_0 = solv_conc_eul(0,res$E*rho*100/1e6,Inf,0,p1,S,Lambda)
    res$C1 = solv_conc_eul(Inf,0,A[i],C1_0,p1,S,Lambda)
    C2_0 = solv_conc_eul(0,res$E*rho*100/1e6,Inf,0,p2,S,Lambda)
    res$C2 = solv_conc_eul(Inf,0,A[i],C2_0,p2,S,Lambda)
    res[nrow(res)+1,]<-NA # for breaks in lines
    if(i==1){res_A=res}else{res_A=rbind(res_A,res)}
  }
  # courbes iso erosion
  for (i in 1:length(E)){
    res = data.frame(A=seq(0,max(A),length.out = n),E=rep(E[i],n))
    C1_0 = solv_conc_eul(0,E[i]*rho*100/1e6,Inf,0,p1,S,Lambda)
    res$C1 = solv_conc_eul(Inf,0,res$A,C1_0,p1,S,Lambda)
    C2_0 = solv_conc_eul(0,E[i]*rho*100/1e6,Inf,0,p2,S,Lambda)
    res$C2 = solv_conc_eul(Inf,0,res$A,C2_0,p2,S,Lambda)
    res[nrow(res)+1,]<-NA # for breaks in lines
    if(i==1){res_E=res}else{res_E=rbind(res_E,res)}
  }
  res = list(res_A,res_E)
  names(res)<-c("cst_age","cst_denudation")
  return(res)
}


