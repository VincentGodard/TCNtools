


#' Monte Carlo exploration for a depth profile
#'
#' Eulerian approach
#'
#' @param C concentrations vector (at/g)
#' @param C_e uncertainties on the concentrations vector (at/g)
#' @param z depth vector (g/cm2)
#' @param prm production and decay parameters (4 elements vector)
#' p[1] -> unscaled spallation production rate (at/g/a)
#' p[2] -> unscaled stopped muons production rate (at/g/a)
#' p[3] -> unscaled fast muons production rate (at/g/a)
#' p[4] -> decay constant (1/a)
#' @param L attenuation lengths (3 elements vector)
#' L[1] -> neutrons
#' L[2] -> stopped muons
#' L[3] -> fast muons
#' @param S scaling factors (2 elements vector)
#' S[1] -> scaling factor for spallation
#' S[2] -> scaling factor for muons
#' @param age age range to explore (a). 2 elements vector. If 1 element this parameter is fixed
#' @param ero erosion rate range to explore (m/Ma). 2 elements vector. If 1 element this parameter is fixed
#' @param inh inherited concentration range to explore (at/g). 2 elements vector. If 1 element this parameter is fixed
#' @param rho density range to explore (g/cm3). 2 elements vector. If 1 element this parameter is fixed
#' @param n1 number of runs with random
#' @param n2
#'
#' @export
#' @examples
depth_profile_mc <- function(C,C_e,depth,prm,L,S,age,ero,inh,rho,n1=0,n2=0){
  if(length(C)!=length(depth)){stop("Vectors C and depth must have same length")}
  if(length(C)!=length(C_e)){stop("Vectors C and C_e must have same length")}
  if( !length(age)%in%c(1,2) |  !length(ero)%in%c(1,2) |
      !length(inh)%in%c(1,2) |  !length(rho)%in%c(1,2)){stop("age, ero, inh and rho must be of length 1 or 2")}
  if(n1==0 & n2==0){stop("at least one of n1 and n2 must be >0")}
  #
  # expand mc
  if(n1>0){
    if(length(age) == 1){age0=rep(age,n1)}else{age0=runif(n1,age[1],age[2])}
    if(length(ero) == 1){ero0=rep(ero,n1)}else{ero0=runif(n1,ero[1],ero[2])}
    if(length(inh) == 1){inh0=rep(inh,n1)}else{inh0=runif(n1,inh[1],inh[2])}
    if(length(rho) == 1){rho0=rep(rho,n1)}else{rho0=runif(n1,rho[1],rho[2])}
    res1 = data.frame(age0,ero0,inh0,rho0)
    colnames(res1)<-c("age","ero","inh","rho")} else {
      res1 = data.frame(ero=double(),age=double(),inh=double(),rho=double())}
  # expand grid
  if(n2>0){
    if(length(age) == 1){age0=age}else{age0=seq(age[1],age[2],length.out=n2)}
    if(length(ero) == 1){ero0=ero}else{ero0=seq(ero[1],ero[2],length.out=n2)}
    if(length(inh) == 1){inh0=inh}else{inh0=seq(inh[1],inh[2],length.out=n2)}
    if(length(rho) == 1){rho0=rho}else{rho0=seq(rho[1],rho[2],length.out=n2)}
    res2 = expand.grid(age0,ero0,inh0,rho0)
    colnames(res2)<-c("age","ero","inh","rho")} else {
      res2 = data.frame(ero=double(),age=double(),inh=double(),rho=double())}
  # build data frame
  res = rbind(res1,res2)
  res$chi2 = NA
  #
  for (i in 1:nrow(res)){
    Cmod = solv_conc_eul(depth*res$rho[i],res$ero[i]*100/1e6*res$rho[i],res$age[i],res$inh[i],prm,S,Lambda)
    res$chi2[i] = sum((C-Cmod)^2/C_e^2)
  }
  #
  return(res)
}

