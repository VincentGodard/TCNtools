

#' Average nuclide production rate in soil
#'
#' Compute average nuclide production rates for vertically mixed soil column
#'
#' @param h soil thickness (in g/cm2)
#' @param p production and decay parameters for the first nuclide (4 elements vector)
#' \itemize{
#'  \item p[1] unscaled spallation production rate (at/g/a)
#'  \item p[2] unscaled stopped muons production rate (at/g/a)
#'  \item p[3] unscaled fast muons production rate (at/g/a)
#'  \item p[4] decay constant (1/a)
#' }
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
#'
#' @return
#' @export
#'
#' @examples
#' #' #' data("prm") # production and decay data
#' p = prm
#' data("Lambda") # attenuation length data
#' L = Lambda
#' altitude = 1000 # site elevation in m
#' latitude = 20 # site latitude in degrees
#' P = atm_pressure(alt=altitude,model="stone2000") # atmospheric pressure at site
#' S = scaling_st(P,latitude) # Stone 2000 scaling parameters
#' rhob = 2.65 # bedrock density (g/cm3)
#' rhos = rhob/2 # soil density (g/cm3)
#' h = seq(0,300,length.out = 100)
#' plot(NA,xlim=range(h)/100,ylim=c(2,6),
#'      xlab="Mixing depth (m)",ylab="Production rate (at/g/a)")
#' prod_soil = depth_averaged_prod(h*rhos,p[,1],L,S)
#' lines(h/100,prod_soil,lwd=3)
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
#' @param p production and decay parameters for the first nuclide (4 elements vector)
#' \itemize{
#'  \item p[1] unscaled spallation production rate (at/g/a)
#'  \item p[2] unscaled stopped muons production rate (at/g/a)
#'  \item p[3] unscaled fast muons production rate (at/g/a)
#'  \item p[4] decay constant (1/a)
#' }
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
#'@param fqz quartz enrichment factor (default 1)
#'
#' @return Nuclide concentration in soil (at/g)
#'
#'
#' @export
#'
#' @examples
#' #' data("prm") # production and decay data
#' p = prm
#' data("Lambda") # attenuation length data
#' L = Lambda
#' altitude = 1000 # site elevation in m
#' latitude = 20 # site latitude in degrees
#' P = atm_pressure(alt=altitude,model="stone2000") # atmospheric pressure at site
#' S = scaling_st(P,latitude) # Stone 2000 scaling parameters
#' rhob = 2.65 # bedrock density (g/cm3)
#' rhos = rhob/2 # soil density  (g/cm3)
#' E = 10^seq(log10(0.1),log10(100),length.out = 100) # denudaton rate (m/Ma)
#' plot(NA,xlim=range(E),ylim=c(0.02,6),log="xy",
#'     xlab="Denudation rate (m/Ma)",ylab="Concentration (10^6 at/g)")
#' h = 100 # soil depth (cm)
#' C = conc_soil_mixing(h,E,rhos,rhob,p[,1],L,S)
#' lines(E,C/1e6)
#' h = 200 # soil depth (cm)
#' C = conc_soil_mixing(h,E,rhos,rhob,p[,1],L,S)
#' lines(E,C/1e6,col="red")
#' h = 1000 # soil depth (cm)
#' C = conc_soil_mixing(h,E,rhos,rhob,p[,1],L,S)
#' lines(E,C/1e6,col="green")
#' legend("topright",c("1 m","2 m","10 m"),lty=1,col=c("black","red","green"),title="Mixing depth")
conc_soil_mixing <- function(h,E,rhos,rhob,p,L,S,fqz=1){
  if((length(h)>1)&(length(E)>1)){stop("h and E can not be both vectors")}
  p = as.numeric(p)
  L = as.numeric(L)
  S = as.numeric(S)
  h = as.numeric(h)
  E = as.numeric(E)
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
  C = ( (Ps*h/E/beta) + Cb*fqz ) / (1 + p[4]*h/E/beta)

  C[E<0]<-NA
  C[h<=0]<-NA
  return(C)
}

#' Compute constant depth and denudation rate curves array for two-nuclides plots
#'
#'
#'
#' @param h soil thickness vector (cm)
#' @param E erosion rate vector (m/Ma)
#' @param rhos soil density (g/cm3)
#' @param rhob bedrock density (g/cm3)
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
#' @param n number of along-curve evaluation points (optional default 100)
#'
#' @return a list with two dataframes containing the concentrations for the two nuclides as a function of soil depth and denudation
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
#' rhob = 2.65 # bedrock density (g/cm3)
#' rhos = rhob/2 # soil density  (g/cm3)
#' N1 = "Be10" # longer half-life
#'N2 = "Al26" # shorter half-life
#'res = tnp_curves(prm[,N1],prm[,N2],Lambda,S,rhob)
#'plot(NA,xlim=c(0.75,5),ylim=c(3,7),log="x",
#'     xlab=paste(N1,"(x10^6 at/g)"),ylab=paste(N2,"/",N1))
#'lines(res[[1]]$C1/1e6,res[[1]]$C2/res[[1]]$C1,lty=2,lwd=2,col="khaki4") # constant exposure
#'lines(res[[2]]$C1/1e6,res[[2]]$C2/res[[2]]$C1,lty=1,lwd=2,col="khaki4") # steady-state erosion
#' h = seq(100,1000,by = 100) # increments in soil depth (cm)
#' E = c(0.1,0.2,0.5,1,2,5,10) # increments in denudation rate (m/Ma)
#' res = tnp_soil_mixing(h,E,rhos,rhob,prm[,N1],prm[,N2],L,S,n=100) # compute array
#' lines(res[[1]]$C1/1e6,res[[1]]$C2/res[[1]]$C1,col="grey") # constant depth
#' lines(res[[2]]$C1/1e6,res[[2]]$C2/res[[2]]$C1,col="black") # constant denudation
tnp_soil_mixing <-function(h,E,rhos,rhob,p1,p2,L,S,n=100){
  p1 = as.numeric(p1)
  p2 = as.numeric(p2)
  L = as.numeric(L)
  S = as.numeric(S)
  h = as.numeric(h)
  E = as.numeric(E)
  for (i in 1:length(h)){
    res = data.frame(h=rep(h[i],n),E=10^seq(log10(min(E)),log10(max(E)),length.out = n))
    res$C1 = conc_soil_mixing(h[i],res$E,rhos,rhob,p1,L,S)
    res$C2 = conc_soil_mixing(h[i],res$E,rhos,rhob,p2,L,S)
    res[nrow(res)+1,]<-NA # for breaks in lines
    if(i==1){res_h=res}else{res_h=rbind(res_h,res)}
  }
  for (i in 1:length(E)){
    res = data.frame(h=10^seq(log10(0.1),log10(max(h)),length.out = n),E=rep(E[i],n))
    res$C1 = conc_soil_mixing(res$h,E[i],rhos,rhob,p1,L,S)
    res$C2 = conc_soil_mixing(res$h,E[i],rhos,rhob,p2,L,S)
    res[nrow(res)+1,]<-NA # for breaks in lines
    if(i==1){res_E=res}else{res_E=rbind(res_E,res)}
  }
  res = list(res_h,res_E)
  names(res)<-c("cst_depth","cst_denudation")
  return(res)
}

