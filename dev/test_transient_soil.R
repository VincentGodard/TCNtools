# todo  : how to handle H=0
# come back to foster

# 0 preparations -----
# libraries
library("TCNtools")
library("pracma")
library("tictoc")
# cosmos
data(Lambda)
data(prm)
p=prm[,1]
L=Lambda
altitude = 1000 # elevation in m
latitude = 45 # latitude in degrees
P = atm_pressure(alt=altitude,model="stone2000") # compute atmospheric pressure at site
S = scaling_st(P,latitude) # compute the scaling parameters according to Stone (2000)
# soil
rho = 2.5
rhos = 1.8
Sl = 0.1 # m/m - basin slope
transition = 0.5e6
lambda = 1e5 # transition time scale
K1 = 1 # initial
K2 = 20  # final
W0 = 1 # m/Ma - bare bedrock soil production rate
h = 1 # m - evolution lengsscale
# function for steady state soil mixing
# H in cm
# E in cm/a
conc_soil_mixing2 <- function(H,W,rho,rhos,p,L,S){
  beta = rho/rhos
  Ps = depth_averaged_prod(H*rhos,p,L,S) # depth averaged production rate
  Cb = solv_conc_eul(H*rhos,W*rho,Inf,0,p,S,L)
  Cs = (Ps+beta*W/H*Cb)/(p[4]+W/H)
  return(Cs)
}
# function transient soil mixing
# H in cm
# W and E in cm/a
transient_soil2 <- function(t,H,W,E,rho,rhos,p,S,L,final=F){
  beta = rho/rhos
  nt =  length(t)
  if(length(H)==1){H=rep(H,nt)}
  Ps = depth_averaged_prod(H*rhos,p,L,S) # depth averaged production rate
  if(length(Ps)==1){Ps=rep(Ps,nt)}
  Cb = solv_conc_eul(H*rhos,W*rho,t,0,p,S,L) # evolution of concentration at the soil bedrock interface
  Cs0 = conc_soil_mixing2(H[1],W[1],rho,rhos,p,L,S) # initial soil concentration (steady state) !different from foster model
  if (!final){
    #
    Cs = rep(NA,nt)
    Cs[1] = Cs0
    for (i in 2:nt){
      dt = t[i] - t[i-1]
      Cs[i] = Cs[i-1] + Ps[i]*dt + beta*W[i]/H[i]*Cb[i]*dt - p[4]*Cs[i-1]*dt - E[i]/H[i]*Cs[i-1]*dt
      Cs[H<=0] <- Cb[H<=0]
    }
    return(Cs)
  }else{
    if(H[nt]<=0){Csf=Cb[nt]}else{
      tt = max(t)-min(t)
      a = rev(pracma::cumtrapz(t,rev(E/H)))
      Csf = pracma::trapz(t, (Ps + beta*W/H*Cb)*exp(-1*p[4]*abs(t[nt]-t))*exp(-a) ) +
        Cs0*exp(-p[4]*tt)*exp(-pracma::trapz(t,E/H))
    }
    return(Csf)
  }
}
# # check granger -----
# ero = 10^seq(log10(0.1),log10(1000),length.out=100)*100/1e6*rho # m/Ma -> g/cm2/a
# C = solv_conc_eul(0,ero,Inf,0,p,S,L)
# Cs1 = conc_soil_mixing(100,ero/100*1e6/rho,rhos,rho,p,L,S) # initial soil concentration (steady state)
# Cs2 = conc_soil_mixing2(100,ero/rho,rho,rhos,p,L,S) # initial soil concentration (steady state) different from foster model
#
# plot(ero,C,type="l",log="xy")
# lines(ero,Cs1,col="red")
# lines(ero,Cs2,col="green")

# scenario 1 -----
data = data.frame(t = seq(0,1e6,length.out=1000))
data$K = 1/(1+exp(-(data$t-transition)/lambda))*(K2 - K1) + K1
#data$K = K1
data$E = data$K*Sl
plot(data$t,data$K*Sl,type="l")
data$H = -h*log(data$K*Sl/W0)
data$H[data$H<=0] <- 0
plot(data$t,data$H,type="l")
data$W = W0*exp(-1*data$H/h)
plot(data$t,data$W,type="l")
#
t = data$t
H = data$H*100 # m -> cm
W = data$W*100/1e6 # m/Ma -> cm/a
E = data$E*100/1e6 # m/Ma -> cm/a
tic()
Cs = transient_soil2(t,H,W,E,rho,rhos,p,S,L,final=F)
toc()
plot(t,Cs,type="l")
tic()
Csf = transient_soil2(t,H,W,E,rho,rhos,p,S,L,final=T)
toc()
points(t[length(t)],Csf)

# scenario 2 (constant soil thickness) -----
data = data.frame(t = seq(0,1e6,length.out=1000))
data$K = 1/(1+exp(-(data$t-transition)/lambda))*(K2 - K1) + K1
data$E = data$K*Sl
data$H = 1
data$W = data$E
#
t = data$t
H = data$H*100 # m -> cm
W = data$W*100/1e6 # m/Ma -> cm/a
E = data$E*100/1e6 # m/Ma -> cm/a
tic()
Cs = transient_soil2(t,H,W,E,rho,rhos,p,S,L,final=F)
toc()
plot(t,Cs,type="l")
tic()
Csf = transient_soil2(t,H,W,E,rho,rhos,p,S,L,final=T)
toc()
points(t[length(t)],Csf)



# # concentrations
# time = data$t
# H = data$H*100 # m -> cm
# W = data$W*100/1e6 # m/Ma -> cm/a
# E = data$E*100/1e6 # m/Ma -> cm/a
# #
# beta = rho/rhos
# nt =  length(time)
# Ps = depth_averaged_prod(H*rhos,p,L,S) # depth averaged production rate
# if(length(Ps)==1){Ps=rep(Ps,nt)}
# Cb = solv_conc_eul(H*rhos,W*rho,time,0,p,S,L) # evolution of concentration at the soil bedrock interface
# plot(time,Cb,type="l")
# #Cs0 = conc_soil_mixing(H[1],data$E[1],rhos,rho,p,L,S) # initial soil concentration (steady state)
# Cs0 = conc_soil_mixing2(H[1],W[1],rho,rhos,p,L,S) # initial soil concentration (steady state) different from foster model
#
# #
# Cs = rep(NA,nt)
# Cs[1] = Cs0
# tic()
# for (i in 2:nt){
#   dt = time[i] - time[i-1]
#   Cs[i] = Cs[i-1] + Ps[i]*dt + beta*W[i]/H[i]*Cb[i]*dt - p[4]*Cs[i-1]*dt - E[i]/H[i]*Cs[i-1]*dt
# }
# toc()
# plot(time,Cs,type="l")
#
# tic()
#   tt = max(time)-min(time)
#   a = rev(pracma::cumtrapz(time,rev(E/H)))
#   Csf = pracma::trapz(time, (Ps + beta*W/H*Cb)*exp(-1*p[4]*abs(time[nt]-time))*exp(-a) ) +
#     Cs0*exp(-p[4]*tt)*exp(-trapz(time,E/H))
# toc()
# abline(h=Csf)






#
#
#  Csf = pracma::trapz(time, (Ps + W*beta/H*Cb)*exp(-1*(p[4]+E/H)*abs(time[nt]-time) ) )
#
# Cs00 =     Cs0*exp(-p[4]*tt)

#abline(h=Csf)
#abline(h=Csf+Cs00)
#*pracma::trapz(time,exp(-E/H*abs(time[nt]-time)))/tt

#       pracma::trapz(time, Cs0*exp(-1*(p[4]+E/H)*abs(time[nt]-time)))/(max(time)-min(time) )


#       Cs0*exp(-1*(p[4]+pracma::trapz(time,E/H))*(max(time)-min(time) ) )
#points(max(time),Csf)
#
#
# ero = ero/100*1e6 #  cm/a -> m/Ma
# Cs1 = conc_soil_mixing(H[1],ero[1],rhos,rho,p,L,S)
# Cs2 = conc_soil_mixing(H[length(H)],ero[length(ero)],rhos,rho,p,L,S)
# Cb2 = solv_conc_eul(0,ero[length(ero)]*100/1e6*rho,Inf,0,p,S,L)
# plot(time,Cs,type="l",ylim=range(Cs,Cs1,Cs2,Cb2))
#
# abline(h=Cs1,lty=2)
# abline(h=Cs2,lty=2)
# abline(h=Cb2,lty=2,col="red")
#
#
# abline(h=2666577)

# # soil evol
# S = 0.1 # m/m
# transition = 0.5e6
# lambda = 1e3
# K1 = 0.5
# K2 = 5
# W0 = 1 # m/Ma - bare bedrock soil production rate
# h = 0.5 # m
# dt = 10
#
# data = data.frame(t=seq(0,1e6,by=dt))
# data$K = 1/(1+exp(-(data$t-transition)/lambda))*(K2 - K1) + K1
# #plot(data$t,data$K,type="l")
# #data$H = NA
# #data$H[1] = -h*log(data$K[1]*S/W0)
# # for (i in 2:nrow(data)){
# #  data$H[i] = data$H[i-1] + W0*exp(-data$H[i-1]/h)*dt -  data$K[i]*S*dt
# # }
# data$H2 = -h*log(data$K*S/W0)
# plot(data$t,data$H,col="grey",lwd=3,type="l")
# lines(data$t,data$H2,col="red")
#
# # library(pracma)
# # fun <-function(t,H){
# #   S = 0.1 # m/m
# #   transition = 0.5e6
# #   lambda = 1e3
# #   K1 = 0.5
# #   K2 = 5
# #   W0 = 1 # m/Ma - bare bedrock soil production rate
# #   h = 0.5 # m
# #   K = 1/(1+exp(-(t-transition)/lambda))*(K2 - K1) + K1
# #   dHdt = W0*exp(-H/h) - K*S
# #   return(dHdt)
# # }
# # sol <- pracma::ode45(fun, 0, 1e6,-h*log(data$K[1]*S/W0))  # ode23, ode23s
# #
# # plot(data$t,data$H2,col="grey",lwd=3,type="l")
# # lines(sol$t,sol$y,col="red")


