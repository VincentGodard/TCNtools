# 0 preparations -----
library("TCNtools")
library("pracma")
data(Lambda)
data(prm)
p=prm[,1]
L=Lambda
rho = 2.5
rhos = 1.8
#
altitude = 1000 # elevation in m
latitude = 45 # latitude in degrees
P = atm_pressure(alt=altitude,model="stone2000") # compute atmospheric pressure at site
S = scaling_st(P,latitude) # compute the scaling parameters according to Stone (2000)


time = seq(0,1e6,length.out=10000)
#ero = seq(100,300,length.out=length(time))
ero = rep(2,length(time))
ero[1] = 1
H = seq(100,100,length.out=length(time))
#ero = runif(length(time),10,20)
#H = runif(length(time),50,200)
#

beta = rho/rhos
nt =  length(time)
Ps = depth_averaged_prod(H*rhos,p,L,S) # depth averaged production rate
if(length(Ps)==1){Ps=rep(Ps,nt)}
Cb = solv_conc_eul(H*rhos,ero*100/1e6*rho,time,0,p,S,L) # evolution of concentration at the soil bedrock interface
Cs0 = conc_soil_mixing(H[1],ero[1],rhos,rho,p,L,S) # initial soil concentration

Cs = rep(NA,nt)
Cs[1] = Cs0
ero = ero*100/1e6 # m/Ma -> cm/a
for (i in 2:nt){
  dt = time[i] - time[i-1]
  Cs[i] = Ps[i]*dt +
    (1-(p[4]+(ero[i]*beta/H[i]))*dt)*Cs[i-1] +
    (ero[i]*beta/H[i])*dt*Cb[i]
}

Csf = pracma::trapz(time, (Ps + ero*beta/H*Cb)*exp(-1*(p[4]+ero*beta/H)*abs(time[nt]-time) ) ) +
  Cs0*exp(-1*(p[4]+ero*beta/H)*(max(time)-min(time) ) )


ero = ero/100*1e6 #  cm/a -> m/Ma
Cs1 = conc_soil_mixing(H[1],ero[1],rhos,rho,p,L,S)
Cs2 = conc_soil_mixing(H[length(H)],ero[length(ero)],rhos,rho,p,L,S)
Cb2 = solv_conc_eul(0,ero[length(ero)]*100/1e6*rho,Inf,0,p,S,L)
plot(time,Cs,type="l",ylim=range(Cs,Cs1,Cs2,Cb2))

abline(h=Cs1,lty=2)
abline(h=Cs2,lty=2)
abline(h=Cb2,lty=2,col="red")


abline(h=2666577)
