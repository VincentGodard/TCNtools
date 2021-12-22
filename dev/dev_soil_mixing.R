library("TCNtools")

data("prm")
p = prm
data("Lambda")
L = Lambda


altitude = 1000 # elevation in m
latitude = 20 # latitude in degrees
P = atm_pressure(alt=altitude,model="stone2000") # atmospheric pressure at site
S = scaling_st(P,latitude) # Stone 2000 scaling parameters


rhob = 2.65 # bedrock density (g/cm3)
rhos = rhob/2 # soil density
L = Lambda

# test average production
h = seq(0,300,length.out = 100)
plot(NA,xlim=range(h)/100,ylim=c(2,6),
     xlab="Mixing depth (m)",ylab="Production rate (at/g/a)")
prod_soil = depth_averaged_prod(h*rhos,p[,1],L,S)
lines(h/100,prod_soil,lwd=3)

# test concentration

E = 10^seq(log10(0.1),log10(100),length.out = 100) # denudaton rate (m/Ma)
plot(NA,xlim=range(E),ylim=c(0.02,6),log="xy",
     xlab="Denudation rate (m/Ma)",ylab="Concentration (10^6 at/g)")
h = 100 # soil depth (cm)
C = conc_soil_mixing(h,E,rhos,rhob,p[,1],L,S)
lines(E,C/1e6)
h = 200 # soil depth (cm)
C = conc_soil_mixing(h,E,rhos,rhob,p[,1],L,S)
lines(E,C/1e6,col="red")
h = 1000 # soil depth (cm)
C = conc_soil_mixing(h,E,rhos,rhob,p[,1],L,S)
lines(E,C/1e6,col="green")
legend("topright",c("1 m","2 m","10 m"),lty=1,col=c("black","red","green"),title="Mixing depth",cex=0.5)

# tnp plot
N1 = "Be10" # longer half-life
N2 = "Al26" # shorter half-life
res = tnp_curves(prm[,N1],prm[,N2],Lambda,S,rhob) # compute constant exposure and steady state denudation curves
plot(NA,xlim=c(0.75,5),ylim=c(3,7),log="x",
     xlab=paste(N1,"(at/g)"),ylab=paste(N2,"/",N1))
lines(res[[1]]$C1/1e6,res[[1]]$C2/res[[1]]$C1,lty=2,lwd=2,col="khaki4") # constant exposure
lines(res[[2]]$C1/1e6,res[[2]]$C2/res[[2]]$C1,lty=1,lwd=2,col="khaki4") # steady-state erosion

h = seq(100,1000,by = 100) # increments in soil depth (cm)
E = c(0.1,0.2,0.5,1,2,5,10) # increments in denudation rate (m/Ma)
res = tnp_soil_mixing(h,E,rhos,rhob,prm[,N1],prm[,N2],L,S,n=100) # compute array
lines(res[[1]]$C1/1e6,res[[1]]$C2/res[[1]]$C1,col="grey") # constant depth
lines(res[[2]]$C1/1e6,res[[2]]$C2/res[[2]]$C1,col="black") # constant denudation

# test foster
p = prm[,1]
h=seq(1,200,length.out = 100)
h = h * rhos
E = 100 # m/Ma
E = E*100/1e6 # m/Ma to cm/a
Es = E * rhos #  g/cm2/a : erosion rate soil
Eb = E * rhob #  g/cm2/a : erosion rate bedrock
beta = rhob/rhos
Cb = solv_conc_eul(h,Eb,Inf,0,p,S,L)
Ps = depth_averaged_prod(h,p,L,S)
h = h/rhos
Cs = (Ps*h/E/beta)
num = (1 + p[4]*h/E/beta)
V1 = Cs/num
V2 = Cb/num
plot(h,V1+V2,type="l",ylim=range(V1+V2,0))
lines(h,V1)
lines(h,V2)


# ll = unique(tmp$h)
# for (i in 1:length(ll)){
#  tmp0 = tmp[tmp$h==ll[i],]
#  lines(tmp0$C1/1e6,tmp0$C2/tmp0$C1,col="grey")
# }
# # cst denud
# tmp = res[[2]]
# ll = unique(tmp$E)
# for (i in 1:length(ll)){
#   tmp0 = tmp[tmp$E==ll[i],]
#   lines(tmp0$C1/1e6,tmp0$C2/tmp0$C1)
# }

