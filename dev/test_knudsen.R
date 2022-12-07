library("TCNtools")
library("pracma")
data(Lambda)
data(prm)
rho = 2.5
prm[1,1] = 4.5
#
altitude = 1000 # elevation in m
latitude = 45 # latitude in degrees
P = atm_pressure(alt=altitude,model="stone2000") # compute atmospheric pressure at site
S = scaling_st(P,latitude) # compute the scaling parameters according to Stone (2000)
S=c(1,0)
tmax = 1e6# duration in a
time = data.frame(t=seq(0,tmax,length.out=100))
time$ero = seq(1,2,length.out=nrow(time))/1.e6*100*rho
#time$ero = 1/1.e6*100*rho
time$z = cumtrapz(time$t,time$ero) # cumulative erosion trough time (g/cm2)
time$z = max(time$z) - time$z # convert into depth
plot(time$t,time$z/rho)

lambda = as.numeric(prm["lambda",'Be10'])
Pspal = as.numeric(prm["Pspal",'Be10'] * S[1])
Lspal = as.numeric(Lambda["Lspal"])

# Lag
C10_0 = solv_conc_eul(max(time$z),time$ero[1],Inf,0,prm[,'Be10'],c(S[1],0),Lambda) # starting 10Be concentration
test = prm[1,'Be10']*exp(-1*max(time$z)/Lspal)/(lambda+time$ero[1]/Lspal)
time$Ssp  = rep(as.numeric(S[1]),nrow(time)) # scaling spallation
#time$Smu  = rep(as.numeric(S[2]),nrow(time)) # scaling muons
time$Smu  = rep(0,nrow(time)) # scaling muons
# 10Be (not scaled)
time$Psp10 = prm["Pspal",'Be10']*exp(-1*time$z/Lambda["Lspal"]) # spallation
#time$Pmu10 = prm["Pstop",'Be10']*exp(-1*time$z/Lambda["Lstop"]) + prm["Pfast",'Be10']*exp(-1*time$z/Lambda["Lfast"]) # muons
time$Pmu10 =0
time$C10 = solv_conc_lag(time$t,time$z,C10_0,time$Psp10,time$Pmu10,prm["lambda",'Be10'],cbind(time$Ssp,time$Smu),final=FALSE)

# Lag 2
time$C10_3 = NA
time$C10_3[1] = C10_0
for(i in 2:nrow(time)){
  dt =  time$t[i] - time$t[i-1]
  ero = (time$ero[i] + time$ero[i-1])/2
  time$C10_3[i] = time$C10_3[i-1]*exp(-1*lambda*dt) + time$Psp10[i]*time$Ssp[i]*dt

}




time$C10_2 = NA
time$C10_2[1] = solv_conc_eul(0,time$ero[1],Inf,0,prm[,'Be10'],c(S[1],0),Lambda) # starting 10Be concentration
for(i in 2:nrow(time)){
  dt =  time$t[i] - time$t[i-1]
  ero = (time$ero[i] + time$ero[i-1])/2
  time$C10_2[i] = time$C10_2[i-1]*exp(-1*(lambda+ero/Lspal)*dt) +
                  Pspal/(lambda+ero/Lspal)*(1-exp(-1*(lambda+ero/Lspal)*dt))
}


plot(time$t,time$C10_2,type="l",ylim=range(time$C10_2,time$C10))
lines(time$t,time$C10,col="red")
points(time$t,time$C10_3,col="red")
