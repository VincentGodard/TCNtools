# 0 preparations -----
library("TCNtools")
library("pracma")
data(Lambda)
data(prm)
rho = 2.5
#
altitude = 1000 # elevation in m
latitude = 45 # latitude in degrees
P = atm_pressure(alt=altitude,model="stone2000") # compute atmospheric pressure at site
S = scaling_st(P,latitude) # compute the scaling parameters according to Stone (2000)


# 1 steady state test ----
ero = 100 /1e6*100*rho
z = 1000*rho
C1a = solv_conc_eul(z,ero,Inf,0,prm[,'Be10'],S,Lambda)
t = seq(0,1e6,length.out=20)
erob = rep(ero,length(t))
C1b = solv_conc_eul(z,erob,t,0,prm[,'Be10'],S,Lambda)
diff = C1b - C1a
print(diff)

# 2 comparison with lagrangian ----
z = 0*rho
t = seq(0,1e6,length.out=4)
#ero = seq(300,100,length.out=length(t))/1e6*100*rho
ero = runif(length(t),10,20)/1e6*100*rho


C2a = solv_conc_eul(z,ero,t,0,prm[,'Be10'],S,Lambda)

df_l = data.frame(t=seq(min(t),max(t),length.out = 1e5)) # data frame to store results
df_l$ero = approx(t,ero,df_l$t,method="constant",f=1)$y
plot(t,ero,type="S")
points(t,ero,cex=0.5)
lines(df_l$t,df_l$ero,col="red")
points(df_l$t,df_l$ero,col="red",cex=0.5)

df_l$z = pracma::cumtrapz(df_l$t,df_l$ero) # cumulative erosion trough time (g/cm2)
df_l$z = max(df_l$z) - df_l$z + z# convert into depth
plot(df_l$t,df_l$z,type="l")

C10_0 = solv_conc_eul(max(df_l$z),ero[1],Inf,0,prm[,'Be10'],S,Lambda) # starting 10Be concentration
df_l$Ssp  = rep(as.numeric(S[1]),nrow(df_l)) # scaling spallation
df_l$Smu  = rep(as.numeric(S[2]),nrow(df_l)) # scaling muons
# 10Be (not scaled)
df_l$Psp10 = prm["Pspal",'Be10']*exp(-1*df_l$z/Lambda["Lspal"]) # spallation
df_l$Pmu10 = prm["Pstop",'Be10']*exp(-1*df_l$z/Lambda["Lstop"]) + prm["Pfast",'Be10']*exp(-1*df_l$z/Lambda["Lfast"]) # muons
df_l$C2b = solv_conc_lag(df_l$t,df_l$z,C10_0,df_l$Psp10,df_l$Pmu10,prm["lambda",'Be10'],cbind(df_l$Ssp,df_l$Smu),final=FALSE)
df_l$C2c = solv_conc_lag(df_l$t,df_l$z,C10_0,df_l$Psp10,df_l$Pmu10,prm["lambda",'Be10'],cbind(df_l$Ssp,df_l$Smu),final=TRUE)

print(((df_l$C2b[nrow(df_l)] - C2a[length(C2a)]) / C2a[length(C2a)])*100)
print(((df_l$C2c[nrow(df_l)] - C2a[length(C2a)]) / C2a[length(C2a)])*100)



