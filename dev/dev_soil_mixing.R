library("TCNtools")

altitude = 1000 # elevation in m
latitude = -19.2 # latitude in degrees
P = atm_pressure(alt=altitude,model="stone2000") # compute atmospheric pressure at site
S = scaling_st(P,latitude) #

data("prm")
data("Lambda")

rhos = 1.7
rhob = 2.7
L = Lambda
p = as.numeric(prm[,1])
S = as.numeric(S)

# test average production
h = 100*rhos
P1 = sum((p[1:3]*c(S[1],S[2],S[2])*L)/(h)*(1 - exp(-1*h/L)))
P2 = depth_averaged_prod(h,p,L,S)


# test concentration
#h = 100 # cm
h = seq(10,200,length.out = 100)
E = 10 # m/Ma

C = conc_soil_mixing(h,E,rhos,rhob,p,L,S)

