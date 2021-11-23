library("rmatio")
library("pracma")

# data internally used by the package ----

## ERA40 atmospheric model (Balco et al. 2008) ----
ERA40=read.mat("dev/data/data_internal/ERA40.mat")
ERA40$ERA40lat=rev(ERA40$ERA40lat)
ERA40$meanP=flipud(ERA40$meanP)
ERA40$meanT=flipud(ERA40$meanT)

## VDM from Crep ----
tmp = read.mat("dev/data/data_internal/GMDB.mat")
GMDB<-list(musch=tmp$GMDB$Musch[[1]], glopis=tmp$GMDB$GLOPIS[[1]], lsd=tmp$GMDB$LSD[[1]])

## paleomagnetic records for LSD scaling scheme (Lifton et al 2014) ----
Pal_LSD=read.mat("dev/data/data_internal/PMag_Sep12.mat")
# convert annually averaged Usoskin et al. (2011)
# solar modulation potential to Sato Force Field Potential due to different
# assumed Local Interstellar Spectrum and other factors
Pal_LSD$SPhi = 1.1381076*Pal_LSD$SPhi - 1.2738468e-4*Pal_LSD$SPhi^2
Pal_LSD$SPhiInf = mean(Pal_LSD$SPhi)

## Cross sections from Reedy for use in LSD scaling scheme (Lifton et al 2014) ----
XSectsReedyAll=read.mat("dev/data/data_internal/XSectsReedyAll.mat")
XSectsReedyAll$Natoms3 = 2.006e22
XSectsReedyAll$Natoms10 = 2.006e22
XSectsReedyAll$Natoms14 = 2.006e22
XSectsReedyAll$Natoms26 = 1.003e22

## reference values for LSD scaling scheme (Lifton et al 2014) ----
Ref_LSD=read.mat("dev/data/data_internal/Reference_LSD.mat")

# export (R/sysdata.Rda) ----
usethis::use_data(ERA40,GMDB,Pal_LSD,Ref_LSD,XSectsReedyAll, internal = TRUE, overwrite = TRUE)



# # Paleomagnetic data  (Balco et al 2008)
# PMag_Mar07=read.mat("/home/vincent/Documents/development/R_packages/data/PMag_Mar07.mat")
# # Table 1 from Marrero et al. (2017)
# eltb=read.table("/home/vincent/Documents/development/R_packages/data/marrero_table1.dat",header=TRUE,row.names=1)
# usethis::use_data(ERA40,PMag_Mar07,XSectsReedyAll,Ref_LSD,Pal_LSD,GMDB,eltb, internal = TRUE, overwrite = TRUE)

