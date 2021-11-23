


#' Calculates angular distance theta between points on a sphere specified by latitude and longitude
#'
#'
#' @param lat1 latitude of point 1 (degrees)
#' @param lon1 longitude of point 1 (degrees)
#' @param lat2 latitude of point 2 (degrees)
#' @param lon2 longitude of point 2 (degrees)
#' @keywords
#' @examples
#'
angdist<-function(lat1,lon1,lat2,lon2){
  lata = d2r(lat1)
  lona = d2r(lon1)
  latb = d2r(lat2)
  lonb = d2r(lon2)
  #
  theta = acos( (cos(lata)*cos(latb)*((cos(lona)*cos(lonb)) + (sin(lona)*sin(lonb)))) + (sin(lata)*sin(latb)))
  return(r2d(theta))
}



#' Degrees to radians conversion
#'
#'
#' @param a angle (degrees)
#' @keywords
#' @examples
#'
d2r<-function(a){
  return(a/180*pi)
}




#' Radians to degrees conversion
#'
#'
#' @param a angle (radians)
#' @keywords
#' @examples
#'
r2d<-function(a){
  return(a/pi*180)
}









#' Compute atmospheric pressure at site
#'
#'
#'
#' Choice of models :
#'
#' (1) "stone2000" (default) : equation 1 of Stone (2000) JGR paper (only function of site altitude)
#' Stone (2000)  https://doi.org/10.1029/2000JB900181
#'
#' (2) "era40" : Looks up mean sea level pressure and mean 1000 mb temp from ERA-40 reanalysis
#' and calculates site atmospheric pressures using these as inputs to the standard atmosphere equation
#' (function of site latitude, longitude and altitude).
#' Uppala et al (2005)  https://doi.org/10.1256/qj.04.176
#'
#' @param alt Altitude (m)
#' @param lon Longitude (degrees) optional (only for "era40" model)
#' @param lat Latitude (degrees) optional (only for "era40" model)
#' @param model Model used for computation, one of "stone2000" or "era40"
#' @return Atmospheric pressure in hPa
#' @examples
#' z = seq(0,5000,by=100)
#' P1=atm_pressure(alt=z,model="stone2000")
#' P2=atm_pressure(alt=z,lat=-60,lon=0,model="era40")
#' plot(P1,z,type="l",xlab="Pressure (hPa)",ylab="Altitude (m)")
#' lines(P2,z,lty=2)
#' @export
atm_pressure<-function(alt,lon=NULL,lat=NULL,model="stone2000"){
  if (model == "stone2000"){
    Ps=1013.25 # sea level pressure (hPa)
    Ts=288.15 # sea level temperature (K)
    ksi=0.0065 # adiabatic lapse rate (K/m)
    gMR=0.03417 # g*M/R (K/m)
    #hydrostatic pressure (equation 1 of Stone 2000)
    P=Ps*exp(-1*gMR/ksi*(log(Ts) - log(Ts - ksi* alt)))
  }
  else if (model == "era40"){
    if(is.null(lat) | is.null(lon)) {stop("must specify lat and lon for model era40")}
    if(length(lat)>1 | length(lon)>1) {stop("lat and lon parameters must be scalars")}
    z = alt
    # get longitude from 0 to 360
    lon_c=lon+(lon<0)*360
    # Some constants and definitions
    gmr = -0.03417
    lr = c(-6.1517E-03,-3.1831E-06,-1.5014E-07,1.8097E-09,1.1791E-10,-6.5359E-14,-9.5209E-15)
    # interpolate
    # ERA40 packaged in internal data (sysdata.rda)
    site_slp=pracma::interp2(ERA40$ERA40lon,ERA40$ERA40lat,ERA40$meanP, lon_c, lat)
    site_T=pracma::interp2(ERA40$ERA40lon,ERA40$ERA40lat,ERA40$meanT, lon_c, lat)
    # Lifton Lapse Rate Fit to COSPAR CIRA-86 <10 km altitude
    dtdz = lr[1] + lr[2]*lat + lr[3]*lat^2 + lr[4]*lat^3 + lr[5]*lat^4 + lr[6]*lat^5 + lr[7]*lat^6
    dtdz = -dtdz
    P = site_slp * exp( (gmr/dtdz) * ( log(site_T) - log(site_T - (z*dtdz)) ) )
  }
  else {
    stop("model must be stone2000 or era40")
  }


  return(P)
}













#' Get Virtual Dipole Moment time series
#'
#' Original data from Crep calculator (See Martin et al., 2017)
#'
#' @param time time (a)
#' @param model model to be used (one of musch, glopis or lsd)
#'
#' @return
#' @export
#'
#' @examples
#' get_vdm
get_vdm <- function(time,model){
  if (!(model %in%  c("musch","glopis","lsd") )) {stop("Argument model must be one of musch, glopis or lsd")}
  val=0
  if (model == "musch")
  {
    if(max(time)>max(GMDB$musch[1,]*1000)){val = mean(GMDB$musch[2,GMDB$musch[1,]>1800])*1e22}
    vdm = approx(GMDB$musch[1,]*1000,GMDB$musch[2,]*1e22,time,yright = val)$y
  } else if (model == "glopis")
  {
    if(max(time)>max(GMDB$glopis[1,]*1000)){val = mean(GMDB$glopis[2,GMDB$glopis[1,]>1800])*1e22}
    vdm = approx(GMDB$glopis[1,]*1000,GMDB$glopis[2,]*1e22,time,yright = val)$y
  } else
  {
    if(max(time)>max(GMDB$lsd[1,]*1000)){val = mean(GMDB$lsd[2,GMDB$lsd[1,]>1800])*1e22}
    vdm = approx(GMDB$lsd[1,]*1000,GMDB$lsd[2,]*1e22,time,yright = val)$y
  }
  return(vdm)
}




#' calculate cutoff rigidity from Virtual Dipole Moment time series
#'
#'
#'
#' @param vdm Virtual Dipole Moment (A.m^2)
#' @param lat Latitude (deg)
#' @param model Model used to compute Rc from the virtual dipole moment value (one of "elsasser54" or "lifton14")
#'
#' @return
#' @export
#'
#' @examples
#' vdm2rc
vdm2rc <- function(vdm,lat,model="elsasser54"){
  M0 = 7.746e22 # A/m2 reference dipole moment (2010)
  # mu0 = 4*pi*10^-7
  # c = 3.0e8
  # r = 6.3712e6
  if (model == "elsasser54"){
    # Rc = (mu0*c)/(16*pi*10^9*r^2) * vdm * (cos(d2r(lat)))^4
    Rc = 14.31187 * vdm/M0 * (cos(d2r(lat)))^4
  } else if (model == "lifton14") {
    dd = c(6.89901,-103.241,522.061,-1152.15,1189.18,-448.004)
    Rc = vdm/M0 * (dd[1]*cos(d2r(lat)) +
                     dd[2]*(cos(d2r(lat)))^2 +
                     dd[3]*(cos(d2r(lat)))^3 +
                     dd[4]*(cos(d2r(lat)))^4 +
                     dd[5]*(cos(d2r(lat)))^5 +
                     dd[6]*(cos(d2r(lat)))^6)
  } else {
    stop("model must be one of elsasser54 or lifton14")
  }
  #
  return(Rc)
}

# TODO check flipud
#' calculate cutoff rigidity from non-dipolar gridded model
#'
#'
#'
#' @param time time (a)
#' @param lat Latitude (deg)
#' @param lon Longitude (deg)
#'
#' @return
#' @export
#'
#' @examples
#' vdm2rc
rc_ndp <- function(time,lat,lon){
    # get longitude from 0 to 360
    lon_c=lon+(lon<0)*360
    tempRc=rep(0,dim(Pal_LSD$TTRc)[3])
    for(i in 1:dim(Pal_LSD$TTRc)[3]){
      tempRc[i]=pracma::interp2(Pal_LSD$lon_Rc,rev(Pal_LSD$lat_Rc),pracma::flipud(Pal_LSD$TTRc[,,i]), lon_c, lat)
    }
    Rc=approx(Pal_LSD$t_Rc,tempRc,time,rule=1)$y
  return(Rc)
}






#' @export
sapp_launch<-function(app){
  app_list = c("2_nuclides","production_rates","profile_modelling","sato_fluxes",
               "simple_accumulation","transient_erosion")
  if (!(app %in% app_list)) { stop(paste("app must be one of :",paste(app_list,collapse=" ")))}
  appDir <- system.file(paste("shiny_apps/",app,sep=""),package="TCNtools")
  if (appDir == "") {
    stop("Could not find app directory. Try re-installing `TCNtools`.", call. = FALSE)
  }
  shiny::shinyAppDir(appDir,options = list(width = "100%", height = 700))
}

