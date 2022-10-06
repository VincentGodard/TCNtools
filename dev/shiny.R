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
