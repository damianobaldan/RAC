#' Function that loads forcings data for Bass spatialized model and performs the interpolation
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees] and feeding rate [g/individual x d]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats utils
#'

Bass_spatial_dataloader<-function(userpath) {

  # Reads point forcing files
  DaF=read.csv(paste0(userpath,"/Bass_spatial/Inputs/Point forcings//Feeding.csv"),sep=",",header=FALSE)

  # Reads starting and ending dates for spatialized forcings
  Spatial_dates=read.csv(paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings//Spatial_dates.csv"),sep=",",header=FALSE)

  # Reads integration extremes
  Param_matrix=read.csv(paste0(userpath,"/Bass_spatial/Inputs/Parameters//Parameters.csv"),sep=",")                    # Reading the matrix containing parameters and their description

  #Extracts vectors from the forcing files
  Dates=Param_matrix[22:23,3]                                # Vector containing the starting and ending date of the simulation

  timeT=as.matrix(Spatial_dates[2,2:3])                            # Vector of the times of Temperature measurements
  timeG=as.matrix(DaF[,1])                      # Vector of the times of feeding dose

  G=as.double(as.matrix(DaF[,2]))               # Vector of the individual feeding dose time series (daily series)

  # Times needed to perform interpolation
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeG[1], "%d/%m/%Y")))
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0   # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0   # End of integration [day]

  # Prepare t data for Temperature and Feeding interpolation
  timeGseries=as.numeric(as.Date(timeG, "%d/%m/%Y"))-t0     # Days at which food measurements are available

  # Interpolation of Feeding forcing
  Gtem=as.vector(matrix(0,nrow=ti-1))              # Initialize vector Gint
  i=ti:tf+1                                              # Interpolation base points
  Gtem2=approx(timeGseries,G,xout=i)               # G interpolation according to base points
  Gint=c(Gtem, Gtem2$y)                            # Interpolated G values starting at t0

  daysG <- seq(as.Date(timeG[1], format = "%d/%m/%Y"), by = "days", length = length(Gint))

  # Check if forcings are Ok with integration extremes
  if ((ti<(as.numeric(as.Date(timeT[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timeG[1], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcings are beginning after the specified integration start\n")
    cat("Impossible to proceed with interpolation\n")
  }

  if ((ti>(as.numeric(as.Date(timeT[length(timeT)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timeG[length(timeG)], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcing are ending before the specified integration end\n")
    cat("Impossible to proceed with interpolation\n")
  }


  # Spatialized forcings upload

  # read .nc files
  sst <- raster::stack(paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings//sst.nc"))

  # transform rasters to points
  pixel_sst <- t(rasterToPoints(sst))
  sst_export <- pixel_sst[-c(1,2),]
  #colnames(pixel_sst) <- gsub("V", "sst", colnames(pixel_sst))
  write.csv(sst_export,paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings//sst.csv"))
  sst <- read.csv(paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings//sst.csv"),head=T)
  sst <- sst[,-(1)]
  colnames(sst) <- gsub("V", "sst", colnames(sst))

  # save coordinates points
  write.csv(pixel_sst[1:2,],paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings//coordinates.csv"))
  coord <- read.csv(paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings//coordinates.csv"))

  forcings=list(timeT,sst,daysG,Gint,coord)

  return(forcings)
}
