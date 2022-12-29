#' Function that loads forcings data for Mussel spatialized model and performs the interpolation
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats utils
#'

Mussel_spatial_dataloader<-function(userpath) {

  # Reads point forcing files
  POCtem=read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Point forcings//POC.csv"),sep=",",header=FALSE)                      # Reading the POC (daily series) data
  POCchtem=read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Point forcings//POC_characterization.csv"),sep=",",header=FALSE)   # Reading the POC composition  data
  POMtem=read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Point forcings//POM.csv"),sep=",",header=FALSE)                      # Reading the individual POM time series (daily series) data
  TSMtem=read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Point forcings//TSM.csv"),sep=",",header=FALSE)                      # Reading the TSM (daily series) data

  Spatial_dates=read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//Spatial_dates.csv"),sep=",",header=FALSE)

  # Reads integration extremes
  Param_matrix=read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Parameters//Parameters.csv"),sep=",")                    # Reading the matrix containing parameters and their description

  #Extracts vectors from the forcing files
  Dates=Param_matrix[37:38,3]                          # Vector containing the starting and ending date of teh simulation

  # Initial times for spatial forcings
  timeT=as.matrix(Spatial_dates[2,2:3])                            # Vector of the times of Temperature measurements
  timeChla=as.matrix(Spatial_dates[3,2:3])                      #Vector of the times of Chlorophyll a measurements

  # Initial time for time series
  timePOC=as.matrix(POCtem[,1])                        # Vector of the times of POC measurements
  POC=as.double(as.matrix(POCtem[,2]))                 # Vector of  POC  time series (daily series)

  timePOCch=as.matrix(POCchtem[,1])                    # Vector of  POC characterization (C:N:P)  time series (daily series)
  CPOC=as.double(as.matrix(POCchtem[,2]))              # Vector containing the C/P ratio of seston (weight ratio)
  NPOC=as.double(as.matrix(POCchtem[,3]))              # Vector containing the N/P ratio of seston (weight ratio)

  timePOM=as.matrix(POMtem[,1])                        # Vector of the times of POC measurements
  POM=as.double(as.matrix(POMtem[,2]))                 # Vector of  POC  time series (daily series)

  timeTSS=as.matrix(TSMtem[,1])                        # Vector of the times of POC measurements
  TSS=as.double(as.matrix(TSMtem[,2]))                 # Vector of  POC  time series (daily series)

  # Times needed to perform interpolation
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeChla[1], "%d/%m/%Y")), as.numeric(as.Date(timePOC[1], "%d/%m/%Y")), as.numeric(as.Date(timePOM[1], "%d/%m/%Y")), as.numeric(as.Date(timeTSS[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0    # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0    # End of integration [day]

  # Prepare t data for Temperature and Feeding interpolation
  timePOCseries=as.numeric(as.Date(timePOC, "%d/%m/%Y"))-t0        # Days at which POC measurements are available
  timePOCchseries=as.numeric(as.Date(timePOCch, "%d/%m/%Y"))-t0    # Days at which POC characterization is available
  timePOMseries=as.numeric(as.Date(timePOM, "%d/%m/%Y"))-t0        # Days at which POM measurements are available
  timeTSSseries=as.numeric(as.Date(timeTSS, "%d/%m/%Y"))-t0        # Days at which TSS measurements are available

  # Interpolation of forcings
  POCtem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector POCint
  CPOCtem=as.vector(matrix(0,nrow=ti-1))           # Initialize vector CPOCint
  NPOCtem=as.vector(matrix(0,nrow=ti-1))           # Initialize vector NPOCint
  POMtem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector POMint
  TSStem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector TSSint

  i=ti:tf                                          # Interpolation base points
  POCtem2=approx(timePOCseries,POC,xout=i)         # POC interpolation according to base points
  CPOCtem2=approx(timePOCchseries,CPOC,xout=i)     # CPOC interpolation according to base points
  NPOCtem2=approx(timePOCchseries,NPOC,xout=i)     # NPOC interpolation according to base points
  POMtem2=approx(timePOMseries,POM,xout=i)         # POM interpolation according to base points
  TSStem2=approx(timeTSSseries,TSS,xout=i)         # TSS interpolation according to base points

  POCint=c(POCtem, POCtem2$y)      # Interpolated POC values starting at t0
  CPOCint=c(CPOCtem, CPOCtem2$y)   # Interpolated CPOC values starting at t0
  NPOCint=c(NPOCtem, NPOCtem2$y)   # Interpolated NPOC values starting at t0
  POMint=c(POMtem, POMtem2$y)      # Interpolated POM values starting at t0
  TSSint=c(TSStem, TSStem2$y)      # Interpolated TSS values starting at t0

  # Prepare dates vector
  daysPOC <- seq(as.Date(timePOC[1], format = "%d/%m/%Y"), by = "days", length = length(POCint))
  daysCPOC <- seq(as.Date(timePOCch[1], format = "%d/%m/%Y"), by = "days", length = length(CPOCint))
  daysNPOC <- seq(as.Date(timePOCch[1], format = "%d/%m/%Y"), by = "days", length = length(NPOCint))
  daysPOM <- seq(as.Date(timePOM[1], format = "%d/%m/%Y"), by = "days", length = length(POMint))
  daysTSS <- seq(as.Date(timeTSS[1], format = "%d/%m/%Y"), by = "days", length = length(TSSint))

  # Check if forcings are Ok with integration extremes
  if ((ti<(as.numeric(as.Date(timeT[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timeChla[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timePOC[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timePOCch[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timePOM[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timeTSS[1], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcings are beginning after the specified integration start\n")
    cat("Impossible to proceed with interpolation\n")
  }

  if ((ti>(as.numeric(as.Date(timeT[2], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timeChla[2], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timePOC[length(timePOC)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timePOCch[length(timePOCch)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timePOM[length(timePOM)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timeTSS[length(timeTSS)], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcing are ending before the specified integration end\n")
    cat("Impossible to proceed with interpolation\n")
  }

  # Spatialized forcings upload

  # read .nc files
  sst <- raster::stack(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//sst.nc"))
  chl <- raster::stack(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings/chl_a.nc"))

  # transform rasters to points
  pixel_sst <- t(raster::rasterToPoints(sst))
  sst_export <- pixel_sst[-c(1,2),]
  #colnames(pixel_sst) <- gsub("V", "sst", colnames(pixel_sst))
  write.csv(sst_export,paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//sst.csv"))
  sst <- read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//sst.csv"),header = TRUE)
  sst <- sst[,-(1)]
  colnames(sst) <- gsub("V", "sst", colnames(sst))

  pixel_chla <- t(raster::rasterToPoints(chl))
  chla_export <- pixel_chla[-c(1,2),]
  #colnames(pixel_chla) <- gsub("V", "chl", colnames(pixel_chla))
  write.csv(chla_export,paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//chl.csv"))
  chl <- read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//chl.csv"),header = TRUE)
  chl <- chl[,-(1)]
  colnames(chl) <- gsub("V", "chl", colnames(chl))

  # save coordinates points
  write.csv(pixel_sst[1:2,],paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//coordinates.csv"))
  coord <- read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//coordinates.csv"))

 forcings=list(timeT,sst,timeChla,chl,daysPOC,POCint,daysCPOC,CPOCint,daysNPOC,NPOCint,daysPOM,POMint,daysTSS,TSSint,coord)

  return(forcings)
}
