#' Function that loads forcings data for Clam population model and performs the interpolation
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l], particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats utils
#'

Clam_pop_dataloader<-function(userpath) {

  # Reads forcing files
  Ttem=read.csv(paste0(userpath,"/Clam_population/Inputs/Forcings//Water_temperature.csv"),sep=",",header=FALSE)          # Reading the temperature time series (daily series) data
  Chlatem=read.csv(paste0(userpath,"/Clam_population/Inputs/Forcings//Chlorophyll_a.csv"),sep=",",header=FALSE)           # Reading the Chlorophyll_a (daily series) data
  POCtem=read.csv(paste0(userpath,"/Clam_population/Inputs/Forcings//POC.csv"),sep=",",header=FALSE)                      # Reading the POC (daily series) data
  POMtem=read.csv(paste0(userpath,"/Clam_population/Inputs/Forcings//POM.csv"),sep=",",header=FALSE)                      # Reading the population POM time series (daily series) data
  TSMtem=read.csv(paste0(userpath,"/Clam_population/Inputs/Forcings//TSM.csv"),sep=",",header=FALSE)                      # Reading the TSM (daily series) data

  # Reads integration extremes
  Param_matrix=read.csv(paste0(userpath,"/Clam_population/Inputs/Parameters//Parameters.csv"),sep=",")                    # Reading the matrix containing parameters and their description

  #Extracts vectors from the forcing files
  Dates=Param_matrix[24:25,3]                         # Vector containing the starting and ending date of teh simulation

  timeT=as.matrix(Ttem[,1])                            # Vector of the times of Temperature measurements
  Temperature=as.double(as.matrix(Ttem[,2]))           # Vector of  water temperature time series (daily series)

  timeChla=as.matrix(Chlatem[,1])                      #Vector of the times of Chlorophyll a measurements
  Chla=as.double(as.matrix(Chlatem[,2]))               # Vector of  Chlorophyll a  time series (daily series)

  timePOC=as.matrix(POCtem[,1])                        # Vector of the times of POC measurements
  POC=as.double(as.matrix(POCtem[,2]))                 # Vector of  POC  time series (daily series)

  timePOM=as.matrix(POMtem[,1])                        # Vector of the times of POC measurements
  POM=as.double(as.matrix(POMtem[,2]))                 # Vector of  POC  time series (daily series)

  timeTSS=as.matrix(TSMtem[,1])                        # Vector of the times of POC measurements
  TSS=as.double(as.matrix(TSMtem[,2]))                 # Vector of  POC  time series (daily series)

  # Times needed to perform interpolation
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeChla[1], "%d/%m/%Y")), as.numeric(as.Date(timePOC[1], "%d/%m/%Y")), as.numeric(as.Date(timePOM[1], "%d/%m/%Y")), as.numeric(as.Date(timeTSS[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0    # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0    # End of integration [day]

  # Prepare t data for Temperature and Feeding interpolation
  timeTseries=as.numeric(as.Date(timeT, "%d/%m/%Y"))-t0            # Days at which temperature measurements are available
  timeChlaseries=as.numeric(as.Date(timeChla, "%d/%m/%Y"))-t0      # Days at which CHL measurements are available
  timePOCseries=as.numeric(as.Date(timePOC, "%d/%m/%Y"))-t0        # Days at which POC measurements are available
  timePOMseries=as.numeric(as.Date(timePOM, "%d/%m/%Y"))-t0        # Days at which POM measurements are available
  timeTSSseries=as.numeric(as.Date(timeTSS, "%d/%m/%Y"))-t0        # Days at which TSS measurements are available

  # Interpolation of forcings
  Ttem=as.vector(matrix(0,nrow=ti-1))              # Initialize vector Tint
  Chltem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector Chlint
  POCtem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector POCint
  POMtem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector POMint
  TSStem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector TSSint

  i=ti:tf                                          # Interpolation base points
  Ttem2=approx(timeTseries,Temperature,xout=i)     # T interpolation according to base points
  Chltem2=approx(timeChlaseries,Chla,xout=i)       # Chl interpolation according to base points
  POCtem2=approx(timePOCseries,POC,xout=i)         # POC interpolation according to base points
  POMtem2=approx(timePOMseries,POM,xout=i)         # POM interpolation according to base points
  TSStem2=approx(timeTSSseries,TSS,xout=i)         # TSS interpolation according to base points

  Tint=c(Ttem, Ttem2$y)            # Interpolated T values starting at t0
  Chlint=c(Chltem, Chltem2$y)      # Interpolated Chl values starting at t0
  POCint=c(POCtem, POCtem2$y)      # Interpolated POC values starting at t0
  POMint=c(POMtem, POMtem2$y)      # Interpolated POM values starting at t0
  TSSint=c(TSStem, TSStem2$y)      # Interpolated TSS values starting at t0

  # Prepare dates vector
  daysT <- seq(as.Date(timeT[1], format = "%d/%m/%Y"), by = "days", length = length(Tint))
  daysChl <- seq(as.Date(timeChla[1], format = "%d/%m/%Y"), by = "days", length = length(Chlint))
  daysPOC <- seq(as.Date(timePOC[1], format = "%d/%m/%Y"), by = "days", length = length(POCint))
  daysPOM <- seq(as.Date(timePOM[1], format = "%d/%m/%Y"), by = "days", length = length(POMint))
  daysTSS <- seq(as.Date(timeTSS[1], format = "%d/%m/%Y"), by = "days", length = length(TSSint))

  # Check if forcings are Ok with integration extremes
  if ((ti<(as.numeric(as.Date(timeT[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timeChla[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timePOC[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timePOM[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timeTSS[1], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcings are beginning after the specified integration start\n")
    cat("Impossible to proceed with interpolation\n")
  }

  if ((ti>(as.numeric(as.Date(timeT[length(timeT)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timeChla[length(timeChla)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timePOC[length(timePOC)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timePOM[length(timePOM)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timeTSS[length(timeTSS)], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcing are ending before the specified integration end\n")
    cat("Impossible to proceed with interpolation\n")
  }

  forcings=list(daysT,Tint,daysChl,Chlint,daysPOC,POCint,daysPOM,POMint,daysTSS,TSSint)
  return(forcings)
}
