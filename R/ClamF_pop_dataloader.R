#' Function that loads forcings data for Clam population model (alternative version) and performs the interpolation
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats utils
#'

ClamF_pop_dataloader<-function(userpath) {

  # Reads forcing files
  Ttem=read.csv(paste0(userpath,"/ClamF_population/Inputs/Forcings//Water_temperature.csv"),sep=",",header=FALSE)          # Reading the temperature time series (daily series) data
  Chlatem=read.csv(paste0(userpath,"/ClamF_population/Inputs/Forcings//Chlorophyll_a.csv"),sep=",",header=FALSE)           # Reading the Chlorophyll_a (daily series) data


  # Reads integration extremes
  Param_matrix=read.csv(paste0(userpath,"/ClamF_population/Inputs/Parameters//Parameters.csv"),sep=",")                    # Reading the matrix containing parameters and their description

  #Extracts vectors from the forcing files
  Dates=Param_matrix[27:28,3]                           # Vector containing the starting and ending date of teh simulation

  timeT=as.matrix(Ttem[,1])                            # Vector of the times of Temperature measurements
  Temperature=as.double(as.matrix(Ttem[,2]))           # Vector of  water temperature time series (daily series)

  timeChl=as.matrix(Chlatem[,1])                      #Vector of the times of Chlorophyll a measurements
  Chla=as.double(as.matrix(Chlatem[,2]))               # Vector of  Chlorophyll a  time series (daily series)

  # Times needed to perform interpolation
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeChl[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0    # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0    # End of integration [day]

  # Prepare t data for Temperature and Feeding interpolation
  timeTseries=as.numeric(as.Date(timeT, "%d/%m/%Y"))-t0            # Days at which temperature measurements are available
  timeChlaseries=as.numeric(as.Date(timeChl, "%d/%m/%Y"))-t0      # Days at which CHL measurements are available

  # Interpolation of forcings
  Ttem=as.vector(matrix(0,nrow=ti-1))              # Initialize vector Tint
  Chltem=as.vector(matrix(0,nrow=ti-1))            # Initialize vector Chlint

  i=ti:tf                                          # Interpolation base points
  Ttem2=approx(timeTseries,Temperature,xout=i)     # T interpolation according to base points
  Chltem2=approx(timeChlaseries,Chla,xout=i)       # Chl interpolation according to base points

  Tint=c(Ttem, Ttem2$y)            # Interpolated T values starting at t0
  Chlint=c(Chltem, Chltem2$y)      # Interpolated Chl values starting at t0

  # Prepare dates vector
  daysT <- seq(as.Date(timeT[1], format = "%d/%m/%Y"), by = "days", length = length(Tint))
  daysChl <- seq(as.Date(timeChl[1], format = "%d/%m/%Y"), by = "days", length = length(Chlint))

  # Check if forcings are Ok with integration extremes
  if ((ti<(as.numeric(as.Date(timeT[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timeChl[1], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcings are beginning after the specified integration start\n")
    cat("Impossible to proceed with interpolation\n")
  }

  if ((ti>(as.numeric(as.Date(timeT[length(timeT)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timeChl[length(timeChl)], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcing are ending before the specified integration end\n")
    cat("Impossible to proceed with interpolation\n")
  }

  forcings=list(daysT,Tint,daysChl,Chlint)
  return(forcings)
}
