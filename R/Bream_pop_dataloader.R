#' Function that loads forcings data for Seabream population model and performs the interpolation
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees] and feeding rate [g/individual x d]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats utils
#'

Bream_pop_dataloader<-function(userpath) {

  # Reads forcing files
  Ttem=read.csv(paste0(userpath,"/Bream_population/Inputs/Forcings//Water_temperature.csv"),sep=",",header=FALSE)        # Reading the temperature time series (daily series) data
  DaF=read.csv(paste0(userpath,"/Bream_population/Inputs/Forcings//Feeding.csv"),sep=",",header=FALSE)                   # Reading the individual feeding dose time series (daily series) data

  # Reads integration extremes
  Param_matrix=read.csv(paste0(userpath,"/Bream_population/Inputs/Parameters//Parameters.csv"),sep=",")                    # Reading the matrix containing parameters and their description

  #Extracts vectors from the forcing files
  timeT=as.matrix(Ttem[,1])                     # Vector of the times of Temperature measurements
  Temperature=as.double(as.matrix(Ttem[,2]))    # Vector of  water temperature time series (daily series)
  timeG=as.matrix(DaF[,1])                      # Vector of the times of feeding dose
  G=as.double(as.matrix(DaF[,2]))               # Vector of the individual feeding dose time series (daily series)
  Dates=Param_matrix[22:23,3]                   # Vector containing the starting and ending date of the simulation

  # Times needed to perform interpolation
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeG[1], "%d/%m/%Y")))
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0   # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0   # End of integration [day]

  # Prepare t data for Temperature and Feeding interpolation
  timeTseries=as.numeric(as.Date(timeT, "%d/%m/%Y"))-t0     # Days at which temperature measurements are available
  timeGseries=as.numeric(as.Date(timeG, "%d/%m/%Y"))-t0     # Days at which food measurements are available

  # Interpolation of Temperature and Feeding forcings
  Ttem=as.vector(matrix(0,nrow=ti-1))              # Initialize vector Tint
  Gtem=as.vector(matrix(0,nrow=ti-1))              # Initialize vector Gint
  i=ti:tf+1                                               # Interpolation base points
  Ttem2=approx(timeTseries,Temperature,xout=i)     # T interpolation according to base points
  Gtem2=approx(timeGseries,G,xout=i)               # G interpolation according to base points
  Tint=c(Ttem, Ttem2$y)                            # Interpolated T values starting at t0
  Gint=c(Gtem, Gtem2$y)                            # Interpolated G values starting at t0

  # Prepare dates vector
  daysT <- seq(as.Date(timeT[1], format = "%d/%m/%Y"), by = "days", length = length(Tint))
  daysG <- seq(as.Date(timeG[1], format = "%d/%m/%Y"), by = "days", length = length(Tint))

  # Check if forcings are Ok with integration extremes
  if ((ti<(as.numeric(as.Date(timeT[1], "%d/%m/%Y"))-t0))|(ti<(as.numeric(as.Date(timeG[1], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcings are beginning after the specified integration start\n")
    cat("Impossible to proceed with interpolation\n")
  }

  if ((ti>(as.numeric(as.Date(timeT[length(timeT)], "%d/%m/%Y"))-t0))|(ti>(as.numeric(as.Date(timeG[length(timeG)], "%d/%m/%Y"))-t0))) {
    cat("ERROR: forcing are ending before the specified integration end\n")
    cat("Impossible to proceed with interpolation\n")
  }

  forcings=list(daysT,Tint,daysG,Gint)

  return(forcings)
}
