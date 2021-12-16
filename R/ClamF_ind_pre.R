#' Clam bioenergetic individual model preprocessor (alternativer version)
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @param forcings a list containing model forcings
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

ClamF_ind_pre<-function(userpath,forcings){

rm(list=ls())       # Clean workspace

  # Extracts forcings values from the list
  timeT=forcings[[1]]
  Tint=forcings[[2]]
  timeChl=forcings[[3]]
  Chlint=forcings[[4]]

cat("Data preprocessing")

# Read forcings and parameters from .csv files
Param_matrix=read.csv(paste0(userpath,"/ClamF_individual/Inputs/Parameters//Parameters.csv"),sep=",")                      # Reading the matrix containing parameters and their description

# Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
Param=as.double(as.matrix(Param_matrix[1:25,3]))     # Vector containing all parameters
Dates=Param_matrix[27:28,3]                          # Vector containing the starting and ending date of teh simulation
IC=as.double(as.matrix(Param_matrix[26,3]))          # Initial weight condition
CS=as.double(as.matrix(Param_matrix[29,3]))                # Commercial size

# Prepare data for ODE solution
t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeChl[1], "%d/%m/%Y")),  as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
timestep=1                                           # Time step for integration [day]
ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0      # Start of integration [day]
tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0      # End of integration [day]
Ww=as.vector(matrix(0,nrow=ti))                      # Initialize vector weight
Ww[ti]=IC                                            # Wet weight initial value [g]
times<-cbind(ti, tf, timestep)                       # Vector with integration data

# Initial conditions
a=Param[8]               # [cm] Wet weight - length conversion coefficient
k=Param[6]               # [-] Wet weight - length conversion exponent
b=Param[9]               # [gWW] Dry weight - wet weight conversion coefficient
p=Param[5]               # [-] Dry weight - wet weight conversion exponent

Wd=as.vector(matrix(0,nrow=ti))   # Initialize vector dry weight
Wd[ti]=b*Ww[ti]^p                 # Dry weight initial value [g]
L=as.vector(matrix(0,nrow=ti))    # Initialize vector length
L[ti]=(Ww[ti]/a)^k                # Length of the mussel initial value [cm]


# Print to screen inserted parameters

cat(" \n")
cat('The model will be executed with the following parameters:\n');
cat(" \n")
for (i in 1:25){
cat(paste0(toString(Param_matrix[i,2]), ": ", toString(Param_matrix[i,3]), " " ,toString(Param_matrix[i,4])),"\n")
}

cat(" \n")
cat('Weight initial condition is: ', toString(IC)," g\n")
cat(" \n")
cat("Integration is performed between ", toString(Dates[1]), " and ", toString(Dates[2]),"\n")
cat(" \n")
cat('Commercial size is ', toString(CS)," mm")
cat(" \n")

# Plot to screen inserted forcing functions
cat(" \n")
cat("Forcings are represented in graphs available at the following folder\n")
cat(paste0(userpath,"/ClamF_individual/Inputs/Forcings_plots\n"))

# plot Temperature forcing
Tintsave=Tint[ti:tf]
currentpath=getwd()
filepath=paste0(userpath,"/ClamF_individual/Inputs/Forcings_plots//Water_temperature.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
plot(days, Tintsave, ylab="Water temperature (Celsius degrees)", xlab="",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot Chlorophyll a forcing
Chlaintsave=Chlint[ti:tf]
currentpath=getwd()
filepath=paste0(userpath,"/ClamF_individual/Inputs/Forcings_plots//Chlorophyll_a.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
plot(days, Chlaintsave, ylab="Chlorophyll_a (mg/m^3)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

output=list(Param, times, Dates, IC, Tint, Chlint,CS)

}
