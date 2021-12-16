#' Seabass bioenergetic individual model preprocessor
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @param forcings a list containing model forcings
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees] and feeding rate [g/individual x d]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Bass_ind_pre<-function(userpath,forcings){

rm(list=ls())       # Clean workspace

cat("Data preprocessing")
cat(" \n")

# Extracts forcings values from the list
timeT=forcings[[1]]
Tint=forcings[[2]]
timeG=forcings[[3]]
Gint=forcings[[4]]

# Read forcings and parameters from .csv files
Param_matrix=read.csv(paste0(userpath,"/Bass_individual/Inputs/parameters//Parameters.csv"),sep=",")                    # Reading the matrix containing parameters and their description
Food=read.csv(paste0(userpath,"/Bass_individual/Inputs/Forcings//Food_characterization.csv"),sep=",",header=FALSE)    # Reading the food composition (Proteins, Lipids, Carbohydrates) data

# Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
Param=as.double(as.matrix(Param_matrix[1:21,3]))          # Vector containing all parameters
Dates=Param_matrix[22:23,3]                               # Vector containing the starting and ending date of the simulation
IC=as.double(as.matrix(Param_matrix[24,3]))               # Initial weight condition
CS=as.double(as.matrix(Param_matrix[25,3]))                # Commercial size
Food=as.double(as.matrix(Food[,1]))                       # Food composition (Proteins, Lipids, Carbohydrates) data

# Prepare data for ODE solution
t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeG[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
timestep=1                                        # Time step for integration [day]
ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0   # Start of integration [day]
tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0   # End of integration [day]
weight=as.vector(matrix(0,nrow=ti))               # Initialize vector weight
weight[ti]=IC                                     # Weight initial value [g]
times<-cbind(ti, tf, timestep)                    # Vector with integration data

# Food composition vector
Pcont=Food[1]       # [-] Percentage of proteins in the food
Lcont=Food[2]       # [-] Percentage of lipids in the food
Ccont=Food[3]       # [-] Percentage of carbohydrates in the food

# Print to screen inserted parameters

cat(" \n")
cat('The model will be executed with the following parameters:\n');
cat(" \n")
for (i in 1:21){
cat(paste0(toString(Param_matrix[i,2]), ": ", toString(Param_matrix[i,3]), " " ,toString(Param_matrix[i,4])),"\n")
}


cat(" \n")
cat('Weight initial condition is: ', toString(IC)," g\n")
cat(" \n")
cat("Integration is performed between ", toString(Dates[1]), " and ", toString(Dates[2]),"\n")
cat(" \n")
cat("The food has the following composition: \n")
cat(toString(Pcont*100),"% proteins\n")
cat(toString(Lcont*100),"% lipids\n")
cat(toString(Ccont*100),"% carbohydrates\n")
cat(" \n")
cat('Commercial size is ', toString(CS)," g")
cat(" \n")

# Plot to screen inserted forcing functions
cat(" \n")
cat("Forcings are represented in graphs available at the following folder\n")
cat(paste0(userpath,"/Bass_individual/Inputs/Forcings_plots\n"))

# Plot Temperature forcing
Tintsave=Tint[(ti+1):tf]
currentpath=getwd()
filepath=paste0(userpath,"/Bass_individual/Inputs/Forcings_plots//Water_temperature.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)
plot(days, Tintsave, ylab="Water temperature (Celsius degrees)", xlab="",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot Feeding rate forcing
Gintsave=Gint[(ti+1):tf]
currentpath=getwd()
filepath=paste0(userpath,"/Bass_individual/Inputs/Forcings_plots//Feeding.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)
plot(days, Gintsave, ylab="Feed (g/d)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

output=list(Param, Tint, Gint, Food, IC, times, Dates, CS)
return(output)
}
