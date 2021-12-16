#' Clam bioenergetic population model preprocessor
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @param forcings a list containing model forcings
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l], particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Clam_pop_pre<-function(userpath,forcings){


rm(list=ls())       # Clean workspace

cat("Data preprocessing")

# Extracts forcings values from the list
timeT=forcings[[1]]
Tint=forcings[[2]]
timeChl=forcings[[3]]
Chlint=forcings[[4]]
timePOC=forcings[[5]]
POCint=forcings[[6]]
timePOM=forcings[[7]]
POMint=forcings[[8]]
timeTSS=forcings[[9]]
TSSint=forcings[[10]]

# Read forcings and parameters from .csv files
Param_matrix=read.csv(paste0(userpath,"/Clam_population/Inputs/Parameters//Parameters.csv"),sep=",")                      # Reading the matrix containing parameters and their description

# Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
Param=as.matrix(Param_matrix[1:22,3])    # Vector containing all parameters
Param=suppressWarnings(as.numeric(Param))
Dates=Param_matrix[24:25,3]                          # Vector containing the starting and ending date of teh simulation
#IC=as.double(as.matrix(Param_matrix[23,3]))          # Initial weight condition
CS=as.double(as.matrix(Param_matrix[26,3]))                # Commercial size

# Prepare data for ODE solution
t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeChl[1], "%d/%m/%Y")), as.numeric(as.Date(timePOC[1], "%d/%m/%Y")), as.numeric(as.Date(timePOM[1], "%d/%m/%Y")), as.numeric(as.Date(timeTSS[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
timestep=1                                           # Time step for integration [day]
ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0      # Start of integration [day]
tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0      # End of integration [day]
#weight=as.vector(matrix(0,nrow=ti))                  # Initialize vector weight
#weight[ti]=IC                                        # Weight initial value [g]
times<-cbind(ti, tf, timestep,t0)                       # Vector with integration data

# Detritus
lambda=Param[22]    # Chla-C proportionality coefficient: relates Chlorophyll-a to phytoplankton concentration
Phyint=as.vector(matrix(0,nrow=ti-1))  # Initialize vector Phyint
DTint=as.vector(matrix(0,nrow=ti-1))   # Initialize vector DTint
for (i in ti:tf) {
  Phyint[i]=lambda*Chlint[i]/1000    # Phytoplankton concentration [mgC/l]
  DTint[i]=POCint[i]-Phyint[i]       # Detritus concentration [mgC/l]
  if (DTint[i] < 0)   {              # Adjusts "wrong" detritus values
    DTint[i]=0
  }  # end if
}  # end for

# Read files with population parameters and management strategies
Pop_matrix=read.csv(paste0(userpath,"/Clam_population/Inputs/Parameters//Population.csv"),sep=",")   # Reading the matrix containing population parameters and their description
Management=read.csv(paste0(userpath,"/Clam_population/Inputs/Population_management//Management.csv"),sep=",")   # Reading the matrix containing seeding and harvesting management

# Extract population parameters
meanWd=as.double(as.matrix(Pop_matrix[1,3]))     # [g] Dry weight average
IC=meanWd
deltaWd=as.double(as.matrix(Pop_matrix[2,3]))    # [g] Dry weight standard deviation
Wdlb=as.double(as.matrix(Pop_matrix[3,3]))       # [g] Dry weight lower bound
meanCRmax=as.double(as.matrix(Pop_matrix[4,3]))  # [l/d gDW] Clearence rate average
deltaCRmax=as.double(as.matrix(Pop_matrix[5,3])) # [l/d gDW] Clearance rate standard deviation
Nseed=as.double(as.matrix(Pop_matrix[6,3]))      # [-] number of seeded individuals
mort=as.double(as.matrix(Pop_matrix[7,3]))       # [1/d] natural mortality rate
nruns=as.double(as.matrix(Pop_matrix[8,3]))      # [-] number of runs for population simulation

# Prepare management values
manag=as.matrix(matrix(0,nrow=length(Management[,1]),ncol=2))
for (i in 1:length(Management[,1])) {
     manag[i,1]=as.numeric(as.Date(Management[i,1], "%d/%m/%Y"))-t0
  if ((Management[i,2])=="h") {
    manag[i,2]=-as.numeric(Management[i,3])
    } else {
    manag[i,2]=as.numeric(Management[i,3])
  }
}

# Population differential equation solution
N<-Pop_fun(Nseed, mort, manag, times)

# Print to screen inserted parameters

cat(" \n")
cat('The model will be executed with the following parameters:\n');
cat(" \n")
for (i in 1:21){
cat(paste0(toString(Param_matrix[i,2]), ": ", toString(Param_matrix[i,3]), " " ,toString(Param_matrix[i,4])),"\n")
}

cat(" \n")
cat("Integration is performed between ", toString(Dates[1]), " and ", toString(Dates[2]),"\n")
cat(" \n")
cat('Commercial size is ', toString(CS)," mm")
cat(" \n")

# Print to screen population characteristics

cat(" \n")
cat('The population is simulated by assuming that initial weight and maximum clearence rate are normally distributed:\n');
cat(" \n")
for (i in 1:5){
  cat(paste0(toString(Pop_matrix[i,2]), ": ", toString(Pop_matrix[i,3]), " " ,toString(Pop_matrix[i,4])),"\n")
}

cat(" \n")
cat("The population is initially composed by ", toString(Pop_matrix[6,3]), " Individuals\n")
cat(" \n")
cat("The mortality rate is:", toString(Pop_matrix[7,3]),'1/d\n' )

# Print to screen management actions
cat(" \n")
cat('The population is managed according with following list (h:harvesting s:seeding):\n');
cat(" \n")
for (i in 1:length(Management[,1])){
  cat(paste0(toString(Management[i,1])," ", toString(Management[i,2]), " " ,toString(Management[i,3])),"individuals\n")
}

cat(" \n")
cat("The individual model will be executed ", toString(nruns), " times in order to simulate a population\n")
cat(" \n")

# Plot to file inserted forcing functions

cat(" \n")
cat("Forcings are represented in graphs available at the following folder:\n")
cat(paste0(userpath,"/Clam_population/Inputs/Forcings_plots\n"))

# Plot Temperature forcing
Tintsave=Tint[ti:tf]
currentpath=getwd()
filepath=paste0(userpath,"/Clam_population/Inputs/Forcings_plots//Water_temperature.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
plot(days, Tintsave, ylab="Water temperature (Celsius degrees)", xlab="",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot Chla-a forcing
Chlaintsave=Chlint[ti:tf]
currentpath=getwd()
filepath=paste0(userpath,"/Clam_population/Inputs/Forcings_plots//Chlorophyll_a.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
plot(days, Chlaintsave, ylab="Chlorophyll_a (mg/m^3)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot POC forcing
POCintsave=POCint[ti:tf]
currentpath=getwd()
filepath=paste0(userpath,"/Clam_population/Inputs/Forcings_plots//POC.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
plot(days, POCintsave, ylab="POC (mgC/l)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot POM forcing
POMintsave=POMint[ti:tf]
currentpath=getwd()
filepath=paste0(userpath,"/Clam_population/Inputs/Forcings_plots//POM.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
plot(days, POMintsave, ylab="POM (mg/l)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot TSS forcing
TSSintsave=TSSint[ti:tf]
currentpath=getwd()
filepath=paste0(userpath,"/Clam_population/Inputs/Forcings_plots//TSM.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
plot(days, TSSintsave, ylab="TSM (mg/l)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot population dynamics
Nsave=N[(ti+1):tf]
currentpath=getwd()
filepath=paste0(userpath,"/Clam_population/Inputs/Forcings_plots//Population.jpeg")
jpeg(filepath,800,600)
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)
plot(days, Nsave, ylab="Individuals", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

output=list(Param, times, Dates, IC, Tint, Phyint, DTint, POCint, POMint, TSSint,N, CS)

return(output)
}
