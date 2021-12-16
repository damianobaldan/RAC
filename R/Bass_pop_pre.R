#' Seabass bioenergetic population model preprocessor
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

Bass_pop_pre<-function(userpath,forcings){

  rm(list=ls())       # Clean workspace

  cat("Data preprocessing")

  # Extracts forcings values from the list
  timeT=forcings[[1]]
  Tint=forcings[[2]]
  timeG=forcings[[3]]
  Gint=forcings[[4]]


  # Read forcings and parameters from .csv files
  Param_matrix=read.csv(paste0(userpath,"/Bass_population/Inputs/Parameters//Parameters.csv"),sep=",")                    # Reading the matrix containing parameters and their description
  Food=read.csv(paste0(userpath,"/Bass_population/Inputs/Forcings//Food_characterization.csv"),sep=",",header=FALSE)    # Reading the food composition (Proteins, Lipids, Carbohydrates) data

  # Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
  Param=as.matrix(Param_matrix[1:21,3])           # Vector containing all parameters
  Param=suppressWarnings(as.numeric(Param))
  Dates=Param_matrix[22:23,3]                                # Vector containing the starting and ending date of the simulation
  #IC=as.double(as.matrix(Param_matrix[24,3]))                # Initial weight condition
  CS=as.double(as.matrix(Param_matrix[25,3]))                # Commercial size
  Food=as.double(as.matrix(Food[,1]))                        # Food composition (Proteins, Lipids, Carbohydrates) data

  # Prepare data for ODE solution
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeG[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) #  Minimum starting date for forcings and observations
  timestep=1                                        # Time step for integration [day]
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0   # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0   # End of integration [day]
  #weight=as.vector(matrix(0,nrow=ti))               # Initialize vector weight
  #weight[ti]=IC                                     # Weight initial value [g]
  times<-cbind(ti, tf, timestep,t0)                    # Vector with integration data

  # Food composition vector
  Pcont=Food[1]       # [-] Percentage of proteins in the food
  Lcont=Food[2]       # [-] Percentage of lipids in the food
  Ccont=Food[3]       # [-] Percentage of carbohydrates in the food

  # Read files with population parameters and management strategies
  Pop_matrix=read.csv(paste0(userpath,"/Bass_population/Inputs/Parameters//Population.csv"),sep=",")   # Reading the matrix containing population parameters and their description
  Management=read.csv(paste0(userpath,"/Bass_population/Inputs/Population_management//Management.csv"),sep=",")   # Reading the matrix containing seeding and harvesting management

  # Extract population parameters
  meanW=as.double(as.matrix(Pop_matrix[1,3]))      # [g] Dry weight average
  IC=meanW
  deltaW=as.double(as.matrix(Pop_matrix[2,3]))     # [g] Dry weight standard deviation
  Wlb=as.double(as.matrix(Pop_matrix[3,3]))        # [g] Dry weight lower bound
  meanImax=as.double(as.matrix(Pop_matrix[4,3]))   # [l/d gDW] Clearence rate average
  deltaImax=as.double(as.matrix(Pop_matrix[5,3]))  # [l/d gDW] Clearance rate standard deviation
  Nseed=as.double(as.matrix(Pop_matrix[6,3]))      # [-] number of seeded individuals
  mortmyt=as.double(as.matrix(Pop_matrix[7,3]))    # [1/d] natural mortality rate
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
  N<-Pop_fun(Nseed, mortmyt, manag, times)

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
  cat("The food has the following composition: \n")
  cat(toString(Pcont*100),"% proteins\n")
  cat(toString(Lcont*100),"% lipids\n")
  cat(toString(Ccont*100),"% carbohydrates\n")
  cat(" \n")
  cat('Commercial size is ', toString(CS)," g")
  cat(" \n")

  # Print to screen population characteristics
  cat(" \n")
  cat('The population is simulated by assuming that initial weight and maximum ingestion rate are normally distributed:\n');
  cat(" \n")
  for (i in 1:5){
    cat(paste0(toString(Pop_matrix[i,2]), ": ", toString(Pop_matrix[i,3]), " " ,toString(Pop_matrix[i,4])),"\n")
  }

  cat(" \n")
  cat("The population is initially composed by ", toString(Pop_matrix[7,3]), " Individuals\n")
  cat(" \n")
  cat("The mortality rate is:", toString(Pop_matrix[8,3]),'1/d\n' )

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
  cat(paste0(userpath,"/Bass_population/Inputs/Forcings_plots\n"))

  # Plot Temperature forcing
  Tintsave=Tint[(ti+1):tf]
  currentpath=getwd()
  filepath=paste0(userpath,"/Bass_population/Inputs/Forcings_plots//Water_temperature.jpeg")
  jpeg(filepath,800,600)
  days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)
  plot(days, Tintsave, ylab="Water temperature (Celsius degrees)", xlab="",xaxt = "n",type="l",cex.lab=1.4)
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()

  # Plot feeding rate forcing
  Gintsave=Gint[(ti+1):tf]
  currentpath=getwd()
  filepath=paste0(userpath,"/Bass_population/Inputs/Forcings_plots//Feeding.jpeg")
  jpeg(filepath,800,600)
  days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)
  plot(days, Gintsave, ylab="Feed (g/d)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()

  # plot population dynamics
  Nsave=N[(ti+1):tf]
  currentpath=getwd()
  filepath=paste0(userpath,"/Bass_population/Inputs/Forcings_plots//Population.jpeg")
  jpeg(filepath,800,600)
  days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)
  plot(days, Nsave, ylab="Individuals", xlab="", xaxt = "n",type="l",cex.lab=1.4)
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()

  output=list(Param, Tint, Gint, Food, IC, times, Dates, N, CS)

}
