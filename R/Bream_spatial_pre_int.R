#' Bream bioenergetic spatialized model preprocessor - used inside spatialization loop
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @param forcings a list containing forcings used by the model
#' @return a list containing the data used by the main script
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Bream_spatial_pre_int<-function(userpath,forcings){

  # Extracts forcings values from the list
  timeT=forcings[[1]]
  Tint=forcings[[2]]
  timeG=forcings[[3]]
  Gint=forcings[[4]]

  #  parameters from .csv files
  Param_matrix=read.csv(paste0(userpath,"/Bream_spatial/Inputs/parameters//Parameters.csv"),sep=",")                       # Reading the matrix containing parameters and their description
  Food=read.csv(paste0(userpath,"/Bream_spatial/Inputs/Point forcings//Food_characterization.csv"),sep=",",header=FALSE)    # Reading the food composition (Proteins, Lipids, Carbohydrates) data

  # Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
  Param=as.double(as.matrix(Param_matrix[1:21,3]))          # Vector containing all parameters
  Dates=Param_matrix[22:23,3]                               # Vector containing the starting and ending date of the simulation
  IC=as.double(as.matrix(Param_matrix[24,3]))               # Initial weight condition
  CS=as.double(as.matrix(Param_matrix[25,3]))               # Commercial size
  Food=as.double(as.matrix(Food[,1]))                       # Food composition (Proteins, Lipids, Carbohydrates) data

  # Prepare data for ODE solution
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeG[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
  timestep=1                                        # Time step for integration [day]
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0   # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0   # End of integration [day]
  weight=as.vector(matrix(0,nrow=ti))               # Initialize vector weight
  weight[ti]=IC                                     # Weight initial value [g]
  times<-cbind(ti, tf, timestep)

  # Food composition vector
  Pcont=Food[1]       # [-] Percentage of proteins in the food
  Lcont=Food[2]       # [-] Percentage of lipids in the food
  Ccont=Food[3]       # [-] Percentage of carbohydrates in the food

  output=list(Param, Tint, Gint, Food, IC, times, Dates, CS)
}
