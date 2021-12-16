#' Mussel bioenergetic spatialized model preprocessor - used inside spatialization loop
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @param forcings a list containing forcings used by the model
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Mussel_spatial_pre_int<-function(userpath,forcings){

 # rm(list=ls())       # Clean workspace

  # Extracts forcings values from the list
  timeT=forcings[[1]]
  Tint=forcings[[2]]
  timeChl=forcings[[3]]
  Chlint=forcings[[4]]
  timePOC=forcings[[5]]
  POCint=forcings[[6]]
  timeCPOC=forcings[[7]]
  CPOCint=forcings[[8]]
  timeNPOC=forcings[[9]]
  NPOCint=forcings[[10]]
  timePOM=forcings[[11]]
  POMint=forcings[[12]]
  timeTSS=forcings[[13]]
  TSSint=forcings[[14]]

  #  parameters from .csv files
  Param_matrix=read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Parameters//Parameters.csv"),sep=",")                      # Reading the matrix containing parameters and their description

  # Extract parameters and forcing values from parameters matrix and convert to type 'double' the vector contents
  Param=as.double(as.matrix(Param_matrix[1:36,3]))     # Vector containing all parameters
  Dates=Param_matrix[37:38,3]                          # Vector containing the starting and ending date of teh simulation
  IC=as.double(as.matrix(Param_matrix[39,3]))          # Initial weight condition
  CS=as.double(as.matrix(Param_matrix[40,3]))                # Commercial size

  # Prepare data for ODE solution
  t0=min(as.numeric(as.Date(timeT[1], "%d/%m/%Y")), as.numeric(as.Date(timeChl[1], "%d/%m/%Y")), as.numeric(as.Date(timePOC[1], "%d/%m/%Y")), as.numeric(as.Date(timePOM[1], "%d/%m/%Y")), as.numeric(as.Date(timeTSS[1], "%d/%m/%Y")), as.numeric(as.Date(Dates[1], "%d/%m/%Y"))) # starting minimum starting date for forcings and observations
  timestep=1                                         # Time step for integration [day]
  ti=as.numeric(as.Date(Dates[1], "%d/%m/%Y"))-t0    # Start of integration [day]
  tf=as.numeric(as.Date(Dates[2], "%d/%m/%Y"))-t0    # End of integration [day]
  times<-cbind(ti, tf, timestep)                    # Vector with integration data

  # Initial conditions
  Wb=as.vector(matrix(0,nrow=ti))            # Initialize vector somatic tissue weight
  Wb[ti]=IC                                  # Somatic tissue initial value [g]
  R=as.vector(matrix(0,nrow=ti))             # Initialize vector gonadic tissue weight
  R[ti]=0                                    # Gonadic weight initial value [g]
  L=as.vector(matrix(0,nrow=ti))             # Initialize vector length
  L[ti]=Param[20]*(Wb[ti]+R[ti])^Param[21]   # Length of the mussel initial value [cm]

  # Detritus
  lambda=Param[33]    # Chla-C proportionality coefficient: relates Chlorophyll-a to phytoplankton concentration
  Phyint=as.vector(matrix(0,nrow=ti-1))  # Initialize vector Phyint
  DTint=as.vector(matrix(0,nrow=ti-1))   # Initialize vector DTint
  for (i in ti:tf) {
    Phyint[i]=lambda*Chlint[i]/1000   # Phytoplankton concentration [mgC/l]
    DTint[i]=POCint[i]-Phyint[i]       # Detritus concentration [mgC/l]
    if (DTint[i] < 0)   {              # Adjusts "wrong" detritus values
      DTint[i]=0
    }  # end if
  }  # end for

  # POC composition
  Ccont=CPOCint/CPOCint                            # Seston C/C content in weight
  Ncont=NPOCint/CPOCint                            # Seston N/C content in weight
  Pcont=array(1,c(1,length(CPOCint)))/CPOCint      # Seston P/C content in weight

  output=list(Param, times, Dates, IC, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint,CS)
}
