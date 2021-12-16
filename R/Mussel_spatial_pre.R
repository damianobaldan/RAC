#' Mussel bioenergetic spatialized model preprocessor
#'
#' @param userpath the path where folder containing model inputs and outputs is located
#' @param forcings a list containing forcings used by the model
#' @return a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Mussel_spatial_pre<-function(userpath,forcings){

  rm(list=ls())       # Clean workspace

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
  Dates=as.Date(Param_matrix[37:38,3], "%d/%m/%Y")                          # Vector containing the starting and ending date of teh simulation
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

  DTint=0
  Phyint=0

  # POC composition
  Ccont=CPOCint/CPOCint                            # Seston C/C content in weight
  Ncont=NPOCint/CPOCint                            # Seston N/C content in weight
  Pcont=array(1,c(1,length(CPOCint)))/CPOCint      # Seston P/C content in weight

  # plot POC forcing
  POCintsave=POCint[ti:tf]
  currentpath=getwd()
  filepath=paste0(userpath,"/Mussel_spatial/Inputs/Forcings_plots//POC.jpeg")
  jpeg(filepath,800,600)
  days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
  plot(days, POCintsave, ylab="POC (mgC/l)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()

  # Plot POM forcing
  POMintsave=POMint[ti:tf]
  currentpath=getwd()
  filepath=paste0(userpath,"/Mussel_spatial/Inputs/Forcings_plots//POM.jpeg")
  jpeg(filepath,800,600)
  days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
  plot(days, POMintsave, ylab="POM (mg/l)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()

  # Plot TSM forcing
  TSSintsave=TSSint[ti:tf]
  currentpath=getwd()
  filepath=paste0(userpath,"/Mussel_spatial/Inputs/Forcings_plots//TSM.jpeg")
  jpeg(filepath,800,600)
  days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
  plot(days, TSSintsave, ylab="TSM (mg/l)", xlab="", xaxt = "n",type="l",cex.lab=1.4)
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  dev.off()

  # Plot POC charactherization (C/P and N/P molar ratios)
  CPOCintsave=CPOCint[ti:tf]
  NPOCintsave=NPOCint[ti:tf]
  currentpath=getwd()
  filepath=paste0(userpath,"/Mussel_spatial/Inputs/Forcings_plots//POC characterization.jpeg")
  jpeg(filepath,800,600)
  lb=0
  ub=max(CPOCintsave,NPOCintsave)
  days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1)
  plot(days, CPOCintsave, ylab="POC characterization (molar ratios)", xlab="", xaxt = "n",type="l",cex.lab=1.4,col="blue",ylim=c(lb,ub+0.1*ub))
  lines(days,NPOCintsave,col="red")
  labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
  axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
  legend("topright",c("C/P","N/P"),fill=c("red","blue"))
  dev.off()

  output=list(Param, times, Dates, IC, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint,CS)
}
