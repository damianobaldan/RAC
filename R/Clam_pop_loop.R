#' Function that runs the Monte Carlo simulation for the Clam population model
#'
#' @param Param a vector containing model parameters
#' @param times integration extremes and integration timestep
#' @param IC initial condition
#' @param Tint the interpolated water temperature time series
#' @param Phyint the interpolated phytoplankton time series
#' @param DTint the interpolated detritus time series
#' @param POCint the interpolated POC time series
#' @param POMint the interpolated POM time series
#' @param TSSint the interpolated TSS time series
#' @param N time series with number of individuals
#' @param userpath the path where the working folder is located
#'
#' @return a list with RK solver outputs
#'
#' @import matrixStats plotrix rstudioapi
#'

Clam_pop_loop<-function(Param, times, IC, Tint, Phyint, DTint, POCint, POMint, TSSint,N,userpath) {

cat("Population processing\n")


  ti=times[1]
  tf=times[2]
  t0=times[4]

  # Read files with population parameters and management strategies
  Pop_matrix=read.csv(paste0(userpath,"/Clam_population/Inputs/Parameters//Population.csv"),sep=",")   # Reading the matrix containing population parameters and their description
  Management=read.csv(paste0(userpath,"/Clam_population/Inputs/Population_management//Management.csv"),sep=",")   # Reading the matrix containing seeding and harvesting management

  # Extract population parameters
  meanWd=as.double(as.matrix(Pop_matrix[1,3]))     # [g] Dry weight average
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

# Vectors initialization
saveIC=as.vector(matrix(0,nrow=nruns))      # initialize perturbed initial conditions vector
saveCRmax=as.vector(matrix(0,nrow=nruns))   # Initialize perturbed maximum clearamce rate vector
Wd=as.matrix(matrix(0,nrow=nruns,ncol=tf))  # Initialize dry weight vector
Ww=as.matrix(matrix(0,nrow=nruns,ncol=tf))  # Initialize wet weight vector
L=as.matrix(matrix(0,nrow=nruns,ncol=tf))   # Initialize length vector
A=as.matrix(matrix(0,nrow=nruns,ncol=tf))   # Initialize Anabolic rate vector
C=as.matrix(matrix(0,nrow=nruns,ncol=tf))   # Initializ<e catabolic rate vector

# Loop for ODE solution

pb <- txtProgressBar(min = 0, max = nruns, style = 3)

for (ii in 1:nruns){

  # Weight initialization
  IC=rnorm(1,meanWd,deltaWd)   # [g] initial weight extracted from a normal distribution
  IC=max(IC, Wdlb)             # Lower bound for weight distribution
  saveIC[ii]=IC                # Saves initial condition values on Wd for each run

  # Maximum clearance rate initialization
  CRmax=rnorm(1,meanCRmax,deltaCRmax) # [l/d gDW] Maximum ingestion rate extracted from a normal distribution
  CRmax=max(CRmax,0)                  # Forces maximum ingestion rate to be positive
  saveCRmax[ii]=CRmax                 # Saves initial condition values on Imax for each run

  # Perturbe the parameters vector
  Param[6]=CRmax

  # Solves ODE with perturbed parameters
  output<-Clam_pop_RKsolver(Param, times, IC, Tint, Phyint, DTint, POCint, POMint, TSSint)

  # Extract outputs
  weight=t(output[[1]])
  Tfun=output[[2]]
  metab=output[[3]]

  # Saves results of each run to compute statistics
  Wd[ii,1:length(weight[,1])]=weight[,1]        # Clam dry weight [g]
  Ww[ii,1:length(weight[,2])]=weight[,2]        # Clam wet weight [g]
  L[ii,1:length(weight[,3])]=weight[,3]         # Clam length [g]

  A[ii,1:length(metab[,1])]=metab[,1]     # Net anabolism [J/d]
  C[ii,1:length(metab[,2])]=metab[,2]     # Fasting catabolism [J/d]

  setTxtProgressBar(pb, ii)

} # Close population loop

close(pb)

# Temperaure limitation functions

fgT=Tfun[,1] # Optimum  dependance from temperature for ingestion
frT=Tfun[,2] # Exponential dependance from temperature for catabolism

# Statistics computation

Wd_stat=rbind(colMeans(Wd), colSds(Wd))
Ww_stat=rbind(colMeans(Ww), colSds(Ww))
L_stat=rbind(colMeans(L), colSds(L))

A_stat=rbind(colMeans(A), colSds(A))
C_stat=rbind(colMeans(C), colSds(C))

output=list(Wd_stat,Ww_stat,L_stat,A_stat,C_stat,fgT,frT)
return(output)

}
