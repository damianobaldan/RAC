#' Function that runs the Monte Carlo simulation for the Clam population model (alternative version)
#'
#' @param Param a vector containing model parameters
#' @param times integration extremes and integration timestep
#' @param IC initial condition
#' @param Tint the interpolated water temperature time series
#' @param Chlint the interpolated chlorophyll a time series
#' @param N time series with number of individuals
#' @param userpath the path where the working folder is located
#'
#' @return a list with RK solver outputs
#'
#' @import matrixStats plotrix rstudioapi
#'

ClamF_pop_loop<-function(Param, times, IC, Tint, Chlint, N,userpath) {

  cat("Population processing\n")

  ti=times[1]
  tf=times[2]
  t0=times[4]
  # Read files with population parameters and management strategies
  Pop_matrix=read.csv(paste0(userpath,"/ClamF_population/Inputs/Parameters//Population.csv"),sep=",")   # Reading the matrix containing population parameters and their description
  Management=read.csv(paste0(userpath,"/ClamF_population/Inputs/Population_management//Management.csv"),sep=",")   # Reading the matrix containing seeding and harvesting management

  # Extract population parameters
  meanWw=as.double(as.matrix(Pop_matrix[1,3]))    # [g] Wet weight average
  deltaWw=as.double(as.matrix(Pop_matrix[2,3]))   # [g] Wet weight standard deviation
  Wwlb=as.double(as.matrix(Pop_matrix[3,3]))      # [g] Wet weight lower bound
  meanGdmax=as.double(as.matrix(Pop_matrix[4,3])) # [l/d gDW] Maximum growth rate on a dry weight average
  deltaGdmax=as.double(as.matrix(Pop_matrix[5,3]))# [l/d gDW] Maximum growth rate on a dry weight standard deviation
  Nseed=as.double(as.matrix(Pop_matrix[6,3]))     # [-] number of seeded individuals
  mort=as.double(as.matrix(Pop_matrix[7,3]))      # [1/d] natural mortality rate
  nruns=as.double(as.matrix(Pop_matrix[8,3]))     # [-] number of runs for population simulation

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
saveIC=as.vector(matrix(0,nrow=nruns))
saveGrmax=as.vector(matrix(0,nrow=nruns))
Wd=as.matrix(matrix(0,nrow=nruns,ncol=tf))
Ww=as.matrix(matrix(0,nrow=nruns,ncol=tf))
L=as.matrix(matrix(0,nrow=nruns,ncol=tf))
A=as.matrix(matrix(0,nrow=nruns,ncol=tf))
C=as.matrix(matrix(0,nrow=nruns,ncol=tf))

# sources Runge-Kutta solver

pb <- txtProgressBar(min = 0, max = nruns, style = 3)

for (ii in 1:nruns){

  # Weight initialization
  IC=rnorm(1,meanWw,deltaWw)     # [g] initial weight extracted from a normal distribution
  IC=max(IC, Wwlb)               # Lower bound for weight distribution
  saveIC[ii]=IC                  # Saves initial condition values on Wd for each run

  # Maximum clearance rate initialization
  Gdmax=rnorm(1,meanGdmax,deltaGdmax)  # [l/d gDW] Maximum ingestion rate extracted from a normal distribution
  Gdmax=max(Gdmax,0)                  # Forces maximum ingestion rate to be positive
  saveGrmax[ii]=Gdmax                 # Saves initial condition values on Imax for each run

  # Perturbe the parameters vector
  Param[1]=Gdmax

  # Solves ODE with perturbed parameters
  output<-ClamF_pop_RKsolver(Param, times, IC, Tint, Chlint)

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
