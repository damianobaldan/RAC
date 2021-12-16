#' Function that runs the Monte Carlo simulation for the Seabream population model
#'
#' @param Param a vector containing model parameters
#' @param Tint the interpolated water temperature time series
#' @param Gint the interpolated feeding rate time series
#' @param Food the food characterization
#' @param IC initial condition
#' @param times integration extremes and integration timestep
#' @param N time series with number of individuals
#' @param userpath the path where the working folder is located
#'
#' @return a list with RK solver outputs
#'
#' @import matrixStats plotrix rstudioapi
#'

Bream_pop_loop<-function(Param, Tint, Gint, Food, IC, times, N, userpath) {

cat("Population processing\n")

  ti=times[1]
  tf=times[2]
  t0=times[4]

# Read files with population parameters and management strategies
Pop_matrix=read.csv(paste0(userpath,"/Bream_population/Inputs/Parameters//Population.csv"),sep=",")   # Reading the matrix containing population parameters and their description
Management=read.csv(paste0(userpath,"/Bream_population/Inputs/Population_management//Management.csv"),sep=",")   # Reading the matrix containing seeding and harvesting management

# Extract population parameters
meanW=as.double(as.matrix(Pop_matrix[1,3]))      # [g] Dry weight average
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

# Vectors initialization
saveIC=as.vector(matrix(0,nrow=nruns))            # Initialize initial conditions records: saves perturbated initial conditions for each run
saveImax=as.vector(matrix(0,nrow=nruns))          # Initialize maximum ingestion rate records: saves perturbated maximum ingestion rate for each run
W=as.matrix(matrix(0,nrow=nruns,ncol=tf))         # initialize weight vector
Pexc=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize excreted proteines vector
Lexc=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize excreted lipids vector
Cexc=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize excreted carbohydrates vector
Pwst=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize proteines to waste vector
Lwst=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize lipids to waste vector
Cwst=as.matrix(matrix(0,nrow=nruns,ncol=tf))      # Initialize carbohydrates to waste vector
ingestion=as.matrix(matrix(0,nrow=nruns,ncol=tf)) # Initialize actual ingestion vector
A=as.matrix(matrix(0,nrow=nruns,ncol=tf))         # Initialize anabolic rate vector
C=as.matrix(matrix(0,nrow=nruns,ncol=tf))         # Initialize catabolic rate vector
O2=as.matrix(matrix(0,nrow=nruns,ncol=tf))        # Initialize oxygen consumption rate vector
NH4=as.matrix(matrix(0,nrow=nruns,ncol=tf))       # Initialize ammonium release rate vector

# Population LOOP

pb <- txtProgressBar(min = 0, max = nruns, style = 3)

for (ii in 1:nruns){

  # Weight initialization
  IC=rnorm(1,meanW,deltaW)     # [g] initial weight extracted from a normal distribution
  IC=max(IC, Wlb)              # Lower bound for weight distribution
  saveIC[ii]=IC                # Saves initial condition values on Wd for each run

  # Maximum clearance rate initialization
  Imax=rnorm(1,meanImax,deltaImax)  # [l/d gDW] Maximum ingestion rate extracted from a normal distribution
  Imax=max(Imax,0)                  # Forces maximum ingestion rate to be positive
  saveImax[ii]=Imax                 # Saves initial condition values on Imax for each run

  # Perturbe the parameters vector
  Param[1]=Imax

  # Solves ODE with perturbed parameters
  output<-Bream_pop_RKsolver(Param, Tint, Gint, Food, IC, times, N)

  # Unlist outputs
  weight=unlist(output[1])
  exc=output[[2]]
  wst=output[[3]]
  ing=unlist(output[4])
  ingvero=unlist(output[5])
  Tfun=output[[6]]
  metab=output[[7]]
  oxygen=output[[8]]
  ammonium=output[[9]]

  # Saves results of each run to compute statistics
  W[ii,1:length(weight)]=weight           # Tissue dry weight [g]

  Pexc[ii,1:length(exc[,1])]=exc[,1]      # Excreted proteins [kg/d]
  Lexc[ii,1:length(exc[,1])]=exc[,2]      # Excreted lipids [kg/d]
  Cexc[ii,1:length(exc[,1])]=exc[,3]      # Excreted carbohydrates [kg/d]

  Pwst[ii,1:length(wst[,1])]=wst[,1]      # Proteins to waste [kg/d]
  Lwst[ii,1:length(wst[,1])]=wst[,2]      # Lipids to waste [kg/d]
  Cwst[ii,1:length(wst[,1])]=wst[,3]      # Carbohydrates to waste [kg/d]

  ingestion[ii,1:length(ingvero)]=t(ingvero)  # Actual ingested food [g/d]

  A[ii,1:length(metab[,1])]=metab[,1]     # Net anabolism [J/d]
  C[ii,1:length(metab[,1])]=metab[,2]     # Fasting catabolism [J/d]

  O2[ii,1:length(oxygen)]=oxygen           # Tissue dry weight [kgO2/d]
  NH4[ii,1:length(ammonium)]=ammonium      # Tissue dry weight [kgN/d]

  setTxtProgressBar(pb, ii)

} # Close population loop

close(pb)

# Temperaure limitation functions

fgT=Tfun[,1] # Optimum  dependance from temperature for ingestion
frT=Tfun[,2] # Exponential dependance from temperature for catabolism

# Statistics computation:
# Statistics vectors contain the mean and the standard deviation of a variable

W_stat=t(rbind(colMeans(W), colSds(W)))

Pexc_stat=t(rbind(colMeans(Pexc), colSds(Pexc)))
Lexc_stat=t(rbind(colMeans(Lexc), colSds(Lexc)))
Cexc_stat=t(rbind(colMeans(Cexc), colSds(Cexc)))

ingestion_stat=t(rbind(colMeans(ingestion), colSds(ingestion)))

Pwst_stat=t(rbind(colMeans(Pwst), colSds(Pwst)))
Lwst_stat=t(rbind(colMeans(Lwst), colSds(Lwst)))
Cwst_stat=t(rbind(colMeans(Cwst), colSds(Cwst)))

A_stat=t(rbind(colMeans(A), colSds(A)))
C_stat=t(rbind(colMeans(C), colSds(C)))

O2_stat=t(rbind(colMeans(O2), colSds(O2)))
NH4_stat=t(rbind(colMeans(NH4), colSds(NH4)))

output=list(W_stat,Pexc_stat,Lexc_stat,Cexc_stat,ingestion_stat,Pwst_stat,Lwst_stat,Cwst_stat,A_stat,C_stat,fgT,frT, O2_stat, NH4_stat)
return(output)

}
