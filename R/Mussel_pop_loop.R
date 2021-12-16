#' Function that runs the Monte Carlo simulation for the Mussel population model
#'
#' @param Param a vector containing model parameters
#' @param times integration extremes and integration timestep
#' @param IC initial condition
#' @param Tint the interpolated water temperature time series
#' @param Phyint the interpolated phytoplankton time series
#' @param DTint the interpolated detritus time series
#' @param POCint the interpolated POC time series
#' @param Ccont the C/C content of the POC
#' @param Ncont the N/C content of POC
#' @param Pcont the P/C content of POC
#' @param POMint the interpolated POM time series
#' @param TSSint the interpolated TSS time series
#' @param N time series with number of individuals
#' @param userpath the path where the working folder is located
#'
#' @return a list with RK solver outputs
#'
#' @import matrixStats plotrix rstudioapi
#'

Mussel_pop_loop<-function(Param, times, IC, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint,N,userpath) {

cat("Population processing\n")

  ti=times[1]
  tf=times[2]
  t0=times[4]

  # Read files with population parameters and management strategies
  Pop_matrix=read.csv(paste0(userpath,"/Mussel_population/Inputs/Parameters//Population.csv"),sep=",")   # Reading the matrix containing population parameters and their description
  Management=read.csv(paste0(userpath,"/Mussel_population/Inputs/Population_management//Management.csv"),sep=",")   # Reading the matrix containing seeding and harvesting management

  # Extract population parameters
  meanWd=as.double(as.matrix(Pop_matrix[1,3]))   # [g] Dry weight average
  deltaWd=as.double(as.matrix(Pop_matrix[2,3]))  # [g] Dry weight standard deviation
  Wdlb=as.double(as.matrix(Pop_matrix[3,3]))     # [g] Dry weight lower bound
  meanCRmax=as.double(as.matrix(Pop_matrix[4,3]))   # [l/d gDW] Clearence rate average
  deltaCRmax=as.double(as.matrix(Pop_matrix[5,3]))  # [l/d gDW] Clearance rate standard deviation
  Nseed=as.double(as.matrix(Pop_matrix[6,3]))    # [-] number of seeded individuals
  mort=as.double(as.matrix(Pop_matrix[7,3]))  # [1/d] natural mortality rate
  nruns=as.double(as.matrix(Pop_matrix[8,3]))    # [-] number of runs for population simulation

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
saveCRmax=as.vector(matrix(0,nrow=nruns))

Wb=as.matrix(matrix(0,nrow=nruns,ncol=tf))
R=as.matrix(matrix(0,nrow=nruns,ncol=tf))
Wd=as.matrix(matrix(0,nrow=nruns,ncol=tf))
W=as.matrix(matrix(0,nrow=nruns,ncol=tf))
L=as.matrix(matrix(0,nrow=nruns,ncol=tf))

fecC=as.matrix(matrix(0,nrow=nruns,ncol=tf))
fecN=as.matrix(matrix(0,nrow=nruns,ncol=tf))
fecP=as.matrix(matrix(0,nrow=nruns,ncol=tf))

psC=as.matrix(matrix(0,nrow=nruns,ncol=tf))
psN=as.matrix(matrix(0,nrow=nruns,ncol=tf))
psP=as.matrix(matrix(0,nrow=nruns,ncol=tf))

Cmyt=as.matrix(matrix(0,nrow=nruns,ncol=tf))
Nmyt=as.matrix(matrix(0,nrow=nruns,ncol=tf))
Pmyt=as.matrix(matrix(0,nrow=nruns,ncol=tf))

A=as.matrix(matrix(0,nrow=nruns,ncol=tf))
C=as.matrix(matrix(0,nrow=nruns,ncol=tf))

O2=as.matrix(matrix(0,nrow=nruns,ncol=tf))
NH4=as.matrix(matrix(0,nrow=nruns,ncol=tf))

# Loop for ODE solution

pb <- txtProgressBar(min = 0, max = nruns, style = 3)

for (ii in 1:nruns){

  # Weight initialization
  IC=rnorm(1,meanWd,deltaWd)     # [g] initial weight extracted from a normal distribution
  IC=max(IC, Wdlb)               # Lower bound for weight distribution
  saveIC[ii]=IC                  # Saves initial condition values on Wd for each run

  # Maximum clearance rate initialization
  CRmax=rnorm(1,meanCRmax,deltaCRmax)  # [l/d gDW] Maximum clearance rate extracted from a normal distribution
  CRmax=max(CRmax,0)                   # Forces maximum clearance rate to be positive
  saveCRmax[ii]=CRmax                  # Saves initial condition values on CRmax for each run

  # Perturbe the parameters vector
  Param[7]=CRmax

  # Solves ODE with perturbed parameters
  output<-Mussel_pop_RKsolver(Param, times, IC, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint, N)

  # Extract outputs
  weight=t(output[[1]])
  pfec=output[[2]]
  fec=output[[3]]
  comp=output[[4]]
  Tfun=output[[5]]
  metab=output[[6]]
  cons=output[[7]]
  amm=output[[8]]

  # Saves results of each run to compute statistics
  Wb[ii,1:length(weight[,1])]=weight[,1]        # Somatic tissue dry weight [g]
  R[ii,1:length(weight[,2])]=weight[,2]         # Gonadic tissue dry weight [g]
  Wd[ii,1:length(weight[,3])]=weight[,3]        # Total dry weight [g]
  W[ii,1:length(weight[,4])]=weight[,4]         # Mussel weight with shell [g]
  L[ii,1:length(weight[,5])]=weight[,5]         # Mussel length [g]

  psC[ii,1:length(fec[,1])]=pfec[,1]        # pseudofecies C production
  psN[ii,1:length(fec[,2])]=pfec[,2]        # pseudofecies N production
  psP[ii,1:length(fec[,3])]=pfec[,3]        # pseudofecies P production

  fecC[ii,1:length(fec[,1])]=fec[,1]        # fecies C production
  fecN[ii,1:length(fec[,1])]=fec[,2]        # fecies N production
  fecP[ii,1:length(fec[,1])]=fec[,3]         # fecies P production

  Cmyt[ii,1:length(comp[,1])]=comp[,1]        # Mussel C content
  Nmyt[ii,1:length(comp[,2])]=comp[,2]        # Mussel N content
  Pmyt[ii,1:length(comp[,3])]=comp[,3]        # Mussel P content

  A[ii,1:length(metab[,1])]=metab[,1]     # Net anabolism [J/d]
  C[ii,1:length(metab[,2])]=metab[,2]     # Fasting catabolism [J/d]

  O2[ii,1:length(cons)]=cons  # Oxygen consumtion rate
  NH4[ii,1:length(amm)]=amm   # ammonium release

  setTxtProgressBar(pb, ii)

} # Close population loop

close(pb)

# Temperaure limitation functions

fgT=Tfun[,1] # Optimum  dependance from temperature for ingestion
frT=Tfun[,2] # Exponential dependance from temperature for catabolism

# Statistics computation

Wb_stat=rbind(colMeans(Wb), colSds(Wb))
R_stat=rbind(colMeans(R), colSds(R))
Wd_stat=rbind(colMeans(Wd), colSds(Wd))
W_stat=rbind(colMeans(W), colSds(W))
L_stat=rbind(colMeans(L), colSds(L))

fecC_stat=rbind(colMeans(fecC), colSds(fecC))
fecN_stat=rbind(colMeans(fecN), colSds(fecN))
fecP_stat=rbind(colMeans(fecP), colSds(fecP))

psC_stat=rbind(colMeans(psC), colSds(psC))
psN_stat=rbind(colMeans(psN), colSds(psN))
psP_stat=rbind(colMeans(psP), colSds(psP))

Cmyt_stat=rbind(colMeans(Cmyt), colSds(Cmyt))
Nmyt_stat=rbind(colMeans(Nmyt), colSds(Nmyt))
Pmyt_stat=rbind(colMeans(Pmyt), colSds(Pmyt))

A_stat=rbind(colMeans(A), colSds(A))
C_stat=rbind(colMeans(C), colSds(C))

O2_stat=rbind(colMeans(O2), colSds(O2))
NH4_stat=rbind(colMeans(NH4), colSds(NH4))

output=list(Wb_stat,R_stat,Wd_stat,W_stat,L_stat,fecC_stat,fecN_stat,fecP_stat,psC_stat,psN_stat,psP_stat,Cmyt_stat,Nmyt_stat,Pmyt_stat,A_stat,C_stat,O2_stat,NH4_stat,fgT,frT)
return(output)

}
