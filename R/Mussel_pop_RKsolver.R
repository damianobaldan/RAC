#' Solves the Mussel population bioenergetic balance with a 4th order Runge Kutta method
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
#' @param N the number of indivduals time series
#'
#' @return a list containing the weights of the mussel, the excreted CNP, the mussel CNP, temperature limitation functions, metabolic rates, oxygen consumption
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats
#'

Mussel_pop_RKsolver <- function(Param, times, IC, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint, N){


  # Parameters definition

  # Spawning times
  tspawn1=Param[26]         # [d] First spawning: 15/12
  tspawn2=Param[27]         # [d] Second spawning:  20/5

  # Resting period duration
  ripi=Param[24]            # [d] Beginning of reproductory resting period
  ripf=Param[25]            # [d] End of reproductory resting period
  trip=cbind(ripi,ripf)

  # allometric parameters
  a=Param[20]               # [m] Weight to length proportionality constant
  b=Param[21]               # [-] Weight to length exponent
  aF=Param[28]              # [-] Dry weight - wet weight conversion coefficient
  atot=Param[29]            # [-] Dry weight - total (with shell) weight conversion coefficient

  # Integration extremes definition
  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end
  timestep=times[3]     # Timestep for integration

  # Initial condition definition
  Wb=as.vector(matrix(0,nrow=ti))            # Initialize vector somatic tissue weight
  Wb[ti]=IC                                  # Somatic tissue initial value [g]
  R=as.vector(matrix(0,nrow=ti))             # Initialize vector gonadic tissue weight
  R[ti]=0                                    # Gonadic weight initial value [g]
  L=as.vector(matrix(0,nrow=ti))             # Initialize vector length
  L[ti]=Param[20]*(Wb[ti]+R[ti])^Param[21]   # Length of the mussel initial value [cm]
  Wd=as.vector(matrix(0,nrow=ti))            # Initialize vector total dry weight
  Wd[ti]=Wb[ti]+R[ti]                        # total dry weight initial value [g]
  Wf=as.vector(matrix(0,nrow=ti))            # Initialize vector mussel wet weight
  Wf[ti]=aF*Wd[ti]                           # Mussel wet weight as a function of dry weight [g]
  Wtot=as.vector(matrix(0,nrow=ti))          # Initialize vector mussel total weight (with shell)
  Wtot[ti]=atot*Wd[ti]                       # Mussel total weight (with shell) as a function of dry weight [g]

  # initialize output
  pfec=as.matrix(matrix(0,nrow=ti,ncol=3))      # Initialize pseudofecies vector
  fec=as.matrix(matrix(0,nrow=ti,ncol=3))      # Initialize fecies vector
  comp=as.matrix(matrix(0,nrow=ti,ncol=3))     # Initialize mytilus composition vector
  tfun=as.matrix(matrix(0,nrow=ti,ncol=2))     # Initialize temperature limitations vector
  metab=as.matrix(matrix(0,nrow=ti,ncol=2))    # Initialize metabolic rates vector
  cons=as.matrix(matrix(0,nrow=ti,ncol=1))     # Initialize oxygen consumption vector
  amm=as.matrix(matrix(0,nrow=ti,ncol=1))     # Initialize oxygen consumption vector

  for (t in ti:(tf-1)) {

  # Compute Runge-Kutta increments

  # 1
  Tapp=Tint[t]
  PHYapp=Phyint[t]
  DTapp=DTint[t]
  POCapp=POCint[t]
  Ccontapp=Ccont[t]
  Ncontapp=Ncont[t]
  Pcontapp=Pcont[t]
  POMapp=POMint[t]
  TSSapp=TSSint[t]
  Napp=N[t]
  output<-Mussel_pop_equations(Param, Napp, Tapp, PHYapp, DTapp, POCapp, Ccontapp, Ncontapp, Pcontapp, POMapp, TSSapp, Wb[t], R[t],t,trip)
  dw=unlist(output[1])
  dr=unlist(output[2])
  k1=timestep*dw
  l1=timestep*dr

  # 2
  Tapp=approx(seq(from=1,to=tf,by=timestep),Tint,xout=(t+timestep/2))
  PHYapp=approx(seq(from=1,to=tf,by=timestep),Phyint,xout=(t+timestep/2))
  DTapp=approx(seq(from=1,to=tf,by=timestep),DTint,xout=(t+timestep/2))
  POCapp=approx(seq(from=1,to=tf,by=timestep),POCint,xout=(t+timestep/2))
  Ccontapp=approx(seq(from=1,to=tf,by=timestep),Ccont,xout=(t+timestep/2))
  Ncontapp=approx(seq(from=1,to=tf,by=timestep),Ncont,xout=(t+timestep/2))
  Pcontapp=approx(seq(from=1,to=tf,by=timestep),Pcont,xout=(t+timestep/2))
  POMapp=approx(seq(from=1,to=tf,by=timestep),POMint,xout=(t+timestep/2))
  TSSapp=approx(seq(from=1,to=tf,by=timestep),TSSint,xout=(t+timestep/2))
  Napp=approx(seq(from=1,to=tf+1,by=timestep),N,xout=(t+timestep/2))
  output<-Mussel_pop_equations(Param, Napp$y, Tapp$y, PHYapp$y, DTapp$y, POCapp$y, Ccontapp$y, Ncontapp$y, Pcontapp$y, POMapp$y, TSSapp$y, Wb[t]+k1/2, R[t]+l1/2,t,trip)
  dw=unlist(output[1])
  dr=unlist(output[2])
  k2=timestep*dw
  l2=timestep*dr

  # 3
  Tapp=approx(seq(from=1,to=tf,by=timestep),Tint,xout=(t+timestep/2))
  PHYapp=approx(seq(from=1,to=tf,by=timestep),Phyint,xout=(t+timestep/2))
  DTapp=approx(seq(from=1,to=tf,by=timestep),DTint,xout=(t+timestep/2))
  POCapp=approx(seq(from=1,to=tf,by=timestep),POCint,xout=(t+timestep/2))
  Ccontapp=approx(seq(from=1,to=tf,by=timestep),Ccont,xout=(t+timestep/2))
  Ncontapp=approx(seq(from=1,to=tf,by=timestep),Ncont,xout=(t+timestep/2))
  Pcontapp=approx(seq(from=1,to=tf,by=timestep),Pcont,xout=(t+timestep/2))
  POMapp=approx(seq(from=1,to=tf,by=timestep),POMint,xout=(t+timestep/2))
  TSSapp=approx(seq(from=1,to=tf,by=timestep),TSSint,xout=(t+timestep/2))
  Napp=approx(seq(from=1,to=tf+1,by=timestep),N,xout=(t+timestep/2))
  output<-Mussel_pop_equations(Param, Napp$y, Tapp$y, PHYapp$y, DTapp$y, POCapp$y, Ccontapp$y, Ncontapp$y, Pcontapp$y, POMapp$y, TSSapp$y, Wb[t]+k2/2, R[t]+l2/2,t,trip)
  dw=unlist(output[1])
  dr=unlist(output[2])
  k3=timestep*dw
  l3=timestep*dr

  # 4
  Tapp=Tint[t+timestep]
  PHYapp=Phyint[t+timestep]
  DTapp=DTint[t+timestep]
  POCapp=POCint[t+timestep]
  Ccontapp=Ccont[t+timestep]
  Ncontapp=Ncont[t+timestep]
  Pcontapp=Pcont[t+timestep]
  POMapp=POMint[t+timestep]
  TSSapp=TSSint[t+timestep]
  Napp=N[t+timestep]
  output<-Mussel_pop_equations(Param, Tapp, Napp, PHYapp, DTapp, POCapp, Ccontapp, Ncontapp, Pcontapp, POMapp, TSSapp, Wb[t]+k3, R[t]+l3,t,trip)
  dw=unlist(output[1])
  dr=unlist(output[2])
  k4=timestep*dw
  l4=timestep*dr

  # Compute weight at t+1 using Runge-Kutta increments

  # Somatic tissue weight
  Wb[t+timestep]=Wb[t]+(k1+2*k2+2*k3+k4)/6

  # Gonadic tissue weight
   R[t+timestep]=R[t]+(l1+2*l2+2*l3+l4)/6

  if (t==tspawn1) {       # If energy for eggs deposition is reached the 15/06
    R[t+timestep]=0
    tspawn1=tspawn1+365  # Updates the reproductory data to be used also the next year
    ripi=ripi+365      # updates the beginning of the resting period to be used the next year
    ripf=ripf+365      # updates the end of the resting period to be used the next year
  }

  if (t==tspawn2)   {      # If energy for eggs deposition is reached the 15/12
    R[t+timestep]=0
    tspawn2=tspawn2+365   # Updates the reproductory data to be used also the next year
  }

  trip=cbind(ripi,ripf)

  # Total dry weight
  Wd[t+timestep]=Wb[t+timestep]+R[t+timestep]

  # Mussel Length as a function of dry weight
  L[t+timestep]=max(L[t], a*Wd[t+timestep]^b, na.rm = FALSE)          # Mussel length at day t [cm]

  # Mussel wet weight as a function of dry weight [g]
  Wf[t+timestep]=aF*Wd[t+timestep]

  # Mussel total weight (with shell) as a function of dry weight [g]
  Wtot[t+timestep]=atot*Wd[t+timestep]

  # Compute the other outputs of the model
  output<-Mussel_pop_equations(Param, Napp, Tapp, PHYapp, DTapp, POCapp, Ccontapp, Ncontapp, Pcontapp, POMapp, TSSapp, Wb[t+timestep], R[t+timestep],t,trip)

  # Extracts outputs from the output list
  pseudofecies=output[[3]]
  fecies=output[[4]]
  composition=output[[5]]
  temperaturefun=output[[6]]
  metabolism=output[[7]]
  consumption=output[[8]]
  ammonium=output[[9]]


  # Outputs creation
  weight=rbind(Wb,R,Wd,Wtot,L)
  pfec=rbind(pfec,pseudofecies)
  fec=rbind(fec,fecies)
  comp=rbind(comp,composition)
  tfun=rbind(tfun,temperaturefun)
  metab=rbind(metab,metabolism)
  cons=rbind(cons,consumption)
  amm=rbind(amm,ammonium)

}  # Close cycle

  output=list(weight,pfec,fec,comp,tfun,metab,cons,amm)
  return(output)

} # Close function
