#' Solves the Clam bioenergetic balance for population with a 4th order Runge Kutta method
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
#'
#' @return a list containing the clam weights, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats
#'

Clam_pop_RKsolver <- function(Param, times, IC, Tint, Phyint, DTint, POCint, POMint, TSSint){

  # Integration extremes definition
  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end
  timestep=times[3]     # Timestep for integration

  # Initial condition definition
  aF=Param[18]             # [gWW] Dry weight - wet weight conversion coefficient
  bF=Param[19]             # [-] Dry weight - wet weight exponent
  a=Param[20]              # [m] Dry weight - length conversion coefficient
  b=Param[21]              # [-] Weight to length exponent

  weight=as.vector(matrix(0,nrow=ti))      # Initialize vector dry weight
  weight[ti]=IC                            # Dry weight initial condition [g]
  Ww=as.vector(matrix(0,nrow=ti))          # Initialize vector total dry weight
  Ww[ti]=(weight[ti]/aF)^(1/bF)            # total dry weight initial value [g]
  L=as.vector(matrix(0,nrow=ti))           # Initialize vector length
  L[ti]=(Ww[ti]/a)^b                   # length initial value [mm]


  # initialize outputs
  tfun=as.matrix(matrix(0,nrow=ti,ncol=2))     # Initialize temperature limitations vector
  metab=as.matrix(matrix(0,nrow=ti,ncol=2))    # Initialize metabolic rates vector

  for (t in ti:(tf-1)) {

  # Compute Runge-Kutta increments

  #1
  Tapp=Tint[t]
  PHYapp=Phyint[t]
  DTapp=DTint[t]
  POCapp=POCint[t]
  POMapp=POMint[t]
  TSSapp=TSSint[t]
  output<-Clam_pop_equations(Param, Tapp, PHYapp, DTapp, POCapp, POMapp, TSSapp, weight[t])
  dw=unlist(output[1])
  k1=timestep*dw

  #2
  Tapp=approx(seq(from=1,to=tf,by=timestep),Tint,xout=(t+timestep/2))
  PHYapp=approx(seq(from=1,to=tf,by=timestep),Phyint,xout=(t+timestep/2))
  DTapp=approx(seq(from=1,to=tf,by=timestep),DTint,xout=(t+timestep/2))
  POCapp=approx(seq(from=1,to=tf,by=timestep),POCint,xout=(t+timestep/2))
  POMapp=approx(seq(from=1,to=tf,by=timestep),POMint,xout=(t+timestep/2))
  TSSapp=approx(seq(from=1,to=tf,by=timestep),TSSint,xout=(t+timestep/2))
  output<-Clam_pop_equations(Param, Tapp$y, PHYapp$y, DTapp$y, POCapp$y, POMapp$y, TSSapp$y, weight[t]+k1/2)
  dw=unlist(output[1])
  k2=timestep*dw

  #3
  Tapp=approx(seq(from=1,to=tf,by=timestep),Tint,xout=(t+timestep/2))
  PHYapp=approx(seq(from=1,to=tf,by=timestep),Phyint,xout=(t+timestep/2))
  DTapp=approx(seq(from=1,to=tf,by=timestep),DTint,xout=(t+timestep/2))
  POCapp=approx(seq(from=1,to=tf,by=timestep),POCint,xout=(t+timestep/2))
  POMapp=approx(seq(from=1,to=tf,by=timestep),POMint,xout=(t+timestep/2))
  TSSapp=approx(seq(from=1,to=tf,by=timestep),TSSint,xout=(t+timestep/2))
  output<-Clam_pop_equations(Param, Tapp$y, PHYapp$y, DTapp$y, POCapp$y, POMapp$y, TSSapp$y, weight[t]+k2/2)
  dw=unlist(output[1])
  k3=timestep*dw

  #4
  Tapp=Tint[t+timestep]
  PHYapp=Phyint[t+timestep]
  DTapp=DTint[t+timestep]
  POMapp=POMint[t+timestep]
  TSSapp=TSSint[t+timestep]
  output<-Clam_pop_equations(Param, Tapp, PHYapp, DTapp, POCapp, POMapp, TSSapp, weight[t]+k3)
  dw=unlist(output[1])
  k4=timestep*dw

  # Compute weight at t+1 using Runge-Kutta increments
  weight[t+timestep]=weight[t]+(k1+2*k2+2*k3+k4)/6   # Dry weight [gDW]
  Ww[t+timestep]=(weight[t+timestep]/aF)^(1/bF)      # Wet weight [gWW]
  L[t+timestep]=max(L[t-timestep],(Ww[t]/a)^b)       # Mussel's length [mm]

  # Compute the other outputs of the model
  output<-Clam_pop_equations(Param, Tint[t+timestep], Phyint[t+timestep], DTint[t+timestep], POCint[t+timestep], POMint[t+timestep], TSSint[t+timestep] , weight[t+timestep])

  # Extracts outputs from the output list
  temperaturefun=output[[2]]
  metabolism=output[[3]]

  # Outputs creation
  w=rbind(weight,Ww,L)
  tfun=rbind(tfun, temperaturefun)
  metab=rbind(metab, metabolism)

}  # Close cycle

  output=list(w,tfun,metab)
  return(output)

} # Close function
