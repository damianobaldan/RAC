#' Solves the Clam bioenergetic balance (alternative version) with a 4th order Runge Kutta method
#'
#' @param Param a vector containing model parameters
#' @param times integration extremes and integration timestep
#' @param IC initial condition
#' @param Tint the interpolated water temperature time series
#' @param Chlint the interpolated chlorophyll a time series
#'
#' @return a list containing the clam weights, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats
#'

ClamF_ind_RKsolver <- function(Param, times, IC, Tint, Chlint){


  cat("ODE solution\n")

  # Integration extremes definition
  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end
  timestep=times[3]     # Timestep for integration

  # Initial condition definition
  p=Param[5]               # [-] Dry weight - wet weight conversion exponent
  k=Param[6]               # [-] Wet weight - length conversion exponent
  a=Param[8]               # [mm] Wet weight - length conversion coefficient
  b=Param[9]               # [gWW] Dry weight - wet weight conversion coefficient

  Ww=as.vector(matrix(0,nrow=ti))          # Initialize vector wet weight
  Ww[ti]=IC                                # Wet weight initial condition [g]
  Wd=as.vector(matrix(0,nrow=ti))          # Initialize vector dry weight
  Wd[ti]=b*(Ww[ti])^p                      # Dry weight initial value [g]
  L=as.vector(matrix(0,nrow=ti))           # Initialize vector length
  L[ti]=(Ww[ti]/a)^k                       # Length of the mussel initial value [cm]


  # initialize outputs
  tfun=as.matrix(matrix(0,nrow=ti,ncol=3))     # Initialize temperature limitations vector
  metab=as.matrix(matrix(0,nrow=ti,ncol=2))    # Initialize metabolic rates vector

  for (t in ti:(tf-1)) {

  # Compute Runge-Kutta increments

  # 1
  Tapp=Tint[t]
  Chlapp=Chlint[t]
  output<-ClamF_ind_equations(Param, Tapp, Chlapp, Ww[t])
  dw=unlist(output[1])
  k1=timestep*dw

  # 2
  Tapp=approx(seq(from=1,to=tf,by=timestep),Tint,xout=(t+timestep/2))
  Chlapp=approx(seq(from=1,to=tf,by=timestep),Chlint,xout=(t+timestep/2))
  output<-ClamF_ind_equations(Param, Tapp$y, Chlapp$y, Ww[t]+k1/2)
  dw=unlist(output[1])
  k2=timestep*dw

  # 3
  Tapp=approx(seq(from=1,to=tf,by=timestep),Tint,xout=(t+timestep/2))
  Chlapp=approx(seq(from=1,to=tf,by=timestep),Chlint,xout=(t+timestep/2))
  output<-ClamF_ind_equations(Param, Tapp$y, Chlapp$y, Ww[t]+k2/2)
  dw=unlist(output[1])
  k3=timestep*dw

  # 4
  Tapp=Tint[t+timestep]
  Chlapp=Chlint[t+timestep]
  output<-ClamF_ind_equations(Param, Tapp, Chlapp, Ww[t]+k3)
  dw=unlist(output[1])
  k4=timestep*dw

  # Compute weight at t+1 using Runge-Kutta increments
  Ww[t+timestep]=Ww[t]+(k1+2*k2+2*k3+k4)/6      # Dry weight [gDW]
  Wd[t+timestep]=b*Ww[t+timestep]^p             # Wet weight [gWW]
  L[t+timestep]=max(L[t-timestep],(Ww[t]/a)^k)  # Mussel's length [mm]

  # Compute the other outputs of the model
  output<-ClamF_ind_equations(Param, Tint[t+timestep], Chlint[t+timestep], Ww[t+timestep])

  # Extracts outputs from the output list
  temperaturefun=output[[2]]
  metabolism=output[[3]]

  # Outputs creation
  w=rbind(Ww,Wd,L)
  tfun=rbind(tfun, temperaturefun)
  metab=rbind(metab, metabolism)

}  # Close cycle

  output=list(w,tfun,metab)
  return(output) # ClamF_ind_RKsolver output

} # Close function
