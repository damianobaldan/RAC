#' Solves the Seabass bioenergetic balance with a 4th order Runge Kutta method - used in spatialized model
#'
#' @param Param vector containing all metabolic parameters
#' @param Temperature water temperature forcing time series
#' @param G food entering the cage time series
#' @param Food food characterization (Proteins, Lipids, Carbohydrates)
#' @param IC initial conditions
#' @param times vector containing integration extremes and integration timestep
#'
#' @return a list containing the fish weight, proteines, lipids and carbohydrates wasted or produced with excretions, potential and actual ingestion rates, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats
#'

Bass_spatial_RKsolver <- function(Param, Temperature, G, Food, IC, times){


  # Integration extremes definition
  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end
  timestep=times[3]     # Timestep for integration

  # Initial condition definition
  weight=as.vector(matrix(0,nrow=ti))      # Initialize vector w
  weight[ti]=IC                            # Define initial condition

  # initialize outputs
  wst=as.matrix(matrix(0,nrow=ti,ncol=3))      # Initialize food to waste vector
  exc=as.matrix(matrix(0,nrow=ti,ncol=3))      # Initialize excretion vector
  ing=as.vector(matrix(0,nrow=ti))             # Initialize potential ingestion vector
  ingvero=as.vector(matrix(0,nrow=ti))         # Initialize actual ingestion vector
  tfun=as.matrix(matrix(0,nrow=ti,ncol=2))     # Initialize temperature limitations vector
  metab=as.matrix(matrix(0,nrow=ti,ncol=2))    # Initialize metabolic rates vector
  O2=as.matrix(matrix(0,nrow=ti,ncol=2))       # Initialize oxygen consumption rates vector
  NH4=as.matrix(matrix(0,nrow=ti,ncol=2))      # Initialize ammonia release rates vector


  for (t in ti:(tf-1)) {

    # Compute Runge-Kutta increments

    # 1
    Tapp=Temperature[t]
    Gapp=G[t]
    output<-Bass_ind_equations(Param, Tapp, Gapp, Food, weight[t])
    dw=unlist(output[1])
    k1=timestep*dw

    # 2
    Tapp=approx(seq(from=1,to=tf,by=timestep),Temperature,xout=(t+timestep/2))
    Gapp=approx(seq(from=1,to=tf,by=timestep),G,xout=(t+timestep/2))
    output<-Bass_ind_equations(Param, Tapp$y, Gapp$y, Food, weight[t]+k1/2)
    dw=unlist(output[1])
    k2=timestep*dw;

    # 3
    Tapp=approx(seq(from=1,to=tf,by=timestep),Temperature,xout=(t+timestep/2))
    Gapp=approx(seq(from=1,to=tf,by=timestep),G,xout=(t+timestep/2))
    output<-Bass_ind_equations(Param, Tapp$y, Gapp$y, Food, weight[t]+k2/2)
    dw=unlist(output[1])
    k3=timestep*dw;

    # 4
    Tapp=Temperature[t+timestep]
    Gapp=G[t+timestep]
    output<-Bass_ind_equations(Param, Tapp, Gapp, Food, weight[t]+k3)
    dw=unlist(output[1])
    k4=timestep*dw;

    # Compute weight at t+1 using Runge-Kutta increments
    weight[t+timestep]=weight[t]+(k1+2*k2+2*k3+k4)/6 # [g]

    # Compute the other outputs of the model
    output<-Bass_ind_equations(Param, Temperature[t+timestep], G[t+timestep], Food, weight[t+timestep])

    # Extracts outputs from the output list
    excretion=output[[2]]
    waste=output[[3]]
    ingestion=unlist(output[4])
    ingestionvero=unlist(output[5])
    temperaturefun=output[[6]]
    metabolism=output[[7]]
    oxygen=output[[8]]
    ammonia=output[[9]]

    # Outputs creation
    wst=rbind(wst, waste)
    exc=rbind(exc, excretion)
    ing=rbind(ing, ingestion)
    ingvero=rbind(ingvero, ingestionvero)
    tfun=rbind(tfun, temperaturefun)
    metab=rbind(metab, metabolism)
    O2=rbind(O2, oxygen)
    NH4=rbind(NH4, ammonia)

  }  # Close cycle

  output=list(weight,exc,wst,ing,ingvero,tfun,metab, O2, NH4)
  return(output) # Bass_ind_RKsolver output

} # Close function

