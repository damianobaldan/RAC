#' Function that solves the population dynamics equations including discontinuities
#'
#' @param Nseed number of seeded individuals
#' @param mort mortality rate
#' @param manag list of management actions (seeded/harvested individuals)
#' @param times vector containing informations on integration times
#' @return a time series with the number of individuals
#'
#' @import matrixStats plotrix rstudioapi
#'

Pop_fun <- function(Nseed, mort, manag, times) {

  # Integration times
  ti=times[1]
  tf=times[2]
  timestep=times[3]

  # Initial condition and vectors initialization
  N=as.vector(matrix(0,nrow=ti))  # Initialize vector N
  N[ti]=Nseed                     # Impose initial condition
  dN=as.vector(matrix(0,nrow=ti)) # Initialize vector dN

  for (t in ti:tf){ # for cycle that solves population ODE with Euler method

  dN[t]=-mort*N[t]             # individuals increment
  N[t+1]=N[t]+dN[t]*timestep   # Individuals at time t+1

  for (i in 1:length(manag[,1])) {  # For cycle that adjusts N according with management strategies
  if (t==manag[i,1]) {              # if statement to check if it is the tiome to adjust N
    N[t+1]=N[t]+manag[i,2]

  } # close if
  } # close for
  } # close for

  output=N
  return(output)
  } # close function
