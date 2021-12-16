#' Mussel bioenergetic spatialized model - spatialization loop
#'
#' Solves the bioenergetic balance for Mussel
#'
#' @param userpath the path where the working folder is located
#' @param forcings a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @return A list containing model outputs that main script saves to .nc; .csv and .asc files
#'
#' @import matrixStats plotrix rstudioapi
#'

Mussel_spatial_loop<-function(userpath,forcings){

 # rm(list=ls())       # Clean workspace

    out_pre<-Mussel_spatial_pre_int(userpath,forcings)

    # Extract preprocessor outputs
    Param=out_pre[[1]]
    times=out_pre[[2]]
    Dates=out_pre[[3]]
    IC=out_pre[[4]]
    Tint=out_pre[[5]]
    Phyint=out_pre[[6]]
    DTint=out_pre[[7]]
    POCint=out_pre[[8]]
    Ccont=out_pre[[9]]
    Ncont=out_pre[[10]]
    Pcont=out_pre[[11]]
    POMint=out_pre[[12]]
    TSSint=out_pre[[13]]
    CS=out_pre[[14]]

    # Solves ODE
    out_RKsolver<-Mussel_spatial_RKsolver(Param, times, IC, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint)

    # Post-process data
    out_post<-Mussel_spatial_post(userpath, out_RKsolver, times, Dates,CS)

  return(out_post)

}
