#' Bass bioenergetic spatialized model - spatialization loop
#'
#' Solves the bioenergetic balance for Bass
#'
#' @param userpath the path where the working folder is located
#' @param forcings a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @return a list containing the outputs that main script saves to .nc; .csv and .asc files
#'
#' @import matrixStats plotrix rstudioapi
#'

Bass_spatial_loop<-function(userpath,forcings){

#rm(list=ls())       # Clean workspace

    out_pre<-Bass_spatial_pre_int(userpath,forcings)

    # Extract preprocessor outputs
    Param=out_pre[[1]]
    Tint=out_pre[[2]]
    Gint=out_pre[[3]]
    Food=out_pre[[4]]
    IC=out_pre[[5]]
    times=out_pre[[6]]
    Dates=out_pre[[7]]
    CS=out_pre[[8]]

    # Solves ODE
    out_RKsolver<-Bass_spatial_RKsolver(Param, Tint, Gint, Food, IC, times)

    # Post-process data
    out_post<-Bass_spatial_post(userpath,out_RKsolver, times, Dates,CS)

  return(out_post)

}
