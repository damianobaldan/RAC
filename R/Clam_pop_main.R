#' Clam bioenergetic population model
#'
#' @param userpath the path where the working folder is located
#' @param forcings a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l], particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @return A list containing model outputs: weights, temperature limitation functions and metabolic rates
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'

Clam_pop_main<-function(userpath,forcings){

rm(list=ls())           # Clean workspace

cat('Clam population bioenergetic model\n')
cat(" \n")

# Run the preprocessor for the first time to print to screen parameters and forcing selected
out_pre<-Clam_pop_pre(userpath,forcings)

# While cycle to repeat the pre-processing until correct inputs are inserted
selector="y"

while (identical(selector,"y")=="TRUE") {
  cat(" \n")
  selector=readline("Do you want to change the inputs? [y/n]")

  if (identical(selector,"n")=="TRUE") {break}

  cat(" \n")
  cat("Insert forcings and parameters in the following folder\n")
  cat(paste0(userpath,"/Clam_population/Inputs\n"))
  cat(" \n")
  cat("Type y if you entered the correct inputs\n")
  cat("The data will be preprocessed again")
  selector=readline(" ")

  out_pre<-Clam_pop_pre(userpath,forcings)
  selector="y"
}

# Extract preprocessor outputs
Param=out_pre[[1]]
times=out_pre[[2]]
Dates=out_pre[[3]]
IC=out_pre[[4]]
Tint=out_pre[[5]]
Phyint=out_pre[[6]]
DTint=out_pre[[7]]
POCint=out_pre[[8]]
POMint=out_pre[[9]]
TSSint=out_pre[[10]]
N=out_pre[[11]]
CS=out_pre[[12]]

# Manages population
out_RKsolver<-Clam_pop_loop(Param, times, IC, Tint, Phyint, DTint, POCint, POMint, TSSint,N,userpath)

# Post-process data
out_post<-Clam_pop_post(userpath, out_RKsolver, times, Dates,N,CS)

cat(" ")
cat("End")

return(out_post)

}
