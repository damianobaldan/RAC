#' Clam bioenergetic individual model (alternative version)
#'
#' @param userpath the path where the working folder is located
#' @param forcings a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3]
#' @return A list containing model outputs: weights, temperature limitation functions and metabolic rates
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'

ClamF_ind_main<-function(userpath,forcings) {

rm(list=ls())       # Clean workspace

cat('Clam individual based bioenergetic model - alternative version\n')
cat(" \n")

# Run the preprocessor for the first time to print to screen parameters and forcing selected
out_pre<-ClamF_ind_pre(userpath,forcings)

# While cycle to repeat the pre-processing until correct inputs are inserted
selector="y"

while (identical(selector,"y")=="TRUE") {
cat(" \n")
selector=readline("Do you want to change the inputs? [y/n]")

if (identical(selector,"n")=="TRUE") {break}

cat(" \n")
cat("Insert forcings and parameters in the following folder\n")
cat(paste0(userpath,"/ClamF_individual/Inputs\n"))
cat(" \n")
cat("Type y if you entered the correct inputs\n")
cat("The data will be preprocessed again")
selector=readline(" ")

out_pre<-ClamF_ind_pre(userpath,forcings)
selector="y"
}

# Extract preprocessor outputs AGGIUSTARE
Param=out_pre[[1]]
times=out_pre[[2]]
Dates=out_pre[[3]]
IC=out_pre[[4]]
Tint=out_pre[[5]]
Chlint=out_pre[[6]]
CS=out_pre[[7]]

# Solves ODE
output<-ClamF_ind_RKsolver(Param, times, IC, Tint, Chlint)

# Post-process data
out_post<-ClamF_ind_post(userpath, output, times, Dates,CS)

cat(" ")
cat("End")

return(out_post)

}
