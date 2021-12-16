#' Seabream bioenergetic population model
#'
#' @param userpath the path where the working folder is located
#' @param forcings a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees] and feeding rate [g/individual x d]
#' @return A list containing model outputs: weight, excreted quantities and quantities to waste, actual and potential ingestion, temperature limitation functions and metabolic rates
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'

Bream_pop_main<-function(userpath,forcings){

rm(list=ls())           # Clean workspace

cat('Sea Bream population bioenergetic model\n')
cat(" \n")

currentpath=getwd() # Save current working directory path

# Run the preprocessor for the first time to print to screen parameters and forcing selected
out_pre<-Bream_pop_pre(userpath,forcings)

# While cycle to repeat the pre-processing until correct inputs are inserted
selector="y"

while (identical(selector,"y")=="TRUE") {
cat(" \n")
selector=readline("Do you want to change the inputs? [y/n]")

if (identical(selector,"n")=="TRUE") {break}

cat(" \n")
cat("Insert forcings and parameters in the following folder\n")
cat(paste0(userpath,"/Bream_population/Inputs\n"))
cat(" \n")
cat("Type y if you entered the correct inputs\n")
cat("The data will be preprocessed again")
selector=readline(" ")

out_pre<-Bream_pop_pre(userpath,forcings)
selector="y"
}

# Extract preprocessor outputs
Param=out_pre[[1]]
Tint=out_pre[[2]]
Gint=out_pre[[3]]
Food=out_pre[[4]]
IC=out_pre[[5]]
times=out_pre[[6]]
Dates=out_pre[[7]]
N=out_pre[[8]]
CS=out_pre[[9]]

# Manages population
out_RKsolver<-Bream_pop_loop(Param, Tint, Gint, Food, IC, times, N, userpath)

# Post-process data
out_post<-Bream_pop_post(userpath, out_RKsolver, times, Dates,N,CS)

cat(" ")
cat("End")

return(out_post)


}
