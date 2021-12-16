#' Mussel bioenergetic individual model
#'
#' Solves the bioenergetic balance for Mussel
#'
#' @param userpath the path where the working folder is located
#' @param forcings a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @return A list containing model outputs: weight, length mussel CNP, pseudofaeces CNP production, faeces CNP production, temperature limitation functions, metabolic rates and oxygen consumption
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'

Mussel_ind_main<-function(userpath,forcings){

rm(list=ls())       # Clean workspace

cat('Mussel bioenergetic individual based model\n')
cat(" \n")

# Run the preprocessor for the first time to print to screen parameters and forcing selected
out_pre<-Mussel_ind_pre(userpath,forcings)

# While cycle to repeat the pre-processing until correct inputs are inserted
selector="y"

while (identical(selector,"y")=="TRUE") {
cat(" \n")
selector=readline("Do you want to change the inputs? [y/n]")

if (identical(selector,"n")=="TRUE") {break}

cat(" \n")
cat("Insert forcings and parameters in the following folder\n")
cat(paste0(userpath,"/Mussel_individual/Inputs\n"))
cat(" \n")
cat("Type y if you entered the correct inputs\n")
cat("The data will be preprocessed again")
selector=readline(" ")

out_pre<-Mussel_ind_pre(userpath,forcings)
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
Ccont=out_pre[[9]]
Ncont=out_pre[[10]]
Pcont=out_pre[[11]]
POMint=out_pre[[12]]
TSSint=out_pre[[13]]
CS=out_pre[[14]]

# Solves ODE
out_RKsolver<-Mussel_ind_RKsolver(Param, times, IC, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint)

# Post-process data
out_post<-Mussel_ind_post(userpath, out_RKsolver, times, Dates,CS)

cat(" ")
cat("End")

return(out_post)

}

