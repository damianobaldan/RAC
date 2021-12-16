#' Postprocess the Mussel spatialized model results
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param CS the commercial size of Mussel
#'
#' @return a list containing the weights of the mussel, the excreted CNP, the mussel CNP, temperature limitation functions, metabolic rates, oxygen consumption
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Mussel_spatial_post<-function(userpath,output,times,Dates,CS) {

ti=times[1]           # Integration beginning
tf=times[2]           # Integration end

# Extracts outputs from the output list
W=output[[1]]
pfec=output[[2]]
fec=output[[3]]
comp=output[[4]]
tfun=output[[5]]
metab=output[[6]]
cons=output[[7]]
amm=output[[8]]

# Adjusts results acoording with integration extremes
# now day 1 coincides with ti
weightSave=t(W[,ti:tf])
pfecSave=pfec[ti:tf,]
fecSave=fec[ti:tf,]
compSave=comp[ti:tf,]
tfunSave=tfun[ti:tf,]
metabSave=metab[ti:tf,]
consSave=cons[ti:tf]
ammSave=amm[ti:tf]

foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(weightSave)
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  daysToSize=99999
}else{  daysToSize <- min(NonNAindex)
}
daysToSize<-as.list(daysToSize)

output=list(weightSave,pfecSave,fecSave,compSave,tfunSave,metabSave,consSave,ammSave,daysToSize)

return(output)

}
