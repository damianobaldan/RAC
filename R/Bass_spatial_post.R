#' Postprocess the Mussel spatialized model results
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param CS the commercial size of Bass
#'
#' @return a list containing the fish weight, proteines, lipids and carbohydrates wasted or produced with excretions, potential and actual ingestion rates, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Bass_spatial_post<-function(userpath,output,times,Dates,CS) {

  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end

  # Extracts outputs from the output list
  weight=unlist(output[1])
  exc=output[[2]]
  wst=output[[3]]
  ing=unlist(output[4])
  ingvero=unlist(output[5])
  Tfun=output[[6]]
  metab=output[[7]]
  O2=output[[8]]
  NH4=output[[9]]

  # Adjusts results acoording with integration extremes
  # now day 1 coincides with ti
  weightSave=weight[ti:tf]
  excSave=exc[ti:tf,]
  wstSave=wst[ti:tf,]
  ingSave=ing[ti:tf]
  ingveroSave=ingvero[ti:tf]
  TfunSave=Tfun[ti:tf,]
  metabSave=metab[ti:tf,]
  O2Save=O2[ti:tf]
  NH4Save=NH4[ti:tf]

  # Days to commercial size
  foo <- function(w,S){which(w>S)[1]}
  arg=as.data.frame(weightSave)
  days <- apply(arg,1,foo,S=CS)
  days_L <- as.data.frame(days)
  NonNAindex <- which(!is.na(days_L))
  if (length(NonNAindex)==0) {
    daysToSize=9999
  }else{  daysToSize <- min(NonNAindex)
  }
  # daysToSize<-as.list(daysToSize)

  output=list(weightSave,excSave,wstSave,ingSave,ingveroSave,metabSave,TfunSave, O2Save, NH4Save, daysToSize)

  return(output)

}
