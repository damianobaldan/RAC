#' Postprocess the Clam indivual bioenergetic model (alternative version) results
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param CS the commercial size of Clam
#'
#' @return a list containing the clam weights, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

ClamF_ind_post<-function(userpath,output,times,Dates,CS) {

cat('Data post-processing\n')
cat('\n')

ti=times[1]           # Integration beginning
tf=times[2]           # Integration end

# Extracts outputs from the output list
weight=t(output[[1]])
Tfun=output[[2]]
metab=output[[3]]

# Adjusts results acoording with integration extremes
# now day 1 coincides with ti
weightSave=weight[ti:tf,]
TfunSave=Tfun[(ti+1):tf,]
metabSave=metab[(ti+1):tf,]

# Days to commercial size
foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(weightSave[,3])
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  daysToSize="Not reaching the commercial size"
}else{  daysToSize <- min(NonNAindex)
}
daysToSize<-as.list(daysToSize)

output=list(weightSave,TfunSave,metabSave,daysToSize)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1) # create a dates vector to plot results
days2 <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti) # create a dates vector to plot results

# Plot weight
filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_plots//wet_weight.jpeg")
jpeg(filepath,800,600)
plot(days,weightSave[,2],ylab="Wet weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot length
filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_plots//length.jpeg")
jpeg(filepath,800,600)
plot(days,weightSave[,3],ylab="length (mm)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_plots//temperature_response.jpeg")
jpeg(filepath,800,600)
ub=max(max(TfunSave[,1]),max(TfunSave[,2]))
plot(days2,TfunSave[,1],ylab="Temperature response function",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,TfunSave[,2],col="blue")
lines(days2,TfunSave[,3],col="green")
legend("topright",c("Growth limitation","Respiration limitation","Filtration limitation"),fill=c("red","blue","green"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
ub=max(max(metabSave[,1]),max(metabSave[,2]))
plot(days2,metabSave[,1],ylab="Metabolic rate (J/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,metabSave[,2],col="blue")
legend("topright",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_csv//biometries.csv")
write.csv(weightSave,filepath)

filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_csv//temperature_response.csv")
write.csv(TfunSave,filepath)

filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_csv//metabolism.csv")
write.csv(metabSave,filepath)

filepath=paste0(userpath,"/ClamF_individual/Outputs/Out_csv//Days_to_commercial_size.csv")
write.csv(daysToSize,filepath)

return(output)
}
