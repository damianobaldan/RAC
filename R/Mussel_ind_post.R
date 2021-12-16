#' Postprocess the Mussel indivual bioenergetic model results
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

Mussel_ind_post<-function(userpath,output,times,Dates,CS) {

cat('Data post-processing\n')
cat('\n')

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
pfecSave=pfec[(ti+1):tf,]
fecSave=fec[(ti+1):tf,]
compSave=comp[(ti+1):tf,]
tfunSave=tfun[(ti+1):tf,]
metabSave=metab[(ti+1):tf,]
consSave=cons[(ti+1):tf]
ammSave=amm[(ti+1):tf]

foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(weightSave)
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  daysToSize="Not reaching the commercial size"
}else{  daysToSize <- min(NonNAindex)
}
daysToSize<-as.list(daysToSize)

output=list(weightSave,pfecSave,fecSave,compSave,tfunSave,metabSave,consSave,ammSave,daysToSize)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1) # create a dates vector to plot results
days2 <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti)

currentpath=getwd()

# Plot weight
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//dry_weight.jpeg")
jpeg(filepath,800,600)
plot(days,weightSave[,1],ylab="Dry weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red")
lines(days,weightSave[,2],col="blue")
lines(days,weightSave[,3],col="black")
legend("topleft",c("Somatic tissue","Gonadic tissue","Total"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot length
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//length.jpeg")
jpeg(filepath,800,600)
plot(days,weightSave[,5],ylab="Length (cm)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot pseudofecies
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//pseudofaeces_production.jpeg")
jpeg(filepath,800,600)
ub=max(max(pfecSave[,1]),max(pfecSave[,2]),max(pfecSave[,3]))
plot(days2,pfecSave[,1],ylab="Pseudofaeces production (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,pfecSave[,2],col="blue")
lines(days2,pfecSave[,3],col="black")
legend("topleft",c("Excreted C","Excreted N","excreted P"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot fecies
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//faeces_production.jpeg")
jpeg(filepath,800,600)
ub=max(max(fecSave[,1]),max(fecSave[,2]),max(fecSave[,3]))
plot(days2,fecSave[,1],ylab="Faeces production (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,fecSave[,2],col="blue")
lines(days2,fecSave[,3],col="black")
legend("topleft",c("Excreted C","Excreted N","Excreted P"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot CNP mytilus
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//CNP_content.jpeg")
jpeg(filepath,800,600)
ub=max(max(compSave[,1]),max(compSave[,2]),max(compSave[,3]))
plot(days2,compSave[,1],ylab="CNP mytilus (g)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,compSave[,2],col="blue")
lines(days2,compSave[,3],col="black")
legend("topleft",c("C","N","P"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//temperature_response.jpeg")
jpeg(filepath,800,600)
ub=max(max(tfunSave[,1]),max(tfunSave[,2]))
plot(days2,tfunSave[,1],ylab="Temperature response function",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,tfunSave[,2],col="blue")
legend("topright",c("Anabolism limitation","Catabolism limitation"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
ub=max(max(metabSave[,1]),max(metabSave[,2]))
plot(days2,metabSave[,1],ylab="Metabolic rate (J/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,metabSave[,2],col="blue")
legend("topright",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot O2 consumption
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//O2_consumption.jpeg")
jpeg(filepath,800,600)
plot(days2,consSave,ylab="O2 consumption (g/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot ammonium release
filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_plots//NH4_release.jpeg")
jpeg(filepath,800,600)
plot(days2,ammSave,ylab="NH4 release (g/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//biometries.csv")
write.csv(weightSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//faeces_production.csv")
write.csv(fecSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//pseudofaeces_production.csv")
write.csv(pfecSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//CNP_content.csv")
write.csv(compSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//O2_consumption.csv")
write.csv(consSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//NH4_release.csv")
write.csv(ammSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//temperature_response.csv")
write.csv(tfunSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//metabolism.csv")
write.csv(metabSave,filepath)

filepath=paste0(userpath,"/Mussel_individual/Outputs/Out_csv//Days_to_commercial_size.csv")
write.csv(daysToSize,filepath)

return(output)

}
