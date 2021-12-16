#' Seabass bioenergetic individual model postprocessor
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param CS the commercial size of Seabass
#'
#' @return a list containing the fish weight, proteines, lipids and carbohydrates wasted or produced with excretions, potential and actual ingestion rates, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Bass_ind_post<-function(userpath,output,times,Dates,CS) {

cat('Data post-processing\n')
cat('\n')

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
weightSave=weight[ti:(tf-2)]
excSave=exc[(ti+1):(tf-2),]
wstSave=wst[(ti+1):(tf-2),]
ingSave=ing[(ti+1):(tf-2)]
ingveroSave=ingvero[(ti+1):(tf-2)]
TfunSave=Tfun[(ti+1):(tf-2),]
metabSave=metab[(ti+1):(tf-2),]
O2Save=O2[(ti+1):(tf-2)]
NH4Save=NH4[(ti+1):(tf-2)]

# Days to commercial size
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

output=list(weightSave,excSave,ingSave,ingveroSave,wstSave,metabSave,TfunSave, O2Save, NH4Save, daysToSize)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti-2) # create a dates vector to plot results
days2 <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti-1) # create a dates vector to plot results

# Plot weight
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//weight.jpeg")
jpeg(filepath,800,600)
plot(days2,weightSave,ylab="Weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days2, 1), by = "months")
axis.Date(side = 1, days2, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot excretion
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//faeces_production.jpeg")
jpeg(filepath,800,600)
ub=max(max(excSave[,1]),max(excSave[,2]),max(excSave[,3]))
plot(days,excSave[,1],ylab="Faeces production (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,excSave[,2],col="blue")
lines(days,excSave[,3],col="black")
legend("topleft",c("Proteins","Lipids","Carbohydrates"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot wasted food
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//wasted_feed.jpeg")
jpeg(filepath,800,600)
ub=max(max(wstSave[,1]),max(wstSave[,2]),max(wstSave[,3]))
plot(days,wstSave[,1],ylab="Wasted feed (g/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,wstSave[,2],col="blue")
lines(days,wstSave[,3],col="black")
legend("topleft",c("Proteins","Lipids","Carbohydrates"),fill=c("red","blue","black"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot ingested food
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//actual_ingestion.jpeg")
jpeg(filepath,800,600)
plot(days,ingveroSave,ylab="Ingested food (g)",xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//temperature_response.jpeg")
jpeg(filepath,800,600)
ub=max(max(TfunSave[,1]),max(TfunSave[,2]))
plot(days,TfunSave[,1],ylab="Temperature response function",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,TfunSave[,2],col="blue")
legend("topright",c("Anabolism","Catabolism"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
ub=max(max(metabSave[,1]),max(metabSave[,2]))
plot(days,metabSave[,1],ylab="Metabolic rate (J/d)",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days,metabSave[,2],col="blue")
legend("topright",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot O2 consumed
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//O2_consumption.jpeg")
jpeg(filepath,800,600)
plot(days,O2Save,ylab="Oxygen consumption (gO2/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot NH4 produced
filepath=paste0(userpath,"/Bass_individual/Outputs/Out_plots//NH4_release.jpeg")
jpeg(filepath,800,600)
plot(days,NH4Save,ylab="NH4 release (gN/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//weight.csv")
write.csv(weightSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//faeces_production.csv")
write.csv(excSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//wasted_feed.csv")
write.csv(wstSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//potential_ingestion.csv")
write.csv(ingSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//actual_ingestion.csv")
write.csv(ingveroSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//temperature_response.csv")
write.csv(TfunSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//metabolism.csv")
write.csv(metabSave,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//O2_consumption.csv")
write.csv(O2Save,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//NH4_release.csv")
write.csv(NH4Save,filepath)

filepath=paste0(userpath,"/Bass_individual/Outputs/Out_csv//Days_to_commercial_size.csv")
write.csv(daysToSize,filepath)

return(output)
}

