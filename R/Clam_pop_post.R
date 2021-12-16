#' Postprocess the Clam population bioenergetic model results
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param N the number of individuals
#' @param CS the commercial size of Clam
#'
#' @return  a list containing the clam weights, temperature limitation functions and metabolic rates
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import grDevices graphics utils stats
#'

Clam_pop_post<-function(userpath,output,times,Dates,N, CS) {

  ti=times[1]           # Integration beginning
  tf=times[2]           # Integration end

  # Extracts results from output list
  Wd_stat=output[[1]]
  Ww_stat=output[[2]]
  L_stat=output[[3]]
  A_stat=output[[4]]
  C_stat=output[[5]]
  fgT=output[[6]]
  frT=output[[7]]

cat('Data post-processing\n')
cat('\n')

# Adjusts results acoording with integration extremes
# now day 1 coincides with ti

WdSave=Wd_stat[,ti:tf]
WwSave=Ww_stat[,ti:tf]
LSave=L_stat[,ti:tf]

ASave=A_stat[,ti:tf]
CSave=C_stat[,ti:tf]

fgT=fgT[(ti+1):tf]
frT=frT[(ti+1):tf]

tfunSave=cbind(fgT,frT)

N=N[ti:tf]

# Days to commercial size

# Lower bound
foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(LSave[1,]-LSave[2,])
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  Lb_daysToSize="Not reaching the commercial size"
}else{  Lb_daysToSize <- min(NonNAindex)
}

# Mean
foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(LSave[1,])
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  Mean_daysToSize="Not reaching the commercial size"
}else{  Mean_daysToSize <- min(NonNAindex)
}

# Upper bound
foo <- function(w,S){which(w>S)[1]}
arg=as.data.frame(LSave[1,]+LSave[2,])
days <- apply(arg,1,foo,S=CS)
days_L <- as.data.frame(days)
NonNAindex <- which(!is.na(days_L))
if (length(NonNAindex)==0) {
  Ub_daysToSize="Not reaching the commercial size"
}else{  Ub_daysToSize <- min(NonNAindex)
}

# List containing days to size
daysToSize<-as.list(cbind(Ub_daysToSize,Mean_daysToSize,Lb_daysToSize))

output=list(WdSave,WwSave,LSave,ASave,CSave,fgT,frT,N,daysToSize)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1) # create a dates vector to plot results
currentpath=getwd()

# Plot dry weight
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//dry_weight.jpeg")
jpeg(filepath,800,600)
ub=WdSave[1,]+WdSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(WdSave[1,]-WdSave[2,])){
lb[i]=max(WdSave[1,i]-WdSave[2,i],0)
}
maxub=max(WdSave[1,]+WdSave[2,])
plot(days,WdSave[1,],ylab="Dry weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,WdSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot wet weight
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//wet_weight.jpeg")
jpeg(filepath,800,600)
ub=WwSave[1,]+WwSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(WwSave[1,]-WwSave[2,])){
  lb[i]=max(WwSave[1,i]-WwSave[2,i],0)
}
maxub=max(WwSave[1,]+WwSave[2,])
plot(days,WwSave[1,],ylab="Wet weight (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,WwSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot length
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//length.jpeg")
jpeg(filepath,800,600)
ub=LSave[1,]+LSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(LSave[1,]-LSave[2,])){
  lb[i]=max(LSave[1,i]-LSave[2,i],0)
}
maxub=max(LSave[1,]+LSave[2,])
plot(days,LSave[1,],ylab="Length (mm)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,LSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
days2 <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti) # create a dates vector to plot results

filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//temperature_response.jpeg")
jpeg(filepath,800,600)
ub=max(max(fgT),max(frT))
plot(days2,fgT,ylab="Temperature response function",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,frT,col="blue")
legend("topright",c("Anabolism limitation","Catabolism limitation"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//metabolism.jpeg")
jpeg(filepath,800,600)
Aub=ASave[1,]+ASave[2,]
Cub=CSave[1,]+CSave[2,]
Alb=as.matrix(matrix(0,nrow=length(Aub),ncol=1))
Clb=as.matrix(matrix(0,nrow=length(Cub),ncol=1))
for (i in 1:length(ASave[1,]-ASave[2,])){
  Alb[i]=max(ASave[1,i]-ASave[2,i],0)
  Clb[i]=max(CSave[1,i]-CSave[2,i],0)
}
maxub=max(Aub,Cub)
plot(days,ASave[1,],ylab="Metabolic rates (J/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(Alb,rev(Aub)),col="grey75",border=FALSE)
lines(days,ASave[1,],lwd=2,col="red")
polygon(c(days,rev(days)),c(Clb,rev(Cub)),col="grey75",border=FALSE)
lines(days,CSave[1,],lwd=2,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
legend("topleft",c("Anabolic rate","Catabolic rate"),fill=c("red","blue"))
dev.off()

# plot population dynamics
filepath=paste0(userpath,"/Clam_population/Outputs/Out_plots//population.jpeg")
jpeg(filepath,800,600)
plot(days, N, ylab="Number of individuals", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//dry_weight.csv")
write.csv(t(WdSave),filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//wet_weight.csv")
write.csv(t(WwSave),filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//length.csv")
write.csv(t(LSave),filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//temperature_response.csv")
write.csv(t(tfunSave),filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//anabolic_rate.csv")
write.csv(ASave,filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//catabolic_rate.csv")
write.csv(CSave,filepath)

filepath=paste0(userpath,"/Clam_population/Outputs/Out_csv//Days_to_commercial_size.csv")
write.csv(daysToSize,filepath)

return(output)

}
