#' Postprocess the Mussel population bioenergetic model results
#'
#' @param userpath the path where the working folder is located
#' @param output output list containing the output of the RK solver
#' @param times the vector containing informations on integration extremes
#' @param Dates the vector containing the date
#' @param N the number of individuals
#' @param CS the commercial size of Seabass
#'
#' @return a list containing the weights of the mussel, the excreted CNP, the mussel CNP, temperature limitation functions, metabolic rates, oxygen consumption
#'
#' @import matrixStats plotrix rstudioapi
#'

Mussel_pop_post<-function(userpath,output,times,Dates,N,CS) {

cat('Data post-processing\n')
cat('\n')

ti=times[1]           # Integration beginning
tf=times[2]           # Integration end

# Extracts results from output list
Wb_stat=output[[1]]
R_stat=output[[2]]
Wd_stat=output[[3]]
W_stat=output[[4]]
L_stat=output[[5]]
fecC_stat=output[[6]]
fecN_stat=output[[7]]
fecP_stat=output[[8]]
psC_stat=output[[9]]
psN_stat=output[[10]]
psP_stat=output[[11]]
Cmyt_stat=output[[12]]
Nmyt_stat=output[[13]]
Pmyt_stat=output[[14]]
A_stat=output[[15]]
C_stat=output[[16]]
O2_stat=output[[17]]
NH4_stat=output[[18]]
fgT=output[[19]]
frT=output[[20]]

# Adjusts results acoording with integration extremes
# now day 1 coincides with ti

WbSave=Wb_stat[,ti:tf]
RSave=R_stat[,ti:tf]
WdSave=Wd_stat[,ti:tf]
WSave=W_stat[,ti:tf]
LSave=L_stat[,ti:tf]

fecCSave=fecC_stat[,ti:tf]
fecNSave=fecN_stat[,ti:tf]
fecPSave=fecP_stat[,ti:tf]

psCSave=psC_stat[,ti:tf]
psNSave=psN_stat[,ti:tf]
psPSave=psP_stat[,ti:tf]

CmytSave=Cmyt_stat[,ti:tf]
NmytSave=Nmyt_stat[,ti:tf]
PmytSave=Pmyt_stat[,ti:tf]

ASave=A_stat[,ti:tf]
CSave=C_stat[,ti:tf]

O2Save=O2_stat[,ti:tf]
NH4Save=NH4_stat[,ti:tf]

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

output=list(WbSave,RSave,WdSave,WSave,LSave,fecCSave,fecNSave,fecPSave,psCSave,psNSave,psPSave,CmytSave,NmytSave,PmytSave,ASave,CSave,fgT,frT,N,daysToSize)

# Plot results
days <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti+1) # create a dates vector to plot results

# Plot dry weight
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//Dry_weight.jpeg")
jpeg(filepath,800,600)
plot(days,WdSave[1,],ylab="Mean dry weight (g)", xlab=" ",xaxt = "n",type="l",lwd=2,cex.lab=1.4,col="red")
lines(days,WbSave[1,],lwd=2,col="green")
lines(days,RSave[1,],lwd=2,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
legend("topleft",c("Total","Somatic tissue","Gonadic tissue"),fill=c("red","green","blue"))
dev.off()

# Plot length
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//Length.jpeg")
jpeg(filepath,800,600)
ub=LSave[1,]+LSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(LSave[1,]-LSave[2,])){
  lb[i]=max(LSave[1,i]-LSave[2,i],0)
}
maxub=max(LSave[1,]+LSave[2,])
plot(days,LSave[1,],ylab="Length (cm)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,LSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot total weight
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//Total_weight.jpeg")
jpeg(filepath,800,600)
ub=WSave[1,]+WSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(WSave[1,]-WSave[2,])){
  lb[i]=max(WSave[1,i]-WSave[2,i],0)
}
maxub=max(WSave[1,]+WSave[2,])
plot(days,WSave[1,],ylab="Total weight - with shell (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,WSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot pseudofaecies C
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//pseudofaeces_C.jpeg")
jpeg(filepath,800,600)
ub=psCSave[1,]+psCSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(psCSave[1,]-psCSave[2,])){
  lb[i]=max(psCSave[1,i]-psCSave[2,i],0)
}
maxub=max(psCSave[1,]+psCSave[2,])
plot(days,psCSave[1,],ylab="C in pseudofaeces (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,psCSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot pseudofaecies N
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//pseudofaeces_N.jpeg")
jpeg(filepath,800,600)
ub=psNSave[1,]+psNSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(psNSave[1,]-psNSave[2,])){
  lb[i]=max(psNSave[1,i]-psNSave[2,i],0)
}
maxub=max(psNSave[1,]+psNSave[2,])
plot(days,psNSave[1,],ylab="N in pseudofaeces (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,psNSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot pseudofaecies P
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//pseudofaeces_P.jpeg")
jpeg(filepath,800,600)
ub=psPSave[1,]+psPSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(psPSave[1,]-psPSave[2,])){
  lb[i]=max(psPSave[1,i]-psPSave[2,i],0)
}
maxub=max(psPSave[1,]+psPSave[2,])
plot(days,psPSave[1,],ylab="P in pseudofaeces (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,psPSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot faeces C
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//faeces_C.jpeg")
jpeg(filepath,800,600)
ub=fecCSave[1,]+fecCSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(fecCSave[1,]-fecCSave[2,])){
  lb[i]=max(fecCSave[1,i]-fecCSave[2,i],0)
}
maxub=max(fecCSave[1,]+fecCSave[2,])
plot(days,fecCSave[1,],ylab="C in faeces (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,fecCSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot faeces N
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//faeces_N.jpeg")
jpeg(filepath,800,600)
ub=fecNSave[1,]+fecNSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(fecNSave[1,]-fecNSave[2,])){
  lb[i]=max(fecNSave[1,i]-fecNSave[2,i],0)
}
maxub=max(fecNSave[1,]+fecNSave[2,])
plot(days,fecNSave[1,],ylab="N in faeces (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,fecNSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot faeces P
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//faeces_P.jpeg")
jpeg(filepath,800,600)
ub=fecPSave[1,]+fecPSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(fecPSave[1,]-fecPSave[2,])){
  lb[i]=max(fecPSave[1,i]-fecPSave[2,i],0)
}
maxub=max(fecPSave[1,]+fecPSave[2,])
plot(days,fecPSave[1,],ylab="P in faeces (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,fecPSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot mussel C content
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//C_content.jpeg")
jpeg(filepath,800,600)
ub=CmytSave[1,]+CmytSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(CmytSave[1,]-CmytSave[2,])){
  lb[i]=max(CmytSave[1,i]-CmytSave[2,i],0)
}
maxub=max(CmytSave[1,]+CmytSave[2,])
plot(days,CmytSave[1,],ylab="Mussel C content (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,CmytSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot mussel N content
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//N_content.jpeg")
jpeg(filepath,800,600)
ub=NmytSave[1,]+NmytSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(NmytSave[1,]-NmytSave[2,])){
  lb[i]=max(NmytSave[1,i]-NmytSave[2,i],0)
}
maxub=max(NmytSave[1,]+NmytSave[2,])
plot(days,NmytSave[1,],ylab="Mussel N content (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,NmytSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot mussel P content
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//P_content.jpeg")
jpeg(filepath,800,600)
ub=PmytSave[1,]+PmytSave[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(PmytSave[1,]-PmytSave[2,])){
  lb[i]=max(PmytSave[1,i]-PmytSave[2,i],0)
}
maxub=max(PmytSave[1,]+PmytSave[2,])
plot(days,PmytSave[1,],ylab="Mussel P content (g)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,PmytSave[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()


# Plot oxygen consumption
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//O2_consumption.jpeg")
jpeg(filepath,800,600)
ub=O2Save[1,]+O2Save[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(O2Save[1,]-O2Save[2,])){
  lb[i]=max(O2Save[1,i]-O2Save[2,i],0)
}
maxub=max(O2Save[1,]+O2Save[2,])
plot(days,O2Save[1,],ylab="Oxygen consumption (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,O2Save[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Plot ammonium release
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//NH4_release.jpeg")
jpeg(filepath,800,600)
ub=NH4Save[1,]+NH4Save[2,]
lb=as.matrix(matrix(0,nrow=length(ub),ncol=1))
for (i in 1:length(NH4Save[1,]-NH4Save[2,])){
  lb[i]=max(NH4Save[1,i]-NH4Save[2,i],0)
}
maxub=max(NH4Save[1,]+NH4Save[2,])
plot(days,NH4Save[1,],ylab="NH4 release (kg/d)", xlab=" ",xaxt = "n",type="l",cex.lab=1.4,col="red",ylim=c(0,maxub+0.05*maxub))
polygon(c(days,rev(days)),c(lb,rev(ub)),col="grey90",border=FALSE)
lines(days,NH4Save[1,],lwd=2,col="red")
lines(days,lb,col="blue")
lines(days,ub,col="blue")
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot limitation functions
days2 <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), by = "days", length = tf-ti) # create a dates vector to plot results

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//temperature_response.jpeg")
jpeg(filepath,800,600)
ub=max(max(fgT),max(frT))
plot(days2,fgT,ylab="Temperature response function",xlab=" ",xaxt = "n",cex.lab=1.4,col="red",type="l",ylim=c(0,ub+0.05*ub))
lines(days2,frT,col="blue")
legend("topright",c("Anabolism limitation","Catabolism limitation"),fill=c("red","blue"))
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# plot metabolic rates
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots//metabolism.jpeg")
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
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_plots/Population.jpeg")
jpeg(filepath,800,600)
plot(days, N, ylab="Number of individuals", xlab="", xaxt = "n",type="l",cex.lab=1.4)
labDates <- seq(as.Date(Dates[1], format = "%d/%m/%Y"), tail(days, 1), by = "months")
axis.Date(side = 1, days, at = labDates, format = "%d %b %y", las = 2)
dev.off()

# Results save

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//Dry_weight.csv")
write.csv(t(WdSave),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//Total_weight.csv")
write.csv(t(WSave),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//Length.csv")
write.csv(t(LSave),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//faeces_C.csv")
write.csv(t(fecCSave),filepath)
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//faeces_N.csv")
write.csv(t(fecNSave),filepath)
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//faeces_P.csv")
write.csv(t(fecPSave),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//pseudofaeces_C.csv")
write.csv(t(psCSave),filepath)
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//pseudofaeces_N.csv")
write.csv(t(psNSave),filepath)
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//pseudofaeces_P.csv")
write.csv(t(psPSave),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//C_content.csv")
write.csv(t(CmytSave),filepath)
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//N_content.csv")
write.csv(t(NmytSave),filepath)
filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//P_content.csv")
write.csv(t(PmytSave),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//O2_consumption.csv")
write.csv(t(O2Save),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//NH4_release.csv")
write.csv(t(NH4Save),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//temperature_response.csv")
write.csv(t(tfunSave),filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//anabolic_rate.csv")
write.csv(ASave,filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//catabolic_rate.csv")
write.csv(CSave,filepath)

filepath=paste0(userpath,"/Mussel_population/Outputs/Out_csv//Days_to_commercial_size.csv")
write.csv(daysToSize,filepath)

return(output)

}
