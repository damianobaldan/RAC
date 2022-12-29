#' Mussel bioenergetic spatialized model - spatialization loop
#'
#' Solves the bioenergetic balance for Mussel
#'
#' @param userpath the path where the working folder is located
#' @param forcings a list containing the time series in the odd positions and realted forcings in the even positions. Forcings returned are: Water temperature [Celsius degrees], Chlorophyll a concentration [mgChl-a/m^3], particulated organic carbon (POC) concentration [mgC/l] and its characterization in terms of C/P and N/P molar ratios, particulated organic matter (POM) concentration [mgC/l], total suspended solids (TSS) concentration [mg/l]
#' @return saves .nc; .csv and .asc outputs in the 'Outputs' folder
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#' @importFrom sp coordinates proj4string gridded CRS
#' @importFrom raster brick setZ writeRaster
#'

Mussel_spatial_main<-function(userpath,forcings){

cat('Mussel bioenergetic individual model spatialized\n')
cat(" \n")

#sst <- read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//sst.csv"),header=T)
#sst <- sst[,-(1)]
#colnames(sst) <- gsub("V", "sst", colnames(sst))

# Chlorophyll
#chl <- read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//chl.csv"),header=T)
#chl <- chl[,-(1)]
#colnames(chl) <- gsub("V", "chl", colnames(chl))

#preprocessor
out_pre<-Mussel_spatial_pre(userpath,forcings)

sst=forcings[[2]]
chl=forcings[[4]]

#integration_times
times=out_pre[[2]]
ti=times[1]
tf=times[2]
Dates=out_pre[[3]]

# Initialize variables to save to maps
weight<-data.frame(matrix(0,ncol(sst),nrow=tf))
length<-data.frame(matrix(0,ncol(sst),nrow=tf))
C_content<-data.frame(matrix(0,ncol(sst),nrow=tf))
N_content<-data.frame(matrix(0,ncol(sst),nrow=tf))
P_content<-data.frame(matrix(0,ncol(sst),nrow=tf))
faecies_C<-data.frame(matrix(0,ncol(sst),nrow=tf))
faecies_N<-data.frame(matrix(0,ncol(sst),nrow=tf))
faecies_P<-data.frame(matrix(0,ncol(sst),nrow=tf))
pseudofaecies_C<-data.frame(matrix(0,ncol(sst),nrow=tf))
pseudofaecies_N<-data.frame(matrix(0,ncol(sst),nrow=tf))
pseudofaecies_P<-data.frame(matrix(0,ncol(sst),nrow=tf))
Tfun_A<-data.frame(matrix(0,ncol(sst),nrow=tf))
Tfun_C<-data.frame(matrix(0,ncol(sst),nrow=tf))
anabolism<-data.frame(matrix(0,ncol(sst),nrow=tf))
catabolism<-data.frame(matrix(0,ncol(sst),nrow=tf))
NH4<-data.frame(matrix(0,ncol(sst),nrow=tf))
O2<-data.frame(matrix(0,ncol(sst),nrow=tf))
days_commercial<-data.frame(matrix(0,ncol(sst),nrow=1))

pb <- txtProgressBar(min = 0, max = ncol(sst), style = 3)

for (i in 1:ncol(sst)) {

  # Solve ODE
  forcings[[2]] = sst[1:tf,i]
  forcings[[4]] = chl[1:tf,i]
  output<- Mussel_spatial_loop(userpath, forcings)

  # Save outputs
  temp=output[[1]]
  weight[ti:tf,i]=temp[,4]
  length[ti:tf,i]=temp[,5]

  temp=output[[2]]
  pseudofaecies_C[ti:tf,i]=temp[,1]
  pseudofaecies_N[ti:tf,i]=temp[,2]
  pseudofaecies_P[ti:tf,i]=temp[,3]

  temp=output[[3]]
  faecies_C[ti:tf,i]=temp[,1]
  faecies_C[ti:tf,i]=temp[,2]
  faecies_C[ti:tf,i]=temp[,3]

  temp=output[[4]]
  C_content[ti:tf,i]=temp[,1]
  N_content[ti:tf,i]=temp[,2]
  P_content[ti:tf,i]=temp[,3]

  temp=output[[5]]
  Tfun_A[ti:tf,i]=temp[,1]
  Tfun_C[ti:tf,i]=temp[,2]

  temp=output[[6]]
  anabolism[ti:tf,i]=temp[,1]
  catabolism[ti:tf,i]=temp[,2]

  temp=output[[7]]
  NH4[ti:tf,i]=temp

  temp=output[[8]]
  O2[ti:tf,i]=temp

  temp=output[[9]]
  days_commercial[i]=temp

  setTxtProgressBar(pb, i)

}
close(pb)

weight<-weight[-(1:(ti-1)),]
length<-length[-(1:(ti-1)),]
C_content<-C_content[-(1:(ti-1)),]
N_content<-N_content[-(1:(ti-1)),]
P_content<-P_content[-(1:(ti-1)),]
faecies_C<-faecies_C[-(1:(ti-1)),]
faecies_N<-faecies_N[-(1:(ti-1)),]
faecies_P<-faecies_P[-(1:(ti-1)),]
pseudofaecies_C<-pseudofaecies_C[-(1:(ti-1)),]
pseudofaecies_N<-pseudofaecies_N[-(1:(ti-1)),]
pseudofaecies_P<-pseudofaecies_P[-(1:(ti-1)),]
Tfun_A<-Tfun_A[-(1:(ti-1)),]
Tfun_C<-Tfun_C[-(1:(ti-1)),]
anabolism<-anabolism[-(1:(ti-1)),]
catabolism<-catabolism[-(1:(ti-1)),]
NH4<-NH4[-(1:(ti-1)),]
O2<-O2[-(1:(ti-1)),]

# Save to .csv
write.table(weight,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/dry_weight.csv"),sep=',')
write.table(length,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/length.csv"),sep=',')
write.table(C_content,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/C_content.csv"),sep=',')
write.table(N_content,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/N_content.csv"),sep=',')
write.table(P_content,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/P_content.csv"),sep=',')
write.table(faecies_C,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/faeces_C.csv"),sep=',')
write.table(faecies_N,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/faeces_N.csv"),sep=',')
write.table(faecies_P,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/faeces_P.csv"),sep=',')
write.table(pseudofaecies_C,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/Pseudofaeces_C.csv"),sep=',')
write.table(pseudofaecies_N,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/Pseudofaeces_N.csv"),sep=',')
write.table(pseudofaecies_P,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/Pseudofaeces_P.csv"),sep=',')
write.table(Tfun_A,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/temperature_response_anabolism.csv"),sep=',')
write.table(Tfun_C,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/temperature_response_catabolism.csv"),sep=',')
write.table(anabolism,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/anabolic_rate.csv"),sep=',')
write.table(catabolism,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/catabolic_rate.csv"),sep=',')
write.table(NH4,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/NH4_release.csv"),sep=',')
write.table(O2,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/O2_consumption.csv"),sep=',')
write.table(days_commercial,paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/days_to_commercial_size.csv"),sep=',')

#days_commercial<- read.csv(paste0(userpath,"/Mussel_spatial/Outputs/Out_csv/days_to_commercial_size.csv"))

# Attach coordinates for maps generation
coord <- read.csv(paste0(userpath,"/Mussel_spatial/Inputs/Spatial forcings//coordinates.csv"))
coord <- coord[,-(1)]

weight_map<-as.data.frame(cbind(t(coord),t(weight)))
length_map<-as.data.frame(cbind(t(coord),t(length)))
C_content_map<-as.data.frame(cbind(t(coord),t(C_content)))
N_content_map<-as.data.frame(cbind(t(coord),t(N_content)))
P_content_map<-as.data.frame(cbind(t(coord),t(P_content)))
faecies_C_map<-as.data.frame(cbind(t(coord),t(faecies_C)))
faecies_N_map<-as.data.frame(cbind(t(coord),t(faecies_N)))
faecies_P_map<-as.data.frame(cbind(t(coord),t(faecies_P)))
pseudofaecies_C_map<-as.data.frame(cbind(t(coord),t(pseudofaecies_C)))
pseudofaecies_N_map<-as.data.frame(cbind(t(coord),t(pseudofaecies_N)))
pseudofaecies_P_map<-as.data.frame(cbind(t(coord),t(pseudofaecies_P)))
Tfun_A_map<-as.data.frame(cbind(t(coord),t(Tfun_A)))
Tfun_C_map<-as.data.frame(cbind(t(coord),t(Tfun_C)))
anabolism_map<-as.data.frame(cbind(t(coord),t(anabolism)))
catabolism_map<-as.data.frame(cbind(t(coord),t(catabolism)))
NH4_map<-as.data.frame(cbind(t(coord),t(NH4)))
O2_map<-as.data.frame(cbind(t(coord),t(O2)))
days_commercial_map<-as.data.frame(cbind(t(coord),t(days_commercial)))

# Weight map
sp::coordinates(weight_map) <- ~V1+V2
sp::proj4string(weight_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(weight_map) = TRUE
weight_brick <- raster::brick(weight_map)
weight_brick <- raster::setZ(weight_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(weight_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/dry_weight.nc"), format="CDF", varname="dry weight", varunit= "cm",
            longname="Dry weight of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# Length map
sp::coordinates(length_map) <- ~V1+V2
sp::proj4string(length_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(length_map) = TRUE
length_brick <- raster::brick(length_map)
length_brick <- raster::setZ(length_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(length_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/length.nc"), format="CDF", varname="length", varunit= "cm",
            longname="Length of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# C content map
sp::coordinates(C_content_map) <- ~V1+V2
sp::proj4string(C_content_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(C_content_map) = TRUE
C_content_brick <- raster::brick(C_content_map)
C_content_brick <- raster::setZ(C_content_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(C_content_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/C_content.nc"), format="CDF", varname="C content", varunit= "gC",
            longname="Carbon content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# N content map
sp::coordinates(N_content_map) <- ~V1+V2
sp::proj4string(N_content_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(N_content_map) = TRUE
N_content_brick <- raster::brick(N_content_map)
N_content_brick <- raster::setZ(N_content_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(N_content_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/N_content.nc"), format="CDF", varname="N content", varunit= "gN",
            longname="Nitrogen content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# P content map
sp::coordinates(P_content_map) <- ~V1+V2
sp::proj4string(P_content_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(P_content_map) = TRUE
P_content_brick <- raster::brick(P_content_map)
P_content_brick <- raster::setZ(P_content_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(P_content_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/P_content.nc"), format="CDF", varname="P content", varunit= "gP",
            longname="Phosphorous content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# faeces C content map
sp::coordinates(faecies_C_map) <- ~V1+V2
sp::proj4string(faecies_C_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(faecies_C_map) = TRUE
faecies_C_brick <- raster::brick(faecies_C_map)
faecies_C_brick <- raster::setZ(faecies_C_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(faecies_C_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/faeces_C.nc"), format="CDF", varname="Faeces C", varunit= "gC",
            longname="Faeces C content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# faeces N content map
sp::coordinates(faecies_N_map) <- ~V1+V2
sp::proj4string(faecies_N_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(faecies_N_map) = TRUE
faecies_N_brick <- raster::brick(faecies_N_map)
faecies_N_brick <- raster::setZ(faecies_N_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(faecies_N_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/faeces_N.nc"), format="CDF", varname="Faeces N", varunit= "gN",
            longname="Faeces N content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# faeces P content map
sp::coordinates(faecies_P_map) <- ~V1+V2
sp::proj4string(faecies_P_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(faecies_P_map) = TRUE
faecies_P_brick <- raster::brick(faecies_P_map)
faecies_P_brick <- raster::setZ(faecies_P_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(faecies_P_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/faeces_P.nc"), format="CDF", varname="Faeces P", varunit= "gP",
            longname="Faeces P content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# pseudofaeces C content map
sp::coordinates(pseudofaecies_C_map) <- ~V1+V2
sp::proj4string(pseudofaecies_C_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(pseudofaecies_C_map) = TRUE
pseudofaecies_C_brick <- raster::brick(pseudofaecies_C_map)
pseudofaecies_C_brick <- raster::setZ(pseudofaecies_C_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(pseudofaecies_C_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/pseudofaeces_C.nc"), format="CDF", varname="Pseudofaeces C", varunit= "gC",
            longname="Pseudofaeces C content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# pseudofaeces N content map
sp::coordinates(pseudofaecies_N_map) <- ~V1+V2
sp::proj4string(pseudofaecies_N_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(pseudofaecies_N_map) = TRUE
pseudofaecies_N_brick <- raster::brick(pseudofaecies_N_map)
pseudofaecies_N_brick <- raster::setZ(pseudofaecies_N_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(pseudofaecies_N_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/pseudofaeces_N.nc"), format="CDF", varname="Pseudofaeces N", varunit= "gN",
            longname="Pseudofaeces N content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# pseudofaeces P content map
sp::coordinates(pseudofaecies_P_map) <- ~V1+V2
sp::proj4string(pseudofaecies_P_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(pseudofaecies_P_map) = TRUE
pseudofaecies_P_brick <- raster::brick(pseudofaecies_P_map)
pseudofaecies_P_brick <- raster::setZ(pseudofaecies_P_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(pseudofaecies_P_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/pseudofaeces_P.nc"), format="CDF", varname="Pseudofaeces P", varunit= "gP",
            longname="Pseudofaeces P content of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# temperature response function for anabolism map
sp::coordinates(Tfun_A_map) <- ~V1+V2
sp::proj4string(Tfun_A_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(Tfun_A_map) = TRUE
Tfun_A_brick <- raster::brick(Tfun_A_map)
Tfun_A_brick <- raster::setZ(Tfun_A_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(Tfun_A_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/temperature_response_A.nc"), format="CDF", varname="Temperature response anabolism", varunit= "-",
            longname="temperature response function for anabolism of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# temperature response function for catabolism map
sp::coordinates(Tfun_C_map) <- ~V1+V2
sp::proj4string(Tfun_C_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(Tfun_C_map) = TRUE
Tfun_C_brick <- raster::brick(Tfun_C_map)
Tfun_C_brick <- raster::setZ(Tfun_C_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(Tfun_C_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/temperature_response_C.nc"), format="CDF", varname="Temperature response catabolism", varunit= "-",
            longname="temperature response function for catabolism of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# anabolic rates map
sp::coordinates(anabolism_map) <- ~V1+V2
sp::proj4string(anabolism_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(anabolism_map) = TRUE
anabolism_brick <- raster::brick(anabolism_map)
anabolism_brick <- raster::setZ(anabolism_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(anabolism_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/anabolism.nc"), format="CDF", varname="Anabolic rate", varunit= "J/d",
            longname="anabolic rate for catabolism of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# catabolic rates map
sp::coordinates(catabolism_map) <- ~V1+V2
sp::proj4string(catabolism_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(catabolism_map) = TRUE
catabolism_brick <- raster::brick(catabolism_map)
catabolism_brick <- raster::setZ(catabolism_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(catabolism_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/catabolism.nc"), format="CDF", varname="Catabolic rate", varunit= "J/d",
            longname="catabolic rate for catabolism of Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# NH4 release map
sp::coordinates(NH4_map) <- ~V1+V2
sp::proj4string(NH4_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(NH4_map) = TRUE
NH4_brick <- raster::brick(NH4_map)
NH4_brick <- raster::setZ(NH4_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(NH4_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/NH4_release.nc"), format="CDF", varname="NH4 release", varunit= "gN",
            longname="Ammonia released by Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# O2 release map
sp::coordinates(O2_map) <- ~V1+V2
sp::proj4string(O2_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(O2_map) = TRUE
O2_brick <- raster::brick(O2_map)
O2_brick <- raster::setZ(O2_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(O2_brick, paste0(userpath,"/Mussel_spatial/Outputs/Out_nc/O2_consumption.nc"), format="CDF", varname="O2 consumption", varunit= "g",
            longname="Oxygen consumed by Mytilus galloprovincialis estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# Days to commercial size
sp::coordinates(days_commercial_map) <- ~V1+V2
sp::proj4string(days_commercial_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(days_commercial_map) = TRUE
days_commercial_raster = raster::raster(days_commercial_map)
raster::projection(days_commercial_raster) = sp::CRS("+proj=longlat +datum=WGS84")
raster::writeRaster(days_commercial_raster,paste0(userpath,"/Mussel_spatial/Outputs/Out_asc/days_to_commercial_size.asc"),format="ascii",overwrite=TRUE)

}
