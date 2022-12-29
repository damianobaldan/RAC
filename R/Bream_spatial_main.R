#' Bream bioenergetic spatialized model - spatialization loop
#'
#' Solves the bioenergetic balance for Bream
#'
#' @param userpath the path where the working folder is located
#' @param forcings list containing the time series in the odd positions and realted forcings in the even positions. Forcings imputted are: Water temperature [Celsius degrees] and feeding rate [g/individual x d]
#' @return saves .nc; .csv and .asc outputs in the 'Outputs' folder
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'

Bream_spatial_main<-function(userpath,forcings){

cat('Bream bioenergetic individual model spatialized\n')
cat(" \n")

#preprocessor
out_pre<-Bream_spatial_pre(userpath,forcings)

sst=forcings[[2]]

#integration_times
times=out_pre[[6]]
ti=times[1]
tf=times[2]
Dates=as.Date(out_pre[[7]],"%d/%m/%Y")

# Initialize variables to save to maps
weight<-data.frame(matrix(0,ncol(sst),nrow=tf))
actual_ingestion<-data.frame(matrix(0,ncol(sst),nrow=tf))
potential_ingestion<-data.frame(matrix(0,ncol(sst),nrow=tf))
faeces_P<-data.frame(matrix(0,ncol(sst),nrow=tf))
faeces_L<-data.frame(matrix(0,ncol(sst),nrow=tf))
faeces_C<-data.frame(matrix(0,ncol(sst),nrow=tf))
waste_P<-data.frame(matrix(0,ncol(sst),nrow=tf))
waste_L<-data.frame(matrix(0,ncol(sst),nrow=tf))
waste_C<-data.frame(matrix(0,ncol(sst),nrow=tf))
Tfun_A<-data.frame(matrix(0,ncol(sst),nrow=tf))
Tfun_C<-data.frame(matrix(0,ncol(sst),nrow=tf))
anabolism<-data.frame(matrix(0,ncol(sst),nrow=tf))
catabolism<-data.frame(matrix(0,ncol(sst),nrow=tf))
NH4<-data.frame(matrix(0,ncol(sst),nrow=tf))
O2<-data.frame(matrix(0,ncol(sst),nrow=tf))
days_commercial<-matrix(0,ncol(sst),nrow=1)

pb <- txtProgressBar(min = 0, max = ncol(sst), style = 3)

for (i in 1:ncol(sst)) {

  # Solve ODE
  forcings[[2]] = sst[1:tf,i]
  output<- Bream_spatial_loop(userpath, forcings)

  # Save outputs
  temp=output[[1]]
  weight[ti:tf,i]=temp

  temp=output[[2]]
  faeces_P[ti:tf,i]=temp[,1]
  faeces_L[ti:tf,i]=temp[,2]
  faeces_C[ti:tf,i]=temp[,3]

  temp=output[[3]]
  waste_P[ti:tf,i]=temp[,1]
  waste_L[ti:tf,i]=temp[,2]
  waste_C[ti:tf,i]=temp[,3]

  temp=output[[4]]
  potential_ingestion[ti:tf,i]=temp

  temp=output[[5]]
  actual_ingestion[ti:tf,i]=temp

  temp=output[[6]]
  Tfun_A[ti:tf,i]=temp[,1]
  Tfun_C[ti:tf,i]=temp[,2]

  temp=output[[7]]
  anabolism[ti:tf,i]=temp[,1]
  catabolism[ti:tf,i]=temp[,2]

  temp=output[[8]]
  NH4[ti:tf,i]=temp

  temp=output[[9]]
  O2[ti:tf,i]=temp

  temp=output[[10]]
  days_commercial[i]=temp

  setTxtProgressBar(pb, i)

}
close(pb)

#
weight<-weight[-(1:(ti-1)),]
actual_ingestion<-actual_ingestion[-(1:(ti-1)),]
potential_ingestion<-potential_ingestion[-(1:(ti-1)),]
faeces_P<-faeces_P[-(1:(ti-1)),]
faeces_L<-faeces_L[-(1:(ti-1)),]
faeces_C<-faeces_C[-(1:(ti-1)),]
waste_P<-waste_P[-(1:(ti-1)),]
waste_L<-waste_L[-(1:(ti-1)),]
waste_C<-waste_C[-(1:(ti-1)),]
Tfun_A<-Tfun_A[-(1:(ti-1)),]
Tfun_C<-Tfun_C[-(1:(ti-1)),]
anabolism<-anabolism[-(1:(ti-1)),]
catabolism<-catabolism[-(1:(ti-1)),]
NH4<-NH4[-(1:(ti-1)),]
O2<-O2[-(1:(ti-1)),]
#days_commercial<-days_commercial[-(1:(ti-1)),]



# Save to .csv
write.table(weight,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/weight.csv"),sep=',')
write.table(potential_ingestion,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/potential_ingestion.csv"),sep=',')
write.table(actual_ingestion,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/actual_ingestion.csv"),sep=',')
write.table(faeces_P,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/faeces_production_Proteins.csv"),sep=',')
write.table(faeces_L,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/faeces_production_Lipids.csv"),sep=',')
write.table(faeces_C,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/faeces_production_Carbohydrates.csv"),sep=',')
write.table(waste_P,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/wasted_feed_Proteins.csv"),sep=',')
write.table(waste_L,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/wasted_feed_Lipids.csv"),sep=',')
write.table(waste_C,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/wasted_feed_Carbohydrates.csv"),sep=',')
write.table(Tfun_A,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/temperature_response_anabolism.csv"),sep=',')
write.table(Tfun_C,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/temperature_response_catabolism.csv"),sep=',')
write.table(anabolism,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/anabolic_rate.csv"),sep=',')
write.table(catabolism,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/catabolic_rate.csv"),sep=',')
write.table(NH4,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/NH4_release.csv"),sep=',')
write.table(O2,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/O2_consumption.csv"),sep=',')
write.table(days_commercial,paste0(userpath,"/Bream_spatial/Outputs/Out_csv/days_to_commercial_size.csv"),sep=',')

# days_commercial<- read.csv(paste0(userpath,"/Bream_spatial/Outputs/Out_csv/days_to_commercial_size.csv"))

# Attach coordinates for maps generation
coord <- read.csv(paste0(userpath,"/Bream_spatial/Inputs/Spatial forcings//coordinates.csv"))
coord <- coord[,-(1)]

weight_map<-as.data.frame(cbind(t(coord),t(weight)))
potential_ingestion_map<-as.data.frame(cbind(t(coord),t(potential_ingestion)))
actual_ingestion_map<-as.data.frame(cbind(t(coord),t(actual_ingestion)))
faeces_P_map<-as.data.frame(cbind(t(coord),t(faeces_P)))
faeces_L_map<-as.data.frame(cbind(t(coord),t(faeces_L)))
faeces_C_map<-as.data.frame(cbind(t(coord),t(faeces_C)))
waste_P_map<-as.data.frame(cbind(t(coord),t(waste_P)))
waste_L_map<-as.data.frame(cbind(t(coord),t(waste_L)))
waste_C_map<-as.data.frame(cbind(t(coord),t(waste_C)))
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

raster::writeRaster(weight_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/weight.nc"), format="CDF", varname="weight", varunit= "g",
            longname="Weight of Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# potential ingestion map
sp::coordinates(potential_ingestion_map) <- ~V1+V2
sp::proj4string(potential_ingestion_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(potential_ingestion_map) = TRUE
potential_ingestion_brick <- raster::brick(potential_ingestion_map)
potential_ingestion_brick <- raster::setZ(potential_ingestion_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(potential_ingestion_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/potential_ingestion.nc"), format="CDF", varname="potential ingestion", varunit= "g/d",
            longname="Potential ingestion of Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# actual ingestion map
sp::coordinates(actual_ingestion_map) <- ~V1+V2
sp::proj4string(actual_ingestion_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(actual_ingestion_map) = TRUE
actual_ingestion_brick <- raster::brick(actual_ingestion_map)
actual_ingestion_brick <- raster::setZ(actual_ingestion_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(actual_ingestion_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/actual_ingestion.nc"), format="CDF", varname="actual ingestion", varunit= "g/d",
            longname="actual ingestion of Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# faeces C map
sp::coordinates(faeces_C_map) <- ~V1+V2
sp::proj4string(faeces_C_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(faeces_C_map) = TRUE
faeces_C_brick <- raster::brick(faeces_C_map)
faeces_C_brick <- raster::setZ(faeces_C_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(faeces_C_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/faeces_carbohydrates.nc"), format="CDF", varname="faeces carbohydrates", varunit= "g/d",
            longname="Carbohydrate in faeces produced by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# faeces P map
sp::coordinates(faeces_P_map) <- ~V1+V2
sp::proj4string(faeces_P_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(faeces_P_map) = TRUE
faeces_P_brick <- raster::brick(faeces_P_map)
faeces_P_brick <- raster::setZ(faeces_P_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(faeces_P_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/faeces_proteins.nc"), format="CDF", varname="faeces proteins", varunit= "g/d",
            longname="Proteins in faeces produced by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# faeces L map
sp::coordinates(faeces_L_map) <- ~V1+V2
sp::proj4string(faeces_L_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(faeces_L_map) = TRUE
faeces_L_brick <- raster::brick(faeces_L_map)
faeces_L_brick <- raster::setZ(faeces_L_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(faeces_L_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/faeces_lipids.nc"), format="CDF", varname="faeces lipids", varunit= "g/d",
            longname="Lipids in faeces produced by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# wasted C map
sp::coordinates(waste_C_map) <- ~V1+V2
sp::proj4string(waste_C_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(waste_C_map) = TRUE
waste_C_brick <- raster::brick(waste_C_map)
waste_C_brick <- raster::setZ(waste_C_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(waste_C_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/wasted_feed_carbohydrates.nc"), format="CDF", varname="wasted carbohydrates", varunit= "g/d",
            longname="Carbohydrates wasted feed by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# wasted L map
sp::coordinates(waste_L_map) <- ~V1+V2
sp::proj4string(waste_L_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(waste_L_map) = TRUE
waste_L_brick <- raster::brick(waste_L_map)
waste_L_brick <- raster::setZ(waste_L_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(waste_L_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/wasted_feed_lipids.nc"), format="CDF", varname="wasted lipids", varunit= "g/d",
            longname="Lipids wasted feed by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# wasted P map
sp::coordinates(waste_P_map) <- ~V1+V2
sp::proj4string(waste_P_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(waste_P_map) = TRUE
waste_P_brick <- raster::brick(waste_P_map)
waste_P_brick <- raster::setZ(waste_P_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(waste_P_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/wasted_feed_proteins.nc"), format="CDF", varname="wasted proteins", varunit= "g/d",
            longname="Proteins wasted by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# temperature response function for anabolism map
sp::coordinates(Tfun_A_map) <- ~V1+V2
sp::proj4string(Tfun_A_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(Tfun_A_map) = TRUE
Tfun_A_brick <- raster::brick(Tfun_A_map)
Tfun_A_brick <- raster::setZ(Tfun_A_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(Tfun_A_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/temperature_response_A.nc"), format="CDF", varname="Temperature response anabolism", varunit= "-",
            longname="temperature response function for anabolism of Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# temperature response function for catabolism map
sp::coordinates(Tfun_C_map) <- ~V1+V2
sp::proj4string(Tfun_C_map)=CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(Tfun_C_map) = TRUE
Tfun_C_brick <- raster::brick(Tfun_C_map)
Tfun_C_brick <- raster::setZ(Tfun_C_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(Tfun_C_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/temperature_response_C.nc"), format="CDF", varname="Temperature response catabolism", varunit= "-",
            longname="temperature response function for catabolism of Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# anabolic rates map
sp::coordinates(anabolism_map) <- ~V1+V2
sp::proj4string(anabolism_map)=CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(anabolism_map) = TRUE
anabolism_brick <- raster::brick(anabolism_map)
anabolism_brick <- raster::setZ(anabolism_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(anabolism_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/anabolic_rate.nc"), format="CDF", varname="Anabolic rate", varunit= "J/d",
            longname="anabolic rate for catabolism of Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# catabolic rates map
sp::coordinates(catabolism_map) <- ~V1+V2
sp::proj4string(catabolism_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(catabolism_map) = TRUE
catabolism_brick <- raster::brick(catabolism_map)
catabolism_brick <- raster::setZ(catabolism_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(catabolism_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/catabolic_rate.nc"), format="CDF", varname="Catabolic rate", varunit= "J/d",
            longname="catabolic rate for catabolism of Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# NH4 release map
sp::coordinates(NH4_map) <- ~V1+V2
sp::proj4string(NH4_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(NH4_map) = TRUE
NH4_brick <- raster::brick(NH4_map)
NH4_brick <- raster::setZ(NH4_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(NH4_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/NH4_release.nc"), format="CDF", varname="NH4 release", varunit= "gN",
            longname="Ammonia released by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# O2 release map
sp::coordinates(O2_map) <- ~V1+V2
sp::proj4string(O2_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(O2_map) = TRUE
O2_brick <- raster::brick(O2_map)
O2_brick <- raster::setZ(O2_brick, as.Date(Dates[1])+ 0:(tf-ti))    ######ADD DATE

raster::writeRaster(O2_brick, paste0(userpath,"/Bream_spatial/Outputs/Out_nc/O2_consumption.nc"), format="CDF", varname="O2 consumption", varunit= "g",
            longname="Oxygen consumed by Sparus aurata estimated through the R RAC package developed by Baldan et al",
            zname="time",
            zunit="day", overwrite=T)

# Days to commercial size
sp::coordinates(days_commercial_map) <- ~V1+V2
sp::proj4string(days_commercial_map)=sp::CRS(("+proj=longlat +datum=WGS84")) # set it to lat-long
sp::gridded(days_commercial_map) = TRUE
days_commercial_raster = raster::raster(days_commercial_map)
raster::projection(days_commercial_raster) = sp::CRS("+proj=longlat +datum=WGS84")
raster::writeRaster(days_commercial_raster,paste0(userpath,"/Bream_spatial/Outputs/Out_asc/days_to_commercial_size.asc"),format="ascii",overwrite=TRUE)

}
