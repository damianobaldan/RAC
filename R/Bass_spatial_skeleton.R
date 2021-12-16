#' Creates the folders structure for Bass spatialized model
#'
#' @param userpath the path where forcing are located
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'


Bass_spatial_skeleton<-function(userpath){

  workingpath=path.package("RAC", quiet = FALSE) # Save current location of R script

  # Create the folders structure in the path set by the user
  dir.create(paste0(userpath,"/Bass_spatial"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Inputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Inputs/Parameters"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Inputs/Point forcings"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Inputs/Forcings_plots"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Outputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Outputs/Out_asc"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Outputs/Out_nc"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bass_spatial/Outputs/Out_csv"),showWarnings=FALSE)

  # Moves the selected data from the package folder to the user folder
  file.copy(paste0(workingpath,"/extdata/Bass_spatial_data//Parameters.csv"),paste0(userpath,"/Bass_spatial/Inputs/Parameters"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Mussel_spatial_data//SST.nc"),paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Bass_spatial_data//Feeding.csv"),paste0(userpath,"/Bass_spatial/Inputs/Point forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Bass_spatial_data//Food_characterization.csv"),paste0(userpath,"/Bass_spatial/Inputs/Point forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Mussel_spatial_data//Spatial_dates.csv"),paste0(userpath,"/Bass_spatial/Inputs/Spatial forcings"), overwrite=TRUE)

  cat("Folder skeleton for Bass spatialized model created at:\n")
  cat(userpath)
  cat("\n")
  cat("ATTENTION: Executing again this function will overwrite the files\n")

}
