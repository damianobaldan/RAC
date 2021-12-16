#' Creates the folders structure for Clam individual bioenergetic model (alternative version)
#'
#' @param userpath the path where forcing are located
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'


ClamF_ind_skeleton<-function(userpath){

  workingpath=path.package("RAC", quiet = FALSE) # Save current location of R script

  # Create the folders structure in the path set by the user
  dir.create(paste0(userpath,"/ClamF_individual"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_individual/Inputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_individual/Inputs/Parameters"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_individual/Inputs/Forcings"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_individual/Inputs/Forcings_plots"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_individual/Outputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_individual/Outputs/Out_csv"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_individual/Outputs/Out_plots"),showWarnings=FALSE)

  # Moves the selected data from the package folder to the user folder
  file.copy(paste0(workingpath,"/extdata/ClamF_ind_data//Parameters.csv"),paste0(userpath,"/ClamF_individual/Inputs/Parameters"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/ClamF_ind_data//Water_temperature.csv"),paste0(userpath,"/ClamF_individual/Inputs/Forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/ClamF_ind_data//Chlorophyll_a.csv"),paste0(userpath,"/ClamF_individual/Inputs/Forcings"), overwrite=TRUE)

  cat("Folder skeleton for Clam individual bioenergetic model (alternative version) created at:\n")
  cat(userpath)
  cat("\n")
  cat("ATTENTION: Executing again this function will overwrite the files\n")

}
