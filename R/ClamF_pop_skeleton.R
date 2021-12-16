#' Creates the folders structure for Clam population model (alternative version)
#'
#' @param userpath the path where forcing are located
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'


ClamF_pop_skeleton<-function(userpath){

  workingpath=path.package("RAC", quiet = FALSE) # Save current location of R script

  # Create the folders structure in the path set by the user
  dir.create(paste0(userpath,"/ClamF_population"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Inputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Inputs/Parameters"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Inputs/Forcings"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Inputs/Forcings_plots"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Inputs/Population_management"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Outputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Outputs/Out_csv"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/ClamF_population/Outputs/Out_plots"),showWarnings=FALSE)

  # Moves the selected data from the package folder to the user folder
  file.copy(paste0(workingpath,"/extdata/ClamF_pop_data//Parameters.csv"),paste0(userpath,"/ClamF_population/Inputs/Parameters"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/ClamF_pop_data//Water_temperature.csv"),paste0(userpath,"/ClamF_population/Inputs/Forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/ClamF_pop_data//Chlorophyll_a.csv"),paste0(userpath,"/ClamF_population/Inputs/Forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/ClamF_pop_data//Population.csv"),paste0(userpath,"/ClamF_population/Inputs/Parameters"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/ClamF_pop_data//Management.csv"),paste0(userpath,"/ClamF_population/Inputs/Population_management"), overwrite=TRUE)


  cat("Folder skeleton for Clam population model (alternative version) created at:\n")
  cat(userpath)
  cat("\n")
  cat("ATTENTION: Executing again this function will overwrite the files\n")

}
