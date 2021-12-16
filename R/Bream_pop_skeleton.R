#' Creates the folders structure for Seabream population model
#'
#' @param userpath the path where forcing are located
#' @export
#'
#' @import matrixStats plotrix rstudioapi
#'


Bream_pop_skeleton<-function(userpath){

  workingpath=path.package("RAC", quiet = FALSE) # Save current location of R script

  # Create the folders structure in the path set by the user
  dir.create(paste0(userpath,"/Bream_population"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Inputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Inputs/Parameters"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Inputs/Forcings"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Inputs/Forcings_plots"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Inputs/Population_management"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Outputs"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Outputs/Out_csv"),showWarnings=FALSE)
  dir.create(paste0(userpath,"/Bream_population/Outputs/Out_plots"),showWarnings=FALSE)

  # Moves the selected data from the package folder to the user folder
  file.copy(paste0(workingpath,"/extdata/Bream_pop_data//Parameters.csv"),paste0(userpath,"/Bream_population/Inputs/Parameters"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Bream_pop_data//Water_temperature.csv"),paste0(userpath,"/Bream_population/Inputs/Forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Bream_pop_data//Feeding.csv"),paste0(userpath,"/Bream_population/Inputs/Forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Bream_pop_data//Food_characterization.csv"),paste0(userpath,"/Bream_population/Inputs/Forcings"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Bream_pop_data//Population.csv"),paste0(userpath,"/Bream_population/Inputs/Parameters"), overwrite=TRUE)
  file.copy(paste0(workingpath,"/extdata/Bream_pop_data//Management.csv"),paste0(userpath,"/Bream_population/Inputs/Population_management"), overwrite=TRUE)

  cat("Folder skeleton for Sea Bream population model created at:\n")
  cat(userpath)
  cat("\n")
  cat("ATTENTION: Executing again this function will overwrite the files\n")

}
