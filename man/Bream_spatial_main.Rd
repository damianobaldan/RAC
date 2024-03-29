% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bream_spatial_main.R
\name{Bream_spatial_main}
\alias{Bream_spatial_main}
\title{Bream bioenergetic spatialized model - spatialization loop}
\usage{
Bream_spatial_main(userpath, forcings)
}
\arguments{
\item{userpath}{the path where the working folder is located}

\item{forcings}{list containing the time series in the odd positions and realted forcings in the even positions. Forcings imputted are: Water temperature [Celsius degrees] and feeding rate [g/individual x d]}
}
\value{
saves .nc; .csv and .asc outputs in the 'Outputs' folder
}
\description{
Solves the bioenergetic balance for Bream
}
