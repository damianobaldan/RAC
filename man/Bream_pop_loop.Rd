% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bream_pop_loop.R
\name{Bream_pop_loop}
\alias{Bream_pop_loop}
\title{Function that runs the Monte Carlo simulation for the Seabream population model}
\usage{
Bream_pop_loop(Param, Tint, Gint, Food, IC, times, N, userpath)
}
\arguments{
\item{Param}{a vector containing model parameters}

\item{Tint}{the interpolated water temperature time series}

\item{Gint}{the interpolated feeding rate time series}

\item{Food}{the food characterization}

\item{IC}{initial condition}

\item{times}{integration extremes and integration timestep}

\item{N}{time series with number of individuals}

\item{userpath}{the path where the working folder is located}
}
\value{
a list with RK solver outputs
}
\description{
Function that runs the Monte Carlo simulation for the Seabream population model
}
