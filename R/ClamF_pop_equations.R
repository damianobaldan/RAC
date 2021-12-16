#' Clam bioenergetic population model differential equations (alternative version)
#'
#' @param Param a vector containing model parameters
#' @param Tint the interpolated water temperature at time t
#' @param Chlint the interpolated chlorophyll at time t
#' @param Ww clam wet weight at time t
#'
#' @return a list containing the clam weights, temperature limitation functions and metabolic rates at time t
#'
#' @import matrixStats plotrix rstudioapi
#'


ClamF_pop_equations <- function(Param, Tint, Chlint, Ww){

  # Parameters definition
  Gdmax=Param[1]           # [gDW^0.265 d^-1] Max growth rate on a dry weight base
  Gwmax=Param[2]           # [gDW^0.333 d^-1] Max growth rate on a wet weight base
  Rdmax=Param[3]           # [d^-1] Max respiration rate on a dry weight base
  Rwmax=Param[4]           # [d^-1] Max respiration rate on a wet weight base
  p=Param[5]               # [-] Dry weight - wet weight conversion exponent
  k=Param[6]               # [-] Wet weight - length conversion exponent
  q=Param[7]               # [-] Coefficient for allometric filter velocity
  a=Param[8]               # [cm] Wet weight - length conversion coefficient
  b=Param[9]               # [gWW] Dry weight - wet weight conversion coefficient
  m=Param[10]              # [-] Weight exponent for catabolism
  n=Param[11]              # [-] Weight exponent for anabolism
  vf=Param[12]             # [J/d g] Maximum ingestion rate for 1g o mussel
  epsT=Param[13]           # [-] Weight exponent for filtration
  epsF=Param[14]           # [-] Half-saturation constant for AE
  TmaxG=Param[15]          # [1/Celsius degree] Temperature exponent fot anabolism
  ToptG=Param[16]          # [1/Celsius degree] Temperature exponent for catabolism
  TmaxR=Param[17]          # [Celsius degree] Maximum temperature for the anabolic process
  ToptR=Param[18]          # [Celsius degree] Optimum temperature for the anabolic process
  TmaxV=Param[19]          # [Celsius degree] Maximum temperature for the catabolic process
  ToptV=Param[20]          # [Celsius degree] Optimum temperature for the catabolic process
  betaG=Param[21]          # [-] Dry weight - wet weight conversion coefficient
  betaV=Param[22]          # [-] Dry weight - wet weight exponent
  betaR=Param[23]          # [-] Dry weight - length conversion coefficient
  niT=Param[24]            # [-] Ratio of energetic demand satisfied by non planktonic material
  kap=Param[25]            # [mgChla/gC] Chlorophyll/Phytoplankton ratio


  Tanab=min(Tint,TmaxG)   # Temperature for anabolism < Max temperature for anabolism
  Tcatab=min(Tint,TmaxR)  # Temperature for catabolism < Max temperature for catabolism
  Tfood=min(Tint,TmaxV)   # Temperature for filtration < Max temperature for filtration

  # Temperature limitation functions
  FgT=(((TmaxG-Tanab)/(TmaxG-ToptG))^(betaG*(TmaxG-ToptG)))*exp(betaG*(Tanab-ToptG))    # Growth limitation [-]
  FrT=(((TmaxR-Tcatab)/(TmaxR-ToptR))^(betaR*(TmaxR-ToptR)))*exp(betaR*(Tcatab-ToptR))  # Respiration limitation [-]
  FvT=(((TmaxV-Tfood)/(TmaxV-ToptV))^(betaV*(TmaxV-ToptV)))*exp(betaV*(Tfood-ToptV))    # Filtration limitation [-]

  # Maximum food ingestion
  Wd=b*Ww^p                  # Dry weight value [g]
  Fstar=Gdmax*FgT*(Wd^(1-0.3333*p))*epsT/(vf*FvT*(Wd^q)*epsF)
  Ff=min(1, Chlint/Fstar+niT)

  # ENERGETIC BALANCE

  # Anabolism [gWW/d]
  if (Chlint>Fstar) {
  A=Gwmax*FgT*Ww^m
  }else{
    A=Ff*Gwmax*FgT*Ww^m
  }

  # Catabolism [gWW/d]
  C=Rwmax*FrT*Ww^n

  # Balance [gWW/d]
  dWw=A-C


tfun=cbind(FgT, FrT, FvT)
metab=cbind(A, C)

# Function outputs
output=list(dWw,tfun,metab)
return(output)
}
