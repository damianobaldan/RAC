#' Clam bioenergetic population model differential equations
#'
#' @param Param a vector containing model parameters
#' @param Tint the interpolated water temperature at time t
#' @param Phy the interpolated phytoplankton at time t
#' @param DT the interpolated detritus at time t
#' @param POCint the interpolated POC at time t
#' @param POMint the interpolated POM at time t
#' @param TSSint the interpolated TSS at time t
#' @param Wd the weight of the clam at time t
#'
#' @return a list containing the clam weights, temperature limitation functions and metabolic rates at time t
#'
#' @import matrixStats plotrix rstudioapi
#'
#' @import stats utils
#'

Clam_pop_equations <- function(Param, Tint, Phy, DT, POCint, POMint, TSSint, Wd){

# Parameters definition
  epsB=Param[1]            # [J/g] Tissue energy content
  epsDT=Param[2]           # [J/mgC] Detritus energy content
  epsPhy=Param[3]          # [J/mgC] Phytoplankton energy content
  epsO2=Param[4]           # [mlO2/h] Energy consumed by the respiration of 1g of oxygen
  alpha=Param[5]           # [-] Feeding catabolism
  CRmax=Param[6]           # [l/d gDM] Maximum filtration rate
  AEmax=Param[7]           # [-] Maximum adsorption efficiency
  Rmax=Param[8]            # [mgO2/d gDM] maximum respiration rate
  Amax=Param[9]            # [J/d  g] Maximum ingestion rate for 1g o mussel
  q=Param[10]              # [-] Weight exponent for filtration
  Ks=Param[11]             # [-] Half-saturation constant for AE
  betaa=Param[12]          # [1/Celsius degree] Temperature exponent fot anabolism
  betac=Param[13]          # [1/Celsius degree] Temperature exponent for catabolism
  Tma=Param[14]            # [Celsius degree] Maximum temperature for the anabolic process
  Toa=Param[15]            # [Celsius degree] Optimum temperature for the anabolic process
  Tmc=Param[16]            # [Celsius degree] Maximum temperature for the catabolic process
  Toc=Param[17]            # [Celsius degree] Optimum temperature for the catabolic process
  aF=Param[18]             # [-] Dry weight - wet weight conversion coefficient
  bF=Param[19]             # [-] Dry weight - wet weight exponent
  a=Param[20]              # [-] Dry weight - length conversion coefficient
  b=Param[21]              # [-] Dry weight - length exponent
  lambda=Param[22]         # [g/mg] Chlorophyll a - Phytoplankton conversion factor


# CATABOLISM
# Optimum temperature dependence for catabolism [dimensionless]
if (Tint >= Tmc) {
  fc=0.0
} else {
  fc=((Tmc-Tint)/(Tmc-Toc))^(betac*(Tmc-Toc))*exp(betac*(Tint-Toc))
}

C=Rmax*epsO2*fc*Wd      # Daily catabolism [J/day]

# ANABOLISM
# Optimum temperature dependence for anabolism [dimensionless]
if (Tint >= Tma){
  fa=0.0
}else {
  fa=((Tma-Tint)/(Tma-Toa))^(betaa*(Tma-Toa))*exp(betaa*(Tint-Toa))
}

I=(CRmax*fa*Wd^q)*(DT*epsDT+Phy*epsPhy) # Daily ingestion [J/day]

encont=(DT*epsDT+Phy*epsPhy)  # Energy content of ingested food [J/l]

Q=((POMint/TSSint)) # POM/TSS ratio [-]
if (Q>=1){
Q=1
}

AE=(Q/(Q+Ks))        # Limitation on ingested energy [-]
E=I*AE               # Total ingested energy [J/d]
Aing=(1-alpha)*E        # Daily anabolism [J/d]

# Daily anabolism limitation
if (Aing<Amax*Wd^q*fa){
A=Aing                  # Anabolic rate [J/d]
filt=CRmax*fa*Wd^q   # [l/d]
A1=Amax*Wd^q*fa      # Maximum anabolic rate [J/d]
}else{
  A=Amax*Wd^q*fa       # Anabolic rate [J/d]
filt=CRmax*fa*Wd^q   # [l/d]
A1=Amax*Wd^q*fa      # Maximum anabolic rate [J/d]
}

# Mass balance
dWd=((A)-C)/epsB     # Weight increment [g/d]

tfun=cbind(fa, fc)
metab=cbind(A, C)

# Function outputs
output=list(dWd,tfun,metab)
return(output)
}
