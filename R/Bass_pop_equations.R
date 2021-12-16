#' Seabass bioenergetic population model differential equations
#'
#' @param Param vector containing all metabolic parameters
#' @param N the number of individuals at time t
#' @param Temp water temperature forcing at time t
#' @param G food entering the cage at time t
#' @param Food food characterization (Proteins, Lipids, Carbohydrates)
#' @param weight individual weight at time t
#'
#' @return model output at time t
#'
#' @import matrixStats plotrix rstudioapi
#'

Bass_pop_equations <- function(Param, N, Temp, G, Food, weight){

  # Parameters definition
  # Parameters definition
  ingmax=Param[1]        # [g/d] Maximum ingestion rate
  alpha=Param[2]         # [-] Feeding catabolism coefficient
  betaprot=Param[3]      # [-] Assimilation coefficient for protein
  betalip=Param[4]       # [-] Assimilation coefficient for lipid
  betacarb=Param[5]      # [-] Assimilation coefficient for carbohydrates
  epsprot=Param[6]       # [J/gprot] Energy content of protein
  epslip=Param[7]        # [J/glip] Energy content of lipid
  epscarb=Param[8]       # [J/gcarb] Energy content of carbohydrate
  epsO2=Param[9]         # [J/gO2] Energy consumed by the respiration of 1g of oxygen
  pk=Param[10]           # [1/day] Temperature coefficient for the fasting catabolism
  k0=Param[11]           # [1/Celsius degree]  Fasting catabolism at 0 Celsius degree
  m=Param[12]            # [-] Weight exponent for the anabolism
  n=Param[13]            # [-] Weight exponent for the catabolism
  betac=Param[14]        # [-]  Shape coefficient for the H(Tw) function
  Tma=Param[15]          # [Celsius degree] Maximum lethal temperature for Dicentrarchus labrax
  Toa=Param[16]          # [Celsius degree] Optimal temperature for Dicentrarchus labrax
  Taa=Param[17]          # [Celsius degree] Lowest feeding temperature for Dicentrarchus labrax
  omega=Param[18]        # [gO2/g] Oxygen consumption - weight loss ratio
  a=Param[19]            # [J/gtissue] Energy content of fish tissue
  k=Param[20]            # [-] Weight exponent for energy content
  eff=Param[21]          # [-] Food ingestion efficiency

  # Food composition definition
  Pcont=Food[1]       # [-] Percentage of proteins in the food
  Lcont=Food[2]       # [-] Percentage of lipids in the food
  Ccont=Food[3]       # [-] Percentage of carbohydrates in the food


  # EQUATIONS

  # Forcing temperature
  fgT=((Tma-Temp)/(Tma-Toa))^(betac*(Tma-Toa))*exp(betac*(Temp-Toa)) # Optimum Temperature dependance for ingestion [-]
  frT= exp(pk*Temp)                                                  # Exponential Temperature dependance for catabolism [-]
  Tfun=cbind(fgT, frT)

  # Ingested mass
  ing=ingmax*(weight^m)*fgT  # Potential ingestion rate  [g/d]

  G=G*eff   # Food ingestion efficiency

  # Lowest feeding temperature threshold
  if (is.na(Temp)) print("Temp")

  if (Temp<Taa) {
    ing=0
  }

  # Available food limitation
  if (ing>G) {
  ingvero=G      # [g/d] Actual ingestion rate
  }
  else {
    ingvero=ing  # [g/d] Actual ingestion rate
  }

  # Energy content of somatic tissue [J/g]. Lupatsch et al. (2003)
  epstiss=a*weight^k

  # Ingested energy
  diet=Pcont*epsprot*betaprot+Lcont*epslip*betalip+Ccont*epscarb*betacarb  # [J/g]
  assE=ingvero*diet   # [J/d]

  # Compute excretion
  Pexc=(1-betaprot)*Pcont*ingvero*N/1e3  # Excreted proteins [kg/d]
  Lexc=(1-betalip)*Lcont*ingvero*N/1e3   # Excreted lipids [kg/d]
  Cexc=(1-betacarb)*Ccont*ingvero*N/1e3  # Excreted carbohydrates [kg/d]
  exc=cbind(Pexc,Lexc,Cexc)

  # Compute waste
  Pwst=((G/eff)-ingvero)*Pcont*N/1e3  # Proteins to waste [kg/d]
  Lwst=((G/eff)-ingvero)*Lcont*N/1e3  # Lipids to waste [kg/d]
  Cwst=((G/eff)-ingvero)*Ccont*N/1e3  # Carbohydrates to waste [kg/d]
  wst=cbind(Pwst,Cwst,Lwst)

  # Metabolism terms
  anab=assE*(1-alpha)                    # Net anabolic rate [J/d]
  catab=epsO2*k0*frT*(weight^n)*omega    # Fasting catabolic rate [J/d]
  metab=cbind(anab,catab)

  # O2 and NH4 produced
  O2=catab/epsO2*N/1e3          # O2 consumed [kg02/d]
  NH4=O2*0.06*N/1e3             # NH4 produced [kgN/d]

  # Mass balance
  dw = (anab-catab)/epstiss   # Weight increment [g/d]

  # Function outputs
  output=list(dw,exc,wst,ing,ingvero,Tfun,metab, O2, NH4)
  return(output)
}
