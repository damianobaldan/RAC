#' Mussel bioenergetic individual model differential equations
#'
#' @param Param a vector containing model parameters
#' @param Tint the interpolated water temperature at time t
#' @param Phyint the interpolated phytoplankton at time t
#' @param DTint the interpolated detritus at time t
#' @param POCint the interpolated POC at time t
#' @param Ccont the C/C content of the POC at time t
#' @param Ncont the N/C content of POC at time t
#' @param Pcont the P/C content of POC at time t
#' @param POMint the interpolated POM at time t
#' @param TSSint the interpolated TSS at time t
#' @param Wb the somatic tissue dry weight at time t
#' @param R the gondadic tissue dry weight at time t
#' @param t the time
#' @param trip vector containing the flags with resting periods
#'
#' @return the outputs at time t
#'
#' @import matrixStats plotrix rstudioapi
#'

Mussel_ind_equations <- function(Param, Tint, Phyint, DTint, POCint, Ccont, Ncont, Pcont, POMint, TSSint, Wb,R,t,trip){

  # Parameters definition
  epsR=Param[1]             # [J/g] Somatic tissue energy content
  epsB=Param[2]             # [J/g] Gonadic tissue energy content
  epsDT=Param[3]            # [J/mgC] Detritus energy content
  epsPhy=Param[4]           # [J/mgC] Phytoplankton energy content
  epsO2=Param[5]            # [J/mgO2] Energy consumed by the respiration of 1g of oxygen
  alpha=Param[6]            # [-] Feeding catabolism
  CRmax=Param[7]            # [l/d gDM] Maximum filtration rate
  AEmax=Param[8]            # [-] Maximum adsorption efficiency
  Rmax=Param[9]             # [mgO2/d gDM] maximum respiration rate
  Amax=Param[10]            # [J/d g] Maximum ingestion rate for 1g o mussel
  q=Param[11]               # [-] Weight exponent for filtration
  n=Param[12]               # [-] Weight exponent for catabolism
  Ks=Param[13]              # [-] Half-saturation constant for AE
  betaa=Param[14]           # [1/Celsius degree] Temperature exponent fot anabolism
  betac=Param[15]           # [1/Celsius degree] Temperature exponent for catabolism
  Tma=Param[16]             # [Celsius degrees] Maximum temperature for the anabolic process
  Toa=Param[17]             # [Celsius degrees] Optimum temperature for the anabolic process
  Tmc=Param[18]             # [Celsius degrees] Maximum temperature for the catabolic process
  Toc=Param[19]             # [Celsius degrees] Optimum temperature for the catabolic process
  a=Param[20]               # [m] Weight to length proportionality constant
  b=Param[21]               # [-] Weight to length exponent
  k=Param[22]               # [-] Energy fraction used for reproduction
  epsNO=Param[23]           # [-] Nexcreted to oxygen consumed ratio
  ripi=Param[24]            # [d] Beginning of reproductory resting period
  ripf=Param[25]            # [d] End of reproductory resting period
  tspawn1=Param[26]         # [d] First spawning
  tspawn2=Param[27]         # [d] Second spawning
  aF=Param[28]              # [-] Dry weight - wet weight conversion coefficient
  atot=Param[29]            # [-] Dry weight - total (with shell) weight conversion coefficient
  CcontMyt=Param[30]        # [gC/gDW] Mussel C content
  NcontMyt=Param[31]        # [gN/gDW] Mussel N content
  PcontMyt=Param[32]        # [gP/gDW] Mussel P content
  lambda=Param[33]          # [-] Chlorophyll to Phytoplankton conversion factor
  NC_Fec=Param[34]          # [-] N/C ratio for feces
  PC_Fec=Param[35]          # [-] P/C ratio for feces
  gamma=Param[36]           # [-] weight exponent for anabolism limitation

  # Update time data
  ripi=trip[1]  # [d] Beginning of reproductory resting period
  ripf=trip[2]  # [d] End of reproductory resting period

  # Soft tissue dry weight [g]
  Wd=Wb+R

    # CATABOLISM
    # Optimum temperature dependence for catabolism [dimensionless]

    if (Tint>=Tmc) {
      fc=0.0
    } else {
      fc=((Tmc-Tint)/(Tmc-Toc))^(betac*(Tmc-Toc))*exp(betac*(Tint-Toc))
    }

    C=Rmax*fc*epsO2*Wb^n       # Daily catabolism [J/d]
    O2=(Rmax*fc*Wb^n)/1e3      # Oxygen consumed [gO2/d]
    NH4=(O2*epsNO)/1e3         # Nitrogen excreted [kgNH4-N/d]


    # ANABOLISM
    # Optimum temperature dependence for anabolism [dimensionless]

    if (Tint>=Tma)  {
      fa=0.0
    } else {
      fa=((Tma-Tint)/(Tma-Toa))^(betaa*(Tma-Toa))*exp(betaa*(Tint-Toa))
    }


    CR=(CRmax*fa*Wb^q)                                     # Actual clearance rate [l/d]
    I=(CRmax*fa*Wb^q)*(DTint*epsDT+Phyint*epsPhy)          # Daily ingestion  [J/d]
    I_P=(CRmax*fa*Wb^q)*(DTint*Pcont+Phyint*Pcont)/1e3     # Daily P ingestion P [g/d]
    I_N=(CRmax*fa*Wb^q)*(DTint*Ncont+Phyint*Ncont)/1e3     # Daily N ingestion [g/d]
    I_C=(CRmax*fa*Wb^q)*(DTint*Ccont+Phyint*Ccont)/1e3     # Daily C ingestion [g/d]

    Q=POMint/TSSint       # POM/TSS ratio [-]
    AE=AEmax*Q/(Q+Ks)     # Actual adsorption efficiency [-]

    E=I*AE                    # Total absorbed energy [J/d]
    Ex_C=I_C*(1-AE)           # C escretion [g/d]
    Ex_P=Ex_C*PC_Fec          # P excretion [g/d]
    Ex_N=Ex_C*NC_Fec          # N excretion [g/d]

    Aing=(1-alpha)*E          # Daily anabolism [J/d]
    Epspsf=((DTint/POCint*epsDT)+(Phyint/POCint*epsPhy))*1e3  # Energy content in pseudofecies [J/g]


    # LIMITATION on daily anabolism

    # Energetic treshold computation and pseudofecies production
    if (Aing<(Amax*fa*Wd^gamma))  {

      A=Aing
      Imass=(CRmax*fa*Wb^q)*(POCint)      # Daily POC ingestion [mg/d]
      ratio=NH4/(Imass*AE*Ncont/1e6)      # Excreted/ingested N [-]
      pseudof = 0.0                       # Pseudofecies production [gC/d]
      psC = 0.0                           # pseudofecies C content [gC/d]
      psN = 0.0                           # pseudofecies N content [gN/d]
      psP = 0.0                           # pseudofecies P content [gP/d]

    } else  {

      A=Amax*fa*Wd^gamma
      Imass=(CRmax*fa*Wb^q)*(POCint)       # Daily POC ingestion [mg/d]
      ratio=NH4/(Imass*AE*Ncont/1e6)       # Excreted/ingested N [-]
      pseudof = (Aing-A)/Epspsf            # Pseudofecies production [gC/d]
      psC = pseudof*Ccont             # pseudofecies C content [gC/d]
      psN = pseudof*Ncont              # pseudofecies N content [gN/d]
      psP = pseudof*Pcont             # pseudofecies P content [gP/d]

    }

    if ((t<=ripi)||(t>=ripf)) {    # If mussel is not in the resting period for gonadic tissue production

      if (A>C)    {              # If daily anabolism is greater than daily catabolism
        dWb=(A-C)*(1-k)/epsB     # Daily somatic tissue weight increment [g/d]
        dR=k*(A-C)/epsR          # Daily gonadic tissue weight increment [g/d]
      } else {                   # If daily anabolism is lower than daily catabolism
        dWb=(A-C)/epsB
        dR=0
      }
    } else {               # If mussel is in the resting period for gonadic tissue production
      dWb=(A-C)/epsB       # Daily somatic tissue weight increment [g/d]
      dR=0                 # Daily gonadic tissue weight increment [g/d]
    }

  # CNP mussel composition
  Cmyt=(Wd*CcontMyt)    # C content of mussel [g]
  Nmyt=(Wd*NcontMyt)    # N content of mussel [g]
  Pmyt=(Wd*PcontMyt)    # P content of mussel [g]

  # Other outputs of the function
  pfec=cbind(psC, psN, psP)
  fec=cbind(Ex_C, Ex_N, Ex_P)
  cont=cbind(Cmyt,Nmyt,Pmyt)
  tfun=cbind(fa, fc)
  metab=cbind(A, C)
  cons=O2
  amm=NH4*1e3

  # Assign outputs to a list
  output=list(dWb,dR,pfec,fec,cont,tfun,metab,cons,amm)
  return(output) # Mussel_ind_equations output

} # end function
