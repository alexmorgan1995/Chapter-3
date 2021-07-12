library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

# Single Model ------------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    
    dIsa = betaAA*Isa*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(1-Sh) - 
      
      betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshH) - 
      (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhH) - 
      
      psi*(betaHD*Isa*Sh) - 
      psi*(1-alpha)*(betaHD*Ira*Sh) - 
      
      (1-psi)*(imp1)*(betaHI1*fracimp1*(1-propres_imp1)*Sh) - 
      (1-psi)*(imp1)*(1-alpha)*(betaHI1*fracimp1*propres_imp1*Sh)-
      
      (1-psi)*(imp2)*(betaHI2*fracimp2*(1-propres_imp2)*Sh) - 
      (1-psi)*(imp2)*(1-alpha)*(betaHI2*fracimp2*propres_imp2*Sh) - 
      
      (1-psi)*(imp3)*(betaHI3*fracimp3*(1-propres_imp3)*Sh) - 
      (1-psi)*(imp3)*(1-alpha)*(betaHI3*fracimp3*propres_imp3*Sh) - uh*Sh

    
    dIshDA = psi*betaHD*Isa*Sh  - rh*IshDA - uh*IshDA 
    dIrhDA = psi*(1-alpha)*betaHD*Ira*Sh - rh*IrhDA - uh*IrhDA  
    
    dIshA1 = (1-psi)*(imp1)*(betaHI1*fracimp1*(1-propres_imp1)*Sh) - rh*IshA1 - uh*IshA1 
    dIrhA1 = (1-psi)*(imp1)*(1-alpha)*(betaHI1*fracimp1*propres_imp1*Sh) - rh*IrhA1 - uh*IrhA1  
    
    dIshA2 = (1-psi)*(imp2)*(betaHI2*fracimp2*(1-propres_imp2)*Sh) - rh*IshA2 - uh*IshA2 
    dIrhA2 = (1-psi)*(imp2)*(1-alpha)*(betaHI2*fracimp2*propres_imp2*Sh)  - rh*IrhA2 - uh*IrhA2  
    
    dIshA3 = (1-psi)*(imp3)*(betaHI3*fracimp3*(1-propres_imp3)*Sh) - rh*IshA3 - uh*IshA3 
    dIrhA3 = (1-psi)*(imp3)*(1-alpha)*(betaHI3*fracimp3*propres_imp3*Sh) - rh*IrhA3 - uh*IrhA3  
    
    dIshH = betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshH) - rh*IshH - uh*IshH 
    dIrhH = (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhH)- rh*IrhH - uh*IrhH 
    
    return(list(c(dSa,dIsa,dIra,
                  dSh,
                  dIshDA,dIrhDA,
                  dIshA1,dIrhA1,
                  dIshA2,dIrhA2,
                  dIshA3,dIrhA3,
                  dIshH, dIrhH)))
  })
}

# Test Run ----------------------------------------------------------------

times <- seq(0,30000, by = 1) 

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshH = 0, IrhH = 0)

parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0001, betaHH = 0.00001, tau = 0.0123,
           betaHD = (0.00001),  
           betaHI1 = (0.00001), betaHI2 = (0.00001), betaHI3 = (0.00001),
           phi = 0.03, theta = 1.1, alpha = 0.2, zeta = 0.0001, psi = 0.7, 
           fracimp1 = 0.5, fracimp2 = 0.5, fracimp3 = 0.5,
           imp1 = 0.333, imp2 = 0.333, imp3 = 0.333,
           propres_imp1 = 0.5, propres_imp2 = 0.5, propres_imp3 = 0.5)


out <- ode(y = init, func = amrimp, times = times, parms = parms)

