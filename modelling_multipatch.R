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
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - 
      psi*(betaHD*Isa*Sh) - psi*(1-alpha)*(betaHD*Ira*Sh) - 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh = betaHH*Ish*Sh + psi*betaHD*Isa*Sh + (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*psi*(betaHD*Ira*Sh) + (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}



# Test Run ----------------------------------------------------------------

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
parmtau <- seq(0,0.035, by = 0.002)

output1 <- data.frame(matrix(ncol = 7, nrow = length(parmtau)))

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(nrow = 1, ncol =7))
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = m_betaAA, betaHH = 0.00001, tau = parmtau[i],
             betaHI = (0.00001), betaHD = (0.00001),  phi = m_phi, theta = m_theta, alpha = m_alpha, zeta = m_zeta, 
             psi = 0.7, fracimp = m_fracimp, propres_imp = m_propres_imp)
  
  out <- ode(y = init, func = amrimp, times = times, parms = parms2)
  temp[,1] <- parmtau[i]
  temp[,2] <- rounding(out[nrow(out),5]) 
  temp[,3] <- rounding(out[nrow(out),6]) 
  temp[,4] <- rounding(out[nrow(out),7])
  temp[,5] <- temp[3] + temp[4]
  temp[,6] <- signif(as.numeric(temp[4]/temp[5]), digits = 3)
  temp[,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  output1[i,] <- temp
}

colnames(output1) <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA")

plotdata <- melt(output1, id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 
