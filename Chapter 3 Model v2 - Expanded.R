setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Mathematical Models")

rm(list=ls())

library("deSolve")
library("fast")
library("sensitivity")
library("ggplot2")
library("plotly")
library("tidyr")
library("nlmeODE")
library("phaseR")
library("reshape2")
library("ggrepel")

#### Model For Domestic, EU and non-EU Livestock-to-human AMR transmission ####

amrhetcountry <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSUKa = - betaUKaa*(IUKas + IUKar*alpha)*SUKa + rA*(IUKas+IUKar) + tauUK*IUKas + eta - eta*SUKa
    dIUKas = betaUKaa*IUKas*SUKa - rA*IUKas - tauUK*IUKas - tauUK*theta*IUKas + phi*IUKar - eta*IUKas
    dIUKar = betaUKaa*IUKar*SUKa*alpha + tauUK*theta*IUKas - phi*IUKar - rA*IUKar - eta*IUKar 
    
    dSNEDa = - betaNEDaa*(INEDas + INEDar*alpha)*SNEDa + rA*(INEDas+INEDar) + tauNED*INEDas + eta - eta*SNEDa
    dINEDas = betaNEDaa*INEDas*SNEDa - rA*INEDas - tauNED*INEDas - tauNED*theta*INEDas + phi*INEDar - eta*INEDas
    dINEDar = betaNEDaa*INEDar*SNEDa*alpha + tauNED*theta*INEDas - phi*INEDar - rA*INEDar - eta*INEDar 
    
    dSPOLa = - betaPOLaa*(IPOLas + IPOLar*alpha)*SPOLa + rA*(IPOLas+IPOLar) + tauPOL*IPOLas + eta - eta*SPOLa
    dIPOLas = betaPOLaa*IPOLas*SPOLa - rA*IPOLas - tauPOL*IPOLas - tauPOL*theta*IPOLas + phi*IPOLar - eta*IPOLas
    dIPOLar = betaPOLaa*IPOLar*SPOLa*alpha + tauPOL*theta*IPOLas - phi*IPOLar - rA*IPOLar - eta*IPOLar 
    
    dSIREa = - betaIREaa*(IIREas + IIREar*alpha)*SIREa + rA*(IIREas+IIREar) + tauIRE*IIREas + eta - eta*SIREa
    dIIREas = betaIREaa*IIREas*SIREa - rA*IIREas - tauIRE*IIREas - tauIRE*theta*IIREas + phi*IIREar - eta*IIREas
    dIIREar = betaIREaa*IIREar*SIREa*alpha + tauIRE*theta*IIREas - phi*IIREar - rA*IIREar - eta*IIREar 
    
    dSBELa = - betaBELaa*(IBELas + IBELar*alpha)*SBELa + rA*(IBELas+IBELar) + tauBEL*IBELas + eta - eta*SBELa
    dIBELas = betaBELaa*IBELas*SBELa - rA*IBELas - tauBEL*IBELas - tauBEL*theta*IBELas + phi*IBELar - eta*IBELas
    dIBELar = betaBELaa*IBELar*SBELa*alpha + tauBEL*theta*IBELas - phi*IBELar - rA*IBELar - eta*IBELar 
    
    dSGERa = - betaGERaa*(IGERas + IGERar*alpha)*SGERa + rA*(IGERas+IGERar) + tauGER*IGERas + eta - eta*SGERa
    dIGERas = betaGERaa*IGERas*SGERa - rA*IGERas - tauGER*IGERas - tauGER*theta*IGERas + phi*IGERar - eta*IGERas
    dIGERar = betaGERaa*IGERar*SGERa*alpha + tauGER*theta*IGERas - phi*IGERar - rA*IGERar - eta*IGERar 
    
    dSITAa = - betaITAaa*(IITAas + IITAar*alpha)*SITAa + rA*(IITAas+IITAar) + tauITA*IITAas + eta - eta*SITAa
    dIITAas = betaITAaa*IITAas*SITAa - rA*IITAas - tauITA*IITAas - tauITA*theta*IITAas + phi*IITAar - eta*IITAas
    dIITAar = betaITAaa*IITAar*SITAa*alpha + tauITA*theta*IITAas - phi*IITAar - rA*IITAar - eta*IITAar 
    
    dSFRAa = - betaFRAaa*(IFRAas + IFRAar*alpha)*SFRAa + rA*(IFRAas+IFRAar) + tauFRA*IFRAas + eta - eta*SFRAa
    dIFRAas = betaFRAaa*IFRAas*SFRAa - rA*IFRAas - tauFRA*IFRAas - tauFRA*theta*IFRAas + phi*IFRAar - eta*IFRAas
    dIFRAar = betaFRAaa*IFRAar*SFRAa*alpha + tauFRA*theta*IFRAas - phi*IFRAar - rA*IFRAar - eta*IFRAar 
    
    dSSPAa = - betaSPAaa*(ISPAas + ISPAar*alpha)*SSPAa + rA*(ISPAas+ISPAar) + tauSPA*ISPAas + eta - eta*SSPAa
    dISPAas = betaSPAaa*ISPAas*SSPAa - rA*ISPAas - tauSPA*ISPAas - tauSPA*theta*ISPAas + phi*ISPAar - eta*ISPAas
    dISPAar = betaSPAaa*ISPAar*SSPAa*alpha + tauSPA*theta*ISPAas - phi*ISPAar - rA*ISPAar - eta*ISPAar 
    
    dSROMa = - betaROMaa*(IROMas + IROMar*alpha)*SROMa + rA*(IROMas+IROMar) + tauROM*IROMas + eta - eta*SROMa
    dIROMas = betaROMaa*IROMas*SROMa - rA*IROMas - tauROM*IROMas - tauROM*theta*IROMas + phi*IROMar - eta*IROMas
    dIROMar = betaROMaa*IROMar*SROMa*alpha + tauROM*theta*IROMas - phi*IROMar - rA*IROMar - eta*IROMar 
    
    dSTHAa = - betaTHAaa*(ITHAas + ITHAar*alpha)*STHAa + rA*(ITHAas+ITHAar) + tauTHA*ITHAas + eta - eta*STHAa
    dITHAas = betaTHAaa*ITHAas*STHAa - rA*ITHAas - tauTHA*ITHAas - tauTHA*theta*ITHAas + phi*ITHAar - eta*ITHAas
    dITHAar = betaTHAaa*ITHAar*STHAa*alpha + tauTHA*theta*ITHAas - phi*ITHAar - rA*ITHAar - eta*ITHAar 
    
    dSBRAa = - betaBRAaa*(IBRAas + IBRAar*alpha)*SBRAa + rA*(IBRAas+IBRAar) + tauBRA*IBRAas + eta - eta*SBRAa
    dIBRAas = betaBRAaa*IBRAas*SBRAa - rA*IBRAas - tauBRA*IBRAas - tauBRA*theta*IBRAas + phi*IBRAar - eta*IBRAas
    dIBRAar = betaBRAaa*IBRAar*SBRAa*alpha + tauBRA*theta*IBRAas - phi*IBRAar - rA*IBRAar - eta*IBRAar 
    
    dSCHIa = - betaCHIaa*(ICHIas + ICHIar*alpha)*SCHIa + rA*(ICHIas+ICHIar) + tauCHI*ICHIas + eta - eta*SCHIa
    dICHIas = betaCHIaa*ICHIas*SCHIa - rA*ICHIas - tauCHI*ICHIas - tauCHI*theta*ICHIas + phi*ICHIar - eta*ICHIas
    dICHIar = betaCHIaa*ICHIar*SCHIa*alpha + tauCHI*theta*ICHIas - phi*ICHIar - rA*ICHIar - eta*ICHIar 
    
    dSh = - betaUKha*(eta*IUKas)*Sh*psiUK - betaUKha*(eta*IUKar)*Sh*psiUK*alpha - 
      betaNEDha*(eta*INEDas)*Sh*psiNED - betaNEDha*(eta*INEDar)*Sh*psiNED*alpha - 
      betaPOLha*(eta*IPOLas)*Sh*psiPOL - betaPOLha*(eta*IPOLar)*Sh*psiPOL*alpha -
      betaIREha*(eta*IIREas)*Sh*psiIRE - betaIREha*(eta*IIREar)*Sh*psiIRE*alpha - 
      betaBELha*(eta*IBELas)*Sh*psiBEL - betaBELha*(eta*IBELar)*Sh*psiBEL*alpha - 
      betaGERha*(eta*IGERas)*Sh*psiGER - betaGERha*(eta*IGERar)*Sh*psiGER*alpha -
      betaITAha*(eta*IITAas)*Sh*psiITA - betaITAha*(eta*IITAar)*Sh*psiITA*alpha -
      betaFRAha*(eta*IFRAas)*Sh*psiFRA - betaFRAha*(eta*IFRAar)*Sh*psiFRA*alpha -
      betaSPAha*(eta*ISPAas)*Sh*psiSPA - betaSPAha*(eta*ISPAar)*Sh*psiSPA*alpha -
      betaROMha*(eta*IROMas)*Sh*psiROM - betaSPAha*(eta*ISPAar)*Sh*psiSPA*alpha - 
      betaTHAha*(eta*ITHAas)*Sh*psiTHA - betaTHAha*(eta*ITHAar)*Sh*psiTHA*alpha -
      betaBRAha*(eta*IBRAas)*Sh*psiBRA - betaBRAha*(eta*IBRAar)*Sh*psiBRA*alpha -
      betaCHIha*(eta*ICHIas)*Sh*psiCHI - betaCHIha*(eta*ICHIar)*Sh*psiCHI*alpha +
      rH*(IUKhs + IUKhr + INEDhs + INEDhr + IPOLhs + IPOLhr + IIREhs + IIREhr + IBELhs + IBELhr +
            IGERhs + IGERhr + IITAhs + IITAhr + IFRAhs + IFRAhr + ISPAhs + ISPAhr + IROMhs + IROMhr + 
            ITHAhs + ITHAhr + IBRAhs + IBRAhr + ICHIhs + ICHIhr) + mu - mu*Sh
      
    dIUKhs = betaUKha*(eta*IUKas)*Sh*psiUK - rH*IUKhs - mu*IUKhs 
    dIUKhr = betaUKha*(eta*IUKar)*Sh*psiUK*alpha - rH*IUKhr - mu*IUKhr 
    
    dINEDhs = betaNEDha*(eta*INEDas)*Sh*psiNED - rH*INEDhs - mu*INEDhs 
    dINEDhr = betaNEDha*(eta*INEDar)*Sh*psiNED*alpha - rH*INEDhr - mu*INEDhr 
    
    dIPOLhs = betaPOLha*(eta*IPOLas)*Sh*psiPOL - rH*IPOLhs - mu*IPOLhs 
    dIPOLhr = betaPOLha*(eta*IPOLar)*Sh*psiPOL*alpha - rH*IPOLhr - mu*IPOLhr 
    
    dIIREhs = betaIREha*(eta*IIREas)*Sh*psiIRE - rH*IIREhs - mu*IIREhs 
    dIIREhr = betaIREha*(eta*IIREar)*Sh*psiIRE*alpha - rH*IIREhr - mu*IIREhr 
    
    dIBELhs = betaBELha*(eta*IBELas)*Sh*psiBEL - rH*IBELhs - mu*IBELhs 
    dIBELhr = betaBELha*(eta*IBELar)*Sh*psiBEL*alpha - rH*IBELhr - mu*IBELhr 
    
    dIGERhs = betaGERha*(eta*IGERas)*Sh*psiGER - rH*IGERhs - mu*IGERhs 
    dIGERhr = betaGERha*(eta*IGERar)*Sh*psiGER*alpha - rH*IGERhr - mu*IGERhr 
    
    dIITAhs = betaITAha*(eta*IITAas)*Sh*psiITA - rH*IITAhs - mu*IITAhs 
    dIITAhr = betaITAha*(eta*IITAar)*Sh*psiITA*alpha - rH*IITAhr - mu*IITAhr 
    
    dIFRAhs = betaFRAha*(eta*IFRAas)*Sh*psiFRA - rH*IFRAhs - mu*IFRAhs 
    dIFRAhr = betaFRAha*(eta*IFRAar)*Sh*psiFRA*alpha - rH*IFRAhr - mu*IFRAhr 
    
    dISPAhs = betaSPAha*(eta*ISPAas)*Sh*psiSPA - rH*ISPAhs - mu*ISPAhs 
    dISPAhr = betaSPAha*(eta*ISPAar)*Sh*psiSPA*alpha - rH*ISPAhr - mu*ISPAhr 
    
    dIROMhs = betaROMha*(eta*IROMas)*Sh*psiROM - rH*IROMhs - mu*IROMhs 
    dIROMhr = betaROMha*(eta*IROMar)*Sh*psiROM*alpha - rH*IROMhr - mu*IROMhr 
    
    dITHAhs = betaTHAha*(eta*ITHAas)*Sh*psiTHA - rH*ITHAhs - mu*ITHAhs 
    dITHAhr = betaTHAha*(eta*ITHAar)*Sh*psiTHA*alpha - rH*ITHAhr - mu*ITHAhr 
    
    dIBRAhs = betaBRAha*(eta*IBRAas)*Sh*psiBRA - rH*IBRAhs - mu*IBRAhs 
    dIBRAhr = betaBRAha*(eta*IBRAar)*Sh*psiBRA*alpha - rH*IBRAhr - mu*IBRAhr 
    
    dICHIhs = betaCHIha*(eta*ICHIas)*Sh*psiCHI - rH*ICHIhs - mu*ICHIhs 
    dICHIhr = betaCHIha*(eta*ICHIar)*Sh*psiCHI*alpha - rH*ICHIhr - mu*ICHIhr 
    
    return(list(c(dSUKa, dIUKas, dIUKar,dSNEDa, dINEDas, dINEDar, dSPOLa, dIPOLas, dIPOLar, dSIREa, dIIREas, dIIREar,
                  dSBELa, dIBELas, dIBELar, dSGERa, dIGERas, dIGERar, dSITAa, dIITAas, dIITAar,
                  dSFRAa, dIFRAas, dIFRAar, dSSPAa, dISPAas, dISPAar, dSROMa, dIROMas, dIROMar,
                  dSTHAa, dITHAas, dITHAar, dSBRAa, dIBRAas, dIBRAar, dSCHIa, dICHIas, dICHIar,
                  dSh, dIUKhs, dIUKhr, dINEDhs, dINEDhr, dIPOLhs, dIPOLhr, dIIREhs, dIIREhr, 
                  dIBELhs, dIBELhr, dIGERhs, dIGERhr, dIITAhs, dIITAhr, dIFRAhs, dIFRAhr, dISPAhs, dISPAhr,
                  dIROMhs, dIROMhr, dITHAhs, dITHAhr, dIBRAhs, dIBRAhr, dICHIhs, dICHIhr
                  )))
  })
}

rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Parameters ####

init <- c(SUKa = 0.98, IUKas = 0.01, IUKar = 0.01, SNEDa = 0.98, INEDas = 0.01, INEDar = 0.01, 
          SPOLa = 0.98, IPOLas = 0.01, IPOLar = 0.01, SIREa = 0.98, IIREas = 0.01, IIREar = 0.01, 
          SBELa = 0.98, IBELas = 0.01, IBELar = 0.01, 
          SGERa = 0.98, IGERas = 0.01, IGERar = 0.01, SITAa = 0.98, IITAas = 0.01, IITAar = 0.01,
          SFRAa = 0.98, IFRAas = 0.01, IFRAar = 0.01, SSPAa = 0.98, ISPAas = 0.01, ISPAar = 0.01, 
          SROMa = 0.98, IROMas = 0.01, IROMar = 0.01, STHAa = 0.98, ITHAas = 0.01, ITHAar = 0.01, 
          SBRAa = 0.98, IBRAas = 0.01, IBRAar = 0.01, SCHIa = 0.98, ICHIas = 0.01, ICHIar = 0.01,
          Sh = 1, IUKhs = 0, IUKhr = 0, INEDhs = 0, INEDhr = 0, IPOLhs = 0, IPOLhr = 0, IIREhs = 0, IIREhr = 0, 
          IBELhs = 0, IBELhr = 0, IGERhs = 0, IGERhr = 0, IITAhs = 0, IITAhr = 0, IFRAhs = 0, IFRAhr = 0, ISPAhs = 0, 
          ISPAhr = 0, IROMhs = 0, IROMhr = 0, ITHAhs = 0, ITHAhr = 0, IBRAhs = 0, IBRAhr = 0, ICHIhs = 0, ICHIhr = 0)

times1 <- seq(0,10000,by=1)

#Need to Specify Model Parameters

parms = c(betaUKaa = 0.1, betaNEDaa = 0.1, betaPOLaa = 0.1,  betaIREaa = 0.1,
          betaBELaa = 0.1, betaGERaa = 0.1, betaITAaa = 0.1, 
          betaFRAaa = 0.1, betaSPAaa = 0.1, betaROMaa = 0.1, 
          betaTHAaa = 0.1, betaBRAaa = 0.1, betaCHIaa = 0.1, 
          
          betaUKha = 0.01, betaNEDha = 0.01, betaPOLha = 0.01, betaIREha = 0.01, betaBELha = 0.01, betaGERha = 0.01, betaITAha = 0.01, 
          betaFRAha = 0.01, betaSPAha = 0.01, betaROMha = 0.01, betaTHAha = 0.01, betaBRAha = 0.01, betaCHIha = 0.01, 
          
          tauUK = 0.039, tauNED = 0.038, tauPOL = 0.067, tauIRE = 0.046, tauBEL = 0.089, tauGER = 0.101, tauITA = 0.157, tauFRA = 0.41, 
          tauSPA = 0.182, tauROM = 0.068, tauTHA = 0.071, tauBRA = 0.053, tauCHI = 0.319,  
          theta = 0.5, phi = 0.05,
          psiUK = 0.593, psiNED = 0.1177, psiPOL = 0.0661, psiIRE = 0.0403, 
          psiBEL = 0.0215, psiGER = 0.0295, psiITA = 0.0065, psiFRA = 0.0180, psiSPA = 0.0051, 
          psiROM = 0.0051, psiTHA = 0.0726, psiBRA = 0.0214, psiCHI = 0.0032, 
          alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)

out <- ode(y = init, func = amrhetcountry, times = times1, parms = parms)

outdata <- data.frame(out)
outdata$UKIComb <- outdata$IUKhs + outdata$IUKhr 
outdata$NEDIComb <- outdata$INEDhs + outdata$INEDhr 
outdata$POLIComb <- outdata$IPOLhs + outdata$IPOLhr
outdata$IREIComb <- outdata$IIREhs + outdata$IIREhr 
outdata$BELIComb <- outdata$IBELhs + outdata$IBELhr 
outdata$GERIComb <- outdata$IGERhs + outdata$IGERhr
outdata$ITAIComb <- outdata$IITAhs + outdata$IITAhr 
outdata$FRAIComb <- outdata$IFRAhs + outdata$IFRAhr 
outdata$SPAIComb <- outdata$ISPAhs + outdata$ISPAhr
outdata$ROMIComb <- outdata$IROMhs + outdata$IROMhr
outdata$THAIComb <- outdata$ITHAhs + outdata$ITHAhr 
outdata$BRAIComb <- outdata$IBRAhs + outdata$IBRAhr 
outdata$CHIIComb <- outdata$ICHIhs + outdata$ICHIhr

#Creating L I Q U I D code - i.e - Melting for GGPLOT
outdata1 <- outdata

outdata1[41:80] <- outdata[41:80]*100000 # So only human compartments are scaled to per 100,000 population 
meltedout <- melt(outdata1, id = "time", variable.name = "Compartment", value.name = "Value")

#Manipulating the Data for a stacked bar plot - to show the relative proportions from each country
#The bar chart plot will need to be when the model is at equilbirum - length(outdata) or something similar 
equidf <- data.frame(matrix(NA, ncol = 4, nrow = 3))
colnames(equidf) <- c("Category","Domestic","European", "nonEuropean")
equidf[1] <- c("ICombH", "Sensitive", "Resistant")
equidf[2] <- c(outdata$UKIComb[length(outdata$UKIComb)], outdata$IUKhs[length(outdata$IUKhs)], outdata$IUKhr[length(outdata$IUKhr)]) 
equidf[3] <- c(outdata$NEDIComb[length(outdata$NEDIComb)] + outdata$POLIComb[length(outdata$POLIComb)] + 
                 outdata$IREIComb[length(outdata$IREIComb)] + outdata$BELIComb[length(outdata$BELIComb)] + 
                 outdata$GERIComb[length(outdata$GERIComb)] + outdata$ITAIComb[length(outdata$ITAIComb)] + 
                 outdata$FRAIComb[length(outdata$FRAIComb)] + outdata$SPAIComb[length(outdata$SPAIComb)] + 
                 outdata$ROMIComb[length(outdata$ROMIComb)], 
               outdata$INEDhs[length(outdata$INEDhs)] + outdata$IPOLhs[length(outdata$IPOLhs)] + outdata$IIREhs[length(outdata$IIREhs)] + 
                 outdata$IBELhs[length(outdata$IBELhs)] + outdata$IGERhs[length(outdata$IGERhs)] + outdata$IITAhs[length(outdata$IITAhs)] +
                 outdata$IFRAhs[length(outdata$IFRAhs)] + outdata$ISPAhs[length(outdata$ISPAhs)] + outdata$IROMhs[length(outdata$IROMhs)], 
               outdata$INEDhr[length(outdata$INEDhr)] + outdata$IPOLhr[length(outdata$IPOLhr)] + outdata$IIREhr[length(outdata$IIREhr)] + 
                 outdata$IBELhr[length(outdata$IBELhr)] + outdata$IGERhr[length(outdata$IGERhr)] + outdata$IITAhr[length(outdata$IITAhr)] +
                 outdata$IFRAhr[length(outdata$IFRAhr)] + outdata$ISPAhr[length(outdata$ISPAhr)] + outdata$IROMhr[length(outdata$IROMhr)])
equidf[4] <- c(outdata$THAIComb[length(outdata$THAIComb)] + outdata$BRAIComb[length(outdata$BRAIComb)] + 
                 outdata$CHIIComb[length(outdata$CHIIComb)], 
               outdata$ITHAhs[length(outdata$ITHAhs)] + outdata$IBRAhs[length(outdata$IBRAhs)] + outdata$ICHIhs[length(outdata$ICHIhs)], 
               outdata$ITHAhr[length(outdata$ITHAhr)] + outdata$IBRAhr[length(outdata$IBRAhr)] + outdata$ICHIhr[length(outdata$ICHIhr)])  

###Plotting Output for Basic Script####

#Bar Chart Plot - to view the relative infecteds 
plot_ly(equidf, x = ~Category, y = ~(Domestic*100000), type = "bar", name = "Domestic") %>%
  add_trace(y = ~(European*100000), name = "European") %>%
  add_trace(y = ~(nonEuropean*100000), name = "nonEuropean") %>%
  layout(yaxis = list(title = "Proportion of Infecteds (per 100,000 Population", range = c(0, 10)), barmode = "stack")


####TO SHOW THE RELATIONSHIP BETWEEN ANTIBIOTIC USAGE AND RESISTANCE ####
#Simple Model to View Tau Relationship
amranimal <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSUKa = - betaUKaa*(IUKas + IUKar*alpha)*SUKa + rA*(IUKas+IUKar) + tauUK*IUKas + eta - eta*SUKa
    dIUKas = betaUKaa*IUKas*SUKa - rA*IUKas - tauUK*IUKas - tauUK*theta*IUKas + phi*IUKar - eta*IUKas
    dIUKar = betaUKaa*IUKar*SUKa*alpha + tauUK*theta*IUKas - phi*IUKar - rA*IUKar - eta*IUKar 
    return(list(c(dSUKa, dIUKas, dIUKar
    )))
  })
}

parmtau <- seq(0,0.35,by=0.01)
init <- c(SUKa = 0.98, IUKas = 0.01, IUKar = 0.01)
parms = c(betaUKaa = 0.1, tauUK = 0.039, theta = 0.5, phi = 0.05,  alpha = 0.8, rA = 25^-1, eta = 240^-1)
out <- ode(y = init, func = amranimal, times = times1, parms = parms)
outputnew <- data.frame()
output1 <- data.frame()
times1 <- seq(0,1000,by=1)

countryantib <- c(0.039, 0.038, 0.067, 0.046, 0.089, 0.101,0.157, 0.041, 0.182, 0.068, 0.071, 0.053, 0.319)

for (i in 1:length(countryantib)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=6))
  parms2 = c(betaUKaa = 0.1, tauUK = countryantib[i],  
             theta = 0.5, phi = 0.05, 
             alpha = 0.8, rA = 25^-1, eta = 240^-1)
  out <- ode(y = init, func = amranimal, times = times, parms = parms2)
  temp[1,1] <- countryantib[i]
  temp[1,2] <- rounding(out[nrow(out),2]) 
  temp[1,3] <- rounding(out[nrow(out),3]) 
  temp[1,4] <- rounding(out[nrow(out),4])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  print(temp[1,1])
  outputnew <- rbind.data.frame(outputnew, temp)
}

colnames(outputnew)[1:6] <- c("tau", "SuscAnimals","InfUKAnimals","ResUKInfAnimals","IUKCombA","IUKResRatA")
outputnew$country <- c("UK", "NED", "POL", "IRE", "BEL", "GER", "ITA", "FRA", 
                       "SPA", "ROM", "THA", "BRA", "CHI")
for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=6))
  parms2 = c(betaUKaa = 0.1, tauUK = parmtau[i],  
             theta = 0.5, phi = 0.05, 
             alpha = 0.8, rA = 25^-1, eta = 240^-1)
  out <- ode(y = init, func = amranimal, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),2]) 
  temp[1,3] <- rounding(out[nrow(out),3]) 
  temp[1,4] <- rounding(out[nrow(out),4])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  print(temp[1,1])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:6] <- c("tau", "SuscAnimals","InfUKAnimals","ResUKInfAnimals","IUKCombA","IUKResRatA")

#Resistance
ggplot(data = outputnew, aes(x = tau, y = IUKResRatA)) + geom_point() + 
  geom_text_repel(data = outputnew, mapping=aes(x=outputnew$tau, y= outputnew$IUKResRatA,label=outputnew$country), 
                  size=4, box.padding = unit(1.5, "lines")) +
  geom_line(data = output1, aes(x = output1$tau, output1$IUKResRatA, col = "red")) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") + 
  scale_y_continuous(limits = c(0,1), expand = c(0,0))

#Foodborne Disease Carriage 
ggplot(data = outputnew, aes(x = tau, y = IUKCombA)) + geom_point() + 
  geom_text_repel(data = outputnew, mapping=aes(x=outputnew$tau, y= outputnew$IUKCombA,label=outputnew$country), 
                  size=4, box.padding = unit(2, "lines")) +
  geom_line(data = output1, aes(x = output1$tau, output1$IUKCombA, col = "red")) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

#Complex Model to View Relationship

parmtau <- seq(0,0.35,by=0.01)
init1 <- init
outputlambda <- data.frame()
times <- seq(0, 100000, by = 10)

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=5))
  parms2 = c(betaUKaa = 0.1, betaNEDaa = 0.1, betaPOLaa = 0.1,  betaIREaa = 0.1,
             betaBELaa = 0.1, betaGERaa = 0.1, betaITAaa = 0.1, 
             betaFRAaa = 0.1, betaSPAaa = 0.1, betaROMaa = 0.1, betaTHAaa = 0.1, betaBRAaa = 0.1, betaCHIaa = 0.1, 
             
             betaUKha = 0.01, betaNEDha = 0.01, betaPOLha = 0.01, betaIREha = 0.01, betaBELha = 0.01, betaGERha = 0.01, betaITAha = 0.01, 
             betaFRAha = 0.01, betaSPAha = 0.01, betaROMha = 0.01, betaTHAha = 0.01, betaBRAha = 0.01, betaCHIha = 0.01, 
             
             tauUK = parmtau[i], tauNED = 0.05, tauPOL = 0.05, tauIRE = 0.05, tauBEL = 0.05, tauGER = 0.05, tauITA = 0.05, tauFRA = 0.05, 
             tauSPA = 0.05, tauROM = 0.05, tauTHA = 0.05, tauBRA = 0.05, tauCHI = 0.05,  
             theta = 0.5, phi = 0.05,
             psiUK = 0.593, psiNED = 0.1177, psiPOL = 0.0661, psiIRE = 0.0403, 
             psiBEL = 0.0215, psiGER = 0.0295, psiITA = 0.0065, psiFRA = 0.0180, psiSPA = 0.0051, 
             psiROM = 0.0051, psiTHA = 0.0726, psiBRA = 0.0214, psiCHI = 0.0032, 
             alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)
  out <- ode(y = init, func = amrhetcountry, times = times, parms = parms2)
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),2]) 
  temp[1,3] <- rounding(out[nrow(out),3]) 
  temp[1,4] <- rounding(out[nrow(out),4])
  temp[1,5] <- temp[1,3] + temp[1,4]
  temp[1,6] <- temp[1,4]/temp[1,5]
  temp[1,7] <- rounding(out[nrow(out),41]) 
  temp[1,8] <- rounding(out[nrow(out),42]) 
  temp[1,9] <- rounding(out[nrow(out),43])
  temp[1,10] <- temp[1,8] + temp[1,9]
  temp[1,11] <- temp[1,9]/temp[1,10]
  print(temp[1,1])
  outputlambda <- rbind.data.frame(outputlambda, temp)
}

colnames(outputlambda)[1:11] <- c("tau", "SuscAnimals","InfUKAnimals","ResUKInfAnimals","IUKCombA","IUKResRatA",
                                   "SuscHumans","InfUKHumans","ResUKInfHumans", "IUKCombH","IUKResRatH")

ggplot(data = outputlambda, aes(x = parmtau, y = IUKResRatA)) + geom_line()
ggplot(data = outputlambda, aes(x = parmtau, y = IUKCombA)) + geom_line()
ggplot(data = outputlambda, aes(x = parmtau, y = IUKResRatH)) + geom_line()
ggplot(data = outputlambda, aes(x = parmtau, y = IUKCombH)) + geom_line()

#### MISC Plotting ####

#High Res - Human Infect  
ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(meltedout$Compartment == "IDhs" | meltedout$Compartment == "IDhr" | meltedout$Compartment == "IEUhs" | 
                                   meltedout$Compartment == "IEUhr" | meltedout$Compartment == "InEUhs" | 
                                   meltedout$Compartment =="InEUhr"),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_y_continuous(limits = c(0,2), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Human Population 
ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(meltedout$Compartment == "IDhs" | meltedout$Compartment == "IDhr" | meltedout$Compartment == "IEUhs" | 
                                   meltedout$Compartment == "IEUhr" | meltedout$Compartment == "InEUhs" | 
                                   meltedout$Compartment =="InEUhr"),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_y_continuous(limits = c(0,2), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Animal Population 

ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment, linetype = Compartment)) + 
  geom_line(data=(y = meltedout[(meltedout$Compartment == "SDa" | meltedout$Compartment == "IDas" | meltedout$Compartment == "IDar" | 
                                   meltedout$Compartment == "SEUa" | meltedout$Compartment == "IEUar" | meltedout$Compartment =="SnEUa" | 
                                   meltedout$Compartment == "InEUas" | meltedout$Compartment == "InEUar"),]), size = 1.1) +
  scale_linetype_manual(values=c("dashed", "solid", "solid", "dashed", "solid", "dashed", "solid", "solid")) +
  labs(x ="Time (Days)", y = "Proportion of Animal Population") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#### Sensitivity Analysis ####


start_time <- Sys.time()

init <- c(SUKa = 0.98, IUKas = 0.01, IUKar = 0.01, SNEDa = 0.98, INEDas = 0.01, INEDar = 0.01, 
          SPOLa = 0.98, IPOLas = 0.01, IPOLar = 0.01, SIREa = 0.98, IIREas = 0.01, IIREar = 0.01, 
          SBELa = 0.98, IBELas = 0.01, IBELar = 0.01, 
          SGERa = 0.98, IGERas = 0.01, IGERar = 0.01, SITAa = 0.98, IITAas = 0.01, IITAar = 0.01,
          SFRAa = 0.98, IFRAas = 0.01, IFRAar = 0.01, SSPAa = 0.98, ISPAas = 0.01, ISPAar = 0.01, 
          SROMa = 0.98, IROMas = 0.01, IROMar = 0.01, STHAa = 0.98, ITHAas = 0.01, ITHAar = 0.01, 
          SBRAa = 0.98, IBRAas = 0.01, IBRAar = 0.01, SCHIa = 0.98, ICHIas = 0.01, ICHIar = 0.01,
          Sh = 1, IUKhs = 0, IUKhr = 0, INEDhs = 0, INEDhr = 0, IPOLhs = 0, IPOLhr = 0, IIREhs = 0, IIREhr = 0, 
          IBELhs = 0, IBELhr = 0, IGERhs = 0, IGERhr = 0, IITAhs = 0, IITAhr = 0, IFRAhs = 0, IFRAhr = 0, ISPAhs = 0, 
          ISPAhr = 0, IROMhs = 0, IROMhr = 0, ITHAhs = 0, ITHAhr = 0, IBRAhs = 0, IBRAhr = 0, ICHIhs = 0, ICHIhr = 0)

times1 <- seq(0,10000,by=1)

#Need to Specify Model Parameters

parms = c(betaUKaa = 0.1, betaNEDaa = 0.1, betaPOLaa = 0.1,  betaIREaa = 0.1,
          betaBELaa = 0.1, betaGERaa = 0.1, betaITAaa = 0.1, 
          betaFRAaa = 0.1, betaSPAaa = 0.1, betaROMaa = 0.1, 
          betaTHAaa = 0.1, betaBRAaa = 0.1, betaCHIaa = 0.1, 
          
          betaUKha = 0.01, betaNEDha = 0.01, betaPOLha = 0.01, betaIREha = 0.01, betaBELha = 0.01, betaGERha = 0.01, betaITAha = 0.01, 
          betaFRAha = 0.01, betaSPAha = 0.01, betaROMha = 0.01, betaTHAha = 0.01, betaBRAha = 0.01, betaCHIha = 0.01, 
          
          tauUK = 0.039, tauNED = 0.038, tauPOL = 0.067, tauIRE = 0.046, tauBEL = 0.089, tauGER = 0.101, tauITA = 0.157, tauFRA = 0.41, 
          tauSPA = 0.182, tauROM = 0.068, tauTHA = 0.071, tauBRA = 0.053, tauCHI = 0.319,  
          theta = 0.5, phi = 0.05,
          psiUK = 0.593, psiNED = 0.1177, psiPOL = 0.0661, psiIRE = 0.0403, 
          psiBEL = 0.0215, psiGER = 0.0295, psiITA = 0.0065, psiFRA = 0.0180, psiSPA = 0.0051, 
          psiROM = 0.0051, psiTHA = 0.0726, psiBRA = 0.0214, psiCHI = 0.0032, 
          alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)

out <- ode(y = init, func = amrhetcountry, times = times1, parms = parms)


parms = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
          0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  
          0.05, 0.005,
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
          0.5, 
          250^-1, 60^-1, 2400^-1, 288350^-1)


parms = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
          0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  
          2, 0.5,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
          1, 
          2.5^-1, 0.6^-1, 24^-1, 2883.5^-1)

#Create the parameter space for the model analysis/sensitivity analysis
parms = fast_parameters(minimum = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
                                    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  
                                    0.05, 0.005,
                                    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                                    0.5, 
                                    250^-1, 60^-1, 2400^-1, 288350^-1),
                        maximum = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  
                                    2, 0.5,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                    1, 
                                    2.5^-1, 0.6^-1, 24^-1, 2883.5^-1), 
                        factor=10, names = c("betaUKaa", "betaNEDaa", "betaPOLaa", "betaIREaa", "betaBELaa", "betaGERaa", "betaITAaa", "betaFRAaa", "betaSPAaa", 
                                             "betaROMaa", "betaTHAaa", "betaBRAaa", "betaCHIaa", 
                                              "betaUKha", "betaNEDha", "betaPOLha", "betaIREha", "betaBELha", "betaGERha", "betaITAha", "betaFRAha", 
                                              "betaSPAha", "betaROMha", "betaTHAha", "betaBRAha", "betaCHIha", 
                                              "tauUK", "tauNED", "tauPOL", "tauIRE", "tauBEL", "tauGER", "tauITA", "tauFRA", "tauSPA", "tauROM", "tauTHA", 
                                              "tauBRA", "tauCHI",  
                                              "theta", "phi",
                                               "psiUK", "psiNED", "psiPOL", "psiIRE", "psiBEL", "psiGER", "psiITA", "psiFRA", "psiSPA", "psiROM", "psiTHA", 
                                               "psiBRA", "psiCHI", 
                                               "alpha", 
                                               "rA", "rH", "eta", "mu"))

parms = fast_parameters(minimum = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
                                    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  
                                    0.05, 0.005),
                        maximum = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  
                                    2, 0.5), 
                        factor=41, names = c("betaUKaa", "betaNEDaa", "betaPOLaa", "betaIREaa", "betaBELaa", "betaGERaa", "betaITAaa", "betaFRAaa", "betaSPAaa", 
                                             "betaROMaa", "betaTHAaa", "betaBRAaa", "betaCHIaa", 
                                             "betaUKha", "betaNEDha", "betaPOLha", "betaIREha", "betaBELha", "betaGERha", "betaITAha", "betaFRAha", 
                                             "betaSPAha", "betaROMha", "betaTHAha", "betaBRAha", "betaCHIha", 
                                             "tauUK", "tauNED", "tauPOL", "tauIRE", "tauBEL", "tauGER", "tauITA", "tauFRA", "tauSPA", "tauROM", "tauTHA", 
                                             "tauBRA", "tauCHI",  
                                             "theta", "phi"))

#parms = c(betaDaa = 0.01, betaEUaa = 0.01, betanEUaa = 0.01, betaDha = 0.001, betaEUha = 0.001, betanEUha = 0.001,
#          tauD = 0, tauEU = 0, taunEU = 0,
#          theta = 0.05, phi = 0.005,
#          psiEU = 0.01, psinEU = 0.01, psiUK = 0.01,
#          alpha = 0, rA = 250^-1, rH = 60^-1, eta = 2400^-1, mu = 288350^-1)
#parms = c(betaDaa = 0.1, betaEUaa = 0.1, betanEUaa = 0.1, betaDha = 0.01, betaEUha = 0.01, betanEUha = 0.01,
#          tauD = 0.1, tauEU = 0.1, taunEU = 0.1, 
#          theta = 0.5, phi = 0.05,
#          psiEU = 1, psinEU = 1, psiUK = 1,
#          alpha = 1, rA = 2.5^-1, rH = 0.6^-1, eta = 24^-1, mu = 2883.5^-1)
#parms = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "tauD", "tauEU", "taunEU", "theta ", "phi" ,
#          "psiEU", "psinEU", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu")

init <- c(SDa = 0.98, IDas = 0.01, IDar = 0.01, 
          SEUa = 0.98, IEUas = 0.01, IEUar = 0.01, 
          SnEUa = 0.98, InEUas = 0.01, InEUar = 0.01,
          Sh = 1, IDhs = 0, IDhr = 0, IEUhs = 0, IEUhr = 0, InEUhs = 0, InEUhr = 0)
times <- seq(0,10000, by = 10) 

output <- data.frame()

for (i in 1:nrow(parms)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=8))
  parms1 = c(betaDaa = parms$betaDaa[i], betaEUaa = parms$betaEUaa[i], betanEUaa = parms$betanEUaa[i],
             betaDha = parms$betaDha[i], betaEUha = parms$betaEUha[i], betanEUha = parms$betanEUha[i], 
             tauD = parms$tauD[i], tauEU = parms$tauEU[i], taunEU = parms$taunEU[i], 
             theta = parms$theta[i], phi = parms$phi[i],
             psiEU = (parms$psiEU[i]/(parms$psiUK[i]+parms$psiEU[i]+parms$psinEU[i])), 
             psinEU = (parms$psinEU[i]/(parms$psiUK[i]+parms$psiEU[i]+parms$psinEU[i])), 
             psiUK = (parms$psiUK[i]/(parms$psiUK[i]+parms$psiEU[i]+parms$psinEU[i])),
             alpha = parms$alpha[i], rA = parms$rA[i], rH = parms$rH[i], eta = parms$eta[i], mu = parms$mu[i])
  out <- ode(y = init, func = amrhet, times = times, parms = parms1)
  temp[1,1] <- (rounding(out[nrow(out),12]) + rounding(out[nrow(out),13])) + (rounding(out[nrow(out),14]) + rounding(out[nrow(out),15])) +
    (rounding(out[nrow(out),16]) + rounding(out[nrow(out),17]))
  temp[1,2] <- (rounding(out[nrow(out),13]) + rounding(out[nrow(out),15]) + rounding(out[nrow(out),17])) /  temp[1,1]
  temp[1,3] <- (rounding(out[nrow(out),12]) + rounding(out[nrow(out),14]) + rounding(out[nrow(out),16]))
  temp[1,4] <- (rounding(out[nrow(out),12]) + rounding(out[nrow(out),13]))
  temp[1,5] <- rounding(out[nrow(out),13]) / temp[1,4]
  temp[1,6] <- (rounding(out[nrow(out),14]) + rounding(out[nrow(out),15])) + (rounding(out[nrow(out),16]) + rounding(out[nrow(out),17]))
  temp[1,7] <- (rounding(out[nrow(out),15]) + rounding(out[nrow(out),17])) / temp[1,6]
  temp[1,8] <- as.numeric(parms1[14])
  print(temp[1,8])
  output <- rbind.data.frame(output, temp)
}

colnames(output)[1:8] <- c("OverallInf","OverallResProp","OverallSensProp","DomOverInf","DomResProp", 
                           "ImpOverInf","ImpResProp","UKUse")

end_time <- Sys.time()
end_time - start_time

#### Sensitivity Analysis Plotting and Analysis #### 
sensit1 <- output$OverallInf #Creating Variable for the output variable of interest
sens1 <- sensitivity(x=sensit1, numberf=19, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                                                   "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                                   "alpha", "rA" , "rH" , "eta" , "mu"))

df.equilibrium <- NULL
df.equilibrium <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                             "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                             "alpha", "rA" , "rH" , "eta" , "mu"), value=sens1)

ggplot(df.equilibrium, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")

sensit2 <- output$OverallResProp #Creating Variable for the output variable of interest
sensit2[is.nan(sensit2)] <- 0

sens2 <- sensitivity(x=sensit2, numberf=19, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                                                   "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                                   "alpha", "rA" , "rH" , "eta" , "mu"))

df.equilibrium1 <- NULL
df.equilibrium1 <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                              "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                              "alpha", "rA" , "rH" , "eta" , "mu"), value=sens2)

ggplot(df.equilibrium1, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")


sensit3 <- output$ImpOverInf #Creating Variable for the output variable of interest
sensit4 <- output$ImpResProp #Creating Variable for the output variable of interest
sensit4[is.nan(sensit4)] <- 0

sens4 <- sensitivity(x=sensit4, numberf=19, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                                                   "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                                   "alpha", "rA" , "rH" , "eta" , "mu"))

df.equilibrium3 <- NULL
df.equilibrium3 <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                              "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                              "alpha", "rA" , "rH" , "eta" , "mu"), value=sens4)

ggplot(df.equilibrium3, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")
