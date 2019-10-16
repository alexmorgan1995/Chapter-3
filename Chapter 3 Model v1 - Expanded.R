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
      rh*(IUKhs + IUKhr + INEDhs + INEDhr + IPOLhs + IPOLhr + IIREhs + IIREhr + IBELhs + IBELhr +
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
    
    return(list(c(dSUKa, dIUKas, dIUKar,dSNEDa, dINEDas, dINEDar, dSPOLa, dIPOLas, dIPOLar,
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
          SPOLa = 0.98, IPOLas = 0.01, IPOLar = 0.01, SBELa = 0.98, IBELas = 0.01, IBELar = 0.01, 
          SGERa = 0.98, IGERas = 0.01, IGERar = 0.01, SITAa = 0.98, IITAas = 0.01, IITAar = 0.01,
          SFRAa = 0.98, IFRAas = 0.01, IFRAar = 0.01, SSPAa = 0.98, ISPAas = 0.01, ISPAar = 0.01, 
          SROMa = 0.98, IROMas = 0.01, IROMar = 0.01, STHAa = 0.98, ITHAas = 0.01, ITHAar = 0.01, 
          SBRAa = 0.98, IBRAas = 0.01, IBRAar = 0.01, SCHIa = 0.98, ICHIas = 0.01, ICHIar = 0.01,
          Sh = 1, IUKhs = 0, IUKhr = 0, INEDhs = 0, INEDhr = 0, IPOLhs = 0, IPOLhr = 0, IIREhs = 0, IIREhr = 0, 
          IBELhs = 0, IBELhr = 0, IGERhs = 0, IGERhr = 0, IITAhs = 0, IITAhr = 0, IFRAhs = 0, IFRAhr = 0, ISPAhs = 0, 
          ISPAhr = 0, IROMhs = 0, IROMhr = 0, ITHAhs = 0, ITHAhr = 0, IBRAhs = 0, IBRAhr = 0, ICHIhs = 0, ICHIhr = 0)

times1 <- seq(0,1000,by=1)

#Need to Specify Model Parameters

parms = c(betaUKaa = 0.1, betaNEDaa = 0.1, betaPOLaa = 0.1, betaBELaa = 0.1, betaGERaa = 0.1, betaITAaa = 0.1, 
          betaFRAaa = 0.1, betaSPAaa = 0.1, betaROMaa = 0.1, betaTHAaa = 0.1, betaBRAaa = 0.1, betaCHIaa = 0.1, 
          
          betaUKha = 0.1, betaNEDha = 0.1, betaPOLha = 0.1, betaBELha = 0.1, betaGERha = 0.1, betaITAha = 0.1, 
          betaFRAha = 0.1, betaSPAha = 0.1, betaROMha = 0.1, betaTHAha = 0.1, betaBRAha = 0.1, betaCHIha = 0.1, 
          
          tauUK = 0.05, tauNED = 0.08, tauPOL = 0.1, tauBEL = 0.05, tauGER = 0.08, tauITA = 0.1, tauFRA = 0.1, 
          tauSPA = 0.05, tauROM = 0.08, tauTHA = 0.1, tauBRA = 0.05, tauCHI = 0.08,  
          theta = 0.5, phi = 0.05,
          psiUK = 0.2, psiNED = 0.1, psiPOL = 0.7, psiBEL = 0.2, psiGER = 0.1, psiITA = 0.7, psiFRA = 0.7, psiSPA = 0.7, 
          psiROM = 0.2, psiTHA = 0.1, psiBRA = 0.7, psiCHI = 0.2, 
          alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)

out <- ode(y = init, func = amrhetcountry, times = times1, parms = parms)

outdata <- data.frame(out)
outdata$DIComb <- outdata$IDhs + outdata$IDhr 
outdata$EUIComb <- outdata$IEUhs + outdata$IEUhr 
outdata$nEUIComb <- outdata$InEUhs + outdata$InEUhr

#Creating L I Q U I D code - i.e - Melting for GGPLOT
outdata1 <- outdata
outdata1[11:20] <- outdata[11:20]*100000 # So only human compartments are scaled to per 100,000 population 
meltedout <- melt(outdata1, id = "time", variable.name = "Compartment", value.name = "Value")

#Manipulating the Data for a stacked bar plot - to show the relative proportions from each country
#The bar chart plot will need to be when the model is at equilbirum - length(outdata) or something similar 
equidf <- data.frame(matrix(NA, ncol = 4, nrow = 3))
colnames(equidf) <- c("Category","Domestic","European", "nonEuropean")
equidf[1] <- c("ICombH", "Sensitive", "Resistant")
equidf[2] <- c(outdata$DIComb[length(outdata$DIComb)], outdata$IDhs[length(outdata$IDhs)], 
               outdata$IDhr[length(outdata$IDhr)]) 
equidf[3] <- c(outdata$EUIComb[length(outdata$DIComb)], outdata$IEUhs[length(outdata$IEUhs)],
               outdata$IEUhr[length(outdata$IEUhr)]) 
equidf[4] <- c(outdata$nEUIComb[length(outdata$nEUIComb)], outdata$InEUhs[length(outdata$InEUhs)], 
               outdata$InEUhr[length(outdata$InEUhr)])  

###Plotting Output for Basic Script####

#Bar Chart Plot - to view the relative infecteds 
plot_ly(equidf, x = ~Category, y = ~Domestic, type = "bar", name = "Domestic") %>%
  add_trace(y = ~European, name = "European") %>%
  add_trace(y = ~nonEuropean, name = "nonEuropean") %>%
  layout(yaxis = list(title = "Proportion of Infecteds"), barmode = "stack")

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
