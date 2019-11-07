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
    
    dSDENa = - betaDENaa*(IDENas + IDENar*alpha)*SDENa + rA*(IDENas+IDENar) + tauDEN*IDENas + eta - eta*SDENa
    dIDENas = betaDENaa*IDENas*SDENa - rA*IDENas - tauDEN*IDENas - tauDEN*theta*IDENas + phi*IDENar - eta*IDENas
    dIDENar = betaDENaa*IDENar*SDENa*alpha + tauDEN*theta*IDENas - phi*IDENar - rA*IDENar - eta*IDENar 
    
    dSGERa = - betaGERaa*(IGERas + IGERar*alpha)*SGERa + rA*(IGERas+IGERar) + tauGER*IGERas + eta - eta*SGERa
    dIGERas = betaGERaa*IGERas*SGERa - rA*IGERas - tauGER*IGERas - tauGER*theta*IGERas + phi*IGERar - eta*IGERas
    dIGERar = betaGERaa*IGERar*SGERa*alpha + tauGER*theta*IGERas - phi*IGERar - rA*IGERar - eta*IGERar 
    
    dSNEDa = - betaNEDaa*(INEDas + INEDar*alpha)*SNEDa + rA*(INEDas+INEDar) + tauNED*INEDas + eta - eta*SNEDa
    dINEDas = betaNEDaa*INEDas*SNEDa - rA*INEDas - tauNED*INEDas - tauNED*theta*INEDas + phi*INEDar - eta*INEDas
    dINEDar = betaNEDaa*INEDar*SNEDa*alpha + tauNED*theta*INEDas - phi*INEDar - rA*INEDar - eta*INEDar 
    
    dSIREa = - betaIREaa*(IIREas + IIREar*alpha)*SIREa + rA*(IIREas+IIREar) + tauIRE*IIREas + eta - eta*SIREa
    dIIREas = betaIREaa*IIREas*SIREa - rA*IIREas - tauIRE*IIREas - tauIRE*theta*IIREas + phi*IIREar - eta*IIREas
    dIIREar = betaIREaa*IIREar*SIREa*alpha + tauIRE*theta*IIREas - phi*IIREar - rA*IIREar - eta*IIREar 
    
    dSPOLa = - betaPOLaa*(IPOLas + IPOLar*alpha)*SPOLa + rA*(IPOLas+IPOLar) + tauPOL*IPOLas + eta - eta*SPOLa
    dIPOLas = betaPOLaa*IPOLas*SPOLa - rA*IPOLas - tauPOL*IPOLas - tauPOL*theta*IPOLas + phi*IPOLar - eta*IPOLas
    dIPOLar = betaPOLaa*IPOLar*SPOLa*alpha + tauPOL*theta*IPOLas - phi*IPOLar - rA*IPOLar - eta*IPOLar 
    
    dSSPAa = - betaSPAaa*(ISPAas + ISPAar*alpha)*SSPAa + rA*(ISPAas+ISPAar) + tauSPA*ISPAas + eta - eta*SSPAa
    dISPAas = betaSPAaa*ISPAas*SSPAa - rA*ISPAas - tauSPA*ISPAas - tauSPA*theta*ISPAas + phi*ISPAar - eta*ISPAas
    dISPAar = betaSPAaa*ISPAar*SSPAa*alpha + tauSPA*theta*ISPAas - phi*ISPAar - rA*ISPAar - eta*ISPAar 
    
    dSBELa = - betaBELaa*(IBELas + IBELar*alpha)*SBELa + rA*(IBELas+IBELar) + tauBEL*IBELas + eta - eta*SBELa
    dIBELas = betaBELaa*IBELas*SBELa - rA*IBELas - tauBEL*IBELas - tauBEL*theta*IBELas + phi*IBELar - eta*IBELas
    dIBELar = betaBELaa*IBELar*SBELa*alpha + tauBEL*theta*IBELas - phi*IBELar - rA*IBELar - eta*IBELar 
    
    dSh = - betaUKha*(eta*IUKas)*Sh*psiUK - betaUKha*(eta*IUKar)*Sh*psiUK*alpha - 
      betaDENha*(eta*IDENas)*Sh*psiDEN - betaDENha*(eta*IDENar)*Sh*psiDEN*alpha - 
      betaGERha*(eta*IGERas)*Sh*psiGER - betaGERha*(eta*IGERar)*Sh*psiGER*alpha -
      betaNEDha*(eta*INEDas)*Sh*psiNED - betaNEDha*(eta*INEDar)*Sh*psiNED*alpha -
      betaIREha*(eta*IIREas)*Sh*psiIRE - betaIREha*(eta*IIREar)*Sh*psiIRE*alpha - 
      betaPOLha*(eta*IPOLas)*Sh*psiPOL - betaPOLha*(eta*IPOLar)*Sh*psiPOL*alpha -
      betaSPAha*(eta*ISPAas)*Sh*psiSPA - betaSPAha*(eta*ISPAar)*Sh*psiSPA*alpha - 
      betaBELha*(eta*IBELas)*Sh*psiBEL - betaBELha*(eta*IBELar)*Sh*psiBEL*alpha - 
     
      betahh*(IUKhs+IDENhs+IGERhs+INEDhs+IIREhs+IPOLhs+ISPAhs+IBELhs +
                (IUKhr+IDENhr+IGERhr+INEDhr+IIREhr+IPOLhr+ISPAhr+IBELhr)*alpha)*Sh + 
      rH*(IUKhs + IUKhr + IDENhs + IDENhr + IGERhs + IGERhr + INEDhs + INEDhr + IIREhs + IIREhr +
            IPOLhs + IPOLhr + ISPAhs + ISPAhr + IBELhs + IBELhr) + mu - mu*Sh
      
    dIUKhs = betaUKha*(eta*IUKas)*Sh*psiUK + betahh*IUKhs*Sh - rH*IUKhs - mu*IUKhs 
    dIUKhr = betaUKha*(eta*IUKar)*Sh*psiUK*alpha + betahh*IUKhr*Sh*alpha - rH*IUKhr - mu*IUKhr 
    
    dIDENhs = betaDENha*(eta*IDENas)*Sh*psiDEN + betahh*IDENhs*Sh - rH*IDENhs - rH*IDENhs - mu*IDENhs 
    dIDENhr = betaDENha*(eta*IDENar)*Sh*psiDEN*alpha + betahh*IDENhr*Sh*alpha - rH*IDENhr - mu*IDENhr 
    
    dIGERhs = betaGERha*(eta*IGERas)*Sh*psiGER + betahh*IGERhs*Sh - rH*IGERhs - mu*IGERhs 
    dIGERhr = betaGERha*(eta*IGERar)*Sh*psiGER*alpha + betahh*IGERhr*Sh*alpha - rH*IGERhr - mu*IGERhr 
    
    dINEDhs = betaNEDha*(eta*INEDas)*Sh*psiNED + betahh*INEDhs*Sh - rH*INEDhs - rH*INEDhs - mu*INEDhs 
    dINEDhr = betaNEDha*(eta*INEDar)*Sh*psiNED*alpha + betahh*INEDhr*Sh*alpha - rH*INEDhr - mu*INEDhr 
    
    dIIREhs = betaIREha*(eta*IIREas)*Sh*psiIRE + betahh*IIREhs*Sh - rH*IIREhs - mu*IIREhs 
    dIIREhr = betaIREha*(eta*IIREar)*Sh*psiIRE*alpha + betahh*IIREhr*Sh*alpha - rH*IIREhr - mu*IIREhr 
    
    dIPOLhs = betaPOLha*(eta*IPOLas)*Sh*psiPOL + betahh*IPOLhs*Sh - rH*IPOLhs - mu*IPOLhs 
    dIPOLhr = betaPOLha*(eta*IPOLar)*Sh*psiPOL*alpha + betahh*IPOLhr*Sh*alpha - rH*IPOLhr - mu*IPOLhr 
    
    dISPAhs = betaSPAha*(eta*ISPAas)*Sh*psiSPA + betahh*ISPAhs*Sh - rH*ISPAhs - mu*ISPAhs 
    dISPAhr = betaSPAha*(eta*ISPAar)*Sh*psiSPA*alpha + betahh*ISPAhr*Sh*alpha - rH*ISPAhr - mu*ISPAhr 
    
    dIBELhs = betaBELha*(eta*IBELas)*Sh*psiBEL + betahh*IBELhs*Sh - rH*IBELhs - mu*IBELhs 
    dIBELhr = betaBELha*(eta*IBELar)*Sh*psiBEL*alpha + betahh*IBELhr*Sh*alpha - rH*IBELhr - mu*IBELhr 
    
    return(list(c(dSUKa, dIUKas, dIUKar, dSDENa, dIDENas, dIDENar, dSGERa, dIGERas, dIGERar,
                  dSNEDa, dINEDas, dINEDar, dSIREa, dIIREas, dIIREar, dSPOLa, dIPOLas, dIPOLar, 
                  dSSPAa, dISPAas, dISPAar, dSBELa, dIBELas, dIBELar,
                  dSh, dIUKhs, dIUKhr, dIDENhs, dIDENhr, dIGERhs, dIGERhr, dINEDhs, dINEDhr, dIIREhs, dIIREhr, dIPOLhs, dIPOLhr,
                  dISPAhs, dISPAhr, dIBELhs, dIBELhr)))
  })
}

rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Parameters ####

init <- c(SUKa = 0.98, IUKas = 0.01, IUKar = 0.01, SDENa = 0.98, IDENas = 0.01, IDENar = 0.01,
          SGERa = 0.98, IGERas = 0.01, IGERar = 0.01, SNEDa = 0.98, INEDas = 0.01, INEDar = 0.01,
          SIREa = 0.98, IIREas = 0.01, IIREar = 0.01, SPOLa = 0.98, IPOLas = 0.01, IPOLar = 0.01,
          SSPAa = 0.98, ISPAas = 0.01, ISPAar = 0.01, SBELa = 0.98, IBELas = 0.01, IBELar = 0.01, 
          Sh = 1, IUKhs = 0, IUKhr = 0, IDENhs = 0, IDENhr = 0, IGERhs = 0, IGERhr = 0, INEDhs = 0, INEDhr = 0, 
          IIREhs = 0, IIREhr = 0, IPOLhs = 0, IPOLhr = 0, ISPAhs = 0, ISPAhr = 0, IBELhs = 0, IBELhr = 0)

times1 <- seq(0,150000,by=100)

#Need to Specify Model Parameters

parms2 <- parms1

parms = c(betaUKaa = 0.1, betaDENaa = 0.1, betaGERaa = 0.1, betaNEDaa = 0.1, betaIREaa = 0.1, betaPOLaa = 0.1,
          betaSPAaa = 0.1, betaBELaa = 0.1,

          betahh = 0.01,
          
          betaUKha = 0.01, betaDENha = 0.01, betaGERha = 0.01, betaNEDha = 0.01, betaIREha = 0.01, betaPOLha = 0.01,
          betaSPAha = 0.01, betaBELha = 0.01,  

          tauUK = 0.0325, tauDEN = 0.0394, tauGER = 0.089, tauNED = 0.0563, tauIRE = 0.0466, tauPOL = 0.1652, tauSPA = 0.2303,
          tauBEL = 0.1313,

          theta = 0.5, phi = 0.05,
          
          psiUK = 0.6, psiDEN = 0.1445, psiGER = 0.1277, psiNED = 0.1056, psiIRE = 0.0798, psiPOL = 0.0437,
          psiSPA = 0.034, psiBEL = 0.0245,
          
          alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)

out <- ode(y = init, func = amrhetcountry, times = times1, parms = parms2)

sum(as.numeric(out[(length(out[,1])),(2:4)]))
sum(as.numeric(out[(length(out[,1])),(26:42)]))

outdata <- data.frame(out)
outdata$UKIComb <- outdata$IUKhs + outdata$IUKhr 
outdata$DENIComb <- outdata$IDENhs + outdata$IDENhr 
outdata$GERIComb <- outdata$IGERhs + outdata$IGERhr
outdata$NEDIComb <- outdata$INEDhs + outdata$INEDhr 
outdata$IREIComb <- outdata$IIREhs + outdata$IIREhr 
outdata$POLIComb <- outdata$IPOLhs + outdata$IPOLhr
outdata$SPAIComb <- outdata$ISPAhs + outdata$ISPAhr
outdata$BELIComb <- outdata$IBELhs + outdata$IBELhr 

#Creating L I Q U I D code - i.e - Melting for GGPLOT
outdata1 <- outdata

outdata1[26:50] <- outdata1[26:50]*100000 # So only human compartments are scaled to per 100,000 population 
meltedout <- melt(outdata1, id = "time", variable.name = "Compartment", value.name = "Value")

#Manipulating the Data for a stacked bar plot - to show the relative proportions from each country
#The bar chart plot will need to be when the model is at equilbirum - length(outdata) or something similar 

equidf <- data.frame(matrix(NA, ncol = 3, nrow = 3))
colnames(equidf) <- c("Category","Domestic","Imported")
equidf[1] <- c("ICombH", "Sensitive", "Resistant")
equidf[2] <- c(outdata$UKIComb[length(outdata$UKIComb)], outdata$IUKhs[length(outdata$IUKhs)], outdata$IUKhr[length(outdata$IUKhr)]) 

equidf[3] <- c(outdata$DENIComb[length(outdata$DENIComb)] + outdata$GERIComb[length(outdata$GERIComb)] + 
                 outdata$NEDIComb[length(outdata$NEDIComb)] + outdata$IREIComb[length(outdata$IREIComb)] + 
                 outdata$POLIComb[length(outdata$POLIComb)] + outdata$SPAIComb[length(outdata$SPAIComb)] +
                 outdata$BELIComb[length(outdata$BELIComb)], 
               
               outdata$IDENhs[length(outdata$IDENhs)] + outdata$IGERhs[length(outdata$IGERhs)] + 
                 outdata$INEDhs[length(outdata$INEDhs)] + outdata$IIREhs[length(outdata$IIREhs)] +
                 outdata$IPOLhs[length(outdata$IPOLhs)] + outdata$ISPAhs[length(outdata$ISPAhs)] + 
                 outdata$IBELhs[length(outdata$IBELhs)],
               
               outdata$IDENhr[length(outdata$IDENhr)] + outdata$IGERhr[length(outdata$IGERhr)] + 
                 outdata$INEDhr[length(outdata$INEDhr)] + outdata$IIREhr[length(outdata$IIREhr)] + 
                 outdata$IPOLhr[length(outdata$IPOLhr)] + outdata$ISPAhr[length(outdata$ISPAhr)] +
                 outdata$IBELhr[length(outdata$IBELhr)])
equidf[4,] <- c("ResProp", signif((as.numeric(equidf[3,2] / equidf[1,2])), digits = 2), 
                signif(as.numeric(equidf[3,3] / equidf[1,3]), digits = 2))
equidf[,4] <- as.numeric(equidf[,2]) / (as.numeric(equidf[,2]) +  as.numeric(equidf[,3])) 

###Plotting Output for Basic Script####

ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(as.numeric(which((meltedout$Compartment == "UKIComb" & meltedout$time == 0))):
                                   length(meltedout$time)),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Bar Chart Plot - to view the relative infecteds 
plot_ly(equidf, x = ~Category, y = ~(Domestic*100000), type = "bar", name = "Domestic") %>%
  add_trace(y = ~(Imported*100000), name = "Imported", text = y, textposition = "auto")  %>%
  layout(yaxis = list(title = "Proportion of Infecteds (per 100,000 Population", range = c(0, 10)), barmode = "stack",
         annotations = list(x = ~Category, y = ~(Domestic*100000), text = ~(Domestic*100000), yanchor = "top", 
                            showarrow = FALSE,  xshift = 3))

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

#### Sensitivity Analysis ####

start_time <- Sys.time()

init <- c(SUKa = 0.98, IUKas = 0.01, IUKar = 0.01, 
          SDENa = 0.98, IDENas = 0.01, IDENar = 0.01,
          SGERa = 0.98, IGERas = 0.01, IGERar = 0.01,
          SNEDa = 0.98, INEDas = 0.01, INEDar = 0.01, 
          SIREa = 0.98, IIREas = 0.01, IIREar = 0.01,
          SPOLa = 0.98, IPOLas = 0.01, IPOLar = 0.01, 
          SSPAa = 0.98, ISPAas = 0.01, ISPAar = 0.01,
          SBELa = 0.98, IBELas = 0.01, IBELar = 0.01, 
          
          Sh = 1, IUKhs = 0, IUKhr = 0, 
          IDENhs = 0, IDENhr = 0, 
          IGERhs = 0, IGERhr = 0, 
          INEDhs = 0, INEDhr = 0, 
          IIREhs = 0, IIREhr = 0, 
          IPOLhs = 0, IPOLhr = 0, 
          ISPAhs = 0, ISPAhr = 0, 
          IBELhs = 0, IBELhr = 0)

times1 <- seq(0,100000,by=100)

#Create the parameter space for the model analysis/sensitivity analysis
parms = fast_parameters(minimum = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
                                    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
                                    0,
                                    0, 0, 0, 0, 0, 0, 0, 0, 
                                    0.05, 0.005,
                                    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
                                    0.5, 
                                    250^-1, 60^-1, 2400^-1, 288350^-1),
                        maximum = c(1, 1, 1, 1, 1, 1, 1, 1,  
                                    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,  
                                    0.1,
                                    1, 1, 1, 1, 1, 1, 1, 1, 
                                    2, 0.5,
                                    1, 1, 1, 1, 1, 1, 1, 1,  
                                    1, 
                                    2.5^-1, 0.6^-1, 24^-1, 2883.5^-1), 
                        factor=40, names = c("betaUKaa", "betaDENaa", "betaGERaa", "betaNEDaa", "betaIREaa", "betaPOLaa",
                                             "betaSPAaa", "betaBELaa", 
                                             
                                             "betaUKha", "betaDENha", "betaGERha", "betaNEDha", "betaIREha", "betaPOLha", 
                                             "betaSPAha", "betaBELha",
                                             "betahh", 
                                             
                                             "tauUK", "tauDEN", "tauGER", "tauNED", "tauIRE", "tauPOL", "tauSPA", "tauBEL", 
                                             "theta", "phi",
                                             "psiUK", "psiDEN", "psiGER", "psiNED", "psiIRE", "psiPOL", "psiSPA", "psiBEL",
                                               "alpha", 
                                               "rA", "rH", "eta", "mu"))

output <- data.frame()

for (i in 1:nrow(parms)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=8))
  parms1 = c(betaUKaa = parms$betaUKaa[i], betaDENaa = parms$betaDENaa[i], betaGERaa = parms$betaGERaa[i], 
             betaNEDaa = parms$betaNEDaa[i], betaIREaa = parms$betaIREaa[i], betaPOLaa = parms$betaPOLaa[i],  
             betaSPAaa = parms$betaSPAaa[i], betaBELaa = parms$betaBELaa[i], 
             
             betaUKha = parms$betaUKha[i], betaDENha = parms$betaDENha[i], betaGERha = parms$betaGERha[i],
             betaNEDha = parms$betaNEDha[i], betaIREha = parms$betaIREha[i], betaPOLha = parms$betaPOLha[i], 
             betaSPAha = parms$betaSPAha[i], betaBELha = parms$betaBELha[i],
             betahh = parms$betahh[i],
             tauUK = parms$tauUK[i], tauDEN = parms$tauDEN[i], tauGER = parms$tauGER[i], tauNED = parms$tauNED[i], 
             tauIRE = parms$tauIRE[i], tauPOL = parms$tauPOL[i], tauSPA = parms$tauSPA[i], tauBEL = parms$tauBEL[i], 
             theta = parms$theta[i], phi = parms$phi[i],
             
             psiUK = parms$psiUK[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                         parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]), 
             psiDEN = parms$psiDEN[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                           parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]), 
             psiGER = parms$psiGER[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                           parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]), 
             psiNED = parms$psiNED[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                           parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]), 
             psiIRE = parms$psiIRE[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                           parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]), 
             psiPOL = parms$psiPOL[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                           parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]), 
             psiSPA = parms$psiSPA[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                           parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]), 
             psiBEL = parms$psiBEL[i] / (parms$psiUK[i] + parms$psiDEN[i] + parms$psiGER[i] + parms$psiNED[i] + 
                                           parms$psiIRE[i] + parms$psiPOL[i] + parms$psiSPA[i] + parms$psiBEL[i]),
             alpha = parms$alpha[i], 
             rA = parms$rA[i], rH = parms$rH[i], eta = parms$eta[i], mu = parms$mu[i])
    
  out <- ode(y = init, func = amrhetcountry, times = times1, parms = parms1)
  temp[1,1] <- rounding(sum(as.numeric(out[(length(out[,1])),(27:42)])))
  temp[1,2] <-  rounding(sum(as.numeric(out[(length(out[,1])), c(28, 30, 32, 34, 36 ,38, 42)]))) /  temp[1,1]
  temp[1,3] <- rounding(sum(as.numeric(out[(length(out[,1])), c(27, 29, 31, 33, 35 ,37, 39, 41)])))
  temp[1,4] <- rounding(sum(as.numeric(out[(length(out[,1])),c(27,28)])))
  temp[1,5] <- rounding(out[nrow(out),28]) / temp[1,4]
  temp[1,6] <- rounding(sum(as.numeric(out[(length(out[,1])),(29:42)])))
  temp[1,7] <- rounding(sum(as.numeric(out[(length(out[,1])), c(30, 32, 34, 36, 38, 40 ,42)]))) / temp[1,6]
  temp[1,8] <- as.numeric(parms1[28])
  print(temp[1,8])
  output <- rbind.data.frame(output, temp)
}

colnames(output)[1:8] <- c("OverallInf","OverallResProp","OverallSensProp","DomOverInf","DomResProp", 
                           "ImpOverInf","ImpResProp","UKUse")
end_time <- Sys.time()
end_time - start_time

#### Sensitivity Analysis Plotting and Analysis #### 
sensit1 <- output$OverallInf #Creating Variable for the output variable of interest
sens1 <- sensitivity(x=sensit1, numberf=40, make.plot=T, names = c("betaUKaa", "betaDENaa", "betaGERaa", "betaNEDaa", "betaIREaa", "betaPOLaa",
                                                                   "betaSPAaa", "betaBELaa", 
                                                                   
                                                                   "betaUKha", "betaDENha", "betaGERha", "betaNEDha", "betaIREha", "betaPOLha", 
                                                                   "betaSPAha", "betaBELha",
                                                                   "betahh", 
                                                                   
                                                                   "tauUK", "tauDEN", "tauGER", "tauNED", "tauIRE", "tauPOL", "tauSPA", "tauBEL", 
                                                                   "theta", "phi",
                                                                   "psiUK", "psiDEN", "psiGER", "psiNED", "psiIRE", "psiPOL", "psiSPA", "psiBEL",
                                                                   "alpha", 
                                                                   "rA", "rH", "eta", "mu"))

df.equilibrium <- NULL
df.equilibrium <- data.frame(parameter=rbind("betaUKaa", "betaDENaa", "betaGERaa", "betaNEDaa", "betaIREaa", "betaPOLaa",
                                             "betaSPAaa", "betaBELaa", 
                                             
                                             "betaUKha", "betaDENha", "betaGERha", "betaNEDha", "betaIREha", "betaPOLha", 
                                             "betaSPAha", "betaBELha",
                                             "betahh", 
                                             
                                             "tauUK", "tauDEN", "tauGER", "tauNED", "tauIRE", "tauPOL", "tauSPA", "tauBEL", 
                                             "theta", "phi",
                                             "psiUK", "psiDEN", "psiGER", "psiNED", "psiIRE", "psiPOL", "psiSPA", "psiBEL",
                                             "alpha", 
                                             "rA", "rH", "eta", "mu"), value=sens1)

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
