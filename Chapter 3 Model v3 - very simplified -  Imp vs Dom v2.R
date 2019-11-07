rm(list=ls())

library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("tidyr")

#### Model For Domestic, EU and non-EU Livestock-to-human AMR transmission ####

amrhet <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSDa = - betaDaa*(IDas + IDar*alpha)*SDa + rA*(IDas+IDar) + tauD*IDas + eta - eta*SDa
    dIDas = betaDaa*IDas*SDa - rA*IDas - tauD*IDas - tauD*theta*IDas + phi*IDar - eta*IDas
    dIDar = betaDaa*IDar*SDa*alpha + tauD*theta*IDas - phi*IDar - rA*IDar - eta*IDar 
    
    dSIMPa = - betaIMPaa*(IIMPas + IIMPar*alpha)*SIMPa + rA*(IIMPas+IIMPar) + tauIMP*IIMPas + eta - eta*SIMPa
    dIIMPas = betaIMPaa*IIMPas*SIMPa - rA*IIMPas - tauIMP*IIMPas - tauIMP*theta*IIMPas + phi*IIMPar - eta*IIMPas
    dIIMPar = betaIMPaa*IIMPar*SIMPa*alpha + tauIMP*theta*IIMPas - phi*IIMPar - rA*IIMPar - eta*IIMPar 
    
    dSh = - betaDha*(eta*IDas)*Sh*psiUK - betaDha*(eta*IDar)*Sh*alpha*psiUK - 
      betaIMPha*(eta*IIMPas)*Sh*psiIMP - betaIMPha*(eta*IIMPar)*Sh*psiIMP*alpha - 
      betahh*(IDhs+IIMPhs+(IDhr+IIMPhr)*alpha)*Sh +
      rH*(IDhs + IDhr + IIMPhs + IIMPhr) + mu - mu*Sh
    
    dIDhs = betahh*IDhs*Sh + betaDha*(eta*IDas)*Sh*psiUK - rH*IDhs - mu*IDhs 
    dIDhr = betahh*IDhr*Sh*alpha + betaDha*(eta*IDar)*Sh*psiUK*alpha - rH*IDhr - mu*IDhr
    
    dIIMPhs = betahh*IIMPhs*Sh + betaIMPha*(eta*IIMPas)*Sh*psiIMP - rH*IIMPhs - mu*IIMPhs
    dIIMPhr = betahh*IIMPhr*Sh*alpha + betaIMPha*(eta*IIMPar)*Sh*psiIMP*alpha - rH*IIMPhr - mu*IIMPhr
    
    return(list(c(dSDa, dIDas, dIDar, dSIMPa, dIIMPas, dIIMPar, dSh, dIDhs, dIDhr, dIIMPhs, dIIMPhr)))
  })
}

rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Parameters ####

init <- c(SDa = 0.98, IDas = 0.01, IDar = 0.01, 
          SIMPa = 0.98, IIMPas = 0.01, IIMPar = 0.01, 
          Sh = 1, IDhs = 0, IDhr = 0, IIMPhs = 0, IIMPhr = 0)

times1 <- seq(0,5000,by=1)

#Need to Specify Model Parameters

parms = c(betaDaa = 0.1, betaIMPaa = 0.1, betaDha = 0.01, betaIMPha = 0.01, 
          betahh = 0.01,
          tauD = 0.039, tauIMP = 0.0876,
          theta = 0.5, phi = 0.05,
          psiIMP = 0.3099, psiUK = 0.593,
          alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)
          
out <- ode(y = init, func = amrhet, times = times1, parms = parms)

outdata <- data.frame(out)
outdata$DIComb <- outdata$IDhs + outdata$IDhr 
outdata$IMPIComb <- outdata$IIMPhs + outdata$IIMPhr 

#Creating L I Q U I D code - i.e - Melting for GGPLOT
outdata1 <- outdata; outdata1[8:12] <- outdata[8:12]*100000 # So only human compartments are scaled to per 100,000 population 
meltedout <- melt(outdata1, id = "time", variable.name = "Compartment", value.name = "Value")

#Manipulating the Data for a stacked bar plot - to show the relative proportions from each country
#The bar chart plot will need to be when the model is at equilbirum - length(outdata) or something similar 
equidf <- data.frame(matrix(NA, ncol = 3, nrow = 3))
colnames(equidf) <- c("Category","Domestic","Imported")
equidf[1] <- c("ICombH", "Sensitive", "Resistant")
equidf[2] <- c(outdata$DIComb[length(outdata$DIComb)], outdata$IDhs[length(outdata$IDhs)], 
               outdata$IDhr[length(outdata$IDhr)]) 
equidf[3] <- c(outdata$EUIComb[length(outdata$EUIComb)], outdata$IEUhs[length(outdata$IEUhs)],
               outdata$IEUhr[length(outdata$IEUhr)]) 

###Plotting Output for Basic Script####

#Bar Chart Plot - to view the relative infecteds 
plot_ly(equidf, x = ~Category, y = ~(Domestic*100000), type = "bar", name = "Domestic") %>%
  add_trace(y = ~(European*100000), name = "European") %>%
  add_trace(y = ~(nonEuropean*100000), name = "nonEuropean") %>%
  layout(yaxis = list(title = "Proportion of Infecteds (per 100,000)", range = c(0, 15)), barmode = "stack")

#High Res - Human Infect  
ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(meltedout$Compartment == "IDhs" | meltedout$Compartment == "IDhr" | meltedout$Compartment == "IEUhs" | 
                                   meltedout$Compartment == "IEUhr" | meltedout$Compartment == "InEUhs" | 
                                   meltedout$Compartment =="InEUhr"),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
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

#Create the parameter space for the model analysis/sensitivity analysis
parms = fast_parameters(minimum = c(0.01, 0.01, 0.01, 0.0001, 0.0001, 0.0001, 0.0001, 0, 0, 0, 0.05, 0.005, 0.01, 0.01, 0.01, 0.5, 250^-1, 60^-1, 2400^-1, 
                                    288350^-1),
                        maximum = c(1, 1, 1, 0.01, 0.01, 0.01, 0.01 , 1, 1, 1, 2, 0.5, 1, 1, 1, 1, 2.5^-1, 0.6^-1, 24^-1, 2883.5^-1), 
                        factor=20, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "betahh", "tauD", "tauEU", "taunEU", 
                                             "theta", "phi", "psiEU", "psinEU", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"))

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
             betahh = parms$betahh[i],
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

#Plotting and Analysis 
sensit1 <- output$OverallInf #Creating Variable for the output variable of interest
sens1 <- sensitivity(x=sensit1, numberf=20, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "betahh",
                                                                   "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                                   "alpha", "rA" , "rH" , "eta" , "mu"))

df.equilibrium <- NULL
df.equilibrium <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "betahh",
                                             "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                             "alpha", "rA" , "rH" , "eta" , "mu"), value=sens1)

ggplot(df.equilibrium, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")

sensit2 <- output$OverallResProp #Creating Variable for the output variable of interest
sensit2[is.nan(sensit2)] <- 0

sens2 <- sensitivity(x=sensit2, numberf=20, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "betahh",
                                                                   "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                                   "alpha", "rA" , "rH" , "eta" , "mu"))

df.equilibrium1 <- NULL
df.equilibrium1 <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "betahh", 
                                             "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                             "alpha", "rA" , "rH" , "eta" , "mu"), value=sens2)

ggplot(df.equilibrium1, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")


sensit3 <- output$ImpOverInf #Creating Variable for the output variable of interest
sensit4 <- output$ImpResProp #Creating Variable for the output variable of interest
sensit4[is.nan(sensit4)] <- 0

sens4 <- sensitivity(x=sensit4, numberf=20, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "betahh",
                                                               "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                               "alpha", "rA" , "rH" , "eta" , "mu"))

df.equilibrium3 <- NULL
df.equilibrium3 <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "betahh",
                                             "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                             "alpha", "rA" , "rH" , "eta" , "mu"), value=sens4)

ggplot(df.equilibrium3, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")


#### Model For Domestic, EU and non-EU Livestock-to-human AMR transmission - seperated alpha values ####

amrhetalpha <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSDa = - betaDaa*(IDas + IDar*alphaD)*SDa + rA*(IDas+IDar) + tauD*IDas + eta - eta*SDa
    dIDas = betaDaa*IDas*SDa - rA*IDas - tauD*IDas - tauD*theta*IDas + phi*IDar - eta*IDas
    dIDar = betaDaa*IDar*SDa*alphaD + tauD*theta*IDas - phi*IDar - rA*IDar - eta*IDar 
    
    dSEUa = - betaEUaa*(IEUas + IEUar*alphaEU)*SEUa + rA*(IEUas+IEUar) + tauEU*IEUas + eta - eta*SEUa
    dIEUas = betaEUaa*IEUas*SEUa - rA*IEUas - tauEU*IEUas - tauEU*theta*IEUas + phi*IEUar - eta*IEUas
    dIEUar = betaEUaa*IEUar*SEUa*alphaEU + tauEU*theta*IEUas - phi*IEUar - rA*IEUar - eta*IEUar 
    
    dSnEUa = - betanEUaa*(InEUas + InEUar*alphanEU)*SnEUa + rA*(InEUas+InEUar) + taunEU*InEUas + eta - eta*SnEUa
    dInEUas = betanEUaa*InEUas*SnEUa - rA*InEUas - taunEU*InEUas - taunEU*theta*InEUas + phi*InEUar - eta*InEUas
    dInEUar = betanEUaa*InEUar*SnEUa*alphanEU + taunEU*theta*InEUas - phi*InEUar - rA*InEUar - eta*InEUar 
    
    dSh = - betaDha*(eta*IDas)*Sh*psiUK - betaDha*(eta*IDar)*Sh*alphaD*psiUK - 
      betaEUha*(eta*IEUas)*Sh*psiEU - betaEUha*(eta*IEUar)*Sh*psiEU*alphaEU - 
      betanEUha*(eta*InEUas)*Sh*psinEU - betanEUha*(eta*InEUar)*Sh*psinEU*alphanEU + 
      rH*(IDhs + IDhr + IEUhs + IEUhr + InEUhs + InEUhr) + mu - mu*Sh
    
    dIDhs = betaDha*(eta*IDas)*Sh*psiUK - rH*IDhs - mu*IDhs 
    dIDhr = betaDha*(eta*IDar)*Sh*psiUK*alphaD - rH*IDhr - mu*IDhr
    
    dIEUhs = betaEUha*(eta*IEUas)*Sh*psiEU - rH*IEUhs - mu*IEUhs
    dIEUhr = betaEUha*(eta*IEUar)*Sh*psiEU*alphaEU - rH*IEUhr - mu*IEUhr
    
    dInEUhs = betanEUha*(eta*InEUas)*Sh*psinEU - rH*InEUhs - mu*InEUhs
    dInEUhr = betanEUha*(eta*InEUar)*Sh*psinEU*alphanEU - rH*InEUhr - mu*InEUhr
    
    return(list(c(dSDa, dIDas, dIDar, dSEUa, dIEUas, dIEUar, dSnEUa, dInEUas, dInEUar,
                  dSh, dIDhs, dIDhr, dIEUhs, dIEUhr, dInEUhs, dInEUhr)))
  })
}

rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#### Parameters ####

init <- c(SDa = 0.98, IDas = 0.01, IDar = 0.01, 
          SEUa = 0.98, IEUas = 0.01, IEUar = 0.01, 
          SnEUa = 0.98, InEUas = 0.01, InEUar = 0.01,
          Sh = 1, IDhs = 0, IDhr = 0, IEUhs = 0, IEUhr = 0, InEUhs = 0, InEUhr = 0)

times1 <- seq(0,1000,by=1)

#Need to Specify Model Parameters

parms = c(betaDaa = 0.1, betaEUaa = 0.1, betanEUaa = 0.1, betaDha = 0.01, betaEUha = 0.01, betanEUha = 0.01,
          tauD = 0.05, tauEU = 0.08, taunEU = 0.1, 
          theta = 0.5, phi = 0.05,
          psiEU = 0.2, psinEU = 0.1, psiUK = 0.7,
          alphaD = 0.8, alphaEU = 0.8, alphanEU = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)

out <- ode(y = init, func = amrhetalpha, times = times1, parms = parms)

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

#### Sensitivity Analysis ####

start_time <- Sys.time()

#Create the parameter space for the model analysis/sensitivity analysis
parms = fast_parameters(minimum = c(0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0, 0, 0, 0.05, 0.005, 0.01, 0.01, 0.01, 0,0,0, 250^-1, 60^-1, 2400^-1, 
                                    288350^-1),
                        maximum = c(1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 2, 0.5, 1, 1, 1, 1,1,1, 2.5^-1, 0.6^-1, 24^-1, 2883.5^-1), 
                        factor=21, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", "tauD", "tauEU", "taunEU", 
                                             "theta", "phi", "psiEU", "psinEU", "psiUK", "alphaD", "alphaEU", "alphanEU", "rA" , "rH" , "eta" , "mu"))

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
             alphaD = parms$alphaD[i], alphaEU = parms$alphaEU[i], alphanEU = parms$alphanEU[i], 
             rA = parms$rA[i], rH = parms$rH[i], eta = parms$eta[i], mu = parms$mu[i])
  out <- ode(y = init, func = amrhetalpha, times = times, parms = parms1)
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

#Plotting and Analysis 
sensit1 <- output$OverallInf #Creating Variable for the output variable of interest
sens1 <- sensitivity(x=sensit1, numberf=21, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                                                   "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                                   "alphaD", "alphaEU", "alphanEU", "rA" , "rH" , "eta" , "mu"))

df.equilibrium <- NULL
df.equilibrium <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                             "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                             "alphaD", "alphaEU", "alphanEU", "rA" , "rH" , "eta" , "mu"), value=sens1)

ggplot(df.equilibrium, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")

sensit2 <- output$OverallResProp #Creating Variable for the output variable of interest
sensit2[is.nan(sensit2)] <- 0

sens2 <- sensitivity(x=sensit2, numberf=21, make.plot=T, names = c("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                                                   "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                                                   "alphaD", "alphaEU", "alphanEU", "rA" , "rH" , "eta" , "mu"))

df.equilibrium1 <- NULL
df.equilibrium1 <- data.frame(parameter=rbind("betaDaa", "betaEUaa", "betanEUaa", "betaDha", "betaEUha", "betanEUha", 
                                              "tauD", "tauEU", "taunEU", "theta", "phi", "psiEU", "psinEU", "psiUK", 
                                              "alphaD", "alphaEU", "alphanEU", "rA" , "rH" , "eta" , "mu"), value=sens2)

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


