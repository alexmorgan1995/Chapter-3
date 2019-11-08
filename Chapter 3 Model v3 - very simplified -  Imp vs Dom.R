rm(list=ls())

library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("reshape2")

#### Model For Domestic, EU and non-EU Livestock-to-human AMR transmission ####

amrhet <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSDa = - betaDaa*(IDas + IDar*alpha)*SDa + rA*(IDas+IDar) + tauD*IDas*kappa1 + eta - eta*SDa
    dIDas = betaDaa*IDas*SDa - rA*IDas - tauD*IDas*kappa1 - tauD*theta*IDas + phi*IDar - eta*IDas
    dIDar = betaDaa*IDar*SDa*alpha + tauD*theta*IDas - phi*IDar - rA*IDar - eta*IDar 
    
    dSIMPa = - betaIMPaa*(IIMPas + IIMPar*alpha)*SIMPa + rA*(IIMPas+IIMPar) + tauIMP*IIMPas*kappa1 + eta - eta*SIMPa
    dIIMPas = betaIMPaa*IIMPas*SIMPa - rA*IIMPas - tauIMP*IIMPas*kappa1 - tauIMP*theta*IIMPas + phi*IIMPar - eta*IIMPas
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

#### Parameters + Model Run ####

init <- c(SDa = 0.98, IDas = 0.01, IDar = 0.01, 
          SIMPa = 0.98, IIMPas = 0.01, IIMPar = 0.01, 
          Sh = 1, IDhs = 0, IDhr = 0, IIMPhs = 0, IIMPhr = 0)
times1 <- seq(0,5000,by=1)

#Need to Specify Model Parameters

parms = c(betaDaa = 0.1, betaIMPaa = 0.1, betaDha = 0.01, betaIMPha = 0.01, 
          betahh = 0.01,
          tauD = 0.039, tauIMP = 0.0876,
          theta = 0.5, phi = 0.05, kappa1 = 1,
          psiIMP = 0.3099, psiUK = 0.593,
          alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)
          
out <- ode(y = init, func = amrhet, times = times1, parms = parms)

outdata <- data.frame(out)
outdata$DIComb <- outdata$IDhs + outdata$IDhr 
outdata$IMPIComb <- outdata$IIMPhs + outdata$IIMPhr 

#Creating L I Q U I D code - i.e - Melting for GGPLOT
outdata1 <- outdata; outdata1[8:14] <- outdata[8:14]*100000 # So only human compartments are scaled to per 100,000 population 
meltedout <- melt(outdata1, id = "time", variable.name = "Compartment", value.name = "Value")

#Plotting for Humans
ggplot(data = meltedout, aes(x = time, y = Value, col = Compartment)) + 
  geom_line(data=(y = meltedout[(as.numeric(which((meltedout$Compartment == "DIComb" & meltedout$time == 0))):
                                   length(meltedout$time)),]), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population (per 100,000)") +
  scale_y_continuous(limits = c(0,5), expand = c(0,0)) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Plotting the relative proportions - Bar Chart
equidf <- data.frame(matrix(NA, ncol = 3, nrow = 3))
colnames(equidf) <- c("Category","Domestic","Imported")
equidf[1] <- c("ICombH", "Sensitive", "Resistant")
equidf[2] <- c(outdata$DIComb[length(outdata$DIComb)], outdata$IDhs[length(outdata$IDhs)], 
               outdata$IDhr[length(outdata$IDhr)]) 
equidf[3] <- c(outdata$IMPIComb[length(outdata$IMPIComb)], outdata$IIMPhs[length(outdata$IIMPhs)],
               outdata$IIMPhr[length(outdata$IIMPhr)]) 

#Bar Chart Plot - to view the relative infecteds 
plot_ly(equidf, x = ~Category, y = ~(Domestic*100000), type = "bar", name = "Domestic") %>%
  add_trace(y = ~(Imported*100000), name = "Imported") %>%
  layout(yaxis = list(title = "Proportion of Infecteds (per 100,000)", range = c(0, 10)), barmode = "stack")

#### Sensitivity Analysis ####

start_time <- Sys.time()

parms = fast_parameters(minimum = c(0.01, 0.01, 
                                    0.0001, 0.0001, 0.0001, 
                                    0, 0, 
                                    0, 0, 0, 
                                    0.01, 0.01, 0, 
                                    250^-1, 60^-1, 2400^-1, 288350^-1),
                        maximum = c(1, 1, 
                                    0.01, 0.01, 0.01, 
                                    1, 1, 
                                    1, 1, 1, 
                                    1, 1, 1, 
                                    2.5^-1, 0.6^-1, 24^-1, 2883.5^-1), 
                        factor= 17, names = c("betaDaa", "betaIMPaa", 
                                              "betaDha", "betaIMPha", "betahh", 
                                              "tauD", "tauIMP", 
                                              "theta", "phi", "kappa1", 
                                              "psiIMP", "psiUK", "alpha", 
                                              "rA" , "rH" , "eta" , "mu"))

init <- c(SDa = 0.98, IDas = 0.01, IDar = 0.01, 
          SIMPa = 0.98, IIMPas = 0.01, IIMPar = 0.01, 
          Sh = 1, IDhs = 0, IDhr = 0, IIMPhs = 0, IIMPhr = 0)

times <- seq(0,10000, by = 10) 

output <- data.frame()

for (i in 1:nrow(parms)) {
  temp <- data.frame(matrix(NA, nrow = 1, ncol=8))
  parms1 = c(betaDaa = parms$betaDaa[i], betaIMPaa = parms$betaIMPaa[i], 
             betaDha = parms$betaDha[i], betaIMPha = parms$betaIMPha[i], betahh = parms$betahh[i],
             tauD = parms$tauD[i], tauIMP = parms$tauIMP[i], 
             theta = parms$theta[i], phi = parms$phi[i], kappa1 = parms$kappa1[i],
             psiIMP = (parms$psiIMP[i]/(parms$psiUK[i]+parms$psiIMP[i])), psiUK = (parms$psiUK[i]/(parms$psiUK[i]+parms$psiIMP[i])),
             alpha = parms$alpha[i], rA = parms$rA[i], rH = parms$rH[i], eta = parms$eta[i], mu = parms$mu[i])
  out <- ode(y = init, func = amrhet, times = times, parms = parms1)
  temp[1,1] <- rounding(sum(as.numeric(out[(length(out[,1])), (9:12)])))
  temp[1,2] <- rounding(sum(as.numeric(out[(length(out[,1])), c(10, 12)]))) /  temp[1,1]
  temp[1,3] <- rounding(sum(as.numeric(out[(length(out[,1])), c(9, 11)]))) /  temp[1,1]
  temp[1,4] <- rounding(sum(as.numeric(out[(length(out[,1])), c(9, 10)])))
  temp[1,5] <- rounding(out[nrow(out),10]) / temp[1,4]
  temp[1,6] <- rounding(sum(as.numeric(out[(length(out[,1])), c(11, 12)])))
  temp[1,7] <- rounding(out[nrow(out),12]) / temp[1,6]
  temp[1,8] <- as.numeric(parms1[6])
  print(temp[1,8])
  output <- rbind.data.frame(output, temp)
}

colnames(output)[1:8] <- c("OverallInf","OverallResProp","OverallSensProp","DomOverInf","DomResProp", 
                           "ImpOverInf","ImpResProp","UKUse")

end_time <- Sys.time(); end_time - start_time

####Plotting Of Sensitivity Analysis ##### 

#OverallInf
sensit1 <- output$OverallInf #Creating Variable for the output variable of interest
sens1 <- sensitivity(x=sensit1, numberf=17, make.plot=T, names = c("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                                                   "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"))
df.equilibrium <- NULL
df.equilibrium <- data.frame(parameter=rbind("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                             "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"), value=sens1)
ggplot(df.equilibrium, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")

#OverallResProp
sensit2 <- output$OverallResProp #Creating Variable for the output variable of interest
sensit2[is.nan(sensit2)] <- 0

sens2 <- sensitivity(x=sensit2, numberf=17, make.plot=T, names = c("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                                                   "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"))
df.equilibrium1 <- NULL
df.equilibrium1 <- data.frame(parameter=rbind("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                              "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"), value=sens2)
ggplot(df.equilibrium1, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")

#ImpOverInf
sensit3 <- output$ImpOverInf #Creating Variable for the output variable of interest

sens3 <- sensitivity(x=sensit3, numberf=17, make.plot=T, names = c("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                                                   "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"))
df.equilibrium3 <- NULL
df.equilibrium3 <- data.frame(parameter=rbind("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                              "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"), value=sens3)
ggplot(df.equilibrium3, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")

#ImpResProp
sensit4 <- output$ImpResProp #Creating Variable for the output variable of interest
sensit4[is.nan(sensit4)] <- 0

sens4 <- sensitivity(x=sensit4, numberf=17, make.plot=T, names = c("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                                                   "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"))
df.equilibrium4 <- NULL
df.equilibrium4 <- data.frame(parameter=rbind("betaDaa", "betaIMPaa", "betaDha", "betaIMPha", "betahh", "tauD", "tauIMP",  
                                              "theta", "phi", "kappa1", "psiIMP", "psiUK", "alpha", "rA" , "rH" , "eta" , "mu"), value=sens4)
ggplot(df.equilibrium4, aes(parameter, value)) + geom_bar(stat="identity", fill="grey23")
