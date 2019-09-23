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

amrhet <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSDa = - betaDaa*(IDas + IDar*alpha) + rA*(IDas+IDar) + tauD*IDas + eta - eta*SDa
    dIDas = betaDaa*IDas*SDa - rA*IDas - tauD*IDas - tauD*theta*IDas + phi*IDar - eta*IDas
    dIDar = betaDaa*IDar*SDa*alpha + tauD*theta*IDas - phi*IDar - rA*IDar - eta*IDar 
    
    dSEUa = - betaEUaa*(IEUas + IEUar*alpha) + rA*(IEUas+IEUar) + tauEU*IEUas + eta - eta*SEUa
    dIEUas = betaEUaa*IEUas*SEUa - rA*IEUas - tauEU*IEUas - tauEU*theta*IEUas + phi*IEUar - eta*IEUas
    dIEUar = betaEUaa*IEUar*SEUa*alpha + tauEU*theta*IEUas - phi*IEUar - rA*IEUar - eta*IEUar 
    
    dSnEUa = - betanEUaa*(InEUas + InEUar*alpha) + rA*(InEUas+InEUar) + taunEU*InEUas + eta - eta*SnEUa
    dInEUas = betanEUaa*InEUas*SnEUa - rA*InEUas - taunEU*InEUas - taunEU*theta*InEUas + phi*InEUar - eta*InEUas
    dInEUar = betanEUaa*InEUar*SnEUa*alpha + taunEU*theta*InEUas - phi*InEUar - rA*InEUar - eta*InEUar 
    
    dSh = - betaDha*IDas*Sh*(1 - psiEU + psinEU)*(eta*IDas) - betaDha*IDar*Sh*alpha*(1 - psiEU + psinEU)*(eta*IDar) - 
      betaEUha*IEUas*Sh*psiEU*(eta*IEUas) - betaEUha*IEUar*Sh*psiEU*alpha*(eta*IEUar) - 
      betanEUha*InEUas*Sh*psinEU*(eta*InEUas) - betanEUha*InEUar*Sh*psinEU*alpha*(eta*InEUar) + 
      rH*(IDhs + IDhr + IEUhs + IEUhr + InEUhs + InEUhr) + mu - mu*Sh
    dIDhs = betaDha*IDas*Sh*(1 - psiEU + psinEU)*(eta*IDas) - rH*IDhs - mu*IDhs 
    dIDhr = betaDha*IDar*Sh*alpha*(1 - psiEU + psinEU)*(eta*IDar) - rH*IDhr - mu*IDhr
    dIEUhs = betaEUha*IEUas*Sh*psiEU*(eta*IEUas) - rH*IEUhs - mu*IEUhs
    dIEUhr = betaEUha*IEUar*Sh*psiEU*alpha*(eta*IEUar) - rH*IEUhr - mu*IEUhr
    dInEUhs = betanEUha*InEUas*Sh*psinEU*(eta*InEUas) - rH*InEUhs - mu*InEUhs
    dInEUhr = betanEUha*InEUar*Sh*psinEU*alpha*(eta*InEUar) - rH*InEUhr - mu*InEUhr
    
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

times1 <- seq(0,10000,by=1)

#Need to Specify Model Parameters

parms = c(betaDaa = 0.1, betaEUaa = 0.1, betanEUaa = 0.1, betaDha = 0.01, betaEUha = 0.01, betanEUha = 0.01,
          tauD = 0.05, tauEU = 0.08, taunEU = 0.1, 
          theta = 0.5, phi = 0.05,
          psiEU = 0.2, psinEU = 0.1,
          alpha = 0.8, rA = 25^-1, rH = 6^-1, eta = 240^-1, mu = 28835^-1)
          
out <- ode(y = init, func = amrhet, times = times1, parms = parms)

outdata <- data.frame(out)
outdata$DIComb <- outdata$IDhs + outdata$IDhr 
outdata$EUIComb <- outdata$IEUhs + outdata$IEUhr 
outdata$nEUIComb <- outdata$InEUhs + outdata$InEUhr

#Manipulating the Data for a stacked bar plot - to show the relative proportions from each country
equidf <- data.frame(matrix(NA, ncol = 4, nrow = 3))
colnames(equidf) <- c("Category","Domestic","European", "nonEuropean")
equidf[1] <- c("ICombH", "Sensitive", "Resistant")
equidf[2] <- c(outdata$DIComb[length(outdata$DIComb)], outdata$IDhs[length(outdata$IDhs)], 
               outdata$IDhr[length(outdata$IDhr)]) 
equidf[3] <- c(outdata$EUIComb[length(outdata$DIComb)], outdata$IEUhs[length(outdata$IEUhs)],
               outdata$IEUhr[length(outdata$IEUhr)]) 
equidf[4] <- c(outdata$nEUIComb[length(outdata$nEUIComb)], outdata$InEUhs[length(outdata$InEUhs)], 
               outdata$InEUhr[length(outdata$InEUhr)])  

#The bar chart plot will need to be when the model is at equilbirum - length(outdata) or something similar 

###Plotting Output for Basic Script####

#Simplified Human Population - For all human populations 
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sh, colour = "Sh"), size = 1.1) +
  geom_line(aes(y = IDhs, colour = "IDhs"), size = 1.1) + 
  geom_line(aes(y = IDhr, colour = "IDhr"), size = 1.1) +
  geom_line(aes(y = IEUhs, colour = "IEUhs"), size = 1.1) +
  geom_line(aes(y = IEUhr, colour = "IEUhr"), size = 1.1) + 
  geom_line(aes(y = InEUhs, colour = "InEUhs"), size = 1.1) +
  geom_line(aes(y = InEUhr, colour = "InEUhr"), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_x_continuous(limits = c(0,10000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000001), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Bar Chart Plot 
plot_ly(equidf, x = ~Category, y = ~Domestic, type = "bar", name = "Domestic") %>%
  add_trace(y = ~European, name = "European") %>%
  add_trace(y = ~nonEuropean, name = "nonEuropean") %>%
  layout(yaxis = list(title = "Proportion of Infecteds"), barmode = "stack")


#Simplified Animal Populations  
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sh, colour = "Sh"), size = 1.1) +
  geom_line(aes(y = IDhs, colour = "IDhs"), size = 1.1) + 
  geom_line(aes(y = IDhr, colour = "IDhr"), size = 1.1) +
  geom_line(aes(y = IEUhs, colour = "IEUhs"), size = 1.1) +
  geom_line(aes(y = IEUhr, colour = "IEUhr"), size = 1.1) + 
  geom_line(aes(y = InEUhs, colour = "InEUhs"), size = 1.1) +
  geom_line(aes(y = InEUhr, colour = "InEUhr"), size = 1.1) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_x_continuous(limits = c(0,10000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.000001), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))

#Human Population 
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sh, colour = "Sh"), size = 1.1) +
  geom_line(aes(y = Ih, colour = "Ih"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Irh, colour = "Irh"), size = 1.1) +
  geom_line(aes(y = combIh, colour = "combIh"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Sh" = "chartreuse3", "Ih" = "red", "Irh" = "blue", "combIh" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombH"],")")),
                                expression(paste("Resistant Infection ","(I"["RH"],")")),
                                expression(paste("Sensitive Infection ","(I"["H"],")")),
                                expression(paste("Susceptible ","(S"["H"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Human Population") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))


#Animal Population 
par(mfrow=c(1,2))
ggplot(outdata, aes(time)) +
  geom_line(aes(y = Sa, colour = "Sa"), size = 1.1) +
  geom_line(aes(y = Ia, colour = "Ia"), size = 2.5, alpha =0.6) + 
  geom_line(aes(y = Ira, colour = "Ira"), size = 1.1) +
  geom_line(aes(y = combIa, colour = "combIa"), size = 1.1, linetype = "dashed") +
  scale_color_manual(values = c("Sa" = "chartreuse3", "Ia" = "red", "Ira" = "blue", "combIa" = "black"), 
                     labels = c(expression(paste("Overall Infection ","(I"["CombA"],")")),
                                expression(paste("Resistant Infection ","(I"["RA"],")")),
                                expression(paste("Sensitive Infection ","(I"["A"],")")),
                                expression(paste("Susceptible ","(S"["A"],")"))), 
                     guide = guide_legend(reverse = TRUE)) +
  labs(x ="Time (Days)", y = "Proportion of Animal Population") +
  scale_x_continuous(limits = c(0,1000), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1.01), expand = c(0,0)) +
  theme(axis.line.x = element_line(color="black", size = 1),
        legend.position="bottom", legend.title = element_blank(),
        legend.spacing.x = unit(0.2, 'cm'), legend.text=element_text(size=11), plot.margin=unit(c(0.7,0.7,0.8,0.8),"cm"))


