library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate")

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
      (1-psi)*(imp3)*(1-alpha)*(betaHI3*fracimp3*propres_imp3*Sh) -
      
      (1-psi)*(imp4)*(betaHI4*fracimp4*(1-propres_imp4)*Sh) - 
      (1-psi)*(imp4)*(1-alpha)*(betaHI4*fracimp4*propres_imp4*Sh)-
      
      (1-psi)*(imp5)*(betaHI5*fracimp5*(1-propres_imp5)*Sh) - 
      (1-psi)*(imp5)*(1-alpha)*(betaHI5*fracimp5*propres_imp5*Sh) - 
      
      (1-psi)*(imp6)*(betaHI6*fracimp6*(1-propres_imp6)*Sh) - 
      (1-psi)*(imp6)*(1-alpha)*(betaHI6*fracimp6*propres_imp6*Sh) - 
      
      (1-psi)*(imp7)*(betaHI7*fracimp7*(1-propres_imp7)*Sh) - 
      (1-psi)*(imp7)*(1-alpha)*(betaHI7*fracimp7*propres_imp7*Sh) - 
      
      (1-psi)*(imp8)*(betaHI8*fracimp8*(1-propres_imp8)*Sh) - 
      (1-psi)*(imp8)*(1-alpha)*(betaHI8*fracimp8*propres_imp8*Sh) -
      
      (1-psi)*(imp9)*(betaHI9*fracimp9*(1-propres_imp9)*Sh) - 
      (1-psi)*(imp9)*(1-alpha)*(betaHI9*fracimp1*propres_imp9*Sh)-
      
      (1-psi)*(imp10)*(betaHI10*fracimp10*(1-propres_imp10)*Sh) - 
      (1-psi)*(imp10)*(1-alpha)*(betaHI10*fracimp10*propres_imp10*Sh) 
      
    dIshDA = psi*betaHD*Isa*Sh  - rh*IshDA - uh*IshDA 
    dIrhDA = psi*(1-alpha)*betaHD*Ira*Sh - rh*IrhDA - uh*IrhDA  
    
    dIshA1 = (1-psi)*(imp1)*(betaHI1*fracimp1*(1-propres_imp1)*Sh) - rh*IshA1 - uh*IshA1 
    dIrhA1 = (1-psi)*(imp1)*(1-alpha)*(betaHI1*fracimp1*propres_imp1*Sh) - rh*IrhA1 - uh*IrhA1  
    
    dIshA2 = (1-psi)*(imp2)*(betaHI2*fracimp2*(1-propres_imp2)*Sh) - rh*IshA2 - uh*IshA2 
    dIrhA2 = (1-psi)*(imp2)*(1-alpha)*(betaHI2*fracimp2*propres_imp2*Sh)  - rh*IrhA2 - uh*IrhA2  
    
    dIshA3 = (1-psi)*(imp3)*(betaHI3*fracimp3*(1-propres_imp3)*Sh) - rh*IshA3 - uh*IshA3 
    dIrhA3 = (1-psi)*(imp3)*(1-alpha)*(betaHI3*fracimp3*propres_imp3*Sh) - rh*IrhA3 - uh*IrhA3  
    
    dIshA4 = (1-psi)*(imp4)*(betaHI4*fracimp4*(1-propres_imp4)*Sh) - rh*IshA4 - uh*IshA4 
    dIrhA4 = (1-psi)*(imp4)*(1-alpha)*(betaHI4*fracimp1*propres_imp4*Sh) - rh*IrhA4 - uh*IrhA4  
    
    dIshA5 = (1-psi)*(imp5)*(betaHI5*fracimp5*(1-propres_imp5)*Sh) - rh*IshA5 - uh*IshA5 
    dIrhA5 = (1-psi)*(imp5)*(1-alpha)*(betaHI5*fracimp5*propres_imp5*Sh)  - rh*IrhA5 - uh*IrhA5  
    
    dIshA6 = (1-psi)*(imp6)*(betaHI6*fracimp6*(1-propres_imp6)*Sh) - rh*IshA6 - uh*IshA6 
    dIrhA6 = (1-psi)*(imp6)*(1-alpha)*(betaHI6*fracimp6*propres_imp6*Sh) - rh*IrhA6 - uh*IrhA6  
    
    dIshA7 = (1-psi)*(imp7)*(betaHI7*fracimp7*(1-propres_imp7)*Sh) - rh*IshA7 - uh*IshA7 
    dIrhA7 = (1-psi)*(imp7)*(1-alpha)*(betaHI7*fracimp7*propres_imp7*Sh) - rh*IrhA7 - uh*IrhA7  
    
    dIshA8 = (1-psi)*(imp8)*(betaHI8*fracimp8*(1-propres_imp8)*Sh) - rh*IshA8 - uh*IshA8 
    dIrhA8 = (1-psi)*(imp8)*(1-alpha)*(betaHI8*fracimp8*propres_imp8*Sh)  - rh*IrhA8 - uh*IrhA8  
    
    dIshA9 = (1-psi)*(imp9)*(betaHI9*fracimp9*(1-propres_imp9)*Sh) - rh*IshA9 - uh*IshA9 
    dIrhA9 = (1-psi)*(imp9)*(1-alpha)*(betaHI9*fracimp9*propres_imp9*Sh) - rh*IrhA9 - uh*IrhA9  
    
    dIshA10 = (1-psi)*(imp10)*(betaHI10*fracimp10*(1-propres_imp10)*Sh) - rh*IshA10 - uh*IshA10 
    dIrhA10 = (1-psi)*(imp10)*(1-alpha)*(betaHI10*fracimp1*propres_imp10*Sh) - rh*IrhA10 - uh*IrhA10 
    
    dIshH = betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshA4+IshA5+IshA6+IshA7+IshA8+IshA9+IshA10+IshH) - rh*IshH - uh*IshH 
    dIrhH = (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhA4+IrhA5+IrhA6+IrhA7+IrhA8+IrhA9+IrhA10+IrhH)- rh*IrhH - uh*IrhH 
    
    return(list(c(dSa,dIsa,dIra,
                  dSh,
                  dIshDA,dIrhDA,
                  dIshA1,dIrhA1,
                  dIshA2,dIrhA2,
                  dIshA3,dIrhA3,
                  dIshA4,dIrhA4,
                  dIshA5,dIrhA5,
                  dIshA6,dIrhA6,
                  dIshA7,dIrhA7,
                  dIshA8,dIrhA8,
                  dIshA9,dIrhA9,
                  dIshA10,dIrhA10,
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
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)

parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = 0.0123,
          betaHD = (0.00001),  
          
          betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
          betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
          
          phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = 0.656, 
          
          fracimp1 = 0.076719577, fracimp2 = 0.046997389, fracimp3 = 0.022156129, fracimp4 = 0.116216976, fracimp5 = 0.083127572, fracimp6 = 0.057626057,
          fracimp7 = 0.04588015, fracimp8 = 0.010497067, fracimp9 = 0.05, fracimp10 = 0.2, 
          
          imp1 = 0.1642031, imp2 = 0.1356062, imp3 = 0.1355515, imp4 = 0.134563, imp5 = 0.104059, imp6 = 0.08818593,
          imp7 = 0.07600404, imp8 = 0.05488995, imp9 = 0.03716985, imp10 = 0.06976744,
          
          propres_imp1 = 0.347826087, propres_imp2 = 0.258064516, propres_imp3 = 0.266666667, propres_imp4 = 0.588235294, propres_imp5 = 0.674698795, 
          propres_imp6 = 0.54822335, propres_imp7 = 0.436893204, propres_imp8 = 0.4, propres_imp9 = 0.456, propres_imp10 = 0.5)


out <- ode(y = init, func = amrimp, times = times, parms = parms)

tail(out)

colnames(output1) <- c("tau", "ICombH","IResRat")
output1$IResRat[output1$tau == 0]

sum(out[nrow(out), seq(7,29,by = 2)])/ sum(out[nrow(out),6:29])

1-(0.166406/0.3129203)

# Modelling Changes in Tau ------------------------------------------------

parmtau <- seq(0,0.035, by = 0.002)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)

output1 <- data.frame(matrix(ncol = 13, nrow = 0))
times <- seq(0, 200000, by = 100)

for (i in 1:length(parmtau)) {
  
  temp <- data.frame(matrix(NA, nrow = 1, ncol=13))
  
  parms =  c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = parmtau[i],
             betaHD = (0.00001),  
             
             betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
             betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
             
             phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = 0.656, 
             
             fracimp1 = 0.076719577, fracimp2 = 0.046997389, fracimp3 = 0.022156129, fracimp4 = 0.116216976, fracimp5 = 0.083127572, fracimp6 = 0.057626057,
             fracimp7 = 0.04588015, fracimp8 = 0.010497067, fracimp9 = 0.05, fracimp10 = 0.2, 
             
             imp1 = 0.1642031, imp2 = 0.1356062, imp3 = 0.1355515, imp4 = 0.134563, imp5 = 0.104059, imp6 = 0.08818593,
             imp7 = 0.07600404, imp8 = 0.05488995, imp9 = 0.03716985, imp10 = 0.06976744,
             
             propres_imp1 = 0.347826087, propres_imp2 = 0.258064516, propres_imp3 = 0.266666667, propres_imp4 = 0.588235294, propres_imp5 = 0.674698795, 
             propres_imp6 = 0.54822335, propres_imp7 = 0.436893204, propres_imp8 = 0.4, propres_imp9 = 0.456, propres_imp10 = 0.5)
  
  
  out <- ode(y = init, func = amrimp, times = times, parms = parms)
  
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),7])
  temp[1,3] <- rounding(out[nrow(out),9])
  temp[1,4] <- rounding(out[nrow(out),11])
  temp[1,5] <- rounding(out[nrow(out),13])
  temp[1,6] <- rounding(out[nrow(out),15])
  temp[1,7] <- rounding(out[nrow(out),17])
  temp[1,8] <- rounding(out[nrow(out),19])
  temp[1,9] <- rounding(out[nrow(out),21]) 
  temp[1,10] <- rounding(out[nrow(out),23])
  temp[1,11] <- rounding(out[nrow(out),25])
  temp[1,12] <- rounding(out[nrow(out),27])
  
  print(temp[1,2])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:12] <- c("tau", "Domestic","Netherlands","Ireland","Germany","France","Spain","Italy",
                             "Belgium","Poland","Denmark","Non EU")


plotdata <- melt(output1,
                 id.vars = c("tau"), measure.vars = c("Domestic","Netherlands","Ireland","Germany","France","Spain","Italy",
                                                      "Belgium","Poland","Denmark","Non EU")) 


res_comb_prev <- ggplot(plotdata, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.002) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Tetracycline Sales in Fattening Pig (g/PCU)", y = "Infected Humans (per 100,000)") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  +
  labs(x ="Domestic Tetracycline Usage in Fattening Pig (g/PCU)", y = "Prevalence of Resistant Human Salmenollosis (per 100,000 pop)", fill = "Resistance Source")  

ggsave(res_comb_prev, filename = "prev_multi_count_res.png", dpi = 300, type = "cairo", width = 8, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")

# Modelling Changes in Tau NORMALISED------------------------------------------------

parmtau <- seq(0,0.035, by = 0.002)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)

output1 <- data.frame(matrix(ncol = 13, nrow = 0))
times <- seq(0, 200000, by = 100)

for (i in 1:length(parmtau)) {
  
  temp <- data.frame(matrix(NA, nrow = 1, ncol=13))
  
  parms =  c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = parmtau[i],
             betaHD = (0.00001),  
             
             betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
             betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
             
             phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = 0.656, 
             
             fracimp1 = 0.076719577, fracimp2 = 0.046997389, fracimp3 = 0.022156129, fracimp4 = 0.116216976, fracimp5 = 0.083127572, fracimp6 = 0.057626057,
             fracimp7 = 0.04588015, fracimp8 = 0.010497067, fracimp9 = 0.05, fracimp10 = 0.2, 
             
             imp1 = 0.1642031, imp2 = 0.1356062, imp3 = 0.1355515, imp4 = 0.134563, imp5 = 0.104059, imp6 = 0.08818593,
             imp7 = 0.07600404, imp8 = 0.05488995, imp9 = 0.03716985, imp10 = 0.06976744,
             
             propres_imp1 = 0.347826087, propres_imp2 = 0.258064516, propres_imp3 = 0.266666667, propres_imp4 = 0.588235294, propres_imp5 = 0.674698795, 
             propres_imp6 = 0.54822335, propres_imp7 = 0.436893204, propres_imp8 = 0.4, propres_imp9 = 0.456, propres_imp10 = 0.5)
  
  
  out <- ode(y = init, func = amrimp, times = times, parms = parms)
  
  totalAMR <- rounding(out[nrow(out),7]) + rounding(out[nrow(out),9]) + rounding(out[nrow(out),11]) + rounding(out[nrow(out),13]) + 
    rounding(out[nrow(out),15]) + rounding(out[nrow(out),17]) + rounding(out[nrow(out),19]) + rounding(out[nrow(out),21]) + 
    rounding(out[nrow(out),23]) + rounding(out[nrow(out),25]) + rounding(out[nrow(out),27]) + rounding(out[nrow(out),29])
  
  temp[1,1] <- parmtau[i]
  temp[1,2] <- rounding(out[nrow(out),7])/ totalAMR
  temp[1,3] <- rounding(out[nrow(out),9])/ totalAMR  
  temp[1,4] <- rounding(out[nrow(out),11])/ totalAMR
  temp[1,5] <- rounding(out[nrow(out),13])/ totalAMR 
  temp[1,6] <- rounding(out[nrow(out),15])/ totalAMR
  temp[1,7] <- rounding(out[nrow(out),17])/ totalAMR  
  temp[1,8] <- rounding(out[nrow(out),19])/ totalAMR
  temp[1,9] <- rounding(out[nrow(out),21])/ totalAMR 
  temp[1,10] <- rounding(out[nrow(out),23])/ totalAMR
  temp[1,11] <- rounding(out[nrow(out),25])/ totalAMR 
  temp[1,12] <- rounding(out[nrow(out),27])/ totalAMR 
  
  print(temp[1,2])
  output1 <- rbind.data.frame(output1, temp)
}

colnames(output1)[1:12] <- c("tau", "Domestic","Netherlands","Ireland","Germany","France","Spain","Italy",
                            "Belgium","Poland","Denmark","Non EU")


plotdata <- melt(output1,
                 id.vars = c("tau"), measure.vars = c("Domestic","Netherlands","Ireland","Germany","France","Spain","Italy",
                                                      "Belgium","Poland","Denmark","Non EU")) 


res_comb <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.002) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Tetracycline Sales in Fattening Pig (g/PCU)", y = "Infected Humans (per 100,000)") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  +
  labs(x ="Domestic Tetracycline Usage in Fattening Pig (g/PCU)", y = "Proportion of Attributable Human Resistance", fill = "Resistance Source")  

ggsave(res_comb, filename = "multi_count_res.png", dpi = 300, type = "cairo", width = 8, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")


# Testing Effect of Import ------------------------------------------------

#Use the combo parameter to create a combination of: 
# - Just resistance 
# - Just contamination
# - all import parameters 
# - all of the above 
# Then create a volcano plots with the means plotted 
# Have a red line showing the baseline change in resistance 

#Contamination
parameterspace_1 <- data.frame("fracimp1" = runif(10000, 0, 1), "fracimp2" = runif(10000, 0, 1), "fracimp3" = runif(10000, 0, 1), "fracimp4" = runif(10000, 0, 1), 
                             "fracimp5" = runif(10000, 0, 1), "fracimp6" = runif(10000, 0, 1), "fracimp7" = runif(10000, 0, 1), "fracimp8" = runif(10000, 0, 1), 
                             "fracimp9" = runif(10000, 0, 1), "fracimp10" = runif(10000, 0, 1))
#Resistance
parameterspace_2 <- data.frame("propres_imp1" = runif(10000, 0, 1), "propres_imp2" = runif(10000, 0, 1), "propres_imp3" = runif(10000, 0, 1), 
                               "propres_imp4" = runif(10000, 0, 1), "propres_imp5" = runif(10000, 0, 1), "propres_imp6" = runif(10000, 0, 1), 
                               "propres_imp7" = runif(10000, 0, 1), "propres_imp8" = runif(10000, 0, 1), "propres_imp9" = runif(10000, 0, 1), 
                               "propres_imp10" = runif(10000, 0, 1))

#Proportion of Import Usage
t <- RandVec(a=0, b=1, s=1, n=10, m=10000, Seed=sample(1:10000, size = 1))
parameterspace_3 <- data.frame(matrix(unlist(t), nrow=10000, byrow=TRUE),stringsAsFactors=FALSE)

#Proportion of Domestic Usage 
parameterspace_4 <- seq(0,1, by = 0.0001)[-10000]

#All 
parameterspace_comb <- cbind(parameterspace_1, parameterspace_2, parameterspace_3, parameterspace_4)
colnames(parameterspace_comb) <- c(colnames(parameterspace_1), colnames(parameterspace_2), paste0("imp", 1:10), "psi")

# Parameter set 1 ---------------------------------------------------------

#Run the model

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)
times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,0.0123)

usage_frame_1 <- data.frame(matrix(nrow = 0, ncol = 2))

for(x in 1:nrow(parameterspace_1)) {
  dump_data <- vector() 
  output1 <- data.frame(matrix(nrow = 2, ncol =3))
    for (i in 1:2) {
      temp <- data.frame(matrix(nrow = 1, ncol =3))
      
      parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = parmtau1[i],
                 betaHD = (0.00001),  
                 
                 betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
                 betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
                 
                 phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = 0.656, 
                 
                 fracimp1 = parameterspace_1[x,1], fracimp2 = parameterspace_1[x,2], fracimp3 = parameterspace_1[x,3], fracimp4 = parameterspace_1[x,4], fracimp5 = parameterspace_1[x,5], fracimp6 = parameterspace_1[x,6],
                 fracimp7 = parameterspace_1[x,7], fracimp8 = parameterspace_1[x,8], fracimp9 = parameterspace_1[x,9], fracimp10 = parameterspace_1[x,10], 
                 
                 imp1 = 0.1642031, imp2 = 0.1356062, imp3 = 0.1355515, imp4 = 0.134563, imp5 = 0.104059, imp6 = 0.08818593,
                 imp7 = 0.07600404, imp8 = 0.05488995, imp9 = 0.03716985, imp10 = 0.06976744,
                 
                 propres_imp1 = 0.347826087, propres_imp2 = 0.258064516, propres_imp3 = 0.266666667, propres_imp4 = 0.588235294, propres_imp5 = 0.674698795, 
                 propres_imp6 = 0.54822335, propres_imp7 = 0.436893204, propres_imp8 = 0.4, propres_imp9 = 0.456, propres_imp10 = 0.5)
      
      
      out <- ode(y = init, func = amrimp, times = times, parms = parms2)
      temp[,1] <- parmtau1[i]
      temp[,2] <- sum(out[nrow(out),6:29])
      temp[,3] <- signif(sum(out[nrow(out), seq(7,29,by = 2)])/ temp[,2], digits = 3)
      output1[i,] <- temp
    }
    print(output1)
    colnames(output1) <- c("tau", "ICombH","IResRat")
    dump_data[1] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == 0.0123]))*100 #% Reduction
    dump_data[2] <- "Parameter Set 1"
    usage_frame_1 <- rbind(usage_frame_1, dump_data)
    
    print(paste0("Parameter Set 1 - ", round(x/10000, digits = 4)*100, "%"))
}
colnames(usage_frame_1) <- c("relRes", "Parameter")

# Parameter set 2 ---------------------------------------------------------

#Run the model

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)
times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,0.0123)

usage_frame_2 <- data.frame(matrix(nrow = 0, ncol = 2))

for(x in 1:nrow(parameterspace_1)) {
  dump_data <- vector() 
  output1 <- data.frame(matrix(nrow = 2, ncol =3))
  for (i in 1:2) {
    temp <- data.frame(matrix(nrow = 1, ncol =3))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = parmtau1[i],
               betaHD = (0.00001),  
               
               betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
               betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
               
               phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = 0.656, 
               
               fracimp1 = 0.076719577, fracimp2 = 0.046997389, fracimp3 = 0.022156129, fracimp4 = 0.116216976, fracimp5 = 0.083127572, fracimp6 = 0.057626057,
               fracimp7 = 0.04588015, fracimp8 = 0.010497067, fracimp9 = 0.05, fracimp10 = 0.2, 
               
               imp1 = 0.1642031, imp2 = 0.1356062, imp3 = 0.1355515, imp4 = 0.134563, imp5 = 0.104059, imp6 = 0.08818593,
               imp7 = 0.07600404, imp8 = 0.05488995, imp9 = 0.03716985, imp10 = 0.06976744,
               
               propres_imp1 = parameterspace_2[x,1], propres_imp2 = parameterspace_2[x,2], propres_imp3 = parameterspace_2[x,3], propres_imp4 = parameterspace_2[x,4], propres_imp5 = parameterspace_2[x,5], 
               propres_imp6 = parameterspace_2[x,6], propres_imp7 = parameterspace_2[x,7], propres_imp8 = parameterspace_2[x,8], propres_imp9 = parameterspace_2[x,9], propres_imp10 = parameterspace_2[x,10])
    
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau1[i]
    temp[,2] <- sum(out[nrow(out),6:29])
    temp[,3] <- signif(sum(out[nrow(out), seq(7,29,by = 2)])/ temp[,2], digits = 3)
    output1[i,] <- temp
  }
  
  colnames(output1) <- c("tau", "ICombH","IResRat")
  dump_data[1] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == 0.0123]))*100 #% Reduction
  dump_data[2] <- "Parameter Set 2"
  
  usage_frame_2 <- rbind(usage_frame_2, dump_data)
  
  print(paste0("Parameter Set 2 - ",round(x/10000, digits = 4)*100, "%"))
}
colnames(usage_frame_2) <- c("relRes", "Parameter")

# Parameter set 3 ---------------------------------------------------------

#Run the model

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)
times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,0.0123)

usage_frame_3 <- data.frame(matrix(nrow = 0, ncol = 2))

for(x in 1:nrow(parameterspace_1)) {
  dump_data <- vector() 
  output1 <- data.frame(matrix(nrow = 2, ncol =3))
  for (i in 1:2) {
    temp <- data.frame(matrix(nrow = 1, ncol =3))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = parmtau1[i],
               betaHD = (0.00001),  
               
               betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
               betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
               
               phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = 0.656, 
               
               fracimp1 = 0.076719577, fracimp2 = 0.046997389, fracimp3 = 0.022156129, fracimp4 = 0.116216976, fracimp5 = 0.083127572, fracimp6 = 0.057626057,
               fracimp7 = 0.04588015, fracimp8 = 0.010497067, fracimp9 = 0.05, fracimp10 = 0.2, 
               
               imp1 = parameterspace_3[x,1], imp2 = parameterspace_3[x,2], imp3 = parameterspace_3[x,3], imp4 = parameterspace_3[x,4], imp5 = parameterspace_3[x,5], imp6 = parameterspace_3[x,6],
               imp7 = parameterspace_3[x,7], imp8 = parameterspace_3[x,8], imp9 = parameterspace_3[x,9], imp10 = parameterspace_3[x,10],
               
               propres_imp1 = 0.347826087, propres_imp2 = 0.258064516, propres_imp3 = 0.266666667, propres_imp4 = 0.588235294, propres_imp5 = 0.674698795, 
               propres_imp6 = 0.54822335, propres_imp7 = 0.436893204, propres_imp8 = 0.4, propres_imp9 = 0.456, propres_imp10 = 0.5)
    
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau1[i]
    temp[,2] <- sum(out[nrow(out),6:29])
    temp[,3] <- signif(sum(out[nrow(out), seq(7,29,by = 2)])/ temp[,2], digits = 3)
    output1[i,] <- temp
  }
  
  colnames(output1) <- c("tau", "ICombH","IResRat")
  dump_data[1] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == 0.0123]))*100 #% Reduction
  dump_data[2] <- "Parameter Set 3"
  
  usage_frame_3 <- rbind(usage_frame_3, dump_data)
  
  print(paste0("Parameter Set 3 - ",round(x/10000, digits = 4)*100, "%"))
}
colnames(usage_frame_3) <- c("relRes", "Parameter")

# Parameter set 4 ---------------------------------------------------------

#Run the model

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)
times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,0.0123)

usage_frame_4 <- data.frame(matrix(nrow = 0, ncol = 2))

for(x in 1:nrow(parameterspace_1)) {
  dump_data <- vector() 
  output1 <- data.frame(matrix(nrow = 2, ncol =3))
  for (i in 1:2) {
    temp <- data.frame(matrix(nrow = 1, ncol =3))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = parmtau1[i],
               betaHD = (0.00001),  
               
               betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
               betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
               
               phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = parameterspace_4[x], 
               
               fracimp1 = 0.076719577, fracimp2 = 0.046997389, fracimp3 = 0.022156129, fracimp4 = 0.116216976, fracimp5 = 0.083127572, fracimp6 = 0.057626057,
               fracimp7 = 0.04588015, fracimp8 = 0.010497067, fracimp9 = 0.05, fracimp10 = 0.2, 
               
               imp1 = 0.1642031, imp2 = 0.1356062, imp3 = 0.1355515, imp4 = 0.134563, imp5 = 0.104059, imp6 = 0.08818593,
               imp7 = 0.07600404, imp8 = 0.05488995, imp9 = 0.03716985, imp10 = 0.06976744,
               
               propres_imp1 = 0.347826087, propres_imp2 = 0.258064516, propres_imp3 = 0.266666667, propres_imp4 = 0.588235294, propres_imp5 = 0.674698795, 
               propres_imp6 = 0.54822335, propres_imp7 = 0.436893204, propres_imp8 = 0.4, propres_imp9 = 0.456, propres_imp10 = 0.5)
    
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau1[i]
    temp[,2] <- sum(out[nrow(out),6:29])
    temp[,3] <- signif(sum(out[nrow(out), seq(7,29,by = 2)])/ temp[,2], digits = 3)
    output1[i,] <- temp
  }
  
  colnames(output1) <- c("tau", "ICombH","IResRat")
  dump_data[1] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == 0.0123]))*100 #% Reduction
  dump_data[2] <- "Parameter Set 4"
  
  usage_frame_4 <- rbind(usage_frame_4, dump_data)
  
  print(paste0("Parameter Set 4 - ", round(x/10000, digits = 4)*100, "%"))
}
colnames(usage_frame_4) <- c("relRes", "Parameter")

# Parameter set 5 ---------------------------------------------------------

#Run the model

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
          Sh = 1,
          IshDA = 0,IrhDA = 0,
          IshA1 = 0,IrhA1 = 0,
          IshA2 = 0,IrhA2 = 0,
          IshA3 = 0,IrhA3 = 0,
          IshA4 = 0,IrhA4 = 0,
          IshA5 = 0,IrhA5 = 0,
          IshA6 = 0,IrhA6 = 0,
          IshA7 = 0,IrhA7 = 0,
          IshA8 = 0,IrhA8 = 0,
          IshA9 = 0,IrhA9 = 0,
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)
times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,0.0123)

usage_frame_5 <- data.frame(matrix(nrow = 0, ncol = 2))

for(x in 1:nrow(parameterspace_1)) {
  dump_data <- vector() 
  output1 <- data.frame(matrix(nrow = 2, ncol =3))
  for (i in 1:2) {
    temp <- data.frame(matrix(nrow = 1, ncol =3))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = parmtau1[i],
               betaHD = (0.00001),  
               
               betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
               betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
               
               phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = parameterspace_comb[x,"psi"], 
               
               fracimp1 = parameterspace_comb[x,"fracimp1"], fracimp2 =  parameterspace_comb[x,"fracimp2"], fracimp3 =  parameterspace_comb[x,"fracimp3"], fracimp4 =  parameterspace_comb[x,"fracimp4"], 
               fracimp5 =  parameterspace_comb[x,"fracimp5"], fracimp6 =  parameterspace_comb[x,"fracimp6"], fracimp7 =  parameterspace_comb[x,"fracimp7"], fracimp8 =  parameterspace_comb[x,"fracimp8"], 
               fracimp9 =  parameterspace_comb[x,"fracimp9"], fracimp10 =  parameterspace_comb[x,"fracimp10"], 
               
               imp1 = parameterspace_comb[x,"imp1"], imp2 = parameterspace_comb[x,"imp2"], imp3 = parameterspace_comb[x,"imp3"], imp4 = parameterspace_comb[x,"imp4"], imp5 = parameterspace_comb[x,"imp5"], 
               imp6 = parameterspace_comb[x,"imp6"], imp7 = parameterspace_comb[x,"imp7"], imp8 = parameterspace_comb[x,"imp8"], imp9 = parameterspace_comb[x,"imp9"], imp10 = parameterspace_comb[x,"imp10"],
               
               propres_imp1 = parameterspace_comb[x,"propres_imp1"], propres_imp2 = parameterspace_comb[x,"propres_imp2"], propres_imp3 = parameterspace_comb[x,"propres_imp3"], 
               propres_imp4 = parameterspace_comb[x,"propres_imp4"], propres_imp5 = parameterspace_comb[x,"propres_imp5"], 
               propres_imp6 = parameterspace_comb[x,"propres_imp6"], propres_imp7 = parameterspace_comb[x,"propres_imp7"], propres_imp8 = parameterspace_comb[x,"propres_imp8"], 
               propres_imp9 = parameterspace_comb[x,"propres_imp9"], propres_imp10 = parameterspace_comb[x,"propres_imp10"])
    
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau1[i]
    temp[,2] <- sum(out[nrow(out),6:29])
    temp[,3] <- signif(sum(out[nrow(out), seq(7,29,by = 2)])/ temp[,2], digits = 3)
    output1[i,] <- temp
  }
  
  colnames(output1) <- c("tau", "ICombH","IResRat")
  dump_data[1] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == 0.0123]))*100 #% Reduction
  dump_data[2] <- "Parameter Set 5"
  
  usage_frame_5 <- rbind(usage_frame_5, dump_data)
  
  print(paste0("Parameter Set 5 - ", round(x/10000, digits = 4)*100, "%"))
}
colnames(usage_frame_5) <- c("relRes", "Parameter")

# Combination Violin Plot ------------------------------------------------

comb_data <- rbind(usage_frame_1,usage_frame_2, usage_frame_3, usage_frame_4, usage_frame_5)

colnames(comb_data) <- c("relRes", "Parameter")

typeof(comb_data$relRes)

comb_data$Parameter <- as.factor(comb_data$Parameter)

parm_viol_plot <- ggplot(comb_data, aes(x=Parameter, y=as.numeric(relRes), fill=Parameter)) + theme_bw() +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1, fill="white") +
  geom_hline(yintercept = 46.8, size = 1.2, col = "red", lty = 2) +
  labs(title="Effect of Import Parameters on Relative Change in Resistance", 
       x="Import Parameters", y = "Relative Change in Human Resistance Upon Domestic Curtailment") +
  theme_classic() + scale_fill_brewer(palette="Dark2") + 
  scale_x_discrete(labels=c("Parameter Set 1" = "Import Contamination", "Parameter Set 2" = "Proportion Imports Resistant", 
                            "Parameter Set 3" = "Distribution of Importation", "Parameter Set 4" = "Domestic Food Usage",
                            "Parameter Set 5" = "All Import Parameters")) + 
  theme(legend.position= "none", legend.text=element_text(size=12), legend.title =element_text(size=14), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) 

ggsave(parm_viol_plot, filename = "uncertainty_viol_plot.png", dpi = 300, type = "cairo", width = 12, height = 8, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")


#I also need to do a uncertainty analysis for the betaHA parameters