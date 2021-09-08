library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity"); library("metR")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate"); library("viridis")

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
    
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(1-Sh) - uh*Sh - 
      
      betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshA4+IshA5+IshA6+IshA7+IshA8+IshA9+IshAnEU+IshH) - 
      (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhA4+IrhA5+IrhA6+IrhA7+IrhA8+IrhA9+IrhAnEU+IrhH) - 
      
      psi*(betaHD*Isa*Sh) - 
      psi*(1-alpha)*(betaHD*Ira*Sh) - 
      
      (1-psi)*(imp1)*(betaHI_EU*fracimp1*(1-propres_imp1)*Sh) - 
      (1-psi)*(imp1)*(1-alpha)*(betaHI_EU*fracimp1*propres_imp1*Sh)-
      
      (1-psi)*(imp2)*(betaHI_EU*fracimp2*(1-propres_imp2)*Sh) - 
      (1-psi)*(imp2)*(1-alpha)*(betaHI_EU*fracimp2*propres_imp2*Sh) - 
      
      (1-psi)*(imp3)*(betaHI_EU*fracimp3*(1-propres_imp3)*Sh) - 
      (1-psi)*(imp3)*(1-alpha)*(betaHI_EU*fracimp3*propres_imp3*Sh) -
      
      (1-psi)*(imp4)*(betaHI_EU*fracimp4*(1-propres_imp4)*Sh) - 
      (1-psi)*(imp4)*(1-alpha)*(betaHI_EU*fracimp4*propres_imp4*Sh)-
      
      (1-psi)*(imp5)*(betaHI_EU*fracimp5*(1-propres_imp5)*Sh) - 
      (1-psi)*(imp5)*(1-alpha)*(betaHI_EU*fracimp5*propres_imp5*Sh) - 
      
      (1-psi)*(imp6)*(betaHI_EU*fracimp6*(1-propres_imp6)*Sh) - 
      (1-psi)*(imp6)*(1-alpha)*(betaHI_EU*fracimp6*propres_imp6*Sh) - 
      
      (1-psi)*(imp7)*(betaHI_EU*fracimp7*(1-propres_imp7)*Sh) - 
      (1-psi)*(imp7)*(1-alpha)*(betaHI_EU*fracimp7*propres_imp7*Sh) - 
      
      (1-psi)*(imp8)*(betaHI_EU*fracimp8*(1-propres_imp8)*Sh) - 
      (1-psi)*(imp8)*(1-alpha)*(betaHI_EU*fracimp8*propres_imp8*Sh) -
      
      (1-psi)*(imp9)*(betaHI_EU*fracimp9*(1-propres_imp9)*Sh) - 
      (1-psi)*(imp9)*(1-alpha)*(betaHI_EU*fracimp9*propres_imp9*Sh) - 
      
      (1-psi)*(imp_nEU)*(betaHI_EU*fracimp_nEU*(1-propres_impnEU)*Sh) - 
      (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_EU*fracimp_nEU*propres_impnEU*Sh)
    
    dIshDA = psi*betaHD*Isa*Sh  - rh*IshDA - uh*IshDA 
    dIrhDA = psi*(1-alpha)*betaHD*Ira*Sh - rh*IrhDA - uh*IrhDA  
    
    dIshA1 = (1-psi)*(imp1)*(betaHI_EU*fracimp1*(1-propres_imp1)*Sh) - rh*IshA1 - uh*IshA1 
    dIrhA1 = (1-psi)*(imp1)*(1-alpha)*(betaHI_EU*fracimp1*propres_imp1*Sh) - rh*IrhA1 - uh*IrhA1  
    
    dIshA2 = (1-psi)*(imp2)*(betaHI_EU*fracimp2*(1-propres_imp2)*Sh) - rh*IshA2 - uh*IshA2 
    dIrhA2 = (1-psi)*(imp2)*(1-alpha)*(betaHI_EU*fracimp2*propres_imp2*Sh)  - rh*IrhA2 - uh*IrhA2  
    
    dIshA3 = (1-psi)*(imp3)*(betaHI_EU*fracimp3*(1-propres_imp3)*Sh) - rh*IshA3 - uh*IshA3 
    dIrhA3 = (1-psi)*(imp3)*(1-alpha)*(betaHI_EU*fracimp3*propres_imp3*Sh) - rh*IrhA3 - uh*IrhA3  
    
    dIshA4 = (1-psi)*(imp4)*(betaHI_EU*fracimp4*(1-propres_imp4)*Sh) - rh*IshA4 - uh*IshA4 
    dIrhA4 = (1-psi)*(imp4)*(1-alpha)*(betaHI_EU*fracimp1*propres_imp4*Sh) - rh*IrhA4 - uh*IrhA4  
    
    dIshA5 = (1-psi)*(imp5)*(betaHI_EU*fracimp5*(1-propres_imp5)*Sh) - rh*IshA5 - uh*IshA5 
    dIrhA5 = (1-psi)*(imp5)*(1-alpha)*(betaHI_EU*fracimp5*propres_imp5*Sh)  - rh*IrhA5 - uh*IrhA5  
    
    dIshA6 = (1-psi)*(imp6)*(betaHI_EU*fracimp6*(1-propres_imp6)*Sh) - rh*IshA6 - uh*IshA6 
    dIrhA6 = (1-psi)*(imp6)*(1-alpha)*(betaHI_EU*fracimp6*propres_imp6*Sh) - rh*IrhA6 - uh*IrhA6  
    
    dIshA7 = (1-psi)*(imp7)*(betaHI_EU*fracimp7*(1-propres_imp7)*Sh) - rh*IshA7 - uh*IshA7 
    dIrhA7 = (1-psi)*(imp7)*(1-alpha)*(betaHI_EU*fracimp7*propres_imp7*Sh) - rh*IrhA7 - uh*IrhA7  
    
    dIshA8 = (1-psi)*(imp8)*(betaHI_EU*fracimp8*(1-propres_imp8)*Sh) - rh*IshA8 - uh*IshA8 
    dIrhA8 = (1-psi)*(imp8)*(1-alpha)*(betaHI_EU*fracimp8*propres_imp8*Sh)  - rh*IrhA8 - uh*IrhA8  
    
    dIshA9 = (1-psi)*(imp9)*(betaHI_EU*fracimp9*(1-propres_imp9)*Sh) - rh*IshA9 - uh*IshA9 
    dIrhA9 = (1-psi)*(imp9)*(1-alpha)*(betaHI_EU*fracimp9*propres_imp9*Sh) - rh*IrhA9 - uh*IrhA9  
    
    dIshAnEU = (1-psi)*(imp_nEU)*(betaHI_EU*fracimp_nEU*(1-propres_impnEU)*Sh) - rh*IshAnEU - uh*IshAnEU 
    dIrhAnEU = (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_EU*fracimp_nEU*propres_impnEU*Sh) - rh*IrhAnEU - uh*IrhAnEU  
    
    dIshH = betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshA4+IshA5+IshA6+IshA7+IshA8+IshA9+IshAnEU+IshH) - rh*IshH - uh*IshH 
    dIrhH = (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhA4+IrhA5+IrhA6+IrhA7+IrhA8+IrhA9+IrhAnEU+IrhH)- rh*IrhH - uh*IrhH 
    
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
                  dIshAnEU,dIrhAnEU,
                  dIshH, dIrhH)))
  })
}

# Data Import -------------------------------------------------------------

country_data_imp <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Data/FullData_2021_v1_trim.csv") #This is data for pigs 
country_data_imp$Foodborne_Carriage_2019 <- country_data_imp$Foodborne_Carriage_2019/100
country_data_imp$Corrected_Usage_18 <- country_data_imp$Corrected_Usage_18/100
country_data_imp[,12:13] <- country_data_imp[,12:13]/1000

country_data_gen <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Data/res_sales_generalfit.csv") #This is data for pigs 
country_data_gen[,13:14] <- country_data_gen[,13:14]/1000

UK_amp <- country_data_gen$scaled_sales_amp[country_data_gen$Country == "United Kingdom"]

country_data_gen <- country_data_gen[country_data_gen$num_test_amp >= 10,]

plot(country_data_gen$scaled_sales_amp, country_data_gen$propres_amp, ylim = c(0,1))

# Import in Parameters and Set Baseline Parms -----------------------------

amp_post <- read.csv(list.files(path = "C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data", pattern = "^complexmodel_ABC_SMC_gen_amp.*?\\.csv")[3])

# Parameterisation  -------------------------------------------------------

#Generate the Means or the MAPs for each parameter
MAP_parms <- data.frame("parms" = colnames(amp_post), "mean" = colMeans(amp_post))
#MAP_parms <- map_estimate(amp_post)

#Initial Conditions
parmtau <- seq(0, 0.01, by = 0.001)

init = c(Sa=0.98, Isa=0.01, Ira=0.01, 
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
         IshAnEU = 0,IrhAnEU = 0,
         IshH = 0, IrhH = 0)

output1 <- data.frame(matrix(ncol = 28, nrow = 0))
times <- seq(0, 10000, by = 1)

parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = parmtau[i], psi = 0.656,
          
          fracimp1 = as.numeric(country_data_imp[2,"Normalised_Usage_2018"]), fracimp2 = as.numeric(country_data_imp[3,"Normalised_Usage_2018"]), fracimp3 = as.numeric(country_data_imp[4,"Normalised_Usage_2018"]), 
          fracimp4 = as.numeric(country_data_imp[5,"Normalised_Usage_2018"]), fracimp5 = as.numeric(country_data_imp[6,"Normalised_Usage_2018"]), fracimp6 = as.numeric(country_data_imp[7,"Normalised_Usage_2018"]), 
          fracimp7 = as.numeric(country_data_imp[8,"Normalised_Usage_2018"]), fracimp8 = as.numeric(country_data_imp[9,"Normalised_Usage_2018"]), 
          fracimp9 = as.numeric(country_data_imp[10,"Normalised_Usage_2018"]), fracimp_nEU = 1 - sum(as.numeric(country_data_imp[2:10,"Normalised_Usage_2018"])),
          
          imp1 = country_data_imp[2,"Foodborne_Carriage_2019"], imp2 = country_data_imp[3,"Foodborne_Carriage_2019"], imp3 = country_data_imp[4,"Foodborne_Carriage_2019"], imp4 = country_data_imp[5,"Foodborne_Carriage_2019"], 
          imp5 = country_data_imp[6,"Foodborne_Carriage_2019"], imp6 = country_data_imp[7,"Foodborne_Carriage_2019"], imp7 = country_data_imp[8,"Foodborne_Carriage_2019"], imp8 = country_data_imp[9,"Foodborne_Carriage_2019"],
          imp8 = country_data_imp[9,"Foodborne_Carriage_2019"], imp9 = country_data_imp[10,"Foodborne_Carriage_2019"],
          
          propres_imp1 = country_data_imp[2,"Prop_Amp_Res"], propres_imp2 = country_data_imp[3,"Prop_Amp_Res"], propres_imp3 = country_data_imp[4,"Prop_Amp_Res"], propres_imp4 = country_data_imp[5,"Prop_Amp_Res"], 
          propres_imp5 = country_data_imp[6,"Prop_Amp_Res"], propres_imp6 = country_data_imp[7,"Prop_Amp_Res"], propres_imp7 = country_data_imp[8,"Prop_Amp_Res"], propres_imp8 = country_data_imp[9,"Prop_Amp_Res"], 
          propres_imp9 = country_data_imp[10,"Prop_Amp_Res"],
          
          betaAA = MAP_parms[which(MAP_parms == "betaAA"),2], betaHH = MAP_parms[which(MAP_parms == "betaHH"),2], betaHD =  MAP_parms[which(MAP_parms == "betaHD"),2],
          betaHI_EU =  MAP_parms[which(MAP_parms == "betaHI_EU"),2],
          
          imp_nEU =  MAP_parms[which(MAP_parms == "imp_nEU"),2], propres_impnEU =  MAP_parms[which(MAP_parms == "propres_impnEU"),2],
          phi =  MAP_parms[which(MAP_parms == "phi"),2], kappa =  MAP_parms[which(MAP_parms == "kappa"),2], alpha =  MAP_parms[which(MAP_parms == "alpha"),2], 
          zeta =  MAP_parms[which(MAP_parms == "zeta"),2])

#Above is the old parm set based on the MAPs 
#Bottom is a random particle that was accepted in the 4th generation - use the above when we get there.

plot_analysis <- list()

for (j in 1:2) {
  
  parmtau <- list(country_data_gen$scaled_sales_amp,
                  seq(0, 0.014, by = 0.001))[[j]]
  
  plot_analysis[[j]] <- local({
    
    for (i in 1:length(parmtau)) {
      
      temp <- data.frame(matrix(NA, nrow = 1, ncol=28))
      
      parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = parmtau[i], psi = 0.656,
                
                fracimp1 = as.numeric(country_data_imp[2,"Normalised_Usage_2018"]), fracimp2 = as.numeric(country_data_imp[3,"Normalised_Usage_2018"]), fracimp3 = as.numeric(country_data_imp[4,"Normalised_Usage_2018"]), 
                fracimp4 = as.numeric(country_data_imp[5,"Normalised_Usage_2018"]), fracimp5 = as.numeric(country_data_imp[6,"Normalised_Usage_2018"]), fracimp6 = as.numeric(country_data_imp[7,"Normalised_Usage_2018"]), 
                fracimp7 = as.numeric(country_data_imp[8,"Normalised_Usage_2018"]), fracimp8 = as.numeric(country_data_imp[9,"Normalised_Usage_2018"]), 
                fracimp9 = as.numeric(country_data_imp[10,"Normalised_Usage_2018"]), fracimp_nEU = 1 - sum(as.numeric(country_data_imp[2:10,"Normalised_Usage_2018"])),
                
                imp1 = country_data_imp[2,"Foodborne_Carriage_2019"], imp2 = country_data_imp[3,"Foodborne_Carriage_2019"], imp3 = country_data_imp[4,"Foodborne_Carriage_2019"], imp4 = country_data_imp[5,"Foodborne_Carriage_2019"], 
                imp5 = country_data_imp[6,"Foodborne_Carriage_2019"], imp6 = country_data_imp[7,"Foodborne_Carriage_2019"], imp7 = country_data_imp[8,"Foodborne_Carriage_2019"], imp8 = country_data_imp[9,"Foodborne_Carriage_2019"],
                imp8 = country_data_imp[9,"Foodborne_Carriage_2019"], imp9 = country_data_imp[10,"Foodborne_Carriage_2019"],
                
                propres_imp1 = country_data_imp[2,"Prop_Amp_Res"], propres_imp2 = country_data_imp[3,"Prop_Amp_Res"], propres_imp3 = country_data_imp[4,"Prop_Amp_Res"], propres_imp4 = country_data_imp[5,"Prop_Amp_Res"], 
                propres_imp5 = country_data_imp[6,"Prop_Amp_Res"], propres_imp6 = country_data_imp[7,"Prop_Amp_Res"], propres_imp7 = country_data_imp[8,"Prop_Amp_Res"], propres_imp8 = country_data_imp[9,"Prop_Amp_Res"], 
                propres_imp9 = country_data_imp[10,"Prop_Amp_Res"],
                
                betaAA = 1.887899e-02, betaHH = 2.882880e-02, betaHD =   9.660318e-05 ,
                betaHI_EU =  9.473187e-05,
                
                imp_nEU =   6.665432e-01, propres_impnEU =  2.032954e-02 ,
                phi =  2.882880e-02, kappa =  4.492843e+00 , alpha =  2.862765e-01, 
                zeta =  9.247680e-05 )
      
      out <<- ode(y = init, func = amrimp, times = times, parms = parms)
      
      temp[1,1] <- parmtau[i]
      temp[1,2:13] <- out[nrow(out),seq(7, 29, by = 2)]/sum(out[nrow(out),seq(6, 29)]) 
      temp[1,14:25] <- sapply(seq(6,28, by = 2), function(x) out[nrow(out),x]) + sapply(seq(7,29, by = 2), function(x) out[nrow(out),x]) 
      temp[1,26] <- sum(out[nrow(out),seq(7, 29, by = 2)])/sum(out[nrow(out),seq(6, 29)]) 
      temp[1,27] <- out[nrow(out), 4] / sum(out[nrow(out),3:4]) 
      temp[1,28] <- sum(out[nrow(out),3:4]) 
      
      print(paste0("Run ", j, ": ",temp[1,2]))
      output1 <- rbind.data.frame(output1, temp)
    }
    
    colnames(output1)[1:28] <- c("tau", 
                                 
                                 "PropRes_Domestic","PropRes_Netherlands","PropRes_Ireland","PropRes_Germany","PropRes_France","PropRes_Spain","PropRes_Italy",
                                 "PropRes_Belgium","PropRes_Poland","PropRes_Denmark","PropRes_NonEU", "PropRes_Human",
                                 
                                 "TotInf_Domestic","TotInf_Netherlands","TotInf_Ireland","TotInf_Germany","TotInf_France","TotInf_Spain","TotInf_Italy",
                                 "TotInf_Belgium","TotInf_Poland","TotInf_Denmark","TotInf_NonEU", "TotInf_Human", 
                                 
                                 "Res_PropHum", "Res_PropAnim", "Anim_Inf")
    return(output1)
  })
}

resistance_data <- plot_analysis[[2]]
norm_names <- sapply(1:12, function(x) paste0(substring(colnames(resistance_data)[14:25], 8)[x], "_res_norm"))
resistance_data[c(norm_names)] <- resistance_data[14:25]/rowSums(resistance_data[,14:25])

#Non-normalized Resistance 
res_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[2:13]) 
trim_names <- as.factor(sapply(strsplit(as.character(res_plotdata$variable), split = "_", fixed = TRUE), function(x) x[2]))
res_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Normalized Resistance 
res_norm_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[29:40])
res_norm_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Infection
inf_plotdata <- melt(plot_analysis[[2]], id.vars = c("tau"), measure.vars = colnames(plot_analysis[[1]])[14:25]) 
inf_plotdata$variable <- factor(trim_names, levels = unique(trim_names))


# Main Attribution Plots --------------------------------------------------

#Un-normalised resistance
res_comb <- ggplot(res_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  +
  labs(x ="Domestic Ampicillin Usage in Fattening Pig (g/PCU)", y = "Proportion of Human Salmnellosis Resistant to Ampicillin", fill = "Resistance Source")  

#Normalised resistance
res_norm_comb <- ggplot(res_norm_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  +
  labs(x ="Domestic Ampicillin Usage in Fattening Pig (g/PCU)", y = "Proportion of Attributable Human Resistance", fill = "Resistance Source")  


inf_comb <- ggplot(inf_plotdata, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Total Human Salmenollosis (per 100,000)", fill = "Resistance Source") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_vline(xintercept = UK_amp, size = 1.2, col = "red", lty = 2)

# Fit to Data Checks -------------------------------------------------------------

plot_check <- data.frame("tau" = plot_analysis[[1]]$tau, "country" = country_data_gen$Country, 
                         "model_estim" = plot_analysis[[1]]$Res_PropAnim, "observ" = country_data_gen$propres_tet, 
                         "model_estim_inf" = plot_analysis[[1]]$Anim_Inf)

country_data_gen <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Data/res_sales_generalfit.csv") #This is data for pigs 
country_data_gen[,13:14] <- country_data_gen[,13:14]/1000
country_data_gen <- country_data_gen[country_data_gen$num_test_amp > 10,]

#Animal Resistance vs Animal Livestock Antibiotic Usage
anim_plot <- ggplot(plot_check, mapping = aes(x = tau)) + geom_point(aes(y = observ), size = 2) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Domestic Livestock Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp, y = country_data_gen$propres_amp[country_data_gen$Country == "United Kingdom"]) , col = "red", size = 5)

#Animal Infection vs Animal Livestock Antibiotic Usage

anim_plot_inf <- ggplot(plot_check, mapping = aes(x = tau)) + geom_line(aes(y = model_estim_inf), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.01)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Domestic Livestock Contaminated") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp, y = country_data_imp$Foodborne_Carriage_2019[country_data_imp$Country_of_Origin == "UK Origin"]) , col = "red", size = 5)

#Animal Resistance vs Animal Livestock Antibiotic Usage
plot_check_hum <- plot_check; plot_check_hum$model_estim <- plot_analysis[[1]]$Res_PropHum

hum_plot <- ggplot(plot_check_hum, mapping = aes(x = tau)) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Human Infections Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp, y = 0.3108) , col = "red", size = 5)

ggplot(plot_check_hum, mapping = aes(x = tau)) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.015)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Human Infections Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + geom_vline(xintercept = UK_amp, col = "red", size = 1.2, lty = 2)

