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
    
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(1-Sh) - uh*Sh - 
      
      betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshH) - 
      (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhH) - 
      
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
      
      (1-psi)*(imp_nEU)*(betaHI_nEU*fracimp_nEU*(1-propres_impnEU)*Sh) - 
      (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_nEU*fracimp_nEU*propres_impnEU*Sh)
    
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
    
    dIshAnEU = (1-psi)*(imp_nEU)*(betaHI_nEU*fracimp_nEU*(1-propres_impnEU)*Sh) - rh*IshAnEU - uh*IshAnEU 
    dIrhAnEU = (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_nEU*fracimp_nEU*propres_impnEU*Sh) - rh*IrhAnEU - uh*IrhAnEU  
    
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

# Import Posterior --------------------------------------------------------

amp_post <- do.call(rbind,
                     lapply(list.files(path = "C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data", pattern = "^complexmodel_ABC_SMC_gen_amp.*?\\.csv")[1:6], read.csv))
amp_post$gen <- as.vector(sapply(1:6, 
                                 function(x) rep(paste0("gen_",x), 1000)))

# Assess Posterior --------------------------------------------------------

l_amp_post <- lapply(1:length((colnames(amp_post))[-1]), function(x) melt(amp_post, id.vars = "gen", measure.vars = colnames(amp_post)[x])[,c(1,3)])

names = c(expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")),
          expression(paste("Rate of Antibiotic-Resistant to Antibiotic-Sensitive Reversion (", phi, ")")),
          expression(paste("Efficacy of Antibiotic-Mediated Animal Recovery (", kappa, ")")),
          expression(paste("Transmission-related Antibiotic Resistant Fitness Cost (", alpha, ")")),
          expression(paste("Background Infection Rate (", zeta, ")")),
          expression(paste("Rate of Domestic Animal-to-Human Transmission (", beta[HD], ")")),
          expression(paste("Rate of Domestic Human-to-Human Transmission (", beta[HH], ")")),
          expression(paste("Rate of EU Animal-to-Human Transmission (", beta[HIEU], ")")),
          expression(paste("Rate of non-EU Animal-to-Human Transmission (", beta[HInonEU], ")")),
          expression(paste("Proportion of Contaminated non-EU Imports Resistant (", imp[nonEU], ")")),
          expression(paste("Proportion of non-EU Imports Contaminated (", PropRes[nonEU], ")")))

amp_p_list <- list()

for(i in 1:length(names)){ 
  amp_p_list[[i]] <- ggplot(l_amp_post[[i]], aes(x=value, fill= gen)) + geom_density(alpha=.5) + 
    scale_x_continuous(expand = c(0, 0), name = names[i]) + 
    scale_y_continuous(expand = c(0, 0), name = "") +
    labs(fill = NULL) + scale_fill_discrete(labels = sapply(1:length(unique(amp_post$gen)), function(x) paste0("Generation ", x)))+
    theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
}

# Pairs Plot --------------------------------------------------------------

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
    scale_fill_gradientn(colours=viridis::viridis(100)) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  p
}

pairs_amp <- GGally::ggpairs(amp_post[-12], lower=list(continuous=my_fn)) + theme_bw()

# ABC-SMC Diagnostic Plots --------------------------------------------------------

diag_data <- readRDS("dist_amp_list.rds")

accept_gen <- data.frame("accept" = (1000/unlist(lapply(diag_data, "[[", 1)))*100, "gen" = seq(1, length(diag_data)))
diag_dataframe <- do.call("rbind",lapply(diag_data, "[[", 2)); diag_dataframe$gen <- as.vector(sapply(1:length(diag_data), function(x) rep(x,1000)))
colnames(diag_dataframe) <- c("dist", "foodH", "AMRH", "foodA", "AMRA", "gen")

diag_data_next <- as.data.frame(sapply(1:5 , function(j) sapply(1:length(diag_data), function(x) mean(diag_dataframe[diag_dataframe$gen == x, c("dist", "foodH", "AMRH", "foodA", "AMRA")[j]]))))
diag_data_next$gen <- seq(1,6); colnames(diag_data_next) <- c("dist", "foodH", "AMRH", "foodA", "AMRA", "gen")
diag_data_next[,2:5] <- as.numeric(mapply('/', diag_data_next[,2:5], c(3.26, 0.197, 0.017173052, 0.3333333))*100)

data <- melt(diag_data_next[,2:6], id.vars = "gen", measure.vars = c("foodH", "AMRH", "foodA", "AMRA"))

diag_data_next <- ggplot(diag_data_next, aes(x = gen, y = dist)) + geom_line(size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Generation", y = "Sum of Squared Difference Between Observed Data and Model") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  

p_diag <- ggplot(data, aes(x = gen, y = value, col = variable)) + geom_line(size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Generation", y = "% Distance from Target Outcome Value", col = "Outcome Measure") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  



# Baseline Change ---------------------------------------------------------

MAP_parms <- data.frame("parms" = colnames(amp_post)[-12], "mean" = colMeans(amp_post[amp_post$gen == tail(unique(amp_post$gen), 1),][-12]))

parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = 0, psi = 0.656,
          
          fracimp1 = country_data_imp[2,"Corrected_Usage_18"], fracimp2 = country_data_imp[3,"Corrected_Usage_18"], fracimp3 = country_data_imp[4,"Corrected_Usage_18"], fracimp4 = country_data_imp[5,"Corrected_Usage_18"], 
          fracimp5 = country_data_imp[6,"Corrected_Usage_18"], fracimp6 = country_data_imp[7,"Corrected_Usage_18"], fracimp7 = country_data_imp[8,"Corrected_Usage_18"], fracimp8 = country_data_imp[9,"Corrected_Usage_18"], 
          fracimp9 = country_data_imp[10,"Corrected_Usage_18"], fracimp_nEU = 1 - sum(country_data_imp[1:10,"Corrected_Usage_18"]),
          
          imp1 = country_data_imp[2,"Foodborne_Carriage_2019"], imp2 = country_data_imp[3,"Foodborne_Carriage_2019"], imp3 = country_data_imp[4,"Foodborne_Carriage_2019"], imp4 = country_data_imp[5,"Foodborne_Carriage_2019"], 
          imp5 = country_data_imp[6,"Foodborne_Carriage_2019"], imp6 = country_data_imp[7,"Foodborne_Carriage_2019"], imp7 = country_data_imp[8,"Foodborne_Carriage_2019"], imp8 = country_data_imp[9,"Foodborne_Carriage_2019"],
          imp8 = country_data_imp[9,"Foodborne_Carriage_2019"], imp9 = country_data_imp[10,"Foodborne_Carriage_2019"],
          
          propres_imp1 = country_data_imp[2,"Prop_Amp_Res"], propres_imp2 = country_data_imp[3,"Prop_Amp_Res"], propres_imp3 = country_data_imp[4,"Prop_Amp_Res"], propres_imp4 = country_data_imp[5,"Prop_Amp_Res"], 
          propres_imp5 = country_data_imp[6,"Prop_Amp_Res"], propres_imp6 = country_data_imp[7,"Prop_Amp_Res"], propres_imp7 = country_data_imp[8,"Prop_Amp_Res"], propres_imp8 = country_data_imp[9,"Prop_Amp_Res"], 
          propres_imp9 = country_data_imp[10,"Prop_Amp_Res"],
          
          betaAA = MAP_parms[which(MAP_parms == "betaAA"),2], betaHH = MAP_parms[which(MAP_parms == "betaHH"),2], betaHD =  MAP_parms[which(MAP_parms == "betaHD"),2],
          betaHI_EU =  MAP_parms[which(MAP_parms == "betaHI_EU"),2], betaHI_nEU = MAP_parms[which(MAP_parms == "betaHI_nEU"),2], 
          
          imp_nEU =  MAP_parms[which(MAP_parms == "imp_nEU"),2], propres_impnEU =  MAP_parms[which(MAP_parms == "propres_impnEU"),2],
          phi =  MAP_parms[which(MAP_parms == "phi"),2], kappa =  MAP_parms[which(MAP_parms == "kappa"),2], alpha =  MAP_parms[which(MAP_parms == "alpha"),2], 
          zeta =  MAP_parms[which(MAP_parms == "zeta"),2])

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

times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,UK_amp)

output_base <-  data.frame(matrix(nrow = 1, ncol =3))

for (i in 1:2) {
  temp <- data.frame(matrix(nrow = 1, ncol =3))
  parms2["tau"] <- parmtau1[i]
  
  out <- ode(y = init, func = amrimp, times = times, parms = parms2)
  temp[,1] <- parmtau1[i]
  temp[,2] <- sum(out[nrow(out),6:29])
  temp[,3] <- signif(sum(out[nrow(out), seq(7,29,by = 2)])/ temp[,2], digits = 3)
  output_base[i,] <- temp
}
colnames(output_base) <- c("tau", "ICombH","IResRat")

base_change <- (1-(output_base$IResRat[output_base$tau == 0] / output_base$IResRat[output_base$tau == UK_amp]))*100

# Uncertainty Analysis ----------------------------------------------------

#THis involves the creation of a violin plot - we want to explore combinations of fitted parameters - which have uncertainty 
#BetaAA, phi, kappa, alpha, zeta, betaHD, betaHI_EU, beta_HI,nEU, imp_neu and propres_impneu
#Use the combo parameter to create a combination of: 
# - Just resistance 
# - Just contamination
# - all import parameters 
# - all of the above 
# Then create a violin plots with the means plotted 
# Have a red line showing the baseline change in resistance 

#Baseline Parameters

MAP_parms <- data.frame("parms" = colnames(amp_post)[-12], "mean" = colMeans(amp_post[amp_post$gen == tail(unique(amp_post$gen), 1),][-12]))

parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = 0, psi = 0.656,
          
          fracimp1 = country_data_imp[2,"Corrected_Usage_18"], fracimp2 = country_data_imp[3,"Corrected_Usage_18"], fracimp3 = country_data_imp[4,"Corrected_Usage_18"], fracimp4 = country_data_imp[5,"Corrected_Usage_18"], 
          fracimp5 = country_data_imp[6,"Corrected_Usage_18"], fracimp6 = country_data_imp[7,"Corrected_Usage_18"], fracimp7 = country_data_imp[8,"Corrected_Usage_18"], fracimp8 = country_data_imp[9,"Corrected_Usage_18"], 
          fracimp9 = country_data_imp[10,"Corrected_Usage_18"], fracimp_nEU = 1 - sum(country_data_imp[1:10,"Corrected_Usage_18"]),
          
          imp1 = country_data_imp[2,"Foodborne_Carriage_2019"], imp2 = country_data_imp[3,"Foodborne_Carriage_2019"], imp3 = country_data_imp[4,"Foodborne_Carriage_2019"], imp4 = country_data_imp[5,"Foodborne_Carriage_2019"], 
          imp5 = country_data_imp[6,"Foodborne_Carriage_2019"], imp6 = country_data_imp[7,"Foodborne_Carriage_2019"], imp7 = country_data_imp[8,"Foodborne_Carriage_2019"], imp8 = country_data_imp[9,"Foodborne_Carriage_2019"],
          imp8 = country_data_imp[9,"Foodborne_Carriage_2019"], imp9 = country_data_imp[10,"Foodborne_Carriage_2019"],
          
          propres_imp1 = country_data_imp[2,"Prop_Amp_Res"], propres_imp2 = country_data_imp[3,"Prop_Amp_Res"], propres_imp3 = country_data_imp[4,"Prop_Amp_Res"], propres_imp4 = country_data_imp[5,"Prop_Amp_Res"], 
          propres_imp5 = country_data_imp[6,"Prop_Amp_Res"], propres_imp6 = country_data_imp[7,"Prop_Amp_Res"], propres_imp7 = country_data_imp[8,"Prop_Amp_Res"], propres_imp8 = country_data_imp[9,"Prop_Amp_Res"], 
          propres_imp9 = country_data_imp[10,"Prop_Amp_Res"],
          
          betaAA = MAP_parms[which(MAP_parms == "betaAA"),2], betaHH = MAP_parms[which(MAP_parms == "betaHH"),2], betaHD =  MAP_parms[which(MAP_parms == "betaHD"),2],
          betaHI_EU =  MAP_parms[which(MAP_parms == "betaHI_EU"),2], betaHI_nEU = MAP_parms[which(MAP_parms == "betaHI_nEU"),2], 
          
          imp_nEU =  MAP_parms[which(MAP_parms == "imp_nEU"),2], propres_impnEU =  MAP_parms[which(MAP_parms == "propres_impnEU"),2],
          phi =  MAP_parms[which(MAP_parms == "phi"),2], kappa =  MAP_parms[which(MAP_parms == "kappa"),2], alpha =  MAP_parms[which(MAP_parms == "alpha"),2], 
          zeta =  MAP_parms[which(MAP_parms == "zeta"),2])

#Ranges for Ampicillin Fit 

#Using the Posterior for ranges
#last_gen_amp <- lapply(1:(ncol(amp_post)-1), function(x) range(amp_post[colnames(amp_post)[x]]))
#amp_runif <- lapply(1:(ncol(amp_post)-1), function(x) runif(1000, last_gen_amp[[x]][[1]], last_gen_amp[[x]][[2]]))
#names(amp_runif) <- colnames(amp_post)[-12]

last_gen_amp <- lapply(1:(ncol(amp_post)-1), function(x) range(amp_post[colnames(amp_post)[x]]))
amp_runif <- lapply(1:(ncol(amp_post)-1), function(x) runif(10000, last_gen_amp[[x]][[1]], last_gen_amp[[x]][[2]]))
names(amp_runif) <- colnames(amp_post)[-12]


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

times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,UK_amp)

amp_uncert <- list()

for(j in 1:length(names(amp_runif))){
  amp_uncert[[j]] <- local ({
    
    usage_frame_1 <- data.frame(matrix(nrow = 0, ncol = 2))
    
    for(x in 1:length(amp_runif[[j]])) {
      dump_data <- vector() 
      parms2 <- parms
      parms2[names(amp_runif)[j]] <- amp_runif[[j]][[x]]
      output1 <- data.frame(matrix(nrow = 2, ncol =3))
      
      for (i in 1:2) {
        temp <- data.frame(matrix(nrow = 1, ncol =3))
        parms2["tau"] <- parmtau1[i]
        
        out <- ode(y = init, func = amrimp, times = times, parms = parms2)
        temp[,1] <- parmtau1[i]
        temp[,2] <- sum(out[nrow(out),6:29])
        temp[,3] <- signif(sum(out[nrow(out), seq(7,29,by = 2)])/ temp[,2], digits = 3)
        output1[i,] <- temp
      }
      print(output1)
      colnames(output1) <- c("tau", "ICombH","IResRat")
      dump_data[1] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == UK_amp]))*100 #% Reduction
      dump_data[2] <- paste0("Parm_Set_", names(amp_runif)[j])
      usage_frame_1 <- rbind(usage_frame_1, dump_data)
      
      print(paste0("Parameter Set ", names(amp_runif)[j]," - ", round(x/length(amp_runif[[j]]), digits = 4)*100, "%"))
    }
    colnames(usage_frame_1) <- c("relRes", "Parameter")
    return(usage_frame_1)
  })
}

# Combination Violin Plot ------------------------------------------------

comb_data <- do.call("rbind", amp_uncert)
comb_data$relRes <- as.numeric(comb_data$relRes)

typeof(comb_data$relRes)

comb_data$Parameter <- as.factor(comb_data$Parameter)

parm_viol_plot <- ggplot(comb_data, aes(x=Parameter, y=as.numeric(relRes), fill=Parameter)) + theme_bw() +
  geom_violin(trim=TRUE) + geom_boxplot(width=0.1, fill="white") +
  geom_hline(yintercept = base_change, size = 1.2, col = "red", lty = 2) +
  labs(title="Effect of Import Parameters on Relative Change in Human Resistance Upon Domestic Curtailment", 
       x="Fitted Parameters", y = "Relative Change in Human Resistance") +
  theme_classic() + scale_fill_brewer(palette="Spectral") + 
  scale_x_discrete(labels=c(names(amp_runif))) + 
  theme(legend.position= "none", legend.text=element_text(size=12), legend.title =element_text(size=14), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) 

ggsave(parm_viol_plot, filename = "AMP_uncertainty_viol_plot.png", dpi = 300, type = "cairo", width = 12, height = 8, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")

# Acceptable Import Analysis ----------------------------------------------

MAP_parms <- data.frame("parms" = colnames(amp_post)[-12], "mean" = colMeans(amp_post[amp_post$gen == tail(unique(amp_post$gen), 1),][-12]))

parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = 0, psi = 0.656,
          
          fracimp1 = country_data_imp[2,"Corrected_Usage_18"], fracimp2 = country_data_imp[3,"Corrected_Usage_18"], fracimp3 = country_data_imp[4,"Corrected_Usage_18"], fracimp4 = country_data_imp[5,"Corrected_Usage_18"], 
          fracimp5 = country_data_imp[6,"Corrected_Usage_18"], fracimp6 = country_data_imp[7,"Corrected_Usage_18"], fracimp7 = country_data_imp[8,"Corrected_Usage_18"], fracimp8 = country_data_imp[9,"Corrected_Usage_18"], 
          fracimp9 = country_data_imp[10,"Corrected_Usage_18"], fracimp_nEU = 1 - sum(country_data_imp[1:10,"Corrected_Usage_18"]),
          
          imp1 = country_data_imp[2,"Foodborne_Carriage_2019"], imp2 = country_data_imp[3,"Foodborne_Carriage_2019"], imp3 = country_data_imp[4,"Foodborne_Carriage_2019"], imp4 = country_data_imp[5,"Foodborne_Carriage_2019"], 
          imp5 = country_data_imp[6,"Foodborne_Carriage_2019"], imp6 = country_data_imp[7,"Foodborne_Carriage_2019"], imp7 = country_data_imp[8,"Foodborne_Carriage_2019"], imp8 = country_data_imp[9,"Foodborne_Carriage_2019"],
          imp8 = country_data_imp[9,"Foodborne_Carriage_2019"], imp9 = country_data_imp[10,"Foodborne_Carriage_2019"],
         
           propres_imp1 = country_data_imp[2,"Prop_Amp_Res"], propres_imp2 = country_data_imp[3,"Prop_Amp_Res"], propres_imp3 = country_data_imp[4,"Prop_Amp_Res"], propres_imp4 = country_data_imp[5,"Prop_Amp_Res"], 
          propres_imp5 = country_data_imp[6,"Prop_Amp_Res"], propres_imp6 = country_data_imp[7,"Prop_Amp_Res"], propres_imp7 = country_data_imp[8,"Prop_Amp_Res"], propres_imp8 = country_data_imp[9,"Prop_Amp_Res"], 
          propres_imp9 = country_data_imp[10,"Prop_Amp_Res"],
          
          
          betaAA = MAP_parms[which(MAP_parms == "betaAA"),2], betaHH = MAP_parms[which(MAP_parms == "betaHH"),2], betaHD =  MAP_parms[which(MAP_parms == "betaHD"),2],
          betaHI_EU =  MAP_parms[which(MAP_parms == "betaHI_EU"),2], betaHI_nEU = MAP_parms[which(MAP_parms == "betaHI_nEU"),2], 
          
          imp_nEU =  MAP_parms[which(MAP_parms == "imp_nEU"),2], propres_impnEU =  MAP_parms[which(MAP_parms == "propres_impnEU"),2],
          phi =  MAP_parms[which(MAP_parms == "phi"),2], kappa =  MAP_parms[which(MAP_parms == "kappa"),2], alpha =  MAP_parms[which(MAP_parms == "alpha"),2], 
          zeta =  MAP_parms[which(MAP_parms == "zeta"),2])

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
times <- seq(0, 10000, by = 100) 

#Generate 5 Parm Sets - Baseline, Usage amongst Imports,  Contamination, Resistance, All 

parms_base <- parms

amp_cont_data <- data.frame(matrix(runif(1000, 0 , 1), nrow = 1000, ncol = 10)); colnames(amp_cont_data) <- grep("fracimp", names(parms_base), value=TRUE)
amp_res_data <- data.frame(matrix(runif(1000, 0 , 1), nrow = 1000, ncol = 10)); colnames(amp_res_data) <- grep("propres_imp", names(parms_base), value=TRUE)
t_data <- RandVec(a=0, b=1, s=1, n=10, m=1000, Seed=sample(1:1000, size = 1))
imp_usage_data <- data.frame(matrix(unlist(t_data), nrow=1000, byrow=TRUE), stringsAsFactors=FALSE); colnames(imp_usage_data) <- c(sapply(1:9, function(x) paste0("imp",x)), "imp_nEU")
comb_data <- cbind(amp_cont_data, amp_res_data, imp_usage_data) 

data_list <- list(amp_cont_data, amp_res_data, imp_usage_data, comb_data)

usage_threshold <- c(seq(0, 0.656, by = 0.02), 0.656)

#Function for the Import %

import_res_func <- function(parms, init, usage_threshold, uk_usage) {
  parmtau1 <- c(0, uk_usage)
  usage_frame <- data.frame(matrix(nrow = length(usage_threshold), ncol = 2))
  for(x in 1:length(usage_threshold)) {
    parms[["psi"]] <- usage_threshold[x]
    output1 <- data.frame(matrix(nrow = 2, ncol =3))
    dump_data <- vector()
    
    for (i in 1:2) {
      
      temp <- vector()
      parms["tau"] <- parmtau1[i]
      out <- ode(y = init, func = amrimp, times = times, parms = parms)
      
      temp[1] <- parmtau1[i]
      temp[2] <- sum(out[nrow(out),6:29])
      temp[3] <- sum(out[nrow(out), seq(7,29,by = 2)])/ temp[2]
      output1[i,] <- temp
    }
    colnames(output1) <- c("tau", "ICombH","IResRat")
    
    dump_data[1] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == uk_usage]))*100 #% Reductio
    dump_data[2] <- usage_threshold[x]
    usage_frame[x,] <- dump_data
  }
  
  colnames(usage_frame) <- c("relchange", "domusage")
  usage_frame$relchange <- 100 - (usage_frame$relchange / usage_frame$relchange[usage_frame$domusage == 0.656]) * 100
  
  return(usage_frame)
}

#Run the Scenario Analysis

data_imp_list <- list()

for(j in 1:5){
  data_imp_list[[j]] <- local ({
    
    if(j == 1) {
      output <- import_res_func(parms_base, init, usage_threshold, UK_amp)
      print(paste0(c("base", "cont", "res", "usage", "all")[j]))
      
    } else {
      output <- list()
      parms_list <- data_list[[j-1]]
      parms <- parms_base
      
      for(i in 1:nrow(parms_list)) {
        parms[colnames(parms_list)] <- parms_list[i,]
        out <- import_res_func(parms, init, usage_threshold, UK_amp)
        output[[i]] <- out 
        print(paste0(c("base", "cont", "res", "usage", "all")[j], " - ", i/nrow(parms_list)))
      }
    }
    return(output)
  })
}

#Plotting Output 

data_base <- data_imp_list[[1]]
cont_minmax_data <- as.data.frame(cbind(sapply(data_imp_list[[2]], "[[", 1)))
res_minmax_data <- as.data.frame(cbind(sapply(data_imp_list[[3]], "[[", 1)))
usage_minmax_data <- as.data.frame(cbind(sapply(data_imp_list[[4]], "[[", 1)))
all_minmax_data <- as.data.frame(cbind(sapply(data_imp_list[[5]], "[[", 1)))

data_base$min_cont <- apply(cont_minmax_data, 1, FUN = min); data_base$max_cont <- apply(cont_minmax_data, 1, FUN = max)
data_base$min_res <- apply(res_minmax_data, 1, FUN = min); data_base$max_res <- apply(res_minmax_data, 1, FUN = max)
data_base$min_usage <- apply(usage_minmax_data, 1, FUN = min); data_base$max_usage <- apply(usage_minmax_data, 1, FUN = max)
data_base$min_all <- apply(all_minmax_data, 1, FUN = min); data_base$max_all <- apply(all_minmax_data, 1, FUN = max)

p_cont <- ggplot(data_base, aes(domusage, relchange)) + geom_ribbon(aes(ymin = min_cont, ymax = max_cont), fill = "lightblue", alpha = 0.7, col = "lightblue") + 
  geom_line(col = "blue" , size = 1.15) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Proportion of Food From Domestic Origins", y = "Relative Reduction in the Efficacy of Livestock Antibiotic Curtailment on Resistance") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))

p_res <- ggplot(data_base, aes(domusage, relchange)) + geom_ribbon(aes(ymin = min_res, ymax = max_res), fill = "lightblue", alpha = 0.7, col = "lightblue") + 
  geom_line(col = "blue" , size = 1.15) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Proportion of Food From Domestic Origins", y = "Relative Reduction in the Efficacy of Livestock Antibiotic Curtailment on Resistance") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  

p_usage <- ggplot(data_base, aes(domusage, relchange)) + geom_ribbon(aes(ymin = min_usage, ymax = max_usage), fill = "lightblue", alpha = 0.7, col = "lightblue") + 
  geom_line(col = "blue" , size = 1.15) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Proportion of Food From Domestic Origins", y = "Relative Reduction in the Efficacy of Livestock Antibiotic Curtailment on Resistance") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  

p_all <- ggplot(data_base, aes(domusage, relchange)) + geom_ribbon(aes(ymin = min_all, ymax = max_all), fill = "lightblue", alpha = 0.7, col = "lightblue") + 
  geom_line(col = "blue" , size = 1.15) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Proportion of Food From Domestic Origins", y = "Relative Reduction in the Efficacy of Livestock Antibiotic Curtailment on Resistance") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  

ggsave(p_cont, filename = "import_scen_cont.png", dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
ggsave(p_res, filename = "import_scen_res.png", dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
ggsave(p_usage, filename = "import_scen_usage.png", dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
ggsave(p_all, filename = "import_scen_all.png", dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")

ggplot(data_base, aes(domusage, relchange))  + 
  geom_line(col = "blue" , size = 1.15) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Proportion of Food From Domestic Origins", y = "Relative Reduction in the Efficacy of Livestock Antibiotic Curtailment on Resistance") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  