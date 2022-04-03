library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity"); library("parallel")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate"); library("rootSolve")

rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")

# Single Model ------------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      0.5*zeta*Sa*(1-alpha) - 0.5*zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + 0.5*zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + 0.5*zeta*Sa*(1-alpha)
    
    
    dSh = uh + rh*(1-Sh) - uh*Sh - 
      
      psi*betaHA*(Isa*eta)*Sh - 
      psi*(1-alpha)*betaHA*(Ira*eta)*Sh - 
      
      (1-psi)*(imp1)*(betaHA*fracimp1*(1-propres_imp1)*Sh) - 
      (1-psi)*(imp1)*(1-alpha)*(betaHA*fracimp1*propres_imp1*Sh)-
      
      (1-psi)*(imp2)*(betaHA*fracimp2*(1-propres_imp2)*Sh) - 
      (1-psi)*(imp2)*(1-alpha)*(betaHA*fracimp2*propres_imp2*Sh) - 
      
      (1-psi)*(imp3)*(betaHA*fracimp3*(1-propres_imp3)*Sh) - 
      (1-psi)*(imp3)*(1-alpha)*(betaHA*fracimp3*propres_imp3*Sh) -
      
      (1-psi)*(imp4)*(betaHA*fracimp4*(1-propres_imp4)*Sh) - 
      (1-psi)*(imp4)*(1-alpha)*(betaHA*fracimp4*propres_imp4*Sh)-
      
      (1-psi)*(imp5)*(betaHA*fracimp5*(1-propres_imp5)*Sh) - 
      (1-psi)*(imp5)*(1-alpha)*(betaHA*fracimp5*propres_imp5*Sh) - 
      
      (1-psi)*(imp6)*(betaHA*fracimp6*(1-propres_imp6)*Sh) - 
      (1-psi)*(imp6)*(1-alpha)*(betaHA*fracimp6*propres_imp6*Sh) - 
      
      (1-psi)*(imp7)*(betaHA*fracimp7*(1-propres_imp7)*Sh) - 
      (1-psi)*(imp7)*(1-alpha)*(betaHA*fracimp7*propres_imp7*Sh) - 
      
      (1-psi)*(imp8)*(betaHA*fracimp8*(1-propres_imp8)*Sh) - 
      (1-psi)*(imp8)*(1-alpha)*(betaHA*fracimp8*propres_imp8*Sh) -
      
      (1-psi)*(imp9)*(betaHA*fracimp9*(1-propres_imp9)*Sh) - 
      (1-psi)*(imp9)*(1-alpha)*(betaHA*fracimp9*propres_imp9*Sh) - 
      
      (1-psi)*(imp_nEU)*(betaHA*fracimp_nEU*(1-propres_impnEU)*Sh) - 
      (1-psi)*(imp_nEU)*(1-alpha)*(betaHA*fracimp_nEU*propres_impnEU*Sh)
    
    dIshDA = psi*betaHA*(Isa*eta)*Sh - rh*IshDA - uh*IshDA 
    dIrhDA = psi*(1-alpha)*betaHA*(Ira*eta)*Sh - rh*IrhDA - uh*IrhDA  
    
    dIshA1 = (1-psi)*(imp1)*(betaHA*fracimp1*(1-propres_imp1)*Sh) - rh*IshA1 - uh*IshA1 
    dIrhA1 = (1-psi)*(imp1)*(1-alpha)*(betaHA*fracimp1*propres_imp1*Sh) - rh*IrhA1 - uh*IrhA1  
    
    dIshA2 = (1-psi)*(imp2)*(betaHA*fracimp2*(1-propres_imp2)*Sh) - rh*IshA2 - uh*IshA2 
    dIrhA2 = (1-psi)*(imp2)*(1-alpha)*(betaHA*fracimp2*propres_imp2*Sh) - rh*IrhA2 - uh*IrhA2  
    
    dIshA3 = (1-psi)*(imp3)*(betaHA*fracimp3*(1-propres_imp3)*Sh) - rh*IshA3 - uh*IshA3 
    dIrhA3 = (1-psi)*(imp3)*(1-alpha)*(betaHA*fracimp3*propres_imp3*Sh) - rh*IrhA3 - uh*IrhA3  
    
    dIshA4 = (1-psi)*(imp4)*(betaHA*fracimp4*(1-propres_imp4)*Sh) - rh*IshA4 - uh*IshA4 
    dIrhA4 = (1-psi)*(imp4)*(1-alpha)*(betaHA*fracimp4*propres_imp4*Sh) - rh*IrhA4 - uh*IrhA4  
    
    dIshA5 = (1-psi)*(imp5)*(betaHA*fracimp5*(1-propres_imp5)*Sh) - rh*IshA5 - uh*IshA5 
    dIrhA5 = (1-psi)*(imp5)*(1-alpha)*(betaHA*fracimp5*propres_imp5*Sh) - rh*IrhA5 - uh*IrhA5  
    
    dIshA6 = (1-psi)*(imp6)*(betaHA*fracimp6*(1-propres_imp6)*Sh) - rh*IshA6 - uh*IshA6 
    dIrhA6 = (1-psi)*(imp6)*(1-alpha)*(betaHA*fracimp6*propres_imp6*Sh) - rh*IrhA6 - uh*IrhA6  
    
    dIshA7 = (1-psi)*(imp7)*(betaHA*fracimp7*(1-propres_imp7)*Sh) - rh*IshA7 - uh*IshA7 
    dIrhA7 = (1-psi)*(imp7)*(1-alpha)*(betaHA*fracimp7*propres_imp7*Sh) - rh*IrhA7 - uh*IrhA7  
    
    dIshA8 = (1-psi)*(imp8)*(betaHA*fracimp8*(1-propres_imp8)*Sh) - rh*IshA8 - uh*IshA8 
    dIrhA8 = (1-psi)*(imp8)*(1-alpha)*(betaHA*fracimp8*propres_imp8*Sh) - rh*IrhA8 - uh*IrhA8  
    
    dIshA9 = (1-psi)*(imp9)*(betaHA*fracimp9*(1-propres_imp9)*Sh) - rh*IshA9 - uh*IshA9 
    dIrhA9 = (1-psi)*(imp9)*(1-alpha)*(betaHA*fracimp9*propres_imp9*Sh) - rh*IrhA9 - uh*IrhA9  
    
    dIshAnEU = (1-psi)*(imp_nEU)*(betaHA*fracimp_nEU*(1-propres_impnEU)*Sh) - rh*IshAnEU - uh*IshAnEU 
    dIrhAnEU = (1-psi)*(imp_nEU)*(1-alpha)*(betaHA*fracimp_nEU*propres_impnEU*Sh) - rh*IrhAnEU - uh*IrhAnEU  
    
    CumS = psi*betaHA*(Isa*eta)*Sh + 
      (1-psi)*(imp1)*(betaHA*fracimp1*(1-propres_imp1)*Sh) + 
      (1-psi)*(imp2)*(betaHA*fracimp2*(1-propres_imp2)*Sh) +
      (1-psi)*(imp3)*(betaHA*fracimp3*(1-propres_imp3)*Sh) + 
      (1-psi)*(imp4)*(betaHA*fracimp4*(1-propres_imp4)*Sh) + 
      (1-psi)*(imp5)*(betaHA*fracimp5*(1-propres_imp5)*Sh) + 
      (1-psi)*(imp6)*(betaHA*fracimp6*(1-propres_imp6)*Sh) + 
      (1-psi)*(imp7)*(betaHA*fracimp7*(1-propres_imp7)*Sh) + 
      (1-psi)*(imp8)*(betaHA*fracimp8*(1-propres_imp8)*Sh) + 
      (1-psi)*(imp9)*(betaHA*fracimp9*(1-propres_imp9)*Sh) + 
      (1-psi)*(imp_nEU)*(betaHA*fracimp_nEU*(1-propres_impnEU)*Sh)
    
    CumR = psi*(1-alpha)*betaHA*(Ira*eta)*Sh + 
      (1-psi)*(imp1)*(1-alpha)*(betaHA*fracimp1*propres_imp1*Sh) + 
      (1-psi)*(imp2)*(1-alpha)*(betaHA*fracimp2*propres_imp2*Sh) +
      (1-psi)*(imp3)*(1-alpha)*(betaHA*fracimp3*propres_imp3*Sh) + 
      (1-psi)*(imp4)*(1-alpha)*(betaHA*fracimp4*propres_imp4*Sh) + 
      (1-psi)*(imp5)*(1-alpha)*(betaHA*fracimp5*propres_imp5*Sh) + 
      (1-psi)*(imp6)*(1-alpha)*(betaHA*fracimp6*propres_imp6*Sh) +
      (1-psi)*(imp7)*(1-alpha)*(betaHA*fracimp7*propres_imp7*Sh) + 
      (1-psi)*(imp8)*(1-alpha)*(betaHA*fracimp8*propres_imp8*Sh) + 
      (1-psi)*(imp9)*(1-alpha)*(betaHA*fracimp9*propres_imp9*Sh) + 
      (1-psi)*(imp_nEU)*(1-alpha)*(betaHA*fracimp_nEU*propres_impnEU*Sh)
    
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
                  dIshAnEU,dIrhAnEU),
                CumS, CumR))
  })
}


# Livestock Dynamics Dataset ----------------------------------------------

#Import Data
dataamp_pigs_raw <- read.csv("Amp_FatPigs_Comb.csv"); dataamp_pigs <- dataamp_pigs_raw
dataamp_hum_raw <- read.csv("Hum_FatPigs.csv"); dataamp_hum <- dataamp_hum_raw

#Cleaning Data - Animals
dataamp_pigs[,(2+5):(6+5)][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the No. of pos isolates
dataamp_pigs[,(2+10):(6+10)][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for the prop of resistant isolates
dataamp_pigs[,2:6][dataamp_pigs[,2:6] < 10] <- NA #If N > 10, replace the particular country/year with NA for N
dataamp_pigs <- dataamp_pigs[!(is.na(dataamp_pigs$N_2015) & is.na(dataamp_pigs$N_2016) & is.na(dataamp_pigs$N_2017) & 
                                 is.na(dataamp_pigs$N_2018) & is.na(dataamp_pigs$N_2019)),]
pig_yrs <- sub("N_", "", grep("N_20",colnames(dataamp_pigs), value = TRUE)) #Find years of the EFSA and ESVAC data in the dataset
colnames(dataamp_pigs)[12:16] <- pig_yrs

#Create dataset where each row is a different observation. 
melt_amp_pigs <- melt(dataamp_pigs, id.vars = "Country", measure.vars = pig_yrs)
melt_amp_pigs$usage <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("scale_ampusage_2015", "scale_ampusage_2016", 
                                                                                "scale_ampusage_2017", "scale_ampusage_2018", "scale_ampusage_2019"))[,3]
melt_amp_pigs$N <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("N_2015", "N_2016", 
                                                                            "N_2017", "N_2018", "N_2019"))[,3]
melt_amp_pigs$IsolPos <- melt(dataamp_pigs, id.vars = "Country", measure.vars = c("PosIsol_2015", "PosIsol_2016", 
                                                                                  "PosIsol_2017", "PosIsol_2018", "PosIsol_2019"))[,3]
colnames(melt_amp_pigs)[c(2,3)] <- c("Year", "Resistance")

#Cleaning Data - Humans
#only include countries/years which are present in the resistance dataset
dataamp_hum <- dataamp_hum[dataamp_hum$Country %in% intersect(dataamp_hum$Country, dataamp_pigs$Country),]
colnames(dataamp_hum)[26:31] <- as.character(2014:2019)
melt_amp_pigs$ResPropHum <- melt(dataamp_hum, id.vars = "Country", measure.vars = pig_yrs)[,3]
melt_amp_pigs <- melt_amp_pigs[!(is.na(melt_amp_pigs$Resistance) | is.na(melt_amp_pigs$usage)),] # Remove all rows with NAs for usage and resistance

#Add 95% CIs for each datapoint
melt_amp_pigs$lower_amp <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[1]]))
melt_amp_pigs$upper_amp <- unlist(lapply(1:nrow(melt_amp_pigs), function(i) prop.test(melt_amp_pigs$IsolPos[i],melt_amp_pigs$N[i])[[6]][[2]]))

#Rename the columns
colnames(melt_amp_pigs) <- c("Country", "Year", "ResPropAnim", "Usage", "N", "IsolPos", "ResPropHum", "Lower_Amp", "Upper_Amp")
melt_amp_pigs$Usage <- melt_amp_pigs$Usage/1000 #Change from mg/PCU to g/PCU

ggplot(melt_amp_pigs, aes(x = Usage, y= ResPropAnim, color = Country)) + geom_point() +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.055)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

# Food Usage Dataset ------------------------------------------------------

#country_data_imp <- read.csv("FullData_2021_v1_trim.csv") #This is data for pigs 
country_data_imp <- read.csv("ImportDat_AmpPigs_update.csv") #This is data for pigs 
country_data_imp[country_data_imp$Country_of_Origin == "UK Origin",23] <- NA
isolamp_hum_raw <- read.csv("UK_parameterisation.csv")

UK_hum_ampres <- rowMeans(isolamp_hum_raw[,25:28], na.rm = T)[2]
UK_amp_res <- as.numeric(rowMeans(isolamp_hum_raw[isolamp_hum_raw$Country_of_Origin == "UK Livestock",25:28], na.rm = T )) 
UK_amp_usage <- as.numeric(rowMeans(isolamp_hum_raw[isolamp_hum_raw$Country_of_Origin == "UK Livestock",29:32]))/1000
UK_cont <- as.numeric(isolamp_hum_raw[isolamp_hum_raw$Country_of_Origin == "UK Livestock",24])
UK_food_usage <- isolamp_hum_raw[isolamp_hum_raw$Country_of_Origin == "UK Livestock",2]

UK_food_usage <- isolamp_hum_raw[isolamp_hum_raw$Country_of_Origin == "UK Livestock",2]
UK_food_pig_usage <- isolamp_hum_raw[isolamp_hum_raw$Country_of_Origin == "UK Livestock",4]

#Use the mean for the EU as the parameters (minus the UK) - only the main importers 

EU_cont <- mean(rowMeans(country_data_imp[,24:27], na.rm = T))
EU_res <- mean(rowMeans(country_data_imp[28:31], na.rm = T))

max(country_data_imp[24:27], na.rm = T)

country_data_imp[1, 25] <- NA

country_data_imp$FBD_gen <- rowMeans(country_data_imp[,24:27], na.rm = T)
country_data_imp$FBD_res <- rowMeans(country_data_imp[,28:31], na.rm = T)

# Import in Parameters and Set Baseline Parms -----------------------------

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha")

post_amp <- read.csv(tail(list.files(path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha", pattern = "complex"), 1))
MAP_parms <- map_estimate(post_amp)
MAP_parms <- data.frame("Parameter" = names(post_amp), 
                        "MAP_Estimate" = colMeans(post_amp))

# Test Posteriors ---------------------------------------------------------

amp_post <- do.call(rbind,
                    lapply(list.files(path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha", pattern = "complex"), read.csv))

amp_post$gen <- as.vector(sapply(1:(nrow(amp_post)/nrow(read.csv(tail(list.files(path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha", pattern = "complex"), 1)))), 
                                 function(x) rep(paste0("gen_",x), nrow(read.csv(tail(list.files(path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha", pattern = "complex"), 1))))))

l_amp_post <- lapply(1:length(colnames(amp_post)[-1]), function(x) melt(amp_post, id.vars = "gen", measure.vars = colnames(amp_post)[x])[,c(1,3)])

names = c(expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")),
          expression(paste("Rate of Antibiotic-Resistant to Antibiotic-Sensitive Reversion (", phi, ")")),
          expression(paste("Efficacy of Antibiotic-Mediated Animal Recovery (", kappa, ")")),
          expression(paste("Transmission-related Antibiotic Resistant Fitness Cost (", alpha, ")")),
          expression(paste("Background Infection Rate (", zeta, ")")),
          expression(paste("Rate of Domestic Animal-to-Human Transmission (", beta[HA], ")")),
          expression(paste("Proportion of non-EU Imports Contaminated (", Imp[nonEU], ")")),
          expression(paste("Proportion of Contaminated non-EU Imports Resistant (", PropRes[nonEU], ")")))

amp_p_list <- list()

for(i in 1:length(names)){ 
  amp_p_list[[i]] <- ggplot(l_amp_post[[i]], aes(x=value, fill= gen)) + geom_density(alpha=.5) + 
    scale_x_continuous(expand = c(0, 0), name = names[i]) + 
    scale_y_continuous(expand = c(0, 0), name = "") +
    labs(fill = NULL) + scale_fill_discrete(labels = sapply(1:length(unique(amp_post$gen)), function(x) paste0("Generation ", x)))+
    theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
}

ggarrange(amp_p_list[[1]], amp_p_list[[2]], amp_p_list[[3]],
          amp_p_list[[4]], amp_p_list[[5]], amp_p_list[[6]],
          amp_p_list[[7]], amp_p_list[[8]], nrow = 3, ncol = 3,
          common.legend = TRUE, legend = "bottom")

# Parameterisation  -------------------------------------------------------

#Generate the Means or the MAPs for each parameter
#MAP_parms <- map_estimate(amp_post[amp_post$gen == "gen_6",][, -ncol(amp_post)])

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
         IshAnEU = 0,IrhAnEU = 0)

output1 <- data.frame(matrix(ncol = 28, nrow = 0))

thetaparm = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = 0.656,
              
              fracimp1 = country_data_imp[2,"Normalised_Usage_2018"], fracimp2 = country_data_imp[3,"Normalised_Usage_2018"], fracimp3 = country_data_imp[4,"Normalised_Usage_2018"], 
              fracimp4 = country_data_imp[5,"Normalised_Usage_2018"], fracimp5 = country_data_imp[6,"Normalised_Usage_2018"], fracimp6 = country_data_imp[7,"Normalised_Usage_2018"], 
              fracimp7 = country_data_imp[8,"Normalised_Usage_2018"], fracimp8 = country_data_imp[9,"Normalised_Usage_2018"], 
              fracimp9 = country_data_imp[10,"Normalised_Usage_2018"], fracimp_nEU = 1 - sum(country_data_imp[2:10,"Normalised_Usage_2018"]),
              
              betaAA = MAP_parms["betaAA", 2], phi = MAP_parms["phi", 2], kappa = MAP_parms["kappa", 2], alpha = MAP_parms["alpha", 2], 
              zeta = MAP_parms["zeta", 2], betaHA = MAP_parms["betaHA", 2], imp_nEU = MAP_parms["imp_nEU", 2], propres_impnEU = MAP_parms["propres_impnEU", 2], 
              
              imp1 = country_data_imp[2,"FBD_gen"], imp2 = country_data_imp[3,"FBD_gen"], imp3 = country_data_imp[4,"FBD_gen"], imp4 = country_data_imp[5,"FBD_gen"], 
              imp5 = country_data_imp[6,"FBD_gen"], imp6 = country_data_imp[7,"FBD_gen"], imp7 = country_data_imp[8,"FBD_gen"], imp8 = country_data_imp[9,"FBD_gen"],
              imp8 = country_data_imp[9,"FBD_gen"], imp9 = country_data_imp[10,"FBD_gen"],
              
              propres_imp1 = country_data_imp[2,"FBD_res"], propres_imp2 = country_data_imp[3,"FBD_res"], propres_imp3 = country_data_imp[4,"FBD_res"], propres_imp4 = country_data_imp[5,"FBD_res"], 
              propres_imp5 = country_data_imp[6,"FBD_res"], propres_imp6 = country_data_imp[7,"FBD_res"], propres_imp7 = country_data_imp[8,"FBD_res"], propres_imp8 = country_data_imp[9,"FBD_res"], 
              propres_imp9 = country_data_imp[10,"FBD_res"],
              
              eta = 0.11016, tau = UK_amp_usage)

plot_analysis <- list()

for (j in 1:2) {
  
  parmtau <- list(melt_amp_pigs$Usage,
                  seq(0, 0.014, by = 0.001))[[j]]
  
  plot_analysis[[j]] <- local({
    
    for (i in 1:length(parmtau)) {
      
      temp <- data.frame(matrix(NA, nrow = 1, ncol=37))
      
      parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = 0.656, tau = parmtau[i],
                
                fracimp1 = country_data_imp[2,"Normalised_Usage_2018"], fracimp2 = country_data_imp[3,"Normalised_Usage_2018"], fracimp3 = country_data_imp[4,"Normalised_Usage_2018"], 
                fracimp4 = country_data_imp[5,"Normalised_Usage_2018"], fracimp5 = country_data_imp[6,"Normalised_Usage_2018"], fracimp6 = country_data_imp[7,"Normalised_Usage_2018"], 
                fracimp7 = country_data_imp[8,"Normalised_Usage_2018"], fracimp8 = country_data_imp[9,"Normalised_Usage_2018"], 
                fracimp9 = country_data_imp[10,"Normalised_Usage_2018"], fracimp_nEU = 1 - sum(country_data_imp[2:10,"Normalised_Usage_2018"]),
                
                betaAA = MAP_parms["betaAA", 2], phi = MAP_parms["phi", 2], kappa = MAP_parms["kappa", 2], alpha = MAP_parms["alpha", 2], 
                zeta = MAP_parms["zeta", 2], betaHA = MAP_parms["betaHA", 2], imp_nEU = MAP_parms["imp_nEU", 2], propres_impnEU = MAP_parms["propres_impnEU", 2], 
                
                imp1 = country_data_imp[2,"FBD_gen"], imp2 = country_data_imp[3,"FBD_gen"], imp3 = country_data_imp[4,"FBD_gen"], imp4 = country_data_imp[5,"FBD_gen"], 
                imp5 = country_data_imp[6,"FBD_gen"], imp6 = country_data_imp[7,"FBD_gen"], imp7 = country_data_imp[8,"FBD_gen"], imp8 = country_data_imp[9,"FBD_gen"],
                imp9 = country_data_imp[10,"FBD_gen"],
                
                propres_imp1 = country_data_imp[2,"FBD_res"], propres_imp2 = country_data_imp[3,"FBD_res"], propres_imp3 = country_data_imp[4,"FBD_res"], propres_imp4 = country_data_imp[5,"FBD_res"], 
                propres_imp5 = country_data_imp[6,"FBD_res"], propres_imp6 = country_data_imp[7,"FBD_res"], propres_imp7 = country_data_imp[8,"FBD_res"], propres_imp8 = country_data_imp[9,"FBD_res"], 
                propres_imp9 = country_data_imp[10,"FBD_res"],
                
                eta = 0.11016)
      
      out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms)
      
      temp[1,1] <- parmtau[i]
      temp[1,2:12] <- out[[1]][seq(6, 26, by = 2)]/ sum(out[[1]][5:26])
      temp[1,13:23] <- sapply(seq(5,25, by = 2), function(x) out[[1]][x]) + sapply(seq(6,26, by = 2), function(x) out[[1]][x]) 

      temp[1,24:34] <- ((sapply(seq(5,25, by = 2), function(x) out[[1]][x]) + sapply(seq(6,26, by = 2), function(x) out[[1]][x])) * 446000000)/100000
      
      temp[1,35] <- sum(out[[1]][seq(6, 26, by = 2)])/ sum(out[[1]][5:26])
      temp[1,36] <- out[[1]][3] / sum(out[[1]][2:3])
      temp[1,37] <- sum(out[[1]][2:3])
      
      print(paste0("Run ", j, ": ",temp[1,2]))
      output1 <- rbind.data.frame(output1, temp)
    }
    
    colnames(output1)[1:37] <- c("tau", 
                                 
                                 "PropRes_Domestic","PropRes_Netherlands","PropRes_Ireland","PropRes_Germany","PropRes_France","PropRes_Spain","PropRes_Italy",
                                 "PropRes_Belgium","PropRes_Poland","PropRes_Denmark","PropRes_NonEU", 
                                 
                                 "TotInf_Domestic","TotInf_Netherlands","TotInf_Ireland","TotInf_Germany","TotInf_France","TotInf_Spain","TotInf_Italy",
                                 "TotInf_Belgium","TotInf_Poland","TotInf_Denmark","TotInf_NonEU", 
                                 
                                 "Inc_Domestic","Inc_Netherlands","Inc_Ireland","Inc_Germany","Inc_France","Inc_Spain","Inc_Italy",
                                 "Inc_Belgium","Inc_Poland","Inc_Denmark","Inc_NonEU", 
                                 
                                 "Res_PropHum", "Res_PropAnim", "Anim_Inf")
    return(output1)
  })
}

resistance_data <- plot_analysis[[2]]
norm_names <- sapply(1:11, function(x) paste0(substring(colnames(resistance_data)[13:23], 8)[x], "_res_norm"))
resistance_data[c(norm_names)] <- (resistance_data[2:12]*(sum(resistance_data[13:23])))/
  rowSums((resistance_data[2:12]*(sum(resistance_data[13:23]))))

#Non-normalized Resistance 
res_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[2:12]) 
trim_names <- as.factor(sapply(strsplit(as.character(res_plotdata$variable), split = "_", fixed = TRUE), function(x) x[2]))
res_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Normalized Resistance 
res_norm_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[38:48])
res_norm_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Infection
inf_plotdata <- melt(plot_analysis[[2]], id.vars = c("tau"), measure.vars = colnames(plot_analysis[[1]])[24:34]) 
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

inf_comb <- ggplot(inf_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0 , 0.75)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Total Human Salmenollosis (per 100,000)", fill = "Resistance Source") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_vline(xintercept = UK_amp_usage, size = 1.2, col = "red", lty = 2)

# Fit to Data Checks -------------------------------------------------------------

plot_check <- data.frame("tau" = plot_analysis[[1]]$tau, "country" = melt_amp_pigs$Country, 
                         "model_estim" = plot_analysis[[1]]$Res_PropAnim, "observ" = melt_amp_pigs$ResPropAnim, 
                         "model_estim_inf" = plot_analysis[[1]]$Anim_Inf)

#Include to get UK data
country_data_gen1 <-  read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/UK_parameterisation.csv") 

#Animal Resistance vs Animal Livestock Antibiotic Usage

anim_plot <- ggplot(plot_check, mapping = aes(x = tau)) + geom_point(aes(y = observ), size = 2) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Domestic Livestock Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp_usage, y = UK_amp_res) , col = "red", size = 5, inherit.aes = FALSE)

#Animal Infection vs Animal Livestock Antibiotic Usage

anim_plot_inf <- ggplot(plot_check, mapping = aes(x = tau)) + geom_line(aes(y = model_estim_inf*thetaparm[["eta"]]), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Domestic Livestock Contaminated") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp_usage, y = UK_cont) , col = "red", size = 5, inherit.aes = FALSE)

#Animal Resistance vs Animal Livestock Antibiotic Usage
plot_check_hum <- plot_check; plot_check_hum$model_estim <- plot_analysis[[1]]$Res_PropHum

hum_plot <- ggplot(plot_check_hum, mapping = aes(x = tau)) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Human Infections Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp_usage, y = UK_hum_ampres) , col = "red", size = 5)

ggplot(plot_check_hum, mapping = aes(x = tau)) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Human Infections Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + geom_vline(xintercept = UK_amp_usage, col = "red", size = 1.2, lty = 2)


# Psi Compare -------------------------------------------------------------

#Generate the Means or the MAPs for each parameter
#MAP_parms <- map_estimate(amp_post[amp_post$gen == "gen_6",][, -ncol(amp_post)])

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
         IshAnEU = 0,IrhAnEU = 0)

output1 <- data.frame(matrix(ncol = 28, nrow = 0))

plot_analysis <- list()

for (j in 1:2) {
  
  parmtau <- seq(0, 0.014, by = 0.001)
  
  plot_analysis[[j]] <- local({
    
    parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = c(UK_food_usage, UK_food_pig_usage)[j], 
              
              fracimp1 = country_data_imp[2,"Normalised_Usage_2018"], fracimp2 = country_data_imp[3,"Normalised_Usage_2018"], fracimp3 = country_data_imp[4,"Normalised_Usage_2018"], 
              fracimp4 = country_data_imp[5,"Normalised_Usage_2018"], fracimp5 = country_data_imp[6,"Normalised_Usage_2018"], fracimp6 = country_data_imp[7,"Normalised_Usage_2018"], 
              fracimp7 = country_data_imp[8,"Normalised_Usage_2018"], fracimp8 = country_data_imp[9,"Normalised_Usage_2018"], 
              fracimp9 = country_data_imp[10,"Normalised_Usage_2018"], fracimp_nEU = 1 - sum(country_data_imp[2:10,"Normalised_Usage_2018"]),
              
              betaAA = MAP_parms["betaAA", 2], phi = MAP_parms["phi", 2], kappa = MAP_parms["kappa", 2], alpha = MAP_parms["alpha", 2], 
              zeta = MAP_parms["zeta", 2], betaHA = MAP_parms["betaHA", 2], imp_nEU = MAP_parms["imp_nEU", 2], propres_impnEU = MAP_parms["propres_impnEU", 2], 
              
              imp1 = country_data_imp[2,"FBD_gen"], imp2 = country_data_imp[3,"FBD_gen"], imp3 = country_data_imp[4,"FBD_gen"], imp4 = country_data_imp[5,"FBD_gen"], 
              imp5 = country_data_imp[6,"FBD_gen"], imp6 = country_data_imp[7,"FBD_gen"], imp7 = country_data_imp[8,"FBD_gen"], imp8 = country_data_imp[9,"FBD_gen"],
              imp8 = country_data_imp[9,"FBD_gen"], imp9 = country_data_imp[10,"FBD_gen"],
              
              propres_imp1 = country_data_imp[2,"FBD_res"], propres_imp2 = country_data_imp[3,"FBD_res"], propres_imp3 = country_data_imp[4,"FBD_res"], propres_imp4 = country_data_imp[5,"FBD_res"], 
              propres_imp5 = country_data_imp[6,"FBD_res"], propres_imp6 = country_data_imp[7,"FBD_res"], propres_imp7 = country_data_imp[8,"FBD_res"], propres_imp8 = country_data_imp[9,"FBD_res"], 
              propres_imp9 = country_data_imp[10,"FBD_res"],
              
              eta = 0.11016)
    
    for (i in 1:length(parmtau)) {
      parms["tau"] = parmtau[i]
      temp <- data.frame(matrix(NA, nrow = 1, ncol=37))
      out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms)
      
      temp[1,1] <- parmtau[i]
      temp[1,2:12] <- out[[1]][seq(6, 26, by = 2)]/ sum(out[[1]][5:26])
      temp[1,13:23] <- sapply(seq(5,25, by = 2), function(x) out[[1]][x]) + sapply(seq(6,26, by = 2), function(x) out[[1]][x]) 
      
      temp[1,24:34] <- ((sapply(seq(5,25, by = 2), function(x) out[[1]][x]) + sapply(seq(6,26, by = 2), function(x) out[[1]][x])) * 446000000)/100000
      
      temp[1,35] <- sum(out[[1]][seq(6, 26, by = 2)])/ sum(out[[1]][5:26])
      temp[1,36] <- out[[1]][3] / sum(out[[1]][2:3])
      temp[1,37] <- sum(out[[1]][2:3])
      
      print(paste0("Run ", j, ": ",temp[1,2]))
      output1 <- rbind.data.frame(output1, temp)
    }
    
    colnames(output1)[1:37] <- c("tau", 
                                 
                                 "PropRes_Domestic","PropRes_Netherlands","PropRes_Ireland","PropRes_Germany","PropRes_France","PropRes_Spain","PropRes_Italy",
                                 "PropRes_Belgium","PropRes_Poland","PropRes_Denmark","PropRes_NonEU", 
                                 
                                 "TotInf_Domestic","TotInf_Netherlands","TotInf_Ireland","TotInf_Germany","TotInf_France","TotInf_Spain","TotInf_Italy",
                                 "TotInf_Belgium","TotInf_Poland","TotInf_Denmark","TotInf_NonEU", 
                                 
                                 "Inc_Domestic","Inc_Netherlands","Inc_Ireland","Inc_Germany","Inc_France","Inc_Spain","Inc_Italy",
                                 "Inc_Belgium","Inc_Poland","Inc_Denmark","Inc_NonEU", 
                                 
                                 "Res_PropHum", "Res_PropAnim", "Anim_Inf")
    return(output1)
  })
}


# Plot Base Psi ---------------------------------------------------------------

resistance_data <- plot_analysis[[1]]
norm_names <- sapply(1:11, function(x) paste0(substring(colnames(resistance_data)[13:23], 8)[x], "_res_norm"))
resistance_data[c(norm_names)] <- (resistance_data[2:12]*(sum(resistance_data[13:23])))/
  rowSums((resistance_data[2:12]*(sum(resistance_data[13:23]))))

#Non-normalized Resistance 
res_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[2:12]) 
trim_names <- as.factor(sapply(strsplit(as.character(res_plotdata$variable), split = "_", fixed = TRUE), function(x) x[2]))
res_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Normalized Resistance 
res_norm_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[38:48])
res_norm_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Infection
inf_plotdata <- melt(plot_analysis[[1]], id.vars = c("tau"), measure.vars = colnames(plot_analysis[[1]])[24:34]) 
inf_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Un-normalised resistance
res_comb_psi <- ggplot(res_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.55)) +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  +
  labs(x ="Domestic Ampicillin Usage in Fattening Pig (g/PCU)", y = "Proportion of Human Salmnellosis Resistant to Ampicillin", fill = "Resistance Source")  

#Normalised resistance
res_norm_comb_psi <- ggplot(res_norm_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  +
  labs(x ="Domestic Ampicillin Usage in Fattening Pig (g/PCU)", y = "Proportion of Attributable Human Resistance", fill = "Resistance Source")  

inf_comb_psi <- ggplot(inf_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0 , 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Total Human Salmenollosis (per 100,000)", fill = "Resistance Source") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_vline(xintercept = UK_amp_usage, size = 1.2, col = "red", lty = 2)

# Plot Base Psi ---------------------------------------------------------------

resistance_data <- plot_analysis[[2]]
norm_names <- sapply(1:11, function(x) paste0(substring(colnames(resistance_data)[13:23], 8)[x], "_res_norm"))
resistance_data[c(norm_names)] <- (resistance_data[2:12]*(sum(resistance_data[13:23])))/
  rowSums((resistance_data[2:12]*(sum(resistance_data[13:23]))))

#Non-normalized Resistance 
res_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[2:12]) 
trim_names <- as.factor(sapply(strsplit(as.character(res_plotdata$variable), split = "_", fixed = TRUE), function(x) x[2]))
res_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Normalized Resistance 
res_norm_plotdata <- melt(resistance_data, id.vars = c("tau"), measure.vars = colnames(resistance_data)[38:48])
res_norm_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Infection
inf_plotdata <- melt(plot_analysis[[2]], id.vars = c("tau"), measure.vars = colnames(plot_analysis[[1]])[24:34]) 
inf_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Un-normalised resistance
res_comb <- ggplot(res_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.55)) +
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

inf_comb <- ggplot(inf_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0 , 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Total Human Salmenollosis (per 100,000)", fill = "Resistance Source") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_vline(xintercept = UK_amp_usage, size = 1.2, col = "red", lty = 2)


# Comb Plot ---------------------------------------------------------------

p_comb_res <- ggarrange(res_comb_psi, res_comb, common.legend = T , legend = "bottom", labels = c("A", "B"), font.label = c(size = 20))
p_comb_inf <- ggarrange(inf_comb_psi, inf_comb, common.legend = T , legend = "bottom", labels = c("A", "B"), font.label = c(size = 20))



ggsave(inf_comb_psi, filename = "base_attri_psi_gen.png", dpi = 300, type = "cairo", width = 7, height = 6, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

ggsave(p_comb_res, filename = "base_attri_psi_res.png", dpi = 300, type = "cairo", width = 10, height = 6, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")
ggsave(p_comb_inf, filename = "base_attri_psi_inf.png", dpi = 300, type = "cairo", width = 10, height = 6, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")
