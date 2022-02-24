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
      
      (1-psi)*(share1)*(betaHA*fracimp1*(1-propres_imp1)*Sh) - 
      (1-psi)*(share1)*(1-alpha)*(betaHA*fracimp1*propres_imp1*Sh)-
      
      (1-psi)*(share2)*(betaHA*fracimp2*(1-propres_imp2)*Sh) - 
      (1-psi)*(share2)*(1-alpha)*(betaHA*fracimp2*propres_imp2*Sh) - 
      
      (1-psi)*(share3)*(betaHA*fracimp3*(1-propres_imp3)*Sh) - 
      (1-psi)*(share3)*(1-alpha)*(betaHA*fracimp3*propres_imp3*Sh) -
      
      (1-psi)*(share4)*(betaHA*fracimp4*(1-propres_imp4)*Sh) - 
      (1-psi)*(share4)*(1-alpha)*(betaHA*fracimp4*propres_imp4*Sh)-
      
      (1-psi)*(share5)*(betaHA*fracimp5*(1-propres_imp5)*Sh) - 
      (1-psi)*(share5)*(1-alpha)*(betaHA*fracimp5*propres_imp5*Sh) - 
      
      (1-psi)*(share6)*(betaHA*fracimp6*(1-propres_imp6)*Sh) - 
      (1-psi)*(share6)*(1-alpha)*(betaHA*fracimp6*propres_imp6*Sh) - 
      
      (1-psi)*(share7)*(betaHA*fracimp7*(1-propres_imp7)*Sh) - 
      (1-psi)*(share7)*(1-alpha)*(betaHA*fracimp7*propres_imp7*Sh) - 
      
      (1-psi)*(share8)*(betaHA*fracimp8*(1-propres_imp8)*Sh) - 
      (1-psi)*(share8)*(1-alpha)*(betaHA*fracimp8*propres_imp8*Sh) -
      
      (1-psi)*(share9)*(betaHA*fracimp9*(1-propres_imp9)*Sh) - 
      (1-psi)*(share9)*(1-alpha)*(betaHA*fracimp9*propres_imp9*Sh) - 
      
      (1-psi)*(share_nEU)*(betaHA*fracimp_nEU*(1-propres_impnEU)*Sh) - 
      (1-psi)*(share_nEU)*(1-alpha)*(betaHA*fracimp_nEU*propres_impnEU*Sh)
    
    dIshDA = psi*betaHA*(Isa*eta)*Sh - rh*IshDA - uh*IshDA 
    dIrhDA = psi*(1-alpha)*betaHA*(Ira*eta)*Sh - rh*IrhDA - uh*IrhDA  
    
    dIshA1 = (1-psi)*(share1)*(betaHA*fracimp1*(1-propres_imp1)*Sh) - rh*IshA1 - uh*IshA1 
    dIrhA1 = (1-psi)*(share1)*(1-alpha)*(betaHA*fracimp1*propres_imp1*Sh) - rh*IrhA1 - uh*IrhA1  
    
    dIshA2 = (1-psi)*(share2)*(betaHA*fracimp2*(1-propres_imp2)*Sh) - rh*IshA2 - uh*IshA2 
    dIrhA2 = (1-psi)*(share2)*(1-alpha)*(betaHA*fracimp2*propres_imp2*Sh) - rh*IrhA2 - uh*IrhA2  
    
    dIshA3 = (1-psi)*(share3)*(betaHA*fracimp3*(1-propres_imp3)*Sh) - rh*IshA3 - uh*IshA3 
    dIrhA3 = (1-psi)*(share3)*(1-alpha)*(betaHA*fracimp3*propres_imp3*Sh) - rh*IrhA3 - uh*IrhA3  
    
    dIshA4 = (1-psi)*(share4)*(betaHA*fracimp4*(1-propres_imp4)*Sh) - rh*IshA4 - uh*IshA4 
    dIrhA4 = (1-psi)*(share4)*(1-alpha)*(betaHA*fracimp4*propres_imp4*Sh) - rh*IrhA4 - uh*IrhA4  
    
    dIshA5 = (1-psi)*(share5)*(betaHA*fracimp5*(1-propres_imp5)*Sh) - rh*IshA5 - uh*IshA5 
    dIrhA5 = (1-psi)*(share5)*(1-alpha)*(betaHA*fracimp5*propres_imp5*Sh) - rh*IrhA5 - uh*IrhA5  
    
    dIshA6 = (1-psi)*(share6)*(betaHA*fracimp6*(1-propres_imp6)*Sh) - rh*IshA6 - uh*IshA6 
    dIrhA6 = (1-psi)*(share6)*(1-alpha)*(betaHA*fracimp6*propres_imp6*Sh) - rh*IrhA6 - uh*IrhA6  
    
    dIshA7 = (1-psi)*(share7)*(betaHA*fracimp7*(1-propres_imp7)*Sh) - rh*IshA7 - uh*IshA7 
    dIrhA7 = (1-psi)*(share7)*(1-alpha)*(betaHA*fracimp7*propres_imp7*Sh) - rh*IrhA7 - uh*IrhA7  
    
    dIshA8 = (1-psi)*(share8)*(betaHA*fracimp8*(1-propres_imp8)*Sh) - rh*IshA8 - uh*IshA8 
    dIrhA8 = (1-psi)*(share8)*(1-alpha)*(betaHA*fracimp8*propres_imp8*Sh) - rh*IrhA8 - uh*IrhA8  
    
    dIshA9 = (1-psi)*(share9)*(betaHA*fracimp9*(1-propres_imp9)*Sh) - rh*IshA9 - uh*IshA9 
    dIrhA9 = (1-psi)*(share9)*(1-alpha)*(betaHA*fracimp9*propres_imp9*Sh) - rh*IrhA9 - uh*IrhA9  
    
    dIshAnEU = (1-psi)*(share_nEU)*(betaHA*fracimp_nEU*(1-propres_impnEU)*Sh) - rh*IshAnEU - uh*IshAnEU 
    dIrhAnEU = (1-psi)*(share_nEU)*(1-alpha)*(betaHA*fracimp_nEU*propres_impnEU*Sh) - rh*IrhAnEU - uh*IrhAnEU  
    
    
    CumS_D = psi*betaHA*(Isa*eta)*Sh
    CumS_1 = (1-psi)*(share1)*(betaHA*fracimp1*(1-propres_imp1)*Sh)
    CumS_2 = (1-psi)*(share2)*(betaHA*fracimp2*(1-propres_imp2)*Sh)
    CumS_3 = (1-psi)*(share3)*(betaHA*fracimp3*(1-propres_imp3)*Sh)
    CumS_4 = (1-psi)*(share4)*(betaHA*fracimp4*(1-propres_imp4)*Sh)
    CumS_5 = (1-psi)*(share5)*(betaHA*fracimp5*(1-propres_imp5)*Sh)
    CumS_6 = (1-psi)*(share6)*(betaHA*fracimp6*(1-propres_imp6)*Sh)
    CumS_7 = (1-psi)*(share7)*(betaHA*fracimp7*(1-propres_imp7)*Sh)
    CumS_8 = (1-psi)*(share8)*(betaHA*fracimp8*(1-propres_imp8)*Sh)
    CumS_9 = (1-psi)*(share9)*(betaHA*fracimp9*(1-propres_imp9)*Sh)
    CumS_nEU = (1-psi)*(share_nEU)*(betaHA*fracimp_nEU*(1-propres_impnEU)*Sh)
    
    CumR_D = psi*(1-alpha)*betaHA*(Ira*eta)*Sh
    CumR_1 = (1-psi)*(share1)*(1-alpha)*(betaHA*fracimp1*propres_imp1*Sh)
    CumR_2 =  (1-psi)*(share2)*(1-alpha)*(betaHA*fracimp2*propres_imp2*Sh)
    CumR_3 = (1-psi)*(share3)*(1-alpha)*(betaHA*fracimp3*propres_imp3*Sh)
    CumR_4 = (1-psi)*(share4)*(1-alpha)*(betaHA*fracimp4*propres_imp4*Sh)
    CumR_5 = (1-psi)*(share5)*(1-alpha)*(betaHA*fracimp5*propres_imp5*Sh)
    CumR_6 = (1-psi)*(share6)*(1-alpha)*(betaHA*fracimp6*propres_imp6*Sh)
    CumR_7 = (1-psi)*(share7)*(1-alpha)*(betaHA*fracimp7*propres_imp7*Sh)
    CumR_8 = (1-psi)*(share8)*(1-alpha)*(betaHA*fracimp8*propres_imp8*Sh)
    CumR_9 =  (1-psi)*(share9)*(1-alpha)*(betaHA*fracimp9*propres_imp9*Sh)
    CumR_nEU = (1-psi)*(share_nEU)*(1-alpha)*(betaHA*fracimp_nEU*propres_impnEU*Sh)
    
    CumS = psi*betaHA*(Isa*eta)*Sh + 
      (1-psi)*(share1)*(betaHA*fracimp1*(1-propres_imp1)*Sh) + 
      (1-psi)*(share2)*(betaHA*fracimp2*(1-propres_imp2)*Sh) +
      (1-psi)*(share3)*(betaHA*fracimp3*(1-propres_imp3)*Sh) + 
      (1-psi)*(share4)*(betaHA*fracimp4*(1-propres_imp4)*Sh) + 
      (1-psi)*(share5)*(betaHA*fracimp5*(1-propres_imp5)*Sh) + 
      (1-psi)*(share6)*(betaHA*fracimp6*(1-propres_imp6)*Sh) + 
      (1-psi)*(share7)*(betaHA*fracimp7*(1-propres_imp7)*Sh) + 
      (1-psi)*(share8)*(betaHA*fracimp8*(1-propres_imp8)*Sh) + 
      (1-psi)*(share9)*(betaHA*fracimp9*(1-propres_imp9)*Sh) + 
      (1-psi)*(share_nEU)*(betaHA*fracimp_nEU*(1-propres_impnEU)*Sh)
    
    CumR = psi*(1-alpha)*betaHA*(Ira*eta)*Sh + 
      (1-psi)*(share1)*(1-alpha)*(betaHA*fracimp1*propres_imp1*Sh) + 
      (1-psi)*(share2)*(1-alpha)*(betaHA*fracimp2*propres_imp2*Sh) +
      (1-psi)*(share3)*(1-alpha)*(betaHA*fracimp3*propres_imp3*Sh) + 
      (1-psi)*(share4)*(1-alpha)*(betaHA*fracimp4*propres_imp4*Sh) + 
      (1-psi)*(share5)*(1-alpha)*(betaHA*fracimp5*propres_imp5*Sh) + 
      (1-psi)*(share6)*(1-alpha)*(betaHA*fracimp6*propres_imp6*Sh) +
      (1-psi)*(share7)*(1-alpha)*(betaHA*fracimp7*propres_imp7*Sh) + 
      (1-psi)*(share8)*(1-alpha)*(betaHA*fracimp8*propres_imp8*Sh) + 
      (1-psi)*(share9)*(1-alpha)*(betaHA*fracimp9*propres_imp9*Sh) + 
      (1-psi)*(share_nEU)*(1-alpha)*(betaHA*fracimp_nEU*propres_impnEU*Sh)
    
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
                CumS_D, CumS_1,CumS_2,CumS_3,CumS_4,CumS_5,
                CumS_6, CumS_7,CumS_8,CumS_9,CumS_nEU,
                CumR_D, CumR_1,CumR_2,CumR_3,CumR_4,CumR_5,
                CumR_6, CumR_7,CumR_8,CumR_9,CumR_nEU,
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

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha")

post_amp <- read.csv(tail(list.files(path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha", pattern = "complex"), 1))
MAP_parms <- map_estimate(post_amp)
MAP_parms <- data.frame("Parameter" = names(post_amp), 
                        "MAP_Estimate" = colMeans(post_amp))

# Test --------------------------------------------------------------------

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

thetaparm = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = 0.656,
              
              share1 = country_data_imp[2,"Normalised_Usage_2018"], share2 = country_data_imp[3,"Normalised_Usage_2018"], share3 = country_data_imp[4,"Normalised_Usage_2018"], 
              share4 = country_data_imp[5,"Normalised_Usage_2018"], share5 = country_data_imp[6,"Normalised_Usage_2018"], share6 = country_data_imp[7,"Normalised_Usage_2018"], 
              share7 = country_data_imp[8,"Normalised_Usage_2018"], share8 = country_data_imp[9,"Normalised_Usage_2018"], 
              share9 = country_data_imp[10,"Normalised_Usage_2018"], share_nEU = 1 - sum(country_data_imp[2:10,"Normalised_Usage_2018"]),
              
              betaAA = MAP_parms["betaAA", 2], phi = MAP_parms["phi", 2], kappa = MAP_parms["kappa", 2], alpha = MAP_parms["alpha", 2], 
              zeta = MAP_parms["zeta", 2], betaHA = MAP_parms["betaHA", 2], imp_nEU = MAP_parms["imp_nEU", 2], propres_impnEU = MAP_parms["propres_impnEU", 2], 
              
              fracimp1 = country_data_imp[2,"FBD_gen"], fracimp2 = country_data_imp[3,"FBD_gen"], fracimp3 = country_data_imp[4,"FBD_gen"], 
              fracimp4 = country_data_imp[5,"FBD_gen"], fracimp5 = country_data_imp[6,"FBD_gen"], fracimp6 = country_data_imp[7,"FBD_gen"], 
              fracimp7 = country_data_imp[8,"FBD_gen"], fracimp8 = country_data_imp[9,"FBD_gen"],
              fracimp9 = country_data_imp[9,"FBD_gen"], fracimp_nEU = country_data_imp[10,"FBD_gen"],
              
              propres_imp1 = country_data_imp[2,"FBD_res"], propres_imp2 = country_data_imp[3,"FBD_res"], propres_imp3 = country_data_imp[4,"FBD_res"], propres_imp4 = country_data_imp[5,"FBD_res"], 
              propres_imp5 = country_data_imp[6,"FBD_res"], propres_imp6 = country_data_imp[7,"FBD_res"], propres_imp7 = country_data_imp[8,"FBD_res"], propres_imp8 = country_data_imp[9,"FBD_res"], 
              propres_imp9 = country_data_imp[10,"FBD_res"],
              
              eta = 0.11016, tau = UK_amp_usage)

tau_range <- seq(0, 0.055, by = 0.004)
#tau_range <- append(tau_range, UK_amp_usage)
tauoutput <- data.frame(matrix(nrow = length(tau_range), ncol = 16))
out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = thetaparm)

for (i in 1:length(tau_range)) {
  
  thetaparm["tau"] = tau_range[i]
  out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = thetaparm)
  
  tauoutput[i,] <- c(tau_range[i],
                     ((out[[2]] + out[[2+11]])*(446000000))/100000,
                     ((out[[3]] + out[[3+11]])*(446000000))/100000,
                     ((out[[4]] + out[[4+11]])*(446000000))/100000,
                     ((out[[5]] + out[[5+11]])*(446000000))/100000,
                     ((out[[6]] + out[[6+11]])*(446000000))/100000,
                     ((out[[7]] + out[[7+11]])*(446000000))/100000,
                     ((out[[8]] + out[[8+11]])*(446000000))/100000,
                     ((out[[9]] + out[[9+11]])*(446000000))/100000,
                     ((out[[10]] + out[[10+11]])*(446000000))/100000,
                     ((out[[11]] + out[[11+11]])*(446000000))/100000,
                     ((out[[12]] + out[[12+11]])*(446000000))/100000,
                     ((out[[24]] + out[[25]])*(446000000))/100000,
                     (out[[1]][["Isa"]] + out[[1]][["Ira"]]), 
                     out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                     sum(out[[1]][seq(6, 26, by = 2)]) / (sum(out[[1]][5:26])))
}

colnames(tauoutput) <- c("tau", "IncH_D","IncH_1","IncH_2",
                         "IncH_3","IncH_4","IncH_5","IncH_6",
                         "IncH_7","IncH_8","IncH_9","IncH_nEU",
                         "IncH",
                         "ICombA", "ResPropAnim", "ResPropHum")


# Plotting ----------------------------------------------------------------


inf_plotdata <- melt(tauoutput, id.vars = c("tau"), measure.vars = names(tauoutput)[2:12]) 
#inf_plotdata$variable <- factor(names(tauoutput)[2:12], levels = unique(names(tauoutput)[2:12]))

inf_comb <- ggplot(inf_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.004) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0 , 0.15)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Total Human Salmenollosis (per 100,000)", fill = "Resistance Source") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_vline(xintercept = UK_amp_usage, size = 1.2, col = "red", lty = 2)
