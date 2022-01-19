library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity"); library("parallel")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate"); library("rootSolve")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")

# Single Model ------------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      0.5*zeta*Sa*(1-alpha) - 0.5*zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + 0.5*zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + 0.5*zeta*Sa*(1-alpha)
    
    
    dSh = uh + rh*(1-Sh) - uh*Sh - 
      
      psi*betaHD*(Isa*eta)*Sh - 
      psi*(1-alpha)*betaHD*(Ira*eta)*Sh - 
      
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
    
    dIshDA = psi*betaHD*(Isa*eta)*Sh - rh*IshDA - uh*IshDA 
    dIrhDA = psi*(1-alpha)*betaHD*(Ira*eta)*Sh - rh*IrhDA - uh*IrhDA  
    
    dIshA1 = (1-psi)*(imp1)*(betaHI_EU*fracimp1*(1-propres_imp1)*Sh) - rh*IshA1 - uh*IshA1 
    dIrhA1 = (1-psi)*(imp1)*(1-alpha)*(betaHI_EU*fracimp1*propres_imp1*Sh) - rh*IrhA1 - uh*IrhA1  
    
    dIshA2 = (1-psi)*(imp2)*(betaHI_EU*fracimp2*(1-propres_imp2)*Sh) - rh*IshA2 - uh*IshA2 
    dIrhA2 = (1-psi)*(imp2)*(1-alpha)*(betaHI_EU*fracimp2*propres_imp2*Sh) - rh*IrhA2 - uh*IrhA2  
    
    dIshA3 = (1-psi)*(imp3)*(betaHI_EU*fracimp3*(1-propres_imp3)*Sh) - rh*IshA3 - uh*IshA3 
    dIrhA3 = (1-psi)*(imp3)*(1-alpha)*(betaHI_EU*fracimp3*propres_imp3*Sh) - rh*IrhA3 - uh*IrhA3  
    
    dIshA4 = (1-psi)*(imp4)*(betaHI_EU*fracimp4*(1-propres_imp4)*Sh) - rh*IshA4 - uh*IshA4 
    dIrhA4 = (1-psi)*(imp4)*(1-alpha)*(betaHI_EU*fracimp4*propres_imp4*Sh) - rh*IrhA4 - uh*IrhA4  
    
    dIshA5 = (1-psi)*(imp5)*(betaHI_EU*fracimp5*(1-propres_imp5)*Sh) - rh*IshA5 - uh*IshA5 
    dIrhA5 = (1-psi)*(imp5)*(1-alpha)*(betaHI_EU*fracimp5*propres_imp5*Sh) - rh*IrhA5 - uh*IrhA5  
    
    dIshA6 = (1-psi)*(imp6)*(betaHI_EU*fracimp6*(1-propres_imp6)*Sh) - rh*IshA6 - uh*IshA6 
    dIrhA6 = (1-psi)*(imp6)*(1-alpha)*(betaHI_EU*fracimp6*propres_imp6*Sh) - rh*IrhA6 - uh*IrhA6  
    
    dIshA7 = (1-psi)*(imp7)*(betaHI_EU*fracimp7*(1-propres_imp7)*Sh) - rh*IshA7 - uh*IshA7 
    dIrhA7 = (1-psi)*(imp7)*(1-alpha)*(betaHI_EU*fracimp7*propres_imp7*Sh) - rh*IrhA7 - uh*IrhA7  
    
    dIshA8 = (1-psi)*(imp8)*(betaHI_EU*fracimp8*(1-propres_imp8)*Sh) - rh*IshA8 - uh*IshA8 
    dIrhA8 = (1-psi)*(imp8)*(1-alpha)*(betaHI_EU*fracimp8*propres_imp8*Sh) - rh*IrhA8 - uh*IrhA8  
    
    dIshA9 = (1-psi)*(imp9)*(betaHI_EU*fracimp9*(1-propres_imp9)*Sh) - rh*IshA9 - uh*IshA9 
    dIrhA9 = (1-psi)*(imp9)*(1-alpha)*(betaHI_EU*fracimp9*propres_imp9*Sh) - rh*IrhA9 - uh*IrhA9  
    
    dIshAnEU = (1-psi)*(imp_nEU)*(betaHI_EU*fracimp_nEU*(1-propres_impnEU)*Sh) - rh*IshAnEU - uh*IshAnEU 
    dIrhAnEU = (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_EU*fracimp_nEU*propres_impnEU*Sh) - rh*IrhAnEU - uh*IrhAnEU  
    
    CumS = psi*betaHD*(Isa*eta)*Sh + 
      (1-psi)*(imp1)*(betaHI_EU*fracimp1*(1-propres_imp1)*Sh) + 
      (1-psi)*(imp2)*(betaHI_EU*fracimp2*(1-propres_imp2)*Sh) +
      (1-psi)*(imp3)*(betaHI_EU*fracimp3*(1-propres_imp3)*Sh) + 
      (1-psi)*(imp4)*(betaHI_EU*fracimp4*(1-propres_imp4)*Sh) + 
      (1-psi)*(imp5)*(betaHI_EU*fracimp5*(1-propres_imp5)*Sh) + 
      (1-psi)*(imp6)*(betaHI_EU*fracimp6*(1-propres_imp6)*Sh) + 
      (1-psi)*(imp7)*(betaHI_EU*fracimp7*(1-propres_imp7)*Sh) + 
      (1-psi)*(imp8)*(betaHI_EU*fracimp8*(1-propres_imp8)*Sh) + 
      (1-psi)*(imp9)*(betaHI_EU*fracimp9*(1-propres_imp9)*Sh) + 
      (1-psi)*(imp_nEU)*(betaHI_EU*fracimp_nEU*(1-propres_impnEU)*Sh)
    
    CumR = psi*(1-alpha)*betaHD*(Ira*eta)*Sh + 
      (1-psi)*(imp1)*(1-alpha)*(betaHI_EU*fracimp1*propres_imp1*Sh) + 
      (1-psi)*(imp2)*(1-alpha)*(betaHI_EU*fracimp2*propres_imp2*Sh) +
      (1-psi)*(imp3)*(1-alpha)*(betaHI_EU*fracimp3*propres_imp3*Sh) + 
      (1-psi)*(imp4)*(1-alpha)*(betaHI_EU*fracimp4*propres_imp4*Sh) + 
      (1-psi)*(imp5)*(1-alpha)*(betaHI_EU*fracimp5*propres_imp5*Sh) + 
      (1-psi)*(imp6)*(1-alpha)*(betaHI_EU*fracimp6*propres_imp6*Sh) +
      (1-psi)*(imp7)*(1-alpha)*(betaHI_EU*fracimp7*propres_imp7*Sh) + 
      (1-psi)*(imp8)*(1-alpha)*(betaHI_EU*fracimp8*propres_imp8*Sh) + 
      (1-psi)*(imp9)*(1-alpha)*(betaHI_EU*fracimp9*propres_imp9*Sh) + 
      (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_EU*fracimp_nEU*propres_impnEU*Sh)
    
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

#setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/dump")
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/dump")

post_amp <- read.csv(tail(list.files(path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/dump", pattern = "complex"), 1))
MAP_parms <- map_estimate(post_amp)
MAP_parms <- data.frame("Parameter" = names(post_amp), 
                        "MAP_Estimate" = colMeans(post_amp))

#New Import Parms 
thetaparm = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = 0.656,
              
              fracimp1 = country_data_imp[2,"Normalised_Usage_2018"], fracimp2 = country_data_imp[3,"Normalised_Usage_2018"], fracimp3 = country_data_imp[4,"Normalised_Usage_2018"], 
              fracimp4 = country_data_imp[5,"Normalised_Usage_2018"], fracimp5 = country_data_imp[6,"Normalised_Usage_2018"], fracimp6 = country_data_imp[7,"Normalised_Usage_2018"], 
              fracimp7 = country_data_imp[8,"Normalised_Usage_2018"], fracimp8 = country_data_imp[9,"Normalised_Usage_2018"], 
              fracimp9 = country_data_imp[10,"Normalised_Usage_2018"], fracimp_nEU = 1 - sum(country_data_imp[2:10,"Normalised_Usage_2018"]),
              
              betaAA = MAP_parms["betaAA", 2], phi = MAP_parms["phi", 2], kappa = MAP_parms["kappa", 2], alpha = MAP_parms["alpha", 2], 
              zeta = MAP_parms["zeta", 2], betaHD = MAP_parms["betaHD", 2], betaHI_EU = MAP_parms["betaHI_EU", 2], imp_nEU = MAP_parms["imp_nEU", 2], propres_impnEU = MAP_parms["propres_impnEU", 2], 
              
              imp1 = country_data_imp[2,"FBD_gen"], imp2 = country_data_imp[3,"FBD_gen"], imp3 = country_data_imp[4,"FBD_gen"], imp4 = country_data_imp[5,"FBD_gen"], 
              imp5 = country_data_imp[6,"FBD_gen"], imp6 = country_data_imp[7,"FBD_gen"], imp7 = country_data_imp[8,"FBD_gen"], imp8 = country_data_imp[9,"FBD_gen"],
              imp8 = country_data_imp[9,"FBD_gen"], imp9 = country_data_imp[10,"FBD_gen"],
              
              propres_imp1 = country_data_imp[2,"FBD_res"], propres_imp2 = country_data_imp[3,"FBD_res"], propres_imp3 = country_data_imp[4,"FBD_res"], propres_imp4 = country_data_imp[5,"FBD_res"], 
              propres_imp5 = country_data_imp[6,"FBD_res"], propres_imp6 = country_data_imp[7,"FBD_res"], propres_imp7 = country_data_imp[8,"FBD_res"], propres_imp8 = country_data_imp[9,"FBD_res"], 
              propres_imp9 = country_data_imp[10,"FBD_res"],
              
              eta = 0.11016, tau = UK_amp_usage)

# Run the Import Analysis Model -----------------------------------------------------------

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

usage_threshold <- c(seq(0, 1, by = 0.02), 0.656)

#Function for the Import %

import_res_func <- function(parms, init, usage_threshold, uk_usage) {
  parmtau1 <- c(0, UK_amp_usage)
  usage_frame <- data.frame(matrix(nrow = length(usage_threshold), ncol = 2))
  
  for(x in 1:length(usage_threshold)) {
    parms[["psi"]] <- usage_threshold[x]
    output1 <- data.frame(matrix(nrow = 2, ncol =3))
    dump_data <- vector()
    
    for (i in 1:2) {
      
      temp <- vector()
      parms["tau"] <- parmtau1[i]
      out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms)
      
      output1[i,] <- c(parmtau1[i],
                       ((sum(out[[1]][5:26]))*(446000000))/100000,
                       sum(out[[1]][seq(6, 26, by = 2)]) / (sum(out[[1]][5:26])))
    }
    
    colnames(output1) <- c("tau", "ICombH","IResRat")
    
    usage_frame[x,] <- c((1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == UK_amp_usage]))*100,    
                         usage_threshold[x])
  }
  
  colnames(usage_frame) <- c("relchange", "domusage")
  usage_frame$normchange <- (usage_frame$relchange / usage_frame$relchange[usage_frame$domusage == 0.656]) * 100
  
  return(usage_frame)
}

# Granular Contamination/Resistance (Incremental) Sensitivity ---------------------------

amp_cont_data_incr <- data.frame(matrix(NA, ncol = 10, nrow = 0))
amp_res_data_incr <- data.frame(matrix(NA, ncol = 10, nrow = 0))

for(i in 1:16) {
  tempcut <- rep(c(seq(0,0.3, by = 0.02))[i], 10)
  amp_cont_data_incr <- rbind(amp_cont_data_incr, tempcut)
  colnames(amp_cont_data_incr) <- grep("fracimp", names(thetaparm), value=TRUE)
}

for(i in 1:11) {
  temp <- rep(c(seq(0,1, by = 0.1))[i], 10)
  amp_res_data_incr <- rbind(amp_res_data_incr, temp)
  colnames(amp_res_data_incr) <- grep("propres_imp", names(thetaparm), value=TRUE)
}

data_list_incr <- list(amp_cont_data_incr, amp_res_data_incr)

#Run the Scenario Analysis
data_imp_list_incr <- list()

for(j in 1:3){
  data_imp_list_incr[[j]] <- local ({
    if(j == 1) {
      output <- import_res_func(thetaparm, init, usage_threshold, UK_amp_usage)
      print("base")
      
    } else {
      output <- list()
      parms_list <- data_list_incr[[j-1]]
      parms1 <- thetaparm
      
      for(i in 1:nrow(parms_list)) {
        parms1[colnames(parms_list)] <- parms_list[i,]
        out <- import_res_func(parms1, init, usage_threshold, UK_amp_usage)
        output[[i]] <- out 
        print(paste0(c("base", "cont", "res")[j], " - ", i/nrow(parms_list)))
      }
    }
    return(output)
  })
}

#[[3]][[1]]
#Res = 0 
#Domestic Usage is 0

#Plot the Incremental Analysis 
p_incr_list <- list()

for(i in 1:2) {
  p_incr_list[[i]] <- local ({
    
    data <- as.data.frame(cbind(sapply(data_imp_list_incr[[i+1]], "[[", 1))) # the 1 at the end refers to either normalised or un-normalised change in FBD
    colnames(data) <- sapply((list(amp_cont_data_incr[,1], amp_res_data_incr[,1])[[i]])*100, function(x) paste0(x, "%"))
    data$Baseline <- data_imp_list_incr[[1]][,1]
    data$usage <- c(seq(0, 1, by = 0.02), 0.656)
    
    plot_cont <- melt(data, measure.vars = (head(colnames(data), -1)), id.vars = c("usage"))
    print(max(plot_cont$value, na.rm = T)*1.2)
    p_incr <- ggplot(plot_cont, aes(x = usage, y = value, col = variable, size = variable, lty = variable)) + geom_line() +
      scale_color_manual(values = c(viridis::viridis((ncol(data)-2)), "red")) + 
      geom_vline(xintercept = 0.656, size = 2, col = "hotpink" , alpha = 0.6) + 
      scale_size_manual(values = c(rep(1, (ncol(data)-2)), 2)) + 
      scale_linetype_manual(values = c(rep(1, (ncol(data)-2)), 2)) + theme_bw() + 
      theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
            axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
            legend.spacing.x = unit(0.3, 'cm')) + labs(x ="Proportion of Food From Domestic Origins", 
                                                       y = "Reduction in Human Resistance due to Curtailment", 
                                                       color = c("Proportion Contaminated", "Proportion Resistant")[i]) +
      scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot_cont$value, na.rm = T)*1.2)) + guides(size= "none", linetype = "none") 
  
  #ggsave(p_incr, filename = paste0("import_sens_incr_", c("cont","res")[i] ,".png"), dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
  #       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures/New_Figures")
  
  return(p_incr)
  })
}


# Testing -----------------------------------------------------------------

output <- data.frame()

for(i in 1:length(data_imp_list_incr[[2]])) {
  temp <- data_imp_list_incr[[2]][[i]]
  temp$normchange <- (temp$relchange / temp$relchange[temp$domusage == 0.656]) * 100
  temp$group <- as.factor(seq(0,0.3, by = 0.02)[i])
  print(temp)
  output <- rbind(output, temp)
}

output[output$domusage == 0.656,]

ggplot(output, aes(y= normchange,x = domusage, col = group)) + geom_line()

ggplot(output, aes(y= relchange ,x = domusage, col = group)) + geom_line()

# Plotting Baseline Model -------------------------------------------------

output_base <- import_res_func(parms, init, usage_threshold, UK_amp)

# Coordinates of the upper and lower areas
trsup <- data.frame(x=c(0,0,0.656), y=c(0,100,100)) 
trinf <- data.frame(x=c(0,0.656,0.656), y=c(0,0,100))

c("More than proportional decrease in efficacy", "Less than proportional decrease in efficacy")

# Use geom_polygon for coloring the two areas
p_base <- ggplot(data=output_base, aes(x=domusage, y=relchange)) +
  geom_polygon(aes(x=x, y=y), data=trsup, fill="#00FF0066") +
  geom_polygon(aes(x=x, y=y), data=trinf, fill= "#FF000066") +
  geom_abline(intercept = 0, slope = 100/0.656, size = 3) +
  geom_line(color = "red", size = 1.5, lty = 3) + theme_bw() +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="Proportion of Food From Domestic Origins", y = "Reduction in the Efficacy of Livestock Antibiotic Stewardship", 
       color = "Density") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  + 
  geom_line(data = data_base, aes(x = domusage, y = relchange ), color = "red", size = 1.5, lty = 3, inherit.aes = F) +
  annotate("text", x = c(0.15, 0.5), y = c(75, 25), label = c('atop("Less than proportional \n decrease in efficacy", bold("(less impact)"))',
                                                               'atop("More than proportional \n decrease in efficacy", bold("(greater impact)"))'),
           size = 4, col = "black", parse = TRUE)

ggsave(p_base, filename = paste0("base_plot_import.png"), dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
