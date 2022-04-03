library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity"); library("parallel")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate"); library("rootSolve")

rm(list=ls())
#setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")
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

# Parameters --------------------------------------------------------------
#New Import Parms 

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

# Identify the Distributions ----------------------------------------------

#Baseline 
base_parms <- thetaparm

#Homogenous
base_homo_shareimp <- data.frame(matrix(rbeta(1000*10, shape1 =  1, shape2 =  1), 
                  nrow = 1000, ncol = 10))
homo_shareimp <- t(apply(base_homo_shareimp, 1, function(x) x/sum(x)))
colnames(homo_shareimp) <- grep("share", names(thetaparm), value = T)

#Skewed
base_skew_shareimp <- data.frame(matrix(rbeta(1000*10, shape1 =  0.5, shape2 =  2), 
                                  nrow = 1000, ncol = 10))
skew_shareimp <- t(apply(base_skew_shareimp,1, function(x) x/sum(x)))
colnames(skew_shareimp) <- grep("share", names(thetaparm), value = T)

# Run the model -----------------------------------------------------------

explore_parms <- list(homo_shareimp, skew_shareimp)

usage_threshold <- c(seq(0, 1, by = 0.05), 0.656)

explore_parms_frame <- list()

#Run the Uncertainty Analysis

for(z in 1:2) {
  explore_parms_frame[[z]] <- local({ 
    
    explore_sublist <- data.frame(matrix(nrow = 0, ncol = 5))
    explore_parms_sub <- explore_parms[[z]]
    
    for(j in 1:nrow(explore_parms_sub)) {
      
      usage_frame <- data.frame(matrix(nrow = length(usage_threshold), ncol = 3))
      
      thetaparm[colnames(explore_parms_sub)] <- explore_parms_sub[j,]
      
      for(x in 1:length(usage_threshold)) {
        thetaparm[["psi"]] <- usage_threshold[x]
        output1 <- data.frame(matrix(nrow = 2, ncol =3))

        for (i in 1:2) {
          thetaparm["tau"] <- c(0, UK_amp_usage)[i]
          out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = thetaparm)
          
          output1[i,] <- c(c(0, UK_amp_usage)[i],
                           ((sum(out[[1]][5:26]))*(446000000))/100000,
                           sum(out[[1]][seq(6, 26, by = 2)]) / (sum(out[[1]][5:26])))
        }
        
        colnames(output1) <- c("tau", "ICombH","IResRat")
        
        usage_frame[x,] <- c((1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == UK_amp_usage]))*100,    
                             usage_threshold[x],
                             j) #run
      }
      colnames(usage_frame) <- c("relchange", "domusage", "run_no")
      usage_frame$normchange <- (usage_frame$relchange / usage_frame$relchange[usage_frame$domusage == 0.656]) * 100
      usage_frame$explore <- c("homo", "skew")[z]
      explore_sublist <- rbind(usage_frame, explore_sublist)
      print(paste0( c("homo", "skew")[z], " | " ,(j / nrow(explore_parms_sub))*100 , "%"))
      
      }
    return(explore_sublist)
  })
}

comb_explore_list <- list()

for(i in 1:2) {
  temp <- aggregate(relchange ~ domusage, explore_parms_frame[[i]], mean)
  temp$parms <- c("homo", "skew")[i]
  temp$max <- aggregate(relchange ~ domusage, explore_parms_frame[[i]], max)[,2]
  temp$min <- aggregate(relchange ~ domusage, explore_parms_frame[[i]], min)[,2]
  comb_explore_list[[i]] <- temp
}


#Run Baseline

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

usage_frame_baseline <- data.frame(matrix(nrow = length(usage_threshold), ncol = 2))

for(x in 1:length(usage_threshold)) {
  thetaparm[["psi"]] <- usage_threshold[x]
  output1 <- data.frame(matrix(nrow = 2, ncol =3))
  
  for (i in 1:2) {

    thetaparm["tau"] <- c(0, UK_amp_usage)[i]
    out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = thetaparm)
    
    output1[i,] <- c(c(0, UK_amp_usage)[i],
                     ((sum(out[[1]][5:26]))*(446000000))/100000,
                     sum(out[[1]][seq(6, 26, by = 2)]) / (sum(out[[1]][5:26])))
  }
  
  colnames(output1) <- c("tau", "ICombH","IResRat")
  
  usage_frame_baseline[x,] <- c(usage_threshold[x],
                                (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == UK_amp_usage]))*100)
}

colnames(usage_frame_baseline) <- c("domusage", "relchange")

#usage_frame_baseline$normchange <- (usage_frame_baseline$relchange / usage_frame_baseline$relchange[usage_frame_baseline$domusage == 0.656]) * 100

usage_frame_baseline$parms <- "baseline"
usage_frame_baseline$max <- 0
usage_frame_baseline$min <- 0

comb_explore_list[[length(comb_explore_list)+1]] <- usage_frame_baseline

# Stratified Plotting -----------------------------------------------------

plot_imp <- list()

for(i in 1:2) {
  comb_plot <- rbind(comb_explore_list[[3]], comb_explore_list[[i]])
  labs <- c("Beta"*"("*alpha*"="*1*", "*beta*"="~1*")", "Beta"*"("*alpha*"="*0.5*", "*beta*"="~2*")")[i]

  dist_plots <- data.frame("y" = dbeta(seq(0, 1, by = 0.01), c(1, 0.5)[i], c(1,2)[i]),
                           "x" = seq(0, 1, by = 0.01))
  
  dist <- ggplot(dist_plots, aes(x=x, y = y)) + geom_line(size = 1.2) + theme_bw() + 
    theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text.x =element_text(size=12), 
          axis.text.y = element_blank(), 
          axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0,0,0,0), "cm"),
          legend.spacing.x = unit(0.3, 'cm')) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
    labs(x = "Probability", y = "") + scale_color_manual(values = c("darkblue", "red"))
    
  p_import_comp <- ggplot(comb_plot, aes(x = domusage, y = relchange, col = parms, lty = parms)) +  theme_bw() + 
    theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
          axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm')) + 
    geom_ribbon(aes(ymin = min, ymax = max), fill = "grey70", col = "black", alpha = 0.7) + geom_line(size = 1.2) +
    labs(x ="Proportion of Food From Domestic Origins (Psi)", 
         y = "Efficacy of Curtailment (EoC)",
         col = "Heterogeneity \n in Import") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
    scale_color_manual(labels = c("Baseline", labs), 
                       values = c("red", "black"))+ 
    scale_linetype_manual(values = c(2,1)) + guides(linetype = "none")
  
  plot.with.inset <- ggdraw() +
    draw_plot(p_import_comp) +
    draw_plot(dist, x = 0.125, y = .65, width = .3, height = .3)
  
  plot_imp[[i]] <- plot.with.inset
}

p_dist <- ggarrange(plot_imp[[1]], plot_imp[[2]], nrow = 2,ncol = 1, labels = c("A", "B"), font.label = c(size = 20))

ggsave(p_dist, filename = "dist_uncert_sample.png", dpi = 300, type = "cairo", width = 7, height = 11, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")
