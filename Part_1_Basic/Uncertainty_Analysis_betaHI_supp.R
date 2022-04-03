library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve"); library("parallel")

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
    
    dSh = uh + rh*(Ish+Irh) - 
      psi*(betaHA*(Isa*eta)*Sh) - 
      psi*(1-alpha)*(betaHA*(Ira*eta)*Sh) - 
      (1-psi)*(betaHA*fracimp*(1-propres_imp)*Sh) - 
      (1-psi)*(1-alpha)*(betaHA*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh =  psi*betaHA*(Isa*eta)*Sh + 
      (1-psi)*(betaHA*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*psi*(betaHA*(Ira*eta)*Sh) + 
      (1-psi)*(1-alpha)*(betaHA*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    CumS = psi*betaHA*(Isa*eta)*Sh + (1-psi)*(betaHA*fracimp*(1-propres_imp)*Sh)
    CumR = (1-alpha)*psi*(betaHA*(Ira*eta)*Sh) + (1-psi)*(1-alpha)*(betaHA*fracimp*propres_imp*Sh)
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh), CumS, CumR))
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

EU_cont <- mean(rowMeans(country_data_imp[-1,24:27], na.rm = T))
EU_res <- mean(rowMeans(country_data_imp[-1,28:31], na.rm = T))

# Import Fitted Parameters ------------------------------------------------
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1/betaha")
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1/betaha")

import <- function(id) {
  data <- data.frame(matrix(ncol = 9, nrow = 0))
  
  for(i in 1:length(grep(paste0("amppigs_",id), list.files(), value = TRUE))) {
    test  <- cbind(read.csv(paste0("ABC_post_amppigs_",substitute(id),"_",i,".csv"), 
                            header = TRUE), "group" = paste0("data",i), "fit" = as.character(substitute(id)))
    data <- rbind(data, test)
  }
  return(data)
}

data <- import("gen")

map_list <- list()

for(i in 1:1) {
  temp <- data
  map_list[[i]] <- data.frame("Parameters" = colnames(temp[temp$group == "data8",][1:6]), 
                              "MAP_Estimate" = colMeans(temp[temp$group == "data8",][1:6]))
}

# Generic Heatmap Sensitivity Analysis - ICombH + ResRatio ------------------------


MAP_amp <- map_list[[1]]

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

comb_parms <- expand.grid("fracimp" = seq(0, 0.30, by = (0.30-0)/25), 
                          "propres_imp" = seq(0, 1, by = (1-0)/25))  

country_data_imp
max(country_data_imp[-1,24:27], na.rm = T)

betaHI_list <- list()

for(k in 1:3) {
  heat_list <- list()
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
             
             betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], tau = 0,
             betaHA = MAP_amp[MAP_amp$Parameter == "betaHA", 2], phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
             kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = 0, zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
             eta = c(0.05, 0.11016, 0.2)[k])
  
  for(z in 1:2) {
    
    MAP_amp <- map_list[[1]]
    tauoutput <- data.frame(matrix(nrow = length(comb_parms), ncol = 4))

    parms2["psi"] <- c(UK_food_usage, UK_food_pig_usage)[z]
    for (i in 1:nrow(comb_parms)) {
      
      parms2["fracimp"] <- comb_parms[i, 1] 
      parms2["propres_imp"] = comb_parms[i, 2]
      
      dump <- matrix(nrow = 2, ncol = 2)
      
      for(j in 1:2){
        parms2["tau"] = c(0, UK_amp_usage)[j]
        out <-  runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms2)
        dump[j, 1] <- out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]])
        dump[j, 2] <- ((out[[2]] + out[[3]])*(446000000))/100000
      }
      
      tauoutput[i,] <- c(parms2["fracimp"],
                         parms2["propres_imp"], 
                         (1 - (dump[1,1]/ dump[2,1]))*100, # relRes
                         ((dump[1,2]/ dump[2,2]) - 1)*100) # relFBD
      
      print(paste0((i/nrow(comb_parms)*100)))
      
    }
    colnames(tauoutput) <- c("FracImp", "Propres_Imp", "DecRes", "IncFBD")
    
    heat_list[[z]] <- tauoutput
  }
  betaHI_list[[k]] <- heat_list
}


# Plotting Res ----------------------------------------------------------------

plot_list1 <- list()

for(i in 1:3) { 
  
  tauoutput <- betaHI_list[[i]][[1]]
  tauoutput1 <- betaHI_list[[i]][[2]]
         
  breaks1 <- seq(0, 20, by = 1)
                   
  heat_gen <- ggplot(tauoutput, aes(FracImp, Propres_Imp, z = DecRes)) + metR::geom_contour_fill(breaks = breaks1, color = "black", size = 0.1)  + 
    labs(x = bquote("Fraction of Imports Contaminated (Cont"["Imp"]~")"), y = bquote("Fraction of Contaminated Imports Resistant"), 
         title = paste0("Dom Food Usage = 0.656, Reduc in Dom Cont = ", c("0.05","0.011016","0.2")[i]), fill = "Efficacy of Curtailment (%)") +
    scale_fill_viridis_b(breaks = breaks1, direction = 1)+
    metR::geom_text_contour(col = "white",nudge_y = 0, fontface = "bold", size = 5, breaks = breaks1, label.placer = metR::label_placer_fraction(frac = 0.5),
                            stroke = 0.05, stroke.color = "black")+ 
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
    theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
          plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
          legend.spacing.x = unit(0.4, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.6, "cm"),
          legend.key.width =  unit(2, "cm")) + geom_point(x = EU_cont, y = EU_res, size = 5, col = "red") +
    metR::geom_text_contour(col = "white",nudge_y = 0, fontface = "bold", size = 5, breaks = breaks1, label.placer = metR::label_placer_fraction(frac = 0.5),
                            stroke = 0.05, stroke.color = "black")
  
  
  heat_pig <- ggplot(tauoutput1, aes(FracImp, Propres_Imp, z = DecRes)) + metR::geom_contour_fill(breaks = breaks1, color = "black", size = 0.1)  + 
    labs(x = bquote("Fraction of Imports Contaminated (Cont"["Imp"]~")"), y = bquote("Fraction of Contaminated Imports Resistant"), 
         title = paste0("Dom Food Usage = 0.4455, Reduc in Dom Cont = ", c("0.05","0.011016","0.2")[i]), fill = "Efficacy of Curtailment (%)")  +
    scale_fill_viridis_b(breaks = breaks1, direction =  1) + 
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
    theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
          plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
          legend.spacing.x = unit(0.4, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.6, "cm"),
          legend.key.width =  unit(2, "cm")) + geom_point(x = EU_cont, y = EU_res, size = 5, col = "red") +
    metR::geom_text_contour(col = "white",nudge_y = 0, fontface = "bold", size = 5, breaks = breaks1, label.placer = metR::label_placer_fraction(frac = 0.5),
                            stroke = 0.05, stroke.color = "black")
  
  plot_list1[[i]] <- list(heat_gen, heat_pig)
}

# Isolating Res ----------------------------------------------------------

p_isol_res_656 <- plot_list1[[2]][[1]] + labs(title = expression(Domestic~Food~Usage~(psi)~"="~0.656)) 
p_isol_res_4455 <- plot_list1[[2]][[2]] + labs(title = expression(Domestic~Food~Usage~(psi)~"="~0.4455)) 
  
p_isol <- ggarrange(p_isol_res_656,p_isol_res_4455,
          ncol = 1, nrow = 2, labels = c("A","B"), font.label = list(size = 20), vjust = 1.2)

ggsave(p_isol, filename = "heat_comb_res_isol_alpha.png", dpi = 300, type = "cairo", width = 7, height = 11, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")
