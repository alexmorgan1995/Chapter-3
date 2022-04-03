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
#Import Parms
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1/betaha")
import_parms <- read.csv(tail(list.files(pattern = "gen"),1))

#No Import Parms 
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Part_1_Basic/no_import_fit/output")
no_import_parms <- read.csv(tail(list.files(pattern = "noimp"),1))

# Find MAPs ---------------------------------------------------------------

maps_import <- data.frame("Parameters" = colnames(import_parms), 
                       "MAP_Estimate" = colMeans(import_parms))

maps_noimport <- data.frame("Parameters" = colnames(no_import_parms), 
                            "MAP_Estimate" = colMeans(no_import_parms))

comb_MAP <- list(maps_import, maps_noimport)

# Model Run ---------------------------------------------------------------

parmtau <- seq(0,0.01, by = 0.0005)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

plot_list <- list()

for(j in 1:2) {
  
  MAP_amp <- comb_MAP[[j]]
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
             
             betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], tau = 0,
             betaHA = MAP_amp[MAP_amp$Parameter == "betaHA", 2],
             phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
             kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], 
             zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
             psi = 0.656, fracimp = EU_cont, propres_imp = EU_res, eta = 0.11016)

  if(j == 2) { 
    parms2["psi"] = 1
    parms2["fracimp"] = 0
    parms2["propres_imp"] = 0
    
    }
  
  tauoutput <- data.frame(matrix(nrow = length(parmtau), ncol = 7))
  
  for (i in 1:length(parmtau)) {
    parms2["tau"] = parmtau[i]
    out <-  runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms2)
    tauoutput[i,] <- c(parmtau[i],
                       (out[[2]]*(446000000))/100000,
                       (out[[3]]*(446000000))/100000,
                       ((out[[2]] + out[[3]])*(446000000))/100000,
                       (out[[1]][["Isa"]] + out[[1]][["Ira"]]), 
                       out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                       out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
  }
  
  colnames(tauoutput) <- c("tau", "InfHumans", "ResInfHumans","ICombH", "ICombA", "IResRatA","IResRat")
  tauoutput$group <-  c("gen", "pig")[j]
  
  plotdata <- melt(tauoutput,
                   id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 
  
  p_base <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
    geom_col(color = "black",position= "stack", width  = 0.0005) + scale_x_continuous(expand = c(0, 0.0005)) + 
    scale_y_continuous(limits = c(0,1), expand = c(0, 0))  + 
    geom_text(label= c(round(tauoutput$IResRat, digits = 2), rep("",length(parmtau))), vjust=-0.5, hjust = 0.05,
              position = "stack", angle = 45) +
    theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
          axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm')) + 
    scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
    labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)") + 
    geom_vline(xintercept = UK_amp_usage, size = 1.2, col = "red", lty = 2) 
  
  plot_list[[j]] <- list(tauoutput, p_base)
}

# Plotting ----------------------------------------------------------------

p_comb_import <- ggarrange(plot_list[[1]][[2]], plot_list[[2]][[2]],
                        ncol = 1, nrow = 2, labels = c("A", "B"), font.label = list(size = 20), vjust = 1.2)

ggsave(p_comb_import, filename = "baseline_run_import.png", dpi = 300, type = "cairo", width = 7, height = 8, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

