library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve"); library("parallel")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")
#setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")

# Single Model ------------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      0.5*zeta*Sa*(1-alpha) - 0.5*zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + 0.5*zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + 0.5*zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - 
      psi*(betaHD*(Isa*eta)*Sh) - 
      psi*(1-alpha)*(betaHD*(Ira*eta)*Sh) - 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh =  psi*betaHD*(Isa*eta)*Sh + 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*psi*(betaHD*(Ira*eta)*Sh) + 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    CumS = psi*betaHD*(Isa*eta)*Sh + (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh)
    CumR = (1-alpha)*psi*(betaHD*(Ira*eta)*Sh) + (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh)
    
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

#Use the mean for the EU as the parameters (minus the UK) - only the main importers 

EU_cont <- mean(rowMeans(country_data_imp[-1,24:27], na.rm = T))
EU_res <- mean(rowMeans(country_data_imp[-1,28:31], na.rm = T))


# Import Fitted Parameters ------------------------------------------------
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1")

post_amp <- read.csv(tail(list.files(path = "C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1", pattern = "ABC_post"), 1))
MAP_amp <- map_estimate(post_amp)
MAP_amp <- data.frame("Parameters" = colnames(post_amp), "MAP_Estimate" = colMeans(post_amp))

# Generic Heatmap Sensitivity Analysis - ICombH + ResRatio ------------------------

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

comb_parms <- expand.grid("fracimp" = seq(0, 0.25, by = (0.25-0)/25), 
                          "propres_imp" = seq(0, 1, by = (1-0)/25))  
comb_parms[1,]
tauoutput <- data.frame(matrix(nrow = length(comb_parms), ncol = 4))

for (i in 1:nrow(comb_parms)) {
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
             
             betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], betaHH = MAP_amp[MAP_amp$Parameter == "betaHH", 2]*0.4, tau = parmtau[i], tau = 0,
             betaHD = MAP_amp[MAP_amp$Parameter == "betaHD", 2], betaHI = MAP_amp[MAP_amp$Parameter == "betaHI", 2], phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
             kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
             psi = 0.656, 
             fracimp = comb_parms[i, 1], 
             propres_imp = comb_parms[i, 2], 
             eta = 0.11016)
  
  dump <- matrix(nrow = 2, ncol = 2)
  
  for(j in 1:2){
    parms2["tau"] = c(0, UK_amp_usage)[j]
    out <-  runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms2)
    dump[j, 1] <- out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]])
    dump[j, 2] <- ((out[[2]] + out[[3]])*(446000000))/100000
  }
  
  tauoutput[i,] <- c(parms2["fracimp"],
                     parms2["propres_imp"], 
                     1 - (dump[1,1]/ dump[2,1]), # relRes
                     (dump[1,2]/ dump[2,2]) - 1) # relFBD
  
  print(paste0((i/nrow(comb_parms)*100)))
  
}

colnames(tauoutput) <- c("FracImp", "Propres_Imp", "DecRes", "IncFBD")

# Plotting ----------------------------------------------------------------


plot1 <- ggplot(tauoutput, aes(FracImp, Propres_Imp, z = DecRes)) + metR::geom_contour_fill(color = "black", size = 0.1)  + 
  labs(x = bquote("Fraction of Imports Contaminated (Frac"["Imp"]~")"), y = bquote("Fraction of Contaminated Imports Resistant (PropRes"["Imp"]~")"), title = "Baseline Import Scenario",
       fill = "% Decrease in Resistance") +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
        plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
        legend.spacing.x = unit(0.4, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.6, "cm"),
        legend.key.width =  unit(2, "cm"))


plot2 <- ggplot(tauoutput, aes(FracImp, Propres_Imp, z = IncFBD)) + metR::geom_contour_fill(color = "black", size = 0.1)  + 
  labs(x = bquote("Fraction of Imports Contaminated (Frac"["Imp"]~")"), y = bquote("Fraction of Contaminated Imports Resistant (PropRes"["Imp"]~")"), title = "Baseline Import Scenario",
       fill = "% Increase in FBD") +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
        axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
        plot.title = element_text(size = 15, vjust = 3, hjust = 0.1, face = "bold"),
        legend.spacing.x = unit(0.4, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.6, "cm"),
        legend.key.width =  unit(2, "cm"))



# Exploring the "Lump" ----------------------------------------------------


parmtau <- seq(0,0.0025, by = 0.0001)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

tauoutput <- data.frame(matrix(nrow = length(parmtau), ncol =7))

for (i in 1:length(parmtau)) {
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
             
             betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], betaHH = MAP_amp[MAP_amp$Parameter == "betaHH", 2]*0.4, tau = parmtau[i],
             betaHD = MAP_amp[MAP_amp$Parameter == "betaHD", 2], betaHI = MAP_amp[MAP_amp$Parameter == "betaHI", 2], phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
             kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
             psi = 0.656, fracimp = 0.5, propres_imp = 1, eta = 0.11016)
  
  out <-  runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms2)
  tauoutput[i,] <- c(parmtau[i],
                     (out[[2]]*(446000000))/100000,
                     (out[[3]]*(446000000))/100000,
                     ((out[[2]] + out[[3]])*(446000000))/100000,
                     (out[[1]][["Isa"]] + out[[1]][["Ira"]]), 
                     out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                     out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
}

colnames(tauoutput) <- c("tau", "InfHumans", "ResInfHumans","ICombH", "ICombA", "IResRat","IResRatA")

plotdata <- melt(tauoutput,
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = UK_amp_usage, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,2), expand = c(0, 0))  + 
  geom_text(label= c(round(tauoutput$IResRat, digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)") 


plotdata <- melt(tauoutput,
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

ggplot(tauoutput, aes(x = tau, y = ICombA)) + geom_line()  + scale_y_continuous(expand = c(0, 0.045)) 
