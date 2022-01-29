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

# Run Model --------------------------------------------------------------

#Baseline 

parmtau <- seq(0,0.01, by = 0.0005)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

plot_list <- list()

for(j in 1:2) {
  
  MAP_amp <- map_list[[1]]
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
             
             betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], tau = 0,
             betaHA = MAP_amp[MAP_amp$Parameter == "betaHA", 2],
             phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
             kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], 
             zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
             psi = 0, fracimp = EU_cont, propres_imp = EU_res, eta = 0.11016)
  
  parms2[["psi"]] = c(UK_food_usage, UK_food_pig_usage)[j]
  
  
  print(parms2)
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
    geom_vline(xintercept = UK_amp_usage, alpha = 0.3, size = 2) + 
    geom_col(color = "black",position= "stack", width  = 0.0005) + scale_x_continuous(expand = c(0, 0.0005)) + 
    scale_y_continuous(limits = c(0,1), expand = c(0, 0))  + 
    geom_text(label= c(round(tauoutput$IResRat, digits = 2), rep("",length(parmtau))), vjust=-0.5, hjust = 0.05,
              position = "stack", angle = 45) +
    theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
          axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
          legend.spacing.x = unit(0.3, 'cm')) + 
    scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
    labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)") 
  
  plot_list[[j]] <- list(tauoutput, p_base)
}

p_base <- plot_list[[1]][[2]]
ggsave(p_base, filename = "baseplot_gen.png", dpi = 300, type = "cairo", width = 10, height = 6, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

p_pig <- plot_list[[2]][[2]]
ggsave(p_pig, filename = "baseplot_pigs.png", dpi = 300, type = "cairo", width = 10, height = 6, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

ggplot(tauoutput, aes(x = tau, y = ICombA*0.11)) + theme_bw() + geom_line() + 
  geom_point(aes(x = UK_amp_usage, y = UK_cont), size = 2, col = "red")
ggplot(tauoutput, aes(x = tau, y = ICombH)) + theme_bw() + geom_line()+ 
  geom_point(aes(x = UK_amp_usage, y = 0.593), size = 2, col = "red")
ggplot(tauoutput, aes(x = tau, y = IResRat)) + theme_bw() + geom_line()+ 
  geom_point(aes(x = UK_amp_usage, y = UK_hum_ampres), size = 2, col = "red")
ggplot(tauoutput, aes(x = tau, y = IResRatA)) + theme_bw() + geom_line() + scale_y_continuous(limits = c(0,1)) + scale_x_continuous(limits = c(0,0.05))+ 
  geom_point(aes(x = UK_amp_usage, y = UK_amp_res), size = 2, col = "red")

# Diagnostics Distance Measures  -------------------------------------------------------------

sum_square_diff_dist <- function(data.obs, model.obs) {
  sumsquare <- (data.obs$ResPropAnim - model.obs$ResPropAnim)^2
  return(sum(sumsquare))
}

#Compute the distances for all 3 summary statistics - this section involves running the model
computeDistanceABC_ALEX <- function(distanceABC, fitmodel, tau_range, thetaparm, init.state, data) {
  tauoutput <- data.frame(matrix(nrow = length(tau_range), ncol = 5))
  tau_range <- append(tau_range, UK_amp_usage)
  parms2 = thetaparm
  
  for (i in 1:length(tau_range)) {
    
    parms2["tau"] = tau_range[i]
    out <- runsteady(y = init.state, func = fitmodel, times = c(0, Inf), parms = parms2)
    
    tauoutput[i,] <- c(tau_range[i],
                       ((out[[2]] + out[[3]])*(446000000))/100000,
                       (out[[1]][["Isa"]] + out[[1]][["Ira"]]), 
                       out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                       out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
  }
  
  colnames(tauoutput) <- c("tau", "IncH", "ICombA", "ResPropAnim", "ResPropHum")
  
  return(c(distanceABC(data, tauoutput[(!tauoutput$tau == UK_amp_usage & !tauoutput$tau == 0),]),
           abs(tauoutput$IncH[tauoutput$tau == UK_amp_usage] - 0.593),
           abs(tauoutput$ResPropHum[tauoutput$tau == UK_amp_usage] - UK_hum_ampres),
           abs((tauoutput$ICombA[tauoutput$tau == UK_amp_usage]*parms2[["eta"]]) - UK_cont),
           abs(tauoutput$ResPropAnim[tauoutput$tau == UK_amp_usage] - UK_amp_res)))
}

MAP_amp <- map_list[[1]]

dist_mod <- computeDistanceABC_ALEX(sum_square_diff_dist, 
                                    amrimp, 
                                    melt_amp_pigs$Usage, 
                                    thetaparm = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = UK_food_usage,
                                                  betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], betaHH = MAP_amp[MAP_amp$Parameter == "betaHH", 2]*0.4, tau = parmtau[i],
                                                  betaHA = MAP_amp[MAP_amp$Parameter == "betaHA", 2], phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
                                                  kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2],
                                                  fracimp = EU_cont, propres_imp = EU_res, eta = 0.11016), 
                                    init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
                                    melt_amp_pigs)


# Other Distance Calculations ---------------------------------------------


dist_mod/c(1, 0.593, UK_hum_ampres, UK_cont, UK_amp_res)*100


thetaparm = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = UK_food_usage,
              betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], betaHH = MAP_amp[MAP_amp$Parameter == "betaHH", 2]*0.4, tau = parmtau[i],
              betaHA = MAP_amp[MAP_amp$Parameter == "betaHA", 2], phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
              kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2],
              fracimp = EU_cont, propres_imp = EU_res, eta = 0.11016)
init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)


tau_range <- append(melt_amp_pigs$Usage, UK_amp_usage)
tauoutput <- data.frame(matrix(nrow = length(tau_range), ncol = 5))
parms2 = thetaparm

for (i in 1:length(tau_range)) {
  
  parms2["tau"] = tau_range[i]
  out <- runsteady(y = init.state, func = amrimp, times = c(0, Inf), parms = parms2)
  
  tauoutput[i,] <- c(tau_range[i],
                     ((out[[2]] + out[[3]])*(446000000))/100000,
                     (out[[1]][["Isa"]] + out[[1]][["Ira"]]), 
                     out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                     out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
}

colnames(tauoutput) <- c("tau", "IncH", "ICombA", "ResPropAnim", "ResPropHum")


# Model Fit ---------------------------------------------------------------

start_time <- Sys.time()

parmtau <- seq(0,0.08, by = 0.004)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 3, nrow = 0))
ribbon_final <- data.frame()

for(j in 1:1) {
  t_data <- data
  t_data_gen8 <- t_data[t_data$group == "data8",][1:6]
  
  map_amp <- map_list[[j]]
  output1 <- data.frame(matrix(NA, nrow = length(parmtau), ncol = 3))
  output_ribbon <- data.frame()
  
  for (i in 1:length(parmtau)) {
    
    temp_ribbon <- data.frame(matrix(NA, nrow = nrow(t_data_gen8), ncol=3))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
               betaAA = map_amp["betaAA",2], tau = parmtau[i],
               betaHA = map_amp["betaHA",2],
               phi = map_amp["phi",2], 
               kappa = map_amp["kappa",2], alpha = map_amp["alpha",2], 
               zeta = map_amp["zeta",2], 
               psi = c(UK_food_usage, UK_food_pig_usage)[j], fracimp = EU_cont, propres_imp = EU_res, eta = 0.11016)
    
    out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms2)
    output1[i,1] <- parmtau[i]
    output1[i,2] <- ((out[[2]] + out[[3]])*(446000000))/100000
    output1[i,3] <- out$y["Ira"] / (out$y["Isa"] + out$y["Ira"])
    output1$group <- c("gen", "pig")[j]
    
    for(z in 1:nrow(t_data_gen8)) {
      
      parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
                 betaAA = t_data_gen8[z, "betaAA"], tau = parmtau[i],
                 betaHA = t_data_gen8[z, "betaHA"],
                 phi = t_data_gen8[z, "phi"], 
                 kappa = t_data_gen8[z, "kappa"], alpha = t_data_gen8[z, "alpha"], 
                 zeta = t_data_gen8[z, "zeta"], 
                 psi = c(UK_food_usage, UK_food_pig_usage)[j], fracimp = EU_cont, propres_imp = EU_res, eta = 0.11016)

      out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms2)
      temp_ribbon[z,1] <- parmtau[i]
      temp_ribbon[z,2] <- z
      temp_ribbon[z,3] <- out$y["Ira"] / (out$y["Isa"] + out$y["Ira"])
      temp_ribbon$group <- c("gen", "pig")[j]
      
      print(paste0(c("gen", "pig")[j], ", tau: ", temp_ribbon[z,1], ", ", (z/nrow(t_data_gen8))*100, "%"))
    }
    output_ribbon <- rbind.data.frame(output_ribbon, temp_ribbon)
  }
  icombhdata <- rbind(icombhdata, output1)
  ribbon_final <- rbind(ribbon_final, output_ribbon)
}

colnames(icombhdata)[1:4] <- c("tau", "ICombH","IResRatA","group")
colnames(ribbon_final)[1:4] <- c("tau","particle","IResRatA", "group")

HDI_ribbon <- data.frame()
for(j in 1:length(unique(ribbon_final$group))) {
  for(i in 1:length(unique(ribbon_final$tau))) {
    HDI_ribbon <- rbind(HDI_ribbon, 
                        data.frame("tau" = unique(ribbon_final$tau)[i],
                                   "lowHDI" = hdi(ribbon_final$IResRatA[ribbon_final$tau == unique(ribbon_final$tau)[i] & ribbon_final$group == unique(ribbon_final$group)[j]], credMass = 0.95)[[2]],
                                   "highHDI" = hdi(ribbon_final$IResRatA[ribbon_final$tau == unique(ribbon_final$tau)[i] & ribbon_final$group == unique(ribbon_final$group)[j]], credMass = 0.95)[[3]],
                                   "scen" = unique(ribbon_final$group)[j]))
  }
}

end_time <- Sys.time(); end_time - start_time

#Plots 

model_fit_gen <- ggplot(melt_amp_pigs, aes(x = Usage, y= ResPropAnim, color = Country))  + geom_point() + theme_bw() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.048)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Ampicillin Usage in Fattening Pigs (g/PCU)", y = "Ampicillin-Resistant Fattening Pig Carriage") +
  geom_line(data = icombhdata[icombhdata$group == "gen",], aes(x = tau, y= IResRatA), col = "purple", size = 1.1) +
  geom_ribbon(data = HDI_ribbon[HDI_ribbon$scen == "gen",],
              aes(x = tau ,ymin = lowHDI, ymax = highHDI), fill = "hotpink", alpha = 0.7, inherit.aes=FALSE) +
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"),
        legend.position = "right") + 
geom_errorbar(aes(ymin=Lower_Amp , ymax=Upper_Amp, color = Country),  size=1.01, inherit.aes =  TRUE) + 
  geom_point(x = UK_amp_usage, y = UK_amp_res, size = 5, col = "red", shape  = 22, fill = "red", alpha = 0.1)



# Combined Plot -----------------------------------------------------------

comb_gen <- ggarrange(plot_list[[1]][[2]], model_fit_gen, labels = c("A", "B"), font.label = c(size = 20),
                      nrow = 1, ncol = 2, align = "h")

ggsave(comb_gen, filename = "baseplot_andfits_gen.png", dpi = 300, type = "cairo", width = 15, height = 6, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")


# Demonstrating Outcome Measure -------------------------------------------

MAP_amp <- map_list[[1]]
parmtau <- c(seq(0,0.0025, by = 0.0005), UK_amp_usage)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
           
           betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], tau = 0,
           betaHA = MAP_amp[MAP_amp$Parameter == "betaHA", 2],
           phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
           kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], 
           zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
           psi = UK_food_usage, fracimp = EU_cont, propres_imp = EU_res, eta = 0.11016)

tauoutput_test <- data.frame(matrix(nrow = length(parmtau), ncol = 7))

for (i in 1:length(parmtau)) {
  parms2["tau"] = parmtau[i]
  out <-  runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parms2)
  tauoutput_test[i,] <- c(parmtau[i],
                     (out[[2]]*(446000000))/100000,
                     (out[[3]]*(446000000))/100000,
                     ((out[[2]] + out[[3]])*(446000000))/100000,
                     (out[[1]][["Isa"]] + out[[1]][["Ira"]]), 
                     out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                     out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
}

colnames(tauoutput_test) <- c("tau", "InfHumans", "ResInfHumans","ICombH", "ICombA", "IResRatA","IResRat")

baseline = tauoutput_test[tauoutput_test$tau == UK_amp_usage,"IResRat"]
curtailment = tauoutput_test[tauoutput_test$tau == 0,"IResRat"]

outcome_example <- ggplot(tauoutput_test, aes(x = tau, y = IResRat)) + geom_line(size = 1.2) + theme_bw() + 
  theme(legend.text=element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm")) +
  labs(x ="Ampicillin Usage in Fattening Pigs (g/PCU)", y = "Ampicillin-Resistant Human Infection") + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0.2,0.3)) +
  geom_hline(yintercept = curtailment, col= "red", linetype = "dotted", size = 1.2) +  geom_hline(yintercept = baseline, col= "red", linetype = "dotted", size = 1.2) +
  geom_point(x = 0, y = curtailment, col = "red", size = 5) + geom_point(x = UK_amp_usage, y = baseline, col = "red", size = 5) + 
  annotate(geom = "text", x =  c(0, UK_amp_usage, UK_amp_usage), y =  c(curtailment, baseline, (curtailment + baseline)/2), label = c("Curtailment", "Baseline Usage", 
                                                                                                                                      paste0(round(1-(curtailment/ baseline), digits = 3)*100, " % Relative Resistance Reduction (EoC)")),
           hjust = -0.1, vjust = 2, size = 4, col = c("red","red","purple")) + 
  annotate("segment", x = UK_amp_usage, xend = UK_amp_usage, y = baseline, yend = curtailment, colour = "purple", size=2, alpha=0.6, arrow=arrow())


ggsave(outcome_example, filename = "rel_resoutcome_example.png", dpi = 300, type = "cairo", width = 7, height = 7, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

