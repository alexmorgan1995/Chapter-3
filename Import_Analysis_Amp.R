library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
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

# Import in Parameters and Set Baseline Parms -----------------------------

amp_post <- read.csv(tail(list.files(path = "C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data", pattern = "^complexmodel_ABC_SMC_gen_amp.*?\\.csv")[1:6], 1))

MAP_parms <- data.frame("parms" = colnames(amp_post), "mean" = colMeans(amp_post))

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
         IshAnEU = 0,IrhAnEU = 0,
         IshH = 0, IrhH = 0)
times <- seq(0, 10000, by = 100) 


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
  usage_frame$relchange <- (usage_frame$relchange / usage_frame$relchange[usage_frame$domusage == 0.656]) * 100
  
  return(usage_frame)
}

# Realistic (Data Min/Max) Sensitivity ---------------------------------------------------
#We need to obtain the lowest and maximum value for the resistance and the contamination for each country 

#Contamination - 1000 for 10 contamination parameters
amp_cont_data <- data.frame(matrix(runif(1000, range(country_data_imp$Foodborne_Carriage_2019)[1],  range(country_data_imp$Foodborne_Carriage_2019)[2]), nrow = 1000, ncol = 10))
colnames(amp_cont_data) <- grep("fracimp", names(parms), value=TRUE)

#Resistance - 1000 for 10 resistance parameters
amp_res_data <- data.frame(matrix(runif(1000, range(country_data_imp$Prop_Amp_Res)[1],  range(country_data_imp$Prop_Amp_Res)[2]), nrow = 1000, ncol = 10))
colnames(amp_res_data) <- grep("propres_imp", names(parms), value=TRUE)

#% Import - 1000 for 10 import
t_data <- RandVec(a=0, b=1, s=1, n=10, m=1000, Seed=sample(1:1000, size = 1))
imp_usage_data <- data.frame(matrix(unlist(t_data), nrow=1000, byrow=TRUE), stringsAsFactors=FALSE); colnames(imp_usage_data) <- c(sapply(1:9, function(x) paste0("imp",x)), "imp_nEU")

#All the parameters varying together - 1000 for 30 parameters
comb_data <- cbind(amp_cont_data, amp_res_data, imp_usage_data) 

#Overall Dataframe
data_list <- list(amp_cont_data, amp_res_data, imp_usage_data, comb_data)

#Analyse the dataset 

data_imp_list <- list()

for(j in 1:5){
  data_imp_list[[j]] <- local ({
    
    if(j == 1) {
      output <- import_res_func(parms, init, usage_threshold, UK_amp)
      print(paste0(c("base", "cont", "res", "usage", "all")[j]))
      
    } else {
      output <- list()
      parms_list <- data_list[[j-1]]
      parms1 <- parms
      
      for(i in 1:nrow(parms_list)) {
        parms1[colnames(parms_list)] <- parms_list[i,]
        out <- import_res_func(parms1, init, usage_threshold, UK_amp)
        output[[i]] <- out 
        print(paste0(c("base", "cont", "res", "usage", "all")[j], " - ", i/nrow(parms_list)))
      }
    }
    return(output)
  })
}

#Save Output 
saveRDS(data_imp_list, file = "import_analysis_1000.rds")

data_imp_list <- readRDS("import_analysis_1000.rds")

# Plotting Realistic Sensitivity -----------------------------------------------------------

p_list_import <- list()

for(i in 1:(length(data_imp_list)-1)) {
  p_list_import[[i]] <- local ({
    
    data_base <- data_imp_list[[1]]
    data <- as.data.frame(cbind(sapply(data_imp_list[[i+1]], "[[", 1)))
    data$usage <- c(seq(0, 0.656, by = 0.02), 0.656)
    dens_plot <- melt(data, id.vars = c("usage"), measure.vars = c(head(colnames(data), -1)))
    
    p_dens <- ggplot(dens_plot, aes(x = usage, y = value, group = variable)) + geom_line(alpha = 0.05, size = 1) + theme_bw() +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x ="Proportion of Food From Domestic Origins", y = "Reduction in the Efficacy of Livestock Antibiotic Stewardship", 
           color = "Density") +
      theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
            axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
            legend.spacing.x = unit(0.3, 'cm'))  + 
      geom_line(data = data_base, aes(x = domusage, y = relchange ), color = "red", size = 1.5, lty = 3, inherit.aes = F)
    
    
    ggsave(p_dens, filename = paste0("import_sens_", c("cont","res","usage","all")[i] ,".png"), dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
           path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
    
    return(p_dens)
  })
}

# All (0,1) Sensitivity ---------------------------------------------------------

#Contamination - 1000 for 10 contamination parameters
amp_cont_data_all <- data.frame(matrix(runif(1000, 0,  1), nrow = 1000, ncol = 10))
colnames(amp_cont_data_all) <- grep("fracimp", names(parms), value=TRUE)

#Resistance - 1000 for 10 resistance parameters
amp_res_data_all <- data.frame(matrix(runif(1000, 0,  1), nrow = 1000, ncol = 10))
colnames(amp_res_data_all) <- grep("propres_imp", names(parms), value=TRUE)

#% Import - 1000 for 10 import
t_data <- RandVec(a=0, b=1, s=1, n=10, m=1000, Seed=sample(1:1000, size = 1))
imp_usage_data <- data.frame(matrix(unlist(t_data), nrow=1000, byrow=TRUE), stringsAsFactors=FALSE)
colnames(imp_usage_data) <- c(sapply(1:9, function(x) paste0("imp",x)), "imp_nEU")

#All the parameters varying together - 1000 for 30 parameters
comb_data_all <- cbind(amp_cont_data_all, amp_res_data_all, imp_usage_data) 

#Overall Dataframe
data_list_all <- list(amp_cont_data_all, amp_res_data_all, imp_usage_data, comb_data_all)

#Analyse the dataset 

data_imp_list_all <- list()

for(j in 1:5){
  data_imp_list_all[[j]] <- local ({
    
    if(j == 1) {
      output <- import_res_func(parms, init, usage_threshold, UK_amp)
      print(paste0(c("base", "cont", "res", "usage", "all")[j]))
      
    } else {
      output <- list()
      parms_list <- data_list_all[[j-1]]
      parms1 <- parms
      
      for(i in 1:nrow(parms_list)) {
        parms1[colnames(parms_list)] <- parms_list[i,]
        out <- import_res_func(parms1, init, usage_threshold, UK_amp)
        output[[i]] <- out 
        print(paste0(c("base", "cont", "res", "usage", "all")[j], " - ", i/nrow(parms_list)))
      }
    }
    return(output)
  })
}

# Plotting ALL (0,1) Sensitivity -----------------------------------------------------------

p_all_list_import <- list()

for(i in 1:(length(data_imp_list_all)-1)) {
  p_all_list_import[[i]] <- local({
    
    data_base <- data_imp_list_all[[1]]
    data <- as.data.frame(cbind(sapply(data_imp_list_all[[i+1]], "[[", 1)))
    data$usage <- c(seq(0, 0.656, by = 0.02), 0.656)
    dens_plot <- melt(data, id.vars = c("usage"), measure.vars = c(head(colnames(data), -1)))
    
    p_dens <- ggplot(dens_plot, aes(x = usage, y = value, group = variable)) + geom_line(alpha = 0.05, size = 1) + theme_bw() +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x ="Proportion of Food From Domestic Origins", y = "Reduction in the Efficacy of Livestock Antibiotic Stewardship", 
           color = "Density") +
      theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
            axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
            legend.spacing.x = unit(0.3, 'cm'))  + 
      geom_line(data = data_base, aes(x = domusage, y = relchange ), color = "red", size = 1.5, lty = 3, inherit.aes = F)
    
    ggsave(p_dens, filename = paste0("import_sens_all_", c("cont","res","usage","all")[i] ,".png"), dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
           path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
    
    return(p_dens)
  })
}

# Granular Contamination/Resistance (Incremental) Sensitivity ---------------------------

amp_cont_data_incr <- data.frame(matrix(NA, ncol = 10, nrow = 0))
amp_res_data_incr <- data.frame(matrix(NA, ncol = 10, nrow = 0))

for(i in 1:21) {
  tempcut <- rep(c(seq(0,0.2, by = 0.01))[i], 10)
  amp_cont_data_incr <- rbind(amp_cont_data_incr, tempcut)
  colnames(amp_cont_data_incr) <- grep("fracimp", names(parms), value=TRUE)
}

for(i in 1:11) {
  temp <- rep(c(seq(0,1, by = 0.1))[i], 10)
  amp_res_data_incr <- rbind(amp_res_data_incr, temp)
  colnames(amp_res_data_incr) <- grep("propres_imp", names(parms), value=TRUE)
}

data_list_incr <- list(amp_cont_data_incr, amp_res_data_incr)

#Run the Scenario Analysis
data_imp_list_incr <- list()

for(j in 1:3){
  data_imp_list_incr[[j]] <- local ({
    if(j == 1) {
      output <- import_res_func(parms, init, usage_threshold, UK_amp)
      print("base")
      
    } else {
      output <- list()
      parms_list <- data_list_incr[[j-1]]
      parms1 <- parms
      
      for(i in 1:nrow(parms_list)) {
        parms1[colnames(parms_list)] <- parms_list[i,]
        out <- import_res_func(parms1, init, usage_threshold, UK_amp)
        output[[i]] <- out 
        print(paste0(c("base", "cont", "res")[j], " - ", i/nrow(parms_list)))
      }
    }
    return(output)
  })
}

#Plot the Incremental Analysis 
p_incr_list <- list()

for(i in 1:2) {
  p_incr_list[[i]] <- local ({
    
    data <- as.data.frame(cbind(sapply(data_imp_list_incr[[i+1]], "[[", 1)))
    colnames(data) <- sapply(seq(0, ncol(data)-1, by = 1), function(x) paste0(x, "%"))
    data$Baseline <- data_imp_list_incr[[1]][,1]
    data$usage <- c(seq(0, 0.656, by = 0.02), 0.656)

    plot_cont <- melt(data, measure.vars = (head(colnames(data), -1)), id.vars = c("usage"))
    
    p_incr <- ggplot(plot_cont, aes(x = usage, y = value, col = variable, size = variable, lty = variable)) + geom_line() +
      scale_color_manual(values = c(viridis((ncol(data)-2)), "red")) + 
      scale_size_manual(values = c(rep(1, (ncol(data)-2)), 2)) + 
      scale_linetype_manual(values = c(rep(1, (ncol(data)-2)), 2)) + 
      theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
            axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
            legend.spacing.x = unit(0.3, 'cm')) + labs(x ="Proportion of Food From Domestic Origins", 
                                                       y = "Reduction in the Efficacy of Livestock Antibiotic Stewardship", 
                                                       color = c("Proportion Contaminated", "Proportion Resistant")[i]) +
      scale_x_continuous(expand = c(0, 0), limits = c(0,0.656)) + scale_y_continuous(expand = c(0, 0), limits = c(0,100)) + guides(size= "none", linetype = "none")
    
    ggsave(p_incr, filename = paste0("import_sens_incr_", c("cont","res")[i] ,".png"), dpi = 300, type = "cairo", width = 7, height = 7, units = "in",
           path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
    
    return(p_incr)
  })
}
