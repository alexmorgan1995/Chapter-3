library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data")

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
      
      betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshA4+IshA5+IshA6+IshA7+IshA8+IshA9+IshAnEU+IshH) - 
      (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhA4+IrhA5+IrhA6+IrhA7+IrhA8+IrhA9+IrhAnEU+IrhH) - 
      
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
      
      (1-psi)*(imp_nEU)*(betaHI_EU*fracimp_nEU*(1-propres_impnEU)*Sh) - 
      (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_EU*fracimp_nEU*propres_impnEU*Sh)
    
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
    
    dIshAnEU = (1-psi)*(imp_nEU)*(betaHI_EU*fracimp_nEU*(1-propres_impnEU)*Sh) - rh*IshAnEU - uh*IshAnEU 
    dIrhAnEU = (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_EU*fracimp_nEU*propres_impnEU*Sh) - rh*IrhAnEU - uh*IrhAnEU  
    
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

country_data_imp <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model Fit Data/FullData_2021_v1_trim.csv")  #This is data for pigs 
country_data_imp$Foodborne_Carriage_2019 <- country_data_imp$Foodborne_Carriage_2019/100
country_data_imp$Corrected_Usage_18 <- country_data_imp$Corrected_Usage_18/100
country_data_imp[,12:13] <- country_data_imp[,12:13]/1000

country_data_gen <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model Fit Data/res_sales_generalfit.csv") #This is data for pigs 
country_data_gen[,13:14] <- country_data_gen[,13:14]/1000

UK_amp <- country_data_gen$scaled_sales_amp[country_data_gen$Country == "United Kingdom"]

country_data_gen <- country_data_gen[country_data_gen$num_test_amp >= 10,]

plot(country_data_gen$scaled_sales_amp, country_data_gen$propres_amp, ylim = c(0,1))

EU_cont <- mean(country_data_imp$Foodborne_Carriage_2019[2:10])
EU_res <- mean(country_data_imp$Prop_Amp_Res[2:10])

amp_post <- read.csv(paste0("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data/",
                            tail(grep("Simple_FIT_AMP_",list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data"), value = TRUE), 1)))

MAP_amp <- data.frame("Parameter" = colnames(amp_post), "MAP_Estimate" = colMeans(amp_post))

#Use the mean for the EU as the parameters (minus the UK) - only the main importers 

EU_cont <- mean(country_data_imp$Foodborne_Carriage_2019[2:10])
EU_res <- mean(country_data_imp$Prop_Amp_Res [2:10])

# PART 1 - Posteriors ------------------------------------------------------

amp_post <- do.call(rbind,
                    lapply(list.files(path ="//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data", pattern = "PART2"), read.csv))
amp_post$gen <- as.vector(sapply(1:(nrow(amp_post)/1000), 
                                 function(x) rep(paste0("gen_",x), 1000)))

l_amp_post <- lapply(1:length(colnames(amp_post)[-1]), function(x) melt(amp_post, id.vars = "gen", measure.vars = colnames(amp_post)[x])[,c(1,3)])

names = c(expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")),
          expression(paste("Rate of Antibiotic-Resistant to Antibiotic-Sensitive Reversion (", phi, ")")),
          expression(paste("Efficacy of Antibiotic-Mediated Animal Recovery (", kappa, ")")),
          expression(paste("Transmission-related Antibiotic Resistant Fitness Cost (", alpha, ")")),
          expression(paste("Background Infection Rate (", zeta, ")")),
          expression(paste("Rate of Domestic Animal-to-Human Transmission (", beta[HD], ")")),
          expression(paste("Rate of Domestic Human-to-Human Transmission (", beta[HH], ")")),
          expression(paste("Rate of EU Animal-to-Human Transmission (", beta[HI], ")")),
          expression(paste("Fraction of Imported non-EU Food Products Infected (", imp_nEU, ")")),
          expression(paste("Resistant Fraction of Imported non-EU Food Products Infected  (", propres_impnEU, ")")))

amp_p_list <- list()

for(i in 1:length(names)){ 
  amp_p_list[[i]] <- ggplot(l_amp_post[[i]], aes(x=value, fill= gen)) + geom_density(alpha=.5) + 
    scale_x_continuous(expand = c(0, 0), name = names[i]) + 
    scale_y_continuous(expand = c(0, 0), name = "") +
    labs(fill = NULL) + scale_fill_discrete(labels = sapply(1:length(unique(amp_post$gen)), function(x) paste0("Generation ", x)))+
    theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
}

complex_post <- ggarrange(amp_p_list[[1]], amp_p_list[[2]], amp_p_list[[3]],
          amp_p_list[[4]], amp_p_list[[5]], amp_p_list[[6]], 
          amp_p_list[[7]], amp_p_list[[8]], amp_p_list[[9]], 
          amp_p_list[[10]],
          nrow = 5, ncol = 2,
          common.legend = TRUE, legend = "bottom")

ggsave(complex_post, filename = "complex_post.png", dpi = 300, type = "cairo", width = 6, height = 8, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")


# Plot Model Output -------------------------------------------------------



