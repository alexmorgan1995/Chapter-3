library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data/prov_data")

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
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - 
      psi*(betaHD*Isa*Sh) - 
      psi*(1-alpha)*(betaHD*Ira*Sh) - 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh = betaHH*Ish*Sh + psi*betaHD*Isa*Sh + 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*psi*(betaHD*Ira*Sh) + 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}

# Data Import -------------------------------------------------------------

country_data_imp <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Data/FullData_2021_v1_trim.csv") #This is data for pigs 
country_data_imp$Foodborne_Carriage_2019 <- country_data_imp$Foodborne_Carriage_2019/100
country_data_imp$Corrected_Usage_18 <- country_data_imp$Corrected_Usage_18/100
country_data_imp[,12:13] <- country_data_imp[,12:13]/1000

country_data_gen <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Data/res_sales_generalfit.csv") #This is data for pigs 
country_data_gen[,13:14] <- country_data_gen[,13:14]/1000

UK_tet <- country_data_gen$scaled_sales_tet[country_data_gen$Country == "United Kingdom"]
UK_amp <- country_data_gen$scaled_sales_amp[country_data_gen$Country == "United Kingdom"]

country_data_gen <- country_data_gen[country_data_gen$num_test_amp >= 10,]

plot(country_data_gen$scaled_sales_tet, country_data_gen$propres_tet, ylim = c(0,1))
plot(country_data_gen$scaled_sales_amp, country_data_gen$propres_amp, ylim = c(0,1))


#Use the mean for the EU as the parameters (minus the UK) - only the main importers 

EU_cont <- mean(country_data_imp$Foodborne_Carriage_2019[2:10])
EU_res <- mean(country_data_imp$Prop_Amp_Res [2:10])

# Import Fitted Parameters ------------------------------------------------

post_amp <- read.csv(tail(list.files(path = "C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data/prov_data", pattern = "PART1"), 1))
MAP_amp <- map_estimate(post_amp)

# Run Model --------------------------------------------------------------

#Baseline 

parmtau <- seq(0,0.02, by = 0.001)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)

for(j in 1:3) {
  output1 <- data.frame(matrix(ncol = 8, nrow = length(parmtau)))
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(nrow = 1, ncol =8))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
               
               betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], betaHH = MAP_amp[MAP_amp$Parameter == "betaHH", 2]*0.4, tau = parmtau[i],
               betaHD = MAP_amp[MAP_amp$Parameter == "betaHD", 2], betaHI = MAP_amp[MAP_amp$Parameter == "betaHI", 2], phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
               kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
               psi = c(0.656, 1, 0.1)[j], fracimp = EU_cont, propres_imp = EU_res)
    
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau[i]
    temp[,2] <- rounding(out[nrow(out),5]) 
    temp[,3] <- rounding(out[nrow(out),6]) 
    temp[,4] <- rounding(out[nrow(out),7])
    temp[,5] <- temp[3] + temp[4]
    temp[,6] <- signif(as.numeric(temp[4]/temp[5]), digits = 3)
    temp[,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
    temp[,8] <- c("baseline","import_none", "import_90")[j]
    output1[i,] <- temp
  }
  icombhdata <- rbind(icombhdata, output1)
}

colnames(icombhdata) <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")

plotdata <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[1],],
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_base <- ggplot(plotdata, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() + 
  geom_vline(xintercept = UK_amp, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,10), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[1]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)", title = "Food from Domestic Sources = 65.6% (Baseline)") 


plotdata_imp_none <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[2],],
                          id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_imp_none <- ggplot(plotdata_imp_none, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() + 
  geom_vline(xintercept = UK_amp, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,10), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[2]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF"))  +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)", title = "Food from Domestic Sources = 90%") 

plotdata_imp_90 <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[3],],
                        id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_imp_90 <- ggplot(plotdata_imp_90, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() + 
  geom_vline(xintercept = UK_amp, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,10), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[3]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF"))  +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)", title = "Food from Domestic Sources = 10%") 

ggarrange(p_base, p_imp_none, p_imp_90, nrow  = 3, ncol = 1, common.legend = TRUE, legend = "bottom")
