library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate"); library("viridis")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

####

amp_post <- read.csv(tail(list.files(path = "C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data", pattern = "^complexmodel_ABC_SMC_gen_amp.*?\\.csv")[1:6], 1))
MAP <- map_estimate(amp_post)
median(amp_post$betaHI_nEU)
beta_comp <- melt(amp_post, measure.vars = c("betaHI_EU", "betaHI_nEU"))

ggplot(beta_comp, aes(x = value, fill =  variable)) + geom_density(alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x ="", y = "Density", fill = "Beta Parameter") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_vline(xintercept = c(mean(amp_post$betaHI_EU), mean(amp_post$betaHI_nEU)), col = c("red", "blue"), lty = 2, size = c(1.2, 1.2)) + 
  geom_vline(xintercept = c(MAP$MAP_Estimate[MAP$Parameter == "betaHI_EU"], MAP$MAP_Estimate[MAP$Parameter == "betaHI_nEU"]), 
             col = c("red", "blue"), lty = 4, size= c(1.2, 1.2))

