library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Figures")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

# Single Model ------------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - 
      usage_dom*(betaHA*Isa*Sh) - usage_dom*(1-alpha)*(betaHA*Ira*Sh) - 
      (1-usage_dom)*(betaHA*fracimp*(1-propres_imp)*Sh) - (1-usage_dom)*(1-alpha)*(betaHA*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh = betaHH*Ish*Sh + usage_dom*betaHA*Isa*Sh + (1-usage_dom)*(betaHA*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*usage_dom*(betaHA*Ira*Sh) + (1-usage_dom)*(1-alpha)*(betaHA*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}


# Running the Basic Model -------------------------------------------------

#I want to do a baseline run and a model run where I test sampling the uniform distribution of %Domestic and Uniform proportion of FBD in importing patches 

#Baseline 

parmtau <- seq(0,0.035, by = 0.002)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)

for(j in 1:3) {
  output1 <- data.frame(matrix(ncol = 8, nrow = length(parmtau)))
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(nrow = 1, ncol =8))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, tau = parmtau[i],
               betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, usage_dom = c(1, 0.5, 0.1)[j], fracimp = 0.5, propres_imp = 0.5)
    
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau[i]
    temp[,2] <- rounding(out[nrow(out),5]) 
    temp[,3] <- rounding(out[nrow(out),6]) 
    temp[,4] <- rounding(out[nrow(out),7])
    temp[,5] <- temp[3] + temp[4]
    temp[,6] <- signif(as.numeric(temp[4]/temp[5]), digits = 3)
    temp[,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
    temp[,8] <- c("baseline","import_50", "import_90")[j]
    output1[i,] <- temp
  }
  icombhdata <- rbind(icombhdata, output1)
}

colnames(icombhdata) <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")


plotdata <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[1],],
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_base <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = 0.0123, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,0.00005), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[1]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)") 


plotdata_imp_50 <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[2],],
                        id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_imp_50 <- ggplot(plotdata_imp_50, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = 0.0123, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,0.00005), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[2]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF"))  +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)") 

plotdata_imp_90 <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[3],],
                        id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_imp_90 <- ggplot(plotdata_imp_90, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = 0.0123, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,0.00005), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[3]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF"))  +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)") 

ggsave(p_base, filename = "baseline_tauplot.png", dpi = 300, type = "cairo", width = 9, height = 4, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
ggsave(p_imp_50, filename = "import_50_tauplot.png", dpi = 300, type = "cairo", width = 9, height = 4, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
ggsave(p_imp_90, filename = "import_90_tauplot.png", dpi = 300, type = "cairo", width = 9, height = 4, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")

# LHS - Run Full Parameters ---------------------------------------------------------------------

#First thing is to select how many times we want to sample (h) - this is based on the number of "partitions" in our distribution 

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, usage_dom = 0.6, fracimp = 0.5, propres_imp = 0.5)
#This is without tau

h <- 500
lhs <- maximinLHS(h, length(parms))

parmdetails <- data.frame("parms" = c("fracimp", "propres_imp" ,"usage_dom","betaAA" ,"betaHA" ,"betaHH" ,"phi" ,"theta",
                                      "alpha", "zeta", "rh", "ra" ,"uh" , "ua" ),
                          "lbound" = c(0, 0, 0, 0.0029, 0.000001, 0.000001,  0.00131, 0.113,
                                       0, 0.00497, 55^-1, 600^-1, 288350^-1, 2400^-1),
                          "ubound" = c(1, 1, 1, 0.29, 0.0001, 0.0001,  0.131, 11.3,
                                       1, 0.497, 0.55^-1, 6^-1, 2883.5^-1, 24^-1))

lhsscaled <- data.frame(matrix(nrow = nrow(lhs), ncol = ncol(lhs)))
colnames(lhsscaled) <- unique(parmdetails[,1])

for(i in 1:length(parmdetails[,1])) {
  lhsscaled[,i] <- lhs[,i]*(parmdetails[i,3] - parmdetails[i,2]) + parmdetails[i,2]
}

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)
modelrunlhs <- data.frame(matrix(nrow = h, ncol = nrow(parmdetails) + 4)) #+4 because I also want to keep track of 4 extra outcome measures
colnames(modelrunlhs) <- c(unique(parmdetails[,1]), "deltaFBD", "deltaRes", "relFBD", "relRes")

for(i in 1:h) {
  temptau <- matrix(nrow = 2, ncol = 2)
  parmslhs <- as.list(lhsscaled[i,])
  
  for(j in 1:2) { 
    parmslhstau <- c(parmslhs, "tau" =  c(0, 0.0123)[j])
    out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
    temp <- c((rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000,
              rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    temptau[j,] <-temp
    print(temptau)
  }
  modelrunlhs[i,] <- c(parmslhs, 
                       temptau[1,1] - temptau[2,1], 
                       temptau[1,2] - temptau[2,2],
                       temptau[1,1]/temptau[2,1],
                       temptau[1,2]/temptau[2,2])  
  
  print(paste0(round((i/h)*100, digits = 2),"%" ))
}

modelrunlhs$group <- "Uncertainty"


# LHS - Just Import Parameters ---------------------------

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497)

h <- 500
lhs <- cbind(t(replicate(h, parms)), maximinLHS(h, 3)); colnames(lhs) <- c(names(parms), "fracimp", "propres_imp" ,"usage_dom")

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)
modelrunlhs_imp <- data.frame(matrix(nrow = h, ncol = ncol(lhs) + 4)) #+4 because I also want to keep track of 4 extra outcome measures
colnames(modelrunlhs_imp) <- c(colnames(lhs), "deltaFBD", "deltaRes", "relFBD", "relRes")

for(i in 1:h) {
  temptau <- matrix(nrow = 2, ncol = 4)
  parmslhs <- as.list(lhs[i,])
  
  for(j in 1:2) { 
    parmslhstau <- c(parmslhs, "tau" =  c(0, 0.0123)[j])
    out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
    temp <- c((rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000,
              rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    temptau[j,] <-temp
  }
  modelrunlhs_imp[i,] <- c(parmslhs, 
                           temptau[1,1] - temptau[2,1], 
                           temptau[1,2] - temptau[2,2],
                           temptau[1,1]/temptau[2,1],
                           temptau[1,2]/temptau[2,2])  
  
  print(paste0(round((i/h)*100, digits = 2),"%" ))
}

modelrunlhs_imp$group <- "Uncertainty_Imp"

# Uncertainty Analysis - Baseline Runs  ---------------------------------------------------
#This section will involve me overlaying the relationship between antibiotic usage and human resistance and foodborne disease
#With the every single model run obtained from the LHS sample 

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, usage_dom = 1, fracimp = 0.5, propres_imp = 0.5)

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)

#Baseline Run 
base_anal_scat <- matrix(nrow = 2, ncol = 3)
colnames(base_anal_scat) <- c("FBD", "Res", "tau")

for(j in 1:2) {
  parmslhstau <- c(parms, "tau" =  c(0, 0.0123)[j])
  out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
  temp <- c((rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000,
            rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])),
            c(0, 0.0123)[j])
  base_anal_scat[j,] <- temp
}

baseline_scat <- c(parms,
                   "delta_fbd" = base_anal_scat[1,1] - base_anal_scat[2,1],
                   "delta_res" = base_anal_scat[1,2] - base_anal_scat[2,2],
                   "rel_fbd" = base_anal_scat[1,1]/base_anal_scat[2,1],
                   "rel_res" = base_anal_scat[1,2]/base_anal_scat[2,2],
                   "group" = "Baseline")

# Uncertainty Analysis - Scatter Plots - delta_FBD delta_res/rel_res  ---------------------------------------------------

scatter <- rbind(modelrunlhs, baseline_scat)
p_dfbd_relfbd <- ggplot(scatter, aes(x = as.numeric(deltaFBD), y = as.numeric(relFBD), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1))+ 
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  labs(x = bquote(""*Delta*"FBD"), y = "relFBD") + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title= element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_dres_relres <- ggplot(scatter, aes(x = as.numeric(deltaRes), y = as.numeric(relRes), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1))+ 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_x_continuous(limits = c(-1, 0), expand = c(0, 0)) + theme_bw() +
  labs(x = bquote(""*Delta*"Res"), y = "relRES") + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title= element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_dfbd_dres <- ggplot(scatter, aes(x = as.numeric(deltaFBD), y = as.numeric(deltaRes), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1)) + 
  scale_y_continuous(limits = c(-1, 0), expand = c(0, 0))  + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  labs(x = bquote(""*Delta*"FBD"), y = bquote(""*Delta*"Res")) + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title= element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_relfbd_relres <- ggplot(scatter, aes(x = as.numeric(relFBD), y = as.numeric(relRes), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1))+ 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))  + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  labs(x = "relFBD", y = "relRES") + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title= element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_scatter <- ggarrange(p_dfbd_relfbd, p_dres_relres, p_dfbd_dres, p_relfbd_relres,
                       nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")

ggsave(p_scatter, filename = "uncert_scatter_all.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in")

# Uncertainty Analysis - Scatter Plots - IMPORT ONLY - delta_FBD delta_res/rel_res  ---------------------------------------------------

scatter_imp <- rbind(modelrunlhs_imp, baseline_scat)

p_imp_dfbd_relfbd <- ggplot(scatter_imp, aes(x = as.numeric(deltaFBD), y = as.numeric(relFBD), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1))+ 
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  labs(x = bquote(""*Delta*"FBD"), y = "relFBD") + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title= element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_imp_dres_relres <- ggplot(scatter_imp, aes(x = as.numeric(deltaRes), y = as.numeric(relRes), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1))+ 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + scale_x_continuous(limits = c(-1, 0), expand = c(0, 0)) + theme_bw() +
  labs(x = bquote(""*Delta*"Res"), y = "relRES") + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title= element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_imp_dfbd_dres <- ggplot(scatter_imp, aes(x = as.numeric(deltaFBD), y = as.numeric(deltaRes), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1)) + 
  scale_y_continuous(limits = c(-1, 0), expand = c(0, 0))  + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  labs(x = bquote(""*Delta*"FBD"), y = bquote(""*Delta*"Res")) + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title = element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_imp_relfbd_relres <- ggplot(scatter_imp, aes(x = as.numeric(relFBD), y = as.numeric(relRes), col = group, size = group)) + geom_point() + scale_size_manual(values = c(5,1))+ 
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))  + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
  labs(x = "relFBD", y = "relRES") + scale_color_manual(values = c("red", "darkblue")) + 
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"), axis.title= element_text(size=18),axis.text = element_text(size=18), legend.text = element_text(size=18),
        legend.title = element_blank(), legend.position = "bottom", title = element_text(size=18))

p_imp_scatter <- ggarrange(p_imp_dfbd_relfbd, p_imp_dres_relres, p_imp_dfbd_dres, p_imp_relfbd_relres,
                           nrow = 2, ncol = 2, common.legend = TRUE, legend = "bottom")

ggsave(p_imp_scatter, filename = "uncert_scatter_imp.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in")


# Uncertainty Analysis - Compare the All Parameter vs Import Param --------

p_relfbd_relres <- p_relfbd_relres + labs(title = "All Parameters") + 
  scale_x_continuous(limits = c(1, 3), expand = c(0, 0))
p_imp_relfbd_relres <- p_imp_relfbd_relres + labs(title = "Only Import Parameters") + 
  scale_x_continuous(limits = c(1, 3), expand = c(0, 0))

p_compare_scatter <- ggarrange(p_relfbd_relres, p_imp_relfbd_relres,
                           nrow = 2, ncol = 1, common.legend = TRUE, legend = "bottom")

ggsave(p_compare_scatter, filename = "compare_uncert_scatter_imp.png", dpi = 300, type = "cairo", width = 8, height = 10, units = "in")

# Uncertainty Scatter - Heat Map --------------------------------------

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)
imp_comb <- expand.grid("impFBD" = seq(0,1 , by = 0.05), "imp_prop_res" = seq(0,1 , by = 0.05))
usage_vec <- seq(0.2,0.8, by = 0.2)

heatmap <- list()

for(j in 1:length(usage_vec)) {
  
  heatmap[[j]] = local({
    
    parametertest <- list()
    i = 0
    
    scendata <- data.frame(matrix(nrow = nrow(imp_comb), ncol = 7))
    
    for(i in 1:nrow(imp_comb)) {
      
      parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, 
                 betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, 
                 fracimp = imp_comb[i,1], propres_imp = imp_comb[i,2], usage_dom = usage_vec[j])
      temptau <- matrix(nrow = 2, ncol = 2)
      
      for(z in 1:2) {
        parms_temp <- c(parms, "tau" = c(0, 0.0123)[z])
        out <- data.frame(ode(y = init, func = amrimp, times = times, parms = parms_temp))
        temp <- c((rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000,
                  rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
        temptau[z,] <- temp
      }
      
      scendata[i,] <- c("deltafbd" = signif(temptau[1,1] - temptau[2,1], 3), 
                        "deltares" =  temptau[1,2] - temptau[2,2], 
                        "relfbd" = temptau[1,1]/temptau[2,1], 
                        "relres" = temptau[1,2]/temptau[2,2],
                        
                        "fracimp" = parms[["fracimp"]],
                        "propres_imp" = parms[["propres_imp"]],
                        "usage_dom" = parms[["usage_dom"]])
    }
    colnames(scendata) <- c("deltafbd", "deltares", "relfbd","relres", "fracimp", "propres_imp", "usage_dom")
    parametertest[[z]] <- scendata
    return(parametertest)
  })
}


plotheat <- list()

for(i in 1:length(usage_vec)) {
  
  plottemp <- list()
  
  for(j in 1:4) {
    scentest <- heatmap[[i]][[2]]
    if(j==1) {
      breaks1 <- seq(0, 1, by = 0.1)
      plot <- ggplot(scentest, aes(fracimp, propres_imp, z = deltafbd))+ metR::geom_contour_fill(breaks = breaks1, color = "black", size = 0.1) +
        labs(x = bquote("Proportion of Imports Infected"~"(Imp"["Inf"]*")"), y =bquote("Prop Res Imports"~"(PropRes"["Imp"]*")"),
             fill = bquote(Delta*"FBD")) + scale_fill_viridis_c(begin = unique(scentest$deltafbd), end = unique(scentest$deltafbd))
    }
    if(j==2) {
      breaks2 <- seq(-0.3, 0, by = 0.02)
      plot <- ggplot(scentest, aes(fracimp, propres_imp, z = deltares)) + metR::geom_contour_fill(breaks = breaks2, color = "black", size = 0.1) +
        scale_fill_viridis_b(breaks = breaks2,  direction = -1,begin = 1-(abs(max(scentest$deltares))/0.3), end = 1-(abs(min(scentest$deltares))/0.3), values = seq(0,1, by =0.1))+
        labs(x = bquote("Proportion of Imports Infected"~"(Imp"["Inf"]*")"), y =bquote("Prop Res Imports"~"(PropRes"["Imp"]*")"),
             fill = bquote(Delta*"RES"))
    }    
    if(j==3) {
      breaks3 <- seq(1, 1.3, by = 0.02) 
      plot <- ggplot(scentest, aes(fracimp, propres_imp, z = relfbd)) + metR::geom_contour_fill(breaks = breaks3, color = "black", size = 0.1) +
        scale_fill_viridis_b(breaks = breaks3, begin =  (min(scentest$relfbd)-1)/0.3, end = (max(scentest$relfbd)-1)/0.3, values = seq(0,1, by =0.1))+
        labs(x = bquote("Proportion of Imports Infected"~"(Imp"["Inf"]*")"), y =bquote("Prop Res Imports"~"(PropRes"["Imp"]*")"),
             fill = "RelFBD")
    }   
    if(j==4) {
      breaks4 <- seq(0.25, 1, by = 0.05)
      plot <- ggplot(scentest, aes(fracimp, propres_imp, z = relres)) + metR::geom_contour_fill(breaks = breaks4, color = "black", size = 0.1) +
        scale_fill_viridis_b(breaks = breaks4, begin = (min(scentest$relres)-0.25)/ 0.75, end = (max(scentest$relres)-0.25)/ 0.75, values = seq(0, 1, by =0.1)) +
        labs(x = bquote("Proportion of Imports Infected"~"(Imp"["Inf"]*")"), y =bquote("Prop Res Imports"~"(PropRes"["Imp"]*")"),
             fill = "RelRES")
    }
    
    plot <- plot + 
      scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
      theme(legend.position = "right", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
            plot.title = element_text(size = 18, vjust = 3, hjust = 0.1, face = "bold"),
            legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.75, "cm"),
            legend.key.width =  unit(2, "cm"))
    plottemp[[j]] <- plot
  }
  
  combplot <- ggarrange(plottemp[[1]], plottemp[[2]], plottemp[[3]], plottemp[[4]], ncol = 1, nrow = 4,
                        legend = "bottom", font.label = list(size = 25), vjust = 1.2)
  
  plotheat[[i]] <- combplot
}

for(i in 1:4) {
  ggsave(plotheat[[i]], filename = paste0("heatmap_imp",c(0.2,0.4,0.6,0.8)[i],".png"), dpi = 300, type = "cairo", width = 7, height = 16, units = "in")
}


# Uncertainty Analysis - Relationship between Tau and FBD/Res ---------------------------------------------------
#This section will involve me overlaying the relationship between antibiotic usage and human resistance and foodborne disease
#With the every single model run obtained from the LHS sample 

#Run model for baseline parameters
#Run LHS parameter set and make sure I track fro both human FDB and Res
#Overlay 

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, usage_dom = 1, fracimp = 0.5, propres_imp = 0.5)

tau <- seq(0,0.02, by = 0.0005)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)

#Baseline Run 
base_anal <- data.frame(matrix(nrow = length(tau), ncol = 3))
colnames(base_anal) <- c("FBD", "Res", "tau")

for(j in 1:length(tau)) {
  parmslhstau <- c(parms, "tau" =  tau[j])
  out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
  temp <- c(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]),
            rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])),
            tau[j])
  base_anal[j,] <- temp
  print(paste0(round((j/length(tau))*100, digits = 2),"%" ))
}

#Uncertainty Analysis


parmdetails <- data.frame("parms" = c("fracimp", "propres_imp" ,"usage_dom","betaAA" ,"betaHA" ,"betaHH" ,"phi" ,"theta",
                                      "alpha", "zeta", "rh", "ra" ,"uh" , "ua" ),
                          "lbound" = c(0, 0, 0, 0.0029, 0.000001, 0.000001,  0.00131, 0.113,
                                       0, 0.00497, 55^-1, 600^-1, 288350^-1, 2400^-1),
                          "ubound" = c(1, 1, 1, 0.29, 0.0001, 0.0001,  0.131, 11.3,
                                       1, 0.497, 0.55^-1, 6^-1, 2883.5^-1, 24^-1))
h <- 500; lhs <- maximinLHS(h, length(parms))
lhsscaled <- data.frame(matrix(nrow = nrow(lhs), ncol = ncol(lhs)))
colnames(lhsscaled) <- unique(parmdetails[,1])
for(i in 1:length(parmdetails[,1])) {
  lhsscaled[,i] <- lhs[,i]*(parmdetails[i,3] - parmdetails[i,2]) + parmdetails[i,2]
}

uncert_anal <- data.frame(matrix(nrow = h*length(tau), ncol = 4)) #+2 because I also want to keep track of 2 extra outcome measures
colnames(uncert_anal) <- c("FBD", "Res", "tau","run_n")
t <- 0

for(i in 1:h) {
  parmslhs <- as.list(lhsscaled[i,])
  
  for(j in 1:length(tau)) { 
    t <- t + 1 
    parmslhstau <- c(parmslhs, "tau" =  tau[j])
    out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
    temp <- c(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]),
              rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])),
              tau[j],
              i)
    uncert_anal[t,] <- temp
  }
  print(paste0(round((i/h)*100, digits = 2),"%" ))
}

#Combining the two 
uncert_anal$group <- "LHS Runs"; base_anal$run_n <- uncert_anal$run_n[nrow(uncert_anal)] +1; base_anal$group <- "Baseline"
combdata <- rbind(uncert_anal, base_anal)

#Plotting

p1_FBD <- ggplot(combdata, aes(x = tau, y = FBD, group = run_n, col = group, size = group, alpha = group)) + geom_line() +
  scale_y_continuous(limits = c(0, 0.0035), expand = c(0, 0)) +  scale_x_continuous(expand = c(0, 0)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
        axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15),
        legend.position = "bottom", legend.text = element_text(size=14), legend.title = element_blank()) +
  labs(title = "Uncertainty Analysis (Human FBD) with LHS runs", x ="Livestock Antibiotic Usage", y = "Human Foodborne Disease") +
  scale_color_manual(values = c("red", "black")) + scale_size_manual(values = c(1.2, 1)) + scale_alpha_manual(values = c(1, 0.2))

p2_FBD <- ggplot(combdata, aes(x = tau, y = FBD, group = run_n, col = group, size = group, alpha = group)) + geom_line() +
  scale_y_continuous(limits = c(0.00003, 0.000045), expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), 
        axis.text=element_text(size=10), axis.title = element_blank(), title = element_blank(),legend.position = "none") +
  labs(title = "Uncertainty Analysis (FBD) with LHS runs", x ="Livestock Antibiotic Usage", y = "Human Foodborne Disease") +
  scale_color_manual(values = c("red", "black")) + scale_size_manual(values = c(1.2, 1)) + scale_alpha_manual(values = c(1, 0.2))

p_FBD_comb <- p1_FBD + annotation_custom(ggplotGrob(p2_FBD), xmin = 0.0075, xmax = 0.02, 
                                         ymin = 0.002, ymax = 0.0035)

p_Res <- ggplot(combdata, aes(x = tau, y = Res, group = run_n, col = group, size = group, alpha = group)) + geom_line() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +  scale_x_continuous(expand = c(0, 0)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
        axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15),
        legend.position = "bottom", legend.text = element_text(size=14), legend.title = element_blank()) +
  labs(title = "Uncertainty Analysis (Human Res) with LHS runs", x ="Livestock Antibiotic Usage", y = "Proportion of Resistance") +
  scale_color_manual(values = c("red", "black")) + scale_size_manual(values = c(1.2, 1)) + scale_alpha_manual(values = c(1, 0.2))

unc_plot <- ggarrange(p_FBD_comb, p_Res, nrow = 2, ncol = 1,
                      align = "v", labels = c("A","B"), font.label = c(size = 20)) 

ggsave(unc_plot, filename = "uncertainty_anal.png", dpi = 300, type = "cairo", width = 10, height = 14, units = "in")
