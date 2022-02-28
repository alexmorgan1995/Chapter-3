library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve"); library("fast"); library("metR"); 
library("grid"); library("gridExtra")

rm(list=ls())
#setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_2/Chapter2_Fit_Data/Final_Data")
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData/NewFit")

# Model Functions ---------------------------------------------------------


#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

#Model
amr <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    
    CumS = betaHH*Ish*Sh + betaHA*Isa*Sh
    CumR = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh)
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh), CumS, CumR))
  })
}

#### Data Import ####
dataamp <- read.csv("resistanceprofAnim_amp.csv")

dataamp$mgpcuuseage <- dataamp$mgpcuuseage / 1000
dataamp$pig_amp_sales <- dataamp$pig_amp_sales / 1000
dataamp <- dataamp[!dataamp$N < 10,]

dataamp$lower <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[1]]))
dataamp$upper <- unlist(lapply(1:nrow(dataamp), function(i) prop.test(dataamp$Positive.Sample[i],dataamp$N[i])[[6]][[2]]))

mean(dataamp$hum_res_amp ,na.rm = TRUE)
mean(dataamp$pig_amp_sales ,na.rm = TRUE)

ggplot(dataamp, aes(x = pig_amp_sales, y= ResPropAnim))  + geom_point() +
  geom_text(aes(x = pig_amp_sales, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage") +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1, inherit.aes =  TRUE)

# Extract Posterior Distributions -----------------------------------------

amp_post <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/INC_RESULTS_9.csv")
mean_parm <- data.frame("Parameter" = colnames(amp_post), "MAP_Estimate" = colMeans(amp_post))
mean_parm <- as.data.frame(map_estimate(amp_post))

mean_parm[mean_parm$Parameter == "zeta",2]

# Run the Model -----------------------------------------------------------

parmtau <- seq(0,0.02, by = 0.001)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

output1 <- data.frame(matrix(ncol = 6, nrow = length(parmtau)))
for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(nrow = 1, ncol =8))
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = mean_parm[mean_parm$Parameter == "betaAA",2], betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00042), phi = mean_parm[mean_parm$Parameter == "phi",2], kappa = mean_parm[mean_parm$Parameter == "kappa",2], 
             alpha = mean_parm[mean_parm$Parameter == "alpha",2], tau = parmtau[i],
             zeta = mean_parm[mean_parm$Parameter == "zeta",2])
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.004025057, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = 0.0003150 , phi = mean_parm[mean_parm$Parameter == "phi",2], kappa = mean_parm[mean_parm$Parameter == "kappa",2], 
             alpha = mean_parm[mean_parm$Parameter == "alpha",2], tau = parmtau[i],
             zeta = mean_parm[mean_parm$Parameter == "zeta",2])
  
      
  
  output1[i,] <-c(parmtau[i],
                  (out[[2]]*(446000000))/100000,
                  (out[[3]]*(446000000))/100000,
                  ((out[[2]] + out[[3]])*(446000000))/100000,
                  out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                  out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
}

colnames(output1) <- c("tau","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA")

plotdata <- melt(output1,
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 



ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = 0.01156391, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(expand = c(0, 0))  + 
  geom_text(label= c(round(output1$IResRat,digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Daily Incidence (per 100,000)") 



# What Parameters Can Compensate? ------------------------------------------


sensparms <- c("phi" = mean_parm[mean_parm$Parameter == "phi",2], 
               "kappa" = mean_parm[mean_parm$Parameter == "kappa",2],
               "betaAA" = mean_parm[mean_parm$Parameter == "betaAA",2],
               "alpha" = mean_parm[mean_parm$Parameter == "alpha",2],
               "zeta" = mean_parm[mean_parm$Parameter == "zeta",2],
               "tau" = 0.01156391)

parms = fast_parameters(minimum = c(600^-1, 55^-1, 2400^-1, 288350^-1, 
                                    sensparms["betaAA"]/10, 0.000001, 0.000001, 0.000042, 
                                    sensparms["phi"]/10 , sensparms["kappa"]/10, 0, sensparms["zeta"]/10), 
                        maximum = c(6^-1, 0.55^-1, 24^-1, 2883.5^-1, 
                                    sensparms["betaAA"]*10, 0.0001, 0.0001, 0.0042, 
                                    sensparms["phi"]*10 , sensparms["kappa"]*10, 1, sensparms["zeta"]*10), 
                        factor=12, names = c("ra", "rh" ,"ua", "uh", 
                                             "betaAA", "betaAH", "betaHH", "betaHA",
                                             "phi", "kappa", "alpha", "zeta"))

tauoutput <- data.frame(matrix(nrow = 0, ncol = 3))

for (j in 1:nrow(parms)) {
  parms2 = c(ra = parms$ra[j], rh = parms$rh[j], ua = parms$ua[j], uh = parms$uh[j], betaAA = parms$betaAA[j],
             betaAH = parms$betaAH[j], betaHH = parms$betaHH[j], betaHA = parms$betaHA[j], phi=parms$phi[j],
             kappa=parms$kappa[j], alpha = parms$alpha[j], tau = 0, zeta = parms$zeta[j])
  out <- runsteady(y = init, func = amr, times = c(0, Inf), parms = parms2)
  
  temp <- ((out[[2]] + out[[3]])*(446000000))/100000
  tauoutput <- rbind(tauoutput, c(temp ,abs(temp - 0.593)))
  print(j/nrow(parms))
}

colnames(tauoutput) <- c("IComb0","diff") 

tauoutput1 <- tauoutput 

tauoutput1$inc <- ((tauoutput1$IComb0 / 0.593) - 1)* 100 # % Increase from the current usage scenario
tauoutput1$inc[is.nan(tauoutput1$inc)] <- 0; neg <- tauoutput1[tauoutput1$inc < 0,] 
tauanalysis2 <- tauoutput1$inc[!is.infinite(tauoutput1$inc)]
tauanalysis2 <- tauanalysis2[tauanalysis2 < quantile(tauanalysis2, 0.99)]

#Increase

#Compensation
sensit1 <- tauanalysis2 #Creating Variable for the output variable of interest
sens1 <- sensitivity(x=sensit1, numberf=12, make.plot=T, names = c("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                   "phi", "kappa", "alpha", "zeta"))
df.equilibrium1 <- NULL; df.equilibrium1 <- data.frame(parameter=rbind("ra", "rh" ,"ua", "uh", "betaAA", "betaAH", "betaHH", "betaHA",
                                                                       "phi", "kappa", "alpha", "zeta"), value=sens1)

#Plotting

ggplot(df.equilibrium1, aes(x = reorder(parameter, -value), y = value)) + geom_bar(stat="identity", fill="lightgrey", col = "black", width  = 0.8)

p2 <- ggplot(df.equilibrium1, aes(x = reorder(parameter, -value), y = value)) + geom_bar(stat="identity", fill="lightgrey", col = "black", width  = 0.8) + theme_bw() + 
  scale_y_continuous(limits = c(0,  max(df.equilibrium1$value)*1.1), expand = c(0, 0), name = "Partial Variance") + 
  scale_x_discrete(expand = c(0, 0.7), name = "Parameter", 
                   labels = c( expression(beta[HA]), expression(r[H]),expression(r[A]), expression(alpha), 
                              expression(zeta),expression(beta[AA]), expression(kappa), expression(beta[AH]),
                              expression(mu[A]), expression(beta[HH]), expression(phi),expression(mu[H]))) +
  labs(fill = NULL, title = bquote(bold(.(Mitigating ~ Increases ~ from ~ Baseline ~ "I*"["H"] ~ "=" ~ 0.596~ per ~ "100,000")))) + 
  theme(legend.text=element_text(size=14), axis.text=element_text(size=14), plot.title = element_text(size = 15, vjust = 1.5, hjust = 0.5),
        axis.title.y=element_text(size=14), axis.title.x= element_text(size=14), plot.margin = unit(c(0.4,0.4,0.4,0.55), "cm"))


# Heat Map ----------------------------------------------------------------

parms = c(ra = 60^-1, rh =  (5.5^-1), ua = 42^-1, uh = 28835^-1, betaAA = mean_parm[mean_parm$Parameter == "betaAA",2], betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = (0.00042), phi = mean_parm[mean_parm$Parameter == "phi",2], kappa = mean_parm[mean_parm$Parameter == "kappa",2], 
           alpha = mean_parm[mean_parm$Parameter == "alpha",2], tau = 0,
           zeta = mean_parm[mean_parm$Parameter == "zeta",2])

parametertest <- list()

for(z in 1:3) {
  
  if(z == 1) {
    parameterspace <- expand.grid("betaAA" = seq(parms[["betaAA"]]*0.75, parms[["betaAA"]], by = (parms[["betaAA"]]*0.25)/25), 
                                  "betaHA" = seq(0.00042*0.75, 0.00042, by = 0.00042*0.25/25))  
  }
  if(z == 2) {
    parameterspace <- expand.grid("zeta" = seq(parms[["zeta"]]*0.75, parms[["zeta"]], by = (parms[["zeta"]]*0.25)/25), 
                                  "betaHA" = seq(0.00042*0.75, 0.00042, by = 0.00042*0.25/25))
  }
  
  if(z == 3) {
    parameterspace <- expand.grid("both" = seq(0.75, 1, by = 0.01), 
                                  "betaHA" = seq(0.00042*0.75, 0.00042, by = 0.00042*0.25/25))
  }
  
  scendata <- data.frame(matrix(nrow = nrow(parameterspace), ncol = 9))
  parms1 <- parms
  
  for(i in 1:nrow(parameterspace)) {

    #print(paste0("HeatMap - ", c("tet_pigs", "tet_broil", "amp_pigs")[j]," | ", round(i/nrow(parameterspace), digits = 3)*100, "%"))
    parms1["betaHA"] <- parameterspace[i,2] 
    
    if(z == 1) {
      parms1["betaAA"] <- parameterspace[i,1]
    }
    if(z == 2) {
      parms1["zeta"] <- parameterspace[i,1]
    }
    if(z == 3) {
      parms1["betaAA"] <- parms["betaAA"]*parameterspace[i,1]
      parms1["zeta"] <- parms["zeta"]*parameterspace[i,1]
    }
    
    out <- runsteady(y = init, func = amr, times = c(0, Inf), parms = parms1)
    
    print(parms1["betaHA"])
    
    scendata[i,] <- c("icombh" =  ((out[[2]] + out[[3]])*(446000000))/100000, 
                      "resrat" =  out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]), 
                      "betaAA" = parms[["betaAA"]], 
                      "percbetaAA" = (parms1[["betaAA"]]/ parms[["betaAA"]])*100,
                      "betaHA" = parms1[["betaHA"]],
                      "percbetaHA" = (parms1[["betaHA"]]/ parms[["betaHA"]])*100,
                      "zeta" = parms1[["zeta"]],
                      "perczeta" = (parms1[["zeta"]]/ parms[["zeta"]])*100,
                      "percdecrease_both" = parameterspace[i,1]*100)
    
    if(z == 1 | 2) {
      scendata[["percdecrease_both"]] == 1 
    }
    
    print(paste0(c("parmstet_pigs 1", "parms_amppigs 2", "parmstet_broil 3")[j], " - " ,
                 c("beta 1", "zeta 2", "both 3")[z], " - ",round(i/nrow(parameterspace), digit = 2)*100, "%"))
  }
  colnames(scendata) <- c("icombh", "resrat", "betaAA","percbetaAA", "betaHA", "percbetaHA", "zeta", "perczeta", "percdecrease")
  parametertest[[z]] <- scendata
}





plottemp <- list()

for(j in 1:3) {
  scentest <- parametertest[[j]]
  
  breaks <- c(0, 0.593, seq(0.593, max(scentest$icombh)+0.1, by = 0.1))
  
  if(j == 1) {
    plot <- ggplot(scentest, aes(percbetaAA, percbetaHA, z = icombh))
    if(i == 1) {
      plot <- plot + labs(x = bquote("% of Baseline"~beta["AA"]), y = bquote("% of Baseline"~beta["HA"]), fill = bquote("I*"["H"]), 
                          title = paste(""))
    }
    if(i == 2) {
      plot <- plot + labs(x = bquote("% of Baseline"~beta["AA"]), y = bquote("% of Baseline"~beta["HA"]), fill = bquote("I*"["H"]), 
                          title = paste(""))
    }
    if(i == 3) {
      plot <- plot + labs(x = bquote("% of Baseline"~beta["AA"]), y = bquote("% of Baseline"~beta["HA"]), fill = bquote("I*"["H"]), 
                          title = paste(""))
    }
  } 
  if(j == 2) {
    plot <- ggplot(scentest, aes(perczeta, percbetaHA, z = icombh)) + 
      labs(x = bquote("% of Baseline"~zeta), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", title = "")
  } 
  if(j == 3) {
    plot <- ggplot(scentest, aes(percdecrease, percbetaHA, z = icombh)) + 
      labs(x = bquote("% of Baseline"~beta["AA"]~"and"~ zeta), y = bquote("% of Baseline"~beta["HA"]), fill = "ICombH", title = "")
  } 
  
  plot <- plot + metR::geom_contour_fill(breaks = breaks, color = "black", size = 0.1)  + 
    geom_contour(color = "red", size = 1, breaks = 0.593, alpha = 0.8) +
    metR::geom_text_contour(col = "white",nudge_y = -0.4, fontface = "bold", size = 5, breaks = breaks, label.placement = label_placement_fraction(frac = 0.5),
                            stroke = 0.05, stroke.color = "black",) +
    scale_fill_viridis_b(breaks = breaks, direction = -1, begin = 0, end = 0.9, values = c(0, seq(0.75,1, by = 0.25/8))) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0, 0)) + theme_bw() +
    theme(legend.position = "right", legend.title = element_text(size=14), legend.text=element_text(size=12),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x = element_text(size=14),  
          plot.title = element_text(size = 18, vjust = 3, hjust = 0.1, face = "bold"),
          legend.spacing.x = unit(0.3, 'cm'), plot.margin=unit(c(0.5,0.4,0.4,0.4),"cm"), legend.key.height =unit(0.75, "cm"),
          legend.key.width =  unit(2, "cm"))
  plottemp[[j]] <- plot
}
