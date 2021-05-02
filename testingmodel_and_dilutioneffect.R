library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs")

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
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - 
      psi*(betaHD*Isa*Sh) - psi*(1-alpha)*(betaHD*Ira*Sh) - 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh = betaHH*Ish*Sh + psi*betaHD*Isa*Sh + (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*psi*(betaHD*Ira*Sh) + (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}


# Import in Data ----------------------------------------------------------

#Posterior Distributions
data1 <- cbind(read.csv("results_ABC_SMC_gen_tet_1.csv", header = TRUE), "group" = "data1")
data2 <- cbind(read.csv("results_ABC_SMC_gen_tet_2.csv", header = TRUE), "group" = "data2")
data3 <- cbind(read.csv("results_ABC_SMC_gen_tet_3.csv", header = TRUE), "group" = "data3")
data4 <- cbind(read.csv("results_ABC_SMC_gen_tet_4.csv", header = TRUE), "group" = "data4") 
data5 <- cbind(read.csv("results_ABC_SMC_gen_tet_5.csv", header = TRUE), "group" = "data5") 

#Import Original Data - Tet in Pigs
datatetraanim <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data/resistanceprofAnim_v2.csv")
datatetrahum <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data/resistanceprofHum_v2.csv")
datatetraanim$mgpcuuseage <- datatetraanim$mgpcuuseage / 1000; datatetrahum$mgpcuuseage <- datatetrahum$mgpcuuseage / 1000
datatetraanim$pig_tetra_sales <- datatetraanim$pig_tetra_sales / 1000; datatetrahum$pig_tetra_sales <- datatetrahum$pig_tetra_sales / 1000

datatetrahum$ResPropHum <- datatetrahum$ResPropHum/ 100 
datatetraanim <- datatetraanim[!datatetraanim$N < 10,]

mean(datatetraanim$ResPropHum ,na.rm = TRUE)
mean(datatetraanim$pig_tetra_sales ,na.rm = TRUE)

ggplot()  + geom_point(data = datatetraanim, aes(x = pig_tetra_sales, y= ResPropAnim)) +
  geom_text(data = datatetraanim, aes(x = pig_tetra_sales, y= ResPropAnim, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

ggplot()  + geom_point(data = datatetrahum, aes(x = pig_tetra_sales, y= ResPropHum)) +
  geom_text(data = datatetrahum, aes(x = pig_tetra_sales, y= ResPropHum, label = Country), vjust = -0.5, hjust = - 0.05) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.035)) + scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
  labs(x ="Livestock Antibiotic Usage (g/PCU)", y = "Antibiotic-Resistant Livestock Carriage")

m_phi <- mean(data5[,"phi"]) 
m_theta <- mean(data5[,"theta"]) 
m_betaAA <- mean(data5[,"betaAA"]) 
m_alpha <- mean(data5[,"alpha"]) 
m_zeta <- mean(data5[,"zeta"]) 
m_fracimp <- mean(data5[,"fracimp"]) 
m_propres_imp <- mean(data5[,"propres_imp"]) 


# Run Test Model ----------------------------------------------------------

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
parmtau <- seq(0,0.035, by = 0.002)

output1 <- data.frame(matrix(ncol = 7, nrow = length(parmtau)))

for (i in 1:length(parmtau)) {
  temp <- data.frame(matrix(nrow = 1, ncol =7))
  
  parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = m_betaAA, betaHH = 0.00001, tau = parmtau[i],
             betaHI = (0.00001), betaHD = (0.00001),  phi = m_phi, theta = m_theta, alpha = m_alpha, zeta = m_zeta, 
             psi = 0.7, fracimp = m_fracimp, propres_imp = m_propres_imp)
  
  out <- ode(y = init, func = amrimp, times = times, parms = parms2)
  temp[,1] <- parmtau[i]
  temp[,2] <- rounding(out[nrow(out),5]) 
  temp[,3] <- rounding(out[nrow(out),6]) 
  temp[,4] <- rounding(out[nrow(out),7])
  temp[,5] <- temp[3] + temp[4]
  temp[,6] <- signif(as.numeric(temp[4]/temp[5]), digits = 3)
  temp[,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
  output1[i,] <- temp
}

colnames(output1) <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA")

plotdata <- melt(output1, id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_base <- ggplot(plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() + 
  geom_vline(xintercept = 0.0123, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0015) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,0.00005), expand = c(0, 0))  + 
  geom_text(label= c(round(output1$IResRat,digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)") 

# Showing Dilution Effect -------------------------------------------------

#Baseline intervention - Run the previous section

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
usage <- seq(0,1, by = 0.05)
usage_frame <- data.frame(matrix(ncol = 3, nrow = length(usage)))
parmtau1 <- c(0,0.0123)

for(z in 1:length(usage)) {
  output1 <- data.frame(matrix(ncol = 3, nrow = 2))
  
  for (i in 1:2) {
    temp <- data.frame(matrix(nrow = 1, ncol =3))
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = m_betaAA, betaHH = 0.00001, tau = parmtau1[i],
               betaHI = (0.00001), betaHD = (0.00001),  phi = m_phi, theta = m_theta, alpha = m_alpha, zeta = m_zeta, 
               psi = usage[z], fracimp = m_fracimp, propres_imp = m_propres_imp)
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau1[i]
    temp[,2] <- rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])
    temp[,3] <- signif(as.numeric(rounding(out[nrow(out),7])/ temp[,2]), digits = 3)
    output1[i,] <- temp
  }
  colnames(output1) <- c("tau", "ICombH","IResRat")
  usage_frame[z,1] <- usage[z]
  usage_frame[z,2] <- ((output1$ICombH[output1$tau == 0] / output1$ICombH[output1$tau == 0.0123])-1)*100 # Increase
  usage_frame[z,3] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == 0.0123]))*100 #% Reduction
}
  
colnames(usage_frame) <- c("Usage", "RelIH","RelRes")

#Plot 

single_plot <- ggplot(usage_frame, aes(x = usage,y = RelRes)) + geom_line(size = 1.3) + theme_bw() +
  geom_vline(xintercept = 0.5, size = 1.2, col = "grey")  + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))  + 
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), title =  element_text(size=12)) + 
  labs(x ="Proportion of Domestic Food Usage", y = "% Decrease in Resistance upon Curtailment",
       title = "Proportion of Domestic Food Usage") 

ggsave(single_plot, filename = "usage_res.png", dpi = 300, type = "cairo", width = 6, height = 5, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")


# Uncertainty Analysis Dilution -------------------------------------------

#First thing is to select how many times we want to sample (h) - this is based on the number of "partitions" in our distribution 

parmdetails <- data.frame("parms" = c("ra", "ua" ,"betaAA","betaHH" ,"betaHI" ,"betaHD" ,"phi" ,"theta",
                                      "alpha", "zeta", "fracimp", "propres_imp"),
                          "lbound" = c(600^-1, 2400^-1, m_betaAA/10, 0.000001, 0.000001, 0.000001, 
                                       m_phi/10, m_theta/10, 0, m_zeta/10, 0, 0),
                          "ubound" = c(6^-1, 24^-1, m_betaAA*10, 0.0001, 0.0001, 0.0001, 
                                       m_phi*10, m_theta*10, 1, m_zeta*10, 1, 1))

h <- 500
lhs <- maximinLHS(h, nrow(parmdetails))

lhsscaled <- data.frame(matrix(nrow = nrow(lhs), ncol = ncol(lhs)))
colnames(lhsscaled) <- unique(parmdetails[,1])

for(i in 1:length(parmdetails[,1])) {
  lhsscaled[,i] <- lhs[,i]*(parmdetails[i,3] - parmdetails[i,2]) + parmdetails[i,2]
}


lhsscaled <- rbind(lhsscaled, c(60^-1, 240^-1 ,m_betaAA,0.00001,0.00001 ,0.00001 ,m_phi ,m_theta,
                                m_alpha, m_zeta, m_fracimp, m_propres_imp))

#Run the model
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0,30000, by = 100) 
parmtau1 <- c(0,0.0123)
usage <- seq(0,1, by = 0.025)

nrow(lhsscaled)*length(usage)

usage_frame2 <- data.frame(matrix(ncol = 4, nrow = 0))

for(x in 1:nrow(lhsscaled)) {
  
  dump_data <- data.frame(matrix(ncol = 4, nrow = length(usage)))
    
  for(z in 1:length(usage)) {
    
    output1 <- data.frame(matrix(ncol = 3, nrow = 2))
    
    for (i in 1:2) {
      temp <- data.frame(matrix(nrow = 1, ncol =3))
      parms2 = c(ra = lhsscaled[x,"ra"], rh = (5.5^-1), ua = lhsscaled[x,"ua"], uh = 28835^-1, betaAA = lhsscaled[x,"betaAA"], 
                 betaHH = lhsscaled[x,"betaHH"], tau = parmtau1[i],
                 betaHI = lhsscaled[x,"betaHI"], betaHD = lhsscaled[x,"betaHD"],  phi = lhsscaled[x,"phi"], 
                 theta = lhsscaled[x,"theta"], alpha = lhsscaled[x,"alpha"], zeta = lhsscaled[x,"zeta"], 
                 psi = usage[z], fracimp = lhsscaled[x,"fracimp"], propres_imp = lhsscaled[x,"propres_imp"])
      
      out <- ode(y = init, func = amrimp, times = times, parms = parms2)
      temp[,1] <- parmtau1[i]
      temp[,2] <- rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])
      temp[,3] <- signif(as.numeric(rounding(out[nrow(out),7])/ temp[,2]), digits = 3)
      output1[i,] <- temp
    }
    
    colnames(output1) <- c("tau", "ICombH","IResRat")
    
    dump_data[z,1] <- usage[z]
    dump_data[z,2] <- ((output1$ICombH[output1$tau == 0] / output1$ICombH[output1$tau == 0.0123])-1)*100 # Increase
    dump_data[z,3] <- (1-(output1$IResRat[output1$tau == 0] / output1$IResRat[output1$tau == 0.0123]))*100 #% Reduction
    dump_data[z,4] <- x 
  }
  usage_frame2 <- rbind(usage_frame2, dump_data)
  print(paste0(x/nrow(lhsscaled)*100, "%"))
}

colnames(usage_frame2) <- c("Usage", "RelIH","RelRes", "lhsnumber")

#Plot 
uncert_plot <- ggplot(usage_frame2, aes(x = Usage, y = RelRes, color = as.factor(lhsnumber), alpha = as.factor(lhsnumber))) + 
  geom_line(size = 1.3) + theme_bw() +  
  scale_size_manual(values = c(1,8)) + scale_color_manual(values = c(rep("darkgrey", length(unique(usage_frame2$lhsnumber))-1),"red")) + 
  scale_alpha_manual(values = c(rep(0.5, length(unique(usage_frame2$lhsnumber))-1),1)) +
  scale_x_continuous(expand = c(0, 0.0005)) + scale_y_continuous(limits = c(0,100),expand = c(0, 0))  + 
  theme(legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), legend.position = "none", title =  element_text(size=12)) + 
  labs(x ="Proportion of Domestic Food Usage", y = "% Decrease in Resistance upon Curtailment",
       title = "Proportion of Domestic Food Usage") 

ggsave(uncert_plot, filename = "uncert_usage.png", dpi = 300, type = "cairo", width = 6, height = 5, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")


# Distribution of Dilution Effect -----------------------------------------

dist_frame <- data.frame(matrix(nrow = nrow(lhsscaled), ncol = 2))
  
for(i in 1:nrow(lhsscaled)) {
  data <- usage_frame2[usage_frame2$lhsnumber == i,]
  dist_frame[i,1] <- i
  dist_frame[i,2] <- tail(data$Usage[data$RelRes == max(data$RelRes)],1) - tail(data$Usage[data$RelRes <= max(data$RelRes)/2],1)
}

hist(dist_frame[,2])

hist_usage <- ggplot(dist_frame, aes(X2)) + geom_histogram(bins = 20, col = "black", fill = "grey") + theme_bw() +
  scale_y_continuous(limits = c(0, 60), expand = c(0, 0)) + scale_x_continuous(expand = c(0.01, 0.01)) + 
  labs(x = "Proportion of Import Food Usage (before 50% Reduction in Intervention)", y = "",
       title = "Distribution of Acceptable Import") +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), title =  element_text(size=13),
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) 

ggsave(hist_usage, filename = "hist_uncert_usage.png", dpi = 300, type = "cairo", width = 6, height = 5, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")

