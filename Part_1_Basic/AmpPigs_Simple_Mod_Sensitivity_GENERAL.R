library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("cowplot"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve"); library("parallel"); library("lhs")
library("Rcpp")

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

post_dist_names <- grep("ABC_post_amppigs_gen",
                        list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1/betaha"), value = TRUE)

setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1/betaha")

post_dist <- lapply(post_dist_names, read.csv)

post_dist <- mapply(cbind, post_dist, "gen" = sapply(1:length(post_dist), function(x) paste0("gen_", x)), 
                    SIMPLIFY=F)
post_dist <- do.call("rbind", post_dist)
post_dist_gen <- post_dist[post_dist$gen == tail(unique(post_dist$gen),1),][,1:6]

maps_est <- map_estimate(post_dist_gen)
maps_est <- data.frame("Parameter" = colnames(post_dist_gen), "MAP_Estimate" = colMeans(post_dist_gen))

# Testing for Monotonicity - Identify Delta and Rel  -----------------------------------------------

#The aim of this section is to look at the relationship of changing each variable on delta_FBD and delta_res
parmdetails <- rbind(data.frame("Parameter" = "fracimp", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "propres_imp", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "psi", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "eta", "Value" = seq(0, 1, by = 1/100)),
                     
                     data.frame("Parameter" = "betaAA", "Value" = seq(0, maps_est["betaAA",2]*10, by = maps_est["betaAA",2]*10/100)),
                     data.frame("Parameter" = "betaHA", "Value" = seq(0, maps_est["betaHA",2]*10, by = maps_est["betaHA",2]*10/100)),
                     data.frame("Parameter" = "phi", "Value" = seq(0, maps_est["phi",2]*10, by = maps_est["phi",2]*10/100)),
                     data.frame("Parameter" = "kappa", "Value" = seq(0, maps_est["kappa",2]*10, by = maps_est["kappa",2]*10/100)),
                     data.frame("Parameter" = "alpha", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "zeta", "Value" = seq(0, maps_est["zeta",2]*10, by = maps_est["zeta",2]*10/100)),
                     data.frame("Parameter" = "rh", "Value" = seq(0.01, 0.55^-1, by = 0.55^-1/100)),
                     data.frame("Parameter" = "ra", "Value" = seq(0, 6^-1, by = 6^-1/100)),
                     data.frame("Parameter" = "uh", "Value" = seq(0, 2883.5^-1, by = 2883.5^-1/100)),
                     data.frame("Parameter" = "ua", "Value" = seq(0, 24^-1, by = 24^-1/100)),
                     data.frame("Parameter" = "tau", "Value" = seq(0, max(dataamp_pigs_raw[,17:20]/1000, na.rm = T), 
                                                                   by = max(dataamp_pigs_raw[,17:20]/1000, na.rm = T)/100)))

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = maps_est["betaAA",2],
           betaHA = maps_est["betaHA",2], phi = maps_est["phi",2], kappa = maps_est["kappa",2], 
           alpha = maps_est["alpha",2], zeta = maps_est["zeta",2], psi = 0.656, fracimp = EU_cont, propres_imp = EU_res, 
           eta = 0.11016, tau = UK_amp_usage)

suppplotlist <- list()

for (j in 1:length(unique(parmdetails$Parameter))) { 
  
  suppplotlist[[j]] <- local ({ 
    output <- data.frame(matrix(NA, nrow = length(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]),2]),
                                              ncol = 4))
    
    for (x in 1:length(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]),2])) { #for the individual parameter values in the sequence
        
        parmstemp <- parms
        parmstemp[as.character(unique(parmdetails[,1])[j])] <- parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]
        out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parmstemp)
        
        output[x,] <- c(((out[[2]] + out[[3]])*(446000000))/100000,
                         out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]),
                         as.numeric(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]),        #what is the parameter value used
                         as.factor(unique(parmdetails[,1])[j])) # what is the parameter explored 
      
      print(paste0("Parameter ", unique(parmdetails[,1])[j], " | ", round(x/101, digits = 2)*100,"%" ))
    }
    
    colnames(output)[1:4] <- c("ICombH","ResRatH", "ParmValue", "Parm")
    
    plotnames <- c(bquote("fracimp"~Parameter), bquote("propresimp"~Parameter), bquote(psi~Parameter), bquote(eta~Parameter),
                   bquote(beta["AA"]~Parameter), bquote(beta["HA"]~Parameter),
                   bquote(phi~Parameter), bquote(kappa~Parameter), bquote(alpha~Parameter), bquote(zeta~Parameter), bquote(r["H"]~Parameter), bquote(r["A"]~Parameter), 
                   bquote(mu["H"]~Parameter), bquote(mu["A"]~Parameter), bquote(tau~Parameter))[[j]]
    
    p1 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(ICombH))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, max(na.omit(output$ICombH)*1.1)), expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    p2 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(ResRatH))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())

    return(list(p1,p2, output))
  })
}



#deltaFBD
pFBD <- plot_grid(plot_grid(suppplotlist[[1]][[1]], suppplotlist[[2]][[1]], suppplotlist[[3]][[1]],suppplotlist[[4]][[1]], suppplotlist[[5]][[1]], 
                                  suppplotlist[[6]][[1]], suppplotlist[[7]][[1]], suppplotlist[[8]][[1]], suppplotlist[[9]][[1]], suppplotlist[[10]][[1]], suppplotlist[[11]][[1]],
                                  suppplotlist[[12]][[1]], suppplotlist[[13]][[1]], suppplotlist[[14]][[1]], suppplotlist[[15]][[1]], 
                            nrow = 5, ncol =3), scale=0.95) + 
  draw_label("ICombH", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

pres<- plot_grid(plot_grid(suppplotlist[[1]][[2]], suppplotlist[[2]][[2]], suppplotlist[[3]][[2]],suppplotlist[[4]][[2]], suppplotlist[[5]][[2]], 
                                  suppplotlist[[6]][[2]], suppplotlist[[7]][[2]], suppplotlist[[8]][[2]], suppplotlist[[9]][[2]], suppplotlist[[10]][[2]], suppplotlist[[11]][[2]],
                                  suppplotlist[[12]][[2]], suppplotlist[[13]][[2]], suppplotlist[[14]][[2]],  suppplotlist[[15]][[2]], 
                           nrow = 5, ncol =3), scale=0.95) + 
  draw_label("ResRatH", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)


ggsave(pFBD, filename = "ICombH_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")
ggsave(pres, filename = "ResRat_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")


# LHS - Run Full Parameters ---------------------------------------------------------------------

#First thing is to select how many times we want to sample (h) - this is based on the number of "partitions" in our distribution 

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = maps_est["betaAA",2],
           betaHA = maps_est["betaHA",2], phi = maps_est["phi",2], kappa = maps_est["kappa",2], 
           alpha = maps_est["alpha",2], zeta = maps_est["zeta",2], psi = 0.656, fracimp = EU_cont, propres_imp = EU_res, 
           eta = 0.11016, tau = UK_amp_usage)

#This is without tau

h <- 500
lhs <- maximinLHS(h, length(parms))

#this generates a scaling factor (0,1) - using a uniform distribution for every parameter -random values are sampled from each subsection (of h sections - going vertically)
#Uniform distribution can be transformed into any distribution using q... function (e.g qnorm) - different columns can have different distributions 
#the qnorm function - you give it a probability (from the LHS matrix) - put in the parameters of the distribution - and then it gives you the  

parmdetails <- data.frame("parms" = c("ra", "rh" ,"ua", "uh",
                                      "betaAA","betaHA" ,
                                      "phi", "kappa", "alpha", "zeta", 
                                      "psi", "fracimp" , "propres_imp", "eta", "tau" ),
                          "lbound" = c(600^-1, 55^-1, 2400^-1, 288350^-1, 
                                       parms[["betaAA"]]/10, parms[["betaHA"]]/10,  
                                       parms[["phi"]]/10, parms[["kappa"]]/10, 0, parms[["zeta"]]/10,
                                       0, 0, 0, 0, 0),
                          "ubound" = c(6^-1, 0.55^-1, 24^-1, 2883.5^-1, 
                                       parms[["betaAA"]]*10, parms[["betaHA"]]*10,  
                                       parms[["phi"]]*10, parms[["kappa"]]*10, 1, parms[["zeta"]]*10,
                                       1, 1, 1 , 1, max(dataamp_pigs_raw[,17:20]/1000, na.rm = T)))

#I want to multiply each column with the difference between the lower and upperbound for the particular parameter of interest 

lhsscaled <- data.frame(matrix(nrow = nrow(lhs), ncol = ncol(lhs)))
colnames(lhsscaled) <- unique(parmdetails[,1])

for(i in 1:length(parmdetails[,1])) {
  lhsscaled[,i] <- lhs[,i]*(parmdetails[i,3] - parmdetails[i,2]) + parmdetails[i,2]
}

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
modelrunlhs <- data.frame(matrix(nrow = h, ncol = nrow(parmdetails) + 2)) #+2 because I also want to keep track of 2 extra outcome measures
colnames(modelrunlhs) <- c(unique(parmdetails[,1]), "ICombH", "ResRat")

for(i in 1:h) {
  temptau <- matrix(nrow = 2, ncol = 2)
  parmslhs <- as.list(lhsscaled[i,])
  
  parmslhstau <- c(parmslhs, "tau" =  c(0, UK_amp_usage)[j])
  
  out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parmslhstau)
  
  modelrunlhs[i,] <- c(parmslhs, 
                       ((out[[2]] + out[[3]])*(446000000))/100000,
                       out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]]))  
  
  print(paste0(round((i/h)*100, digits = 2), "%"))
}

prcc_fbd <- pcc(modelrunlhs[,1:15], modelrunlhs[,16], nboot = 100, rank=TRUE)
prcc_res <- pcc(modelrunlhs[,1:15], modelrunlhs[,17], nboot = 100, rank=TRUE)

modelrunlhs$group <- "Uncertainty" # For later plotting


# Plotting the PRCC -------------------------------------------------------

#Plotting Delta_FBD PRCC
plotdf_fbd <- data.frame("parm" = rownames(prcc_fbd[[7]]), as.data.frame(prcc_fbd[[7]]))
#plotdf_fbd$parm <- c("ImpInf","ResPropInf", "psi", plotdf_fbd$parm[4:15])
plotdf_fbd$parm <- factor(plotdf_fbd$parm, levels = plotdf_fbd$parm) 

p_fbd <- ggplot(plotdf_fbd, aes(x = parm, y = original)) + geom_hline(yintercept = 0, col ="red", size = 1.05) + geom_point(stat = "identity", size = 3) + 
  theme_bw() + geom_errorbar(aes(ymin=min..c.i., ymax=max..c.i.), width=.1) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
                                                                  axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15)) +
  labs(title = bquote(bold("Sensitivity Analysis of " ~ "Incidence")), x ="Model Parameters", y = "PRCC") + 
  scale_x_discrete(expand = c(0, 0.7),
                   labels = c(expression(r[A]), expression(r[H]), expression(mu[H]), expression(mu[A]), 
                              expression(beta[AA]), expression(beta[HA]), 
                              expression(phi), expression(kappa), expression(alpha), 
                              expression(zeta), expression(psi), expression(Frac[Imp]),
                              expression(PropRes[Imp]), expression(eta), expression(tau)))
  
#Plotting Delta_RES PRCC
plotdf_res <- data.frame("parm" = rownames(prcc_res[[7]]), as.data.frame(prcc_res[[7]]))
#plotdf_res$parm <- c("ImpInf","ResPropInf", "psi", plotdf_res$parm[4:15])
plotdf_res$parm <- factor(plotdf_res$parm, levels = plotdf_res$parm) 

p_res <- ggplot(plotdf_res, aes(x = parm, y = original)) + geom_hline(yintercept = 0, col ="red", size = 1.05) + geom_point(stat = "identity", size = 3) + 
  theme_bw() + geom_errorbar(aes(ymin=min..c.i., ymax=max..c.i.), width=.1) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
                                                                  axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15)) +
  labs(title = bquote(bold("Sensitivity Analysis of " ~ "Human Resistance")), x ="Model Parameters", y = "PRCC") + 
  scale_x_discrete(expand = c(0, 0.7),
                   labels = c(expression(r[A]), expression(r[H]), expression(mu[H]), expression(mu[A]), 
                              expression(beta[AA]), expression(beta[HA]), 
                              expression(phi), expression(kappa), expression(alpha), 
                              expression(zeta), expression(psi), expression(Frac[Imp]),
                              expression(PropRes[Imp]), expression(eta), expression(tau)))

#Plotting Deltas
PRCC_plot <- ggarrange(p_fbd, p_res, nrow = 2, ncol = 1,
                       align = "v", labels = c("A","B"), font.label = c(size = 20)) 

ggsave(PRCC_plot, filename = "LHS_GENERAL_PRCC.png", dpi = 300, width = 10, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

# AMR Wrapper Function - for eFAST ----------------------------------------------------

ode_function_relfbd <- function(x) {
  init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms = c(ra = x[z,1], rh = x[z,2], ua = x[z,3], uh = x[z,4], 
                       betaAA = x[z,5], 
                       betaHA = x[z,6], phi = x[z,7], 
                       kappa = x[z,8], alpha = x[z,9], zeta = x[z,10], 
                       psi = x[z,11], fracimp = x[z,12], propres_imp = x[z,13], eta = x[z,14], tau = x[z,15])
    print(paste0(round(z/nrow(x), digits  = 3)*100,"%"))
    temp1 <- matrix(nrow = 2, ncol = 2)
    
      parmstemp <- parms
      out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parmstemp)

      return_vec[z]  <-  ((out[[2]] + out[[3]])*(446000000))/100000
  }
  return(return_vec)
}

ode_function_relres <- function(x) {
  init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms = c(ra = x[z,1], rh = x[z,2], ua = x[z,3], uh = x[z,4], 
              betaAA = x[z,5], 
              betaHA = x[z,6],  phi = x[z,7], 
              kappa = x[z,8], alpha = x[z,9], zeta = x[z,10], 
              psi = x[z,11], fracimp = x[z,12], propres_imp = x[z,13], eta = x[z,14], tau = x[z,15])
    print(paste0(round(z/nrow(x), digits  = 3)*100,"%"))
    temp1 <- matrix(nrow = 2, ncol = 2)
    
    parmstemp <- parms
    out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parmstemp)
    
    return_vec[z] <- out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]])
  }
  return(return_vec)
}

# eFAST - Run the Model -------------------------------------------------------------------

factors <- c("ra", "rh", "ua", "uh", 
             "betaAA",
             "betaHA",  "phi", 
             "kappa", "alpha", "zeta",
             "psi", "fracimp", "propres_imp", "eta", "tau")

start_time <- Sys.time()

testfbd <- fast99(model = ode_function_relfbd, factors = factors, n = 1000, 
                  q.arg = list(list(min=600^-1, max=6^-1), 
                               list(min=55^-1, max=0.55^-1), 
                               list(min=2400^-1, max=24^-1),
                               list(min=288350^-1, max=2883.5^-1),
                               list(min=parms[["betaAA"]]/10, max=parms[["betaAA"]]*10),
                               list(min=parms[["betaHA"]]/10, max=parms[["betaHA"]]*10),
                               list(min=parms[["phi"]]/10, max=parms[["phi"]]*10),
                               list(min=parms[["kappa"]]/10, max=parms[["kappa"]]*10),
                               list(min=0.001, max=1),
                               list(min=parms[["zeta"]]/10, max=parms[["zeta"]]*10),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max= max(dataamp_pigs_raw[,17:20]/1000, na.rm = T))))

testres <- fast99(model = ode_function_relres, factors = factors, n = 1000, 
                  q.arg = list(list(min=600^-1, max=6^-1), 
                               list(min=55^-1, max=0.55^-1), 
                               list(min=2400^-1, max=24^-1),
                               list(min=288350^-1, max=2883.5^-1),
                               list(min=parms[["betaAA"]]/10, max=parms[["betaAA"]]*10),
                               list(min=parms[["betaHA"]]/10, max=parms[["betaHA"]]*10),
                               list(min=parms[["phi"]]/10, max=parms[["phi"]]*10),
                               list(min=parms[["kappa"]]/10, max=parms[["kappa"]]*10),
                               list(min=0.001, max=1),
                               list(min=parms[["zeta"]]/10, max=parms[["zeta"]]*10),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max= max(dataamp_pigs_raw[,17:20]/1000, na.rm = T))))

end_time <- Sys.time()
end_time - start_time

par(mfrow = c(2,1))

plot(testfbd)
plot(testres)

f99_dataframe_fbd <- data.frame(X=colnames(testfbd$X), testfbd$D1/testfbd$V, 1 - testfbd$Dt/testfbd$V - testfbd$D1/testfbd$V , 1 - testfbd$Dt/testfbd$V)
colnames(f99_dataframe_fbd)[-1] <- c("first.order","higher.order", "total.order")
f99_dataframe_res <- data.frame(X=colnames(testres$X), testres$D1/testres$V, 1 - testres$Dt/testres$V - testres$D1/testres$V , 1 - testres$Dt/testres$V)
colnames(f99_dataframe_res)[-1] <- c("first.order","higher.order", "total.order")

plotdata_FBD <- melt(f99_dataframe_fbd, id.vars = "X", measure.vars = c("higher.order","first.order" ))
plotdata_res <- melt(f99_dataframe_res, id.vars = "X", measure.vars = c("higher.order","first.order"))

plotdata_FBD$X <- factor(plotdata_FBD$X, levels = reorder(unique(plotdata_FBD$X), - f99_dataframe_fbd$total.order))
plotdata_res$X <- factor(plotdata_res$X, levels = reorder(unique(plotdata_res$X), - f99_dataframe_res$total.order))

# Plotting the eFAST ------------------------------------------------------

p_efast_fbd <- ggplot(plotdata_FBD, aes(fill = variable, x =  X, y = value)) + theme_bw() + 
  geom_col(color = "black",position= "stack") + scale_y_continuous(limits = c(0,0.5), expand = c(0, 0)) +
  theme(legend.position=c(0.1, 0.7), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5)) + 
  scale_fill_manual(labels = c("Second Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = bquote(bold("eFAST first/second order sensitivity indices for Incidence")),
       x ="", y = "Sensitivity Index")  + 
  scale_x_discrete(expand = c(0, 0.6),
                   labels = c(expression(r[A]), expression(r[H]), expression(mu[H]), expression(mu[A]), 
                              expression(beta[AA]), expression(beta[HA]), 
                              expression(phi), expression(kappa), expression(alpha), 
                              expression(zeta), expression(psi), expression(Frac[Imp]),
                              expression(PropRes[Imp]), expression(eta), expression(tau)))

p_efast_res <- ggplot(plotdata_res, aes(fill = variable, x =  X, y = value)) + theme_bw() + 
  geom_col(color = "black",position= "stack") + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(legend.position=c(0.1, 0.7), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5)) + 
  scale_fill_manual(labels = c("Second Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = bquote(bold("eFAST first/second order sensitivity indices for Human Resistance")),
       x ="", y = "Sensitivity Index")  + 
  scale_x_discrete(expand = c(0, 0.6),
                   labels = c(expression(r[A]), expression(r[H]), expression(mu[H]), expression(mu[A]), 
                              expression(beta[AA]), expression(beta[HA]), 
                              expression(phi), expression(kappa), expression(alpha), 
                              expression(zeta), expression(psi), expression(Frac[Imp]),
                              expression(PropRes[Imp]), expression(eta), expression(tau)))

#Plotting rel's

PRCC_plot <- ggarrange(p_efast_fbd, p_efast_res, nrow = 2, ncol = 1,
                       align = "v", labels = c("A","B"), font.label = c(size = 20), common.legend = TRUE,
                       legend = "bottom") 

ggsave(PRCC_plot, filename = "eFAST_GENERAL.png", dpi = 300,  width = 8, height = 8, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")


tot_plot <- ggarrange(p_fbd, p_res,
                      p_efast_fbd, p_efast_res, nrow = 4, ncol = 1,
                       align = "v", labels = c("A","B","C","D"), font.label = c(size = 20), common.legend = TRUE,
                       legend = "bottom") 

ggsave(tot_plot, filename = "ses_GENERAL_ALL.png", dpi = 300,  width = 8, height = 14, units = "in", 
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

