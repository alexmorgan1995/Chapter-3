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
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - eta*ua*Isa + 0.5*zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - eta*ua*Ira + 0.5*zeta*Sa*(1-alpha)
    
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

# Import General Data -------------------------------------------------------------
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

# Food Usage Data ------------------------------------------------------

country_data_imp <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/FullData_2021_v1_trim.csv") #This is data for pigs 
country_data_imp$Corrected_Usage_18 <- country_data_imp$Corrected_Usage_18/100
country_data_imp$Foodborne_Carriage_2019 <- country_data_imp$Foodborne_Carriage_2019/100
country_data_imp[,12:13] <- country_data_imp[,12:13]/1000

UK_amp_res <- as.numeric(rowMeans(dataamp_pigs_raw[dataamp_pigs_raw$Country == "United Kingdom",12:15], na.rm = T))
UK_amp_usage <- rowMeans(dataamp_pigs_raw[dataamp_pigs_raw$Country == "United Kingdom",17:19], na.rm = T)/1000
UK_cont <- country_data_imp$Foodborne_Carriage_2019[country_data_imp$Country_of_Origin == "UK Origin"]
UK_food_usage <- country_data_imp$Corrected_Usage_18[country_data_imp$Country_of_Origin == "UK Origin"]

#Use the mean for the EU as the parameters (minus the UK) - only the main importers 

EU_cont <- mean(country_data_imp$Foodborne_Carriage_2019[2:10])
EU_res <- mean(country_data_imp$Prop_Amp_Res[2:10])

# Import in Parameter MAPs ------------------------------------------------

amp_post <- read.csv(paste0("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data/",
                            tail(grep("Simple_FIT_AMP_",list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data"), value = TRUE), 1)))

MAP_amp <- data.frame("Parameter" = colnames(amp_post), "MAP_Estimate" = colMeans(amp_post))

# Testing for Monotonicity - Identify Delta and Rel  -----------------------------------------------

#The aim of this section is to look at the relationship of changing each variable on delta_FBD and delta_res
parmdetails <- rbind(data.frame("Parameter" = "fracimp", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "propres_imp", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "psi", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "eta", "Value" = seq(0, 1, by = 1/100)),
                     
                     data.frame("Parameter" = "betaAA", "Value" = seq(0, MAP_amp["betaAA",2]*10, by = MAP_amp["betaAA",2]*10/100)),
                     data.frame("Parameter" = "betaHI", "Value" = seq(0, MAP_amp["betaHI",2]*10, by = MAP_amp["betaHI",2]*10/100)),
                     data.frame("Parameter" = "betaHD", "Value" = seq(0, MAP_amp["betaHD",2]*10, by = MAP_amp["betaHD",2]*10/100)),
                     data.frame("Parameter" = "phi", "Value" = seq(0, MAP_amp["phi",2]*10, by = MAP_amp["phi",2]*10/100)),
                     data.frame("Parameter" = "kappa", "Value" = seq(0, MAP_amp["kappa",2]*10, by = MAP_amp["kappa",2]*10/100)),
                     data.frame("Parameter" = "alpha", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "zeta", "Value" = seq(0, MAP_amp["zeta",2]*10, by = MAP_amp["zeta",2]*10/100)),
                     data.frame("Parameter" = "rh", "Value" = seq(0.01, 0.55^-1, by = 0.55^-1/100)),
                     data.frame("Parameter" = "ra", "Value" = seq(0, 6^-1, by = 6^-1/100)),
                     data.frame("Parameter" = "uh", "Value" = seq(0, 2883.5^-1, by = 2883.5^-1/100)),
                     data.frame("Parameter" = "ua", "Value" = seq(0, 24^-1, by = 24^-1/100)))

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
tau_range <- c(0, UK_amp_usage)

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = MAP_amp["betaAA",2],
           betaHI = MAP_amp["betaHI",2], betaHD = MAP_amp["betaHD",2], phi = MAP_amp["phi",2], kappa = MAP_amp["kappa",2], alpha = 0.43, 
           zeta = MAP_amp["zeta",2], psi = 0.656, fracimp = EU_cont, propres_imp = EU_res, eta = 0.0534)

suppplotlist <- list()

for (j in 1:length(unique(parmdetails$Parameter))) { 
  
  suppplotlist[[j]] <- local ({ 
    output <- data.frame()
    
    for (x in 1:length(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]),2])) { #for the individual parameter values in the sequence

      temp <- data.frame(matrix(nrow = length(tau_range), ncol = 4))
      
      for (i in 1:length(tau_range)) {
        
        parmstemp <- c(parms, tau = tau_range[i])
        parmstemp[as.character(unique(parmdetails[,1])[j])] <- parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]
        out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parmstemp)
        
        temp[i,1] <- ((out[[2]] + out[[3]])*(446000000))/100000
        temp[i,2] <- out[[1]][["Irh"]] / (out[[1]][["Ish"]] + out[[1]][["Irh"]])
        temp[i,3] <- as.character(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]) #what is the parameter value used
        temp[i,4] <- as.character(unique(parmdetails[,1])[j]) # what is the parameter explored 
      }
      
      output <- rbind(output, stringsAsFactors = FALSE,
                      c(as.numeric(temp[1,1]), #FBD at tau = 0
                        as.numeric(temp[2,1]), #FBD at tau = 0.0123
                        as.numeric(temp[1,2]), #ResH at tau = 0
                        as.numeric(temp[2,2]), #ResH at tau = 0.0123 
                        as.numeric(temp[1,1] - temp[2,1]), #delta_FBD
                        as.numeric(temp[1,2] - temp[2,2]), #delta_Res
                        as.numeric(temp[1,1]/temp[2,1]), #rel_FBD
                        1-as.numeric(temp[1,2]/temp[2,2]), #rel_Res
                        as.numeric(temp[i,3]), #Parm Value
                        as.factor(temp[i,4]))) #Parm Name 
      
      print(paste0("Parameter ",unique(parmdetails[,1])[j], " | ", round(x/101, digits = 2)*100,"%" ))
    }
    
    colnames(output)[1:10] <- c("FBD_H_0", "FBD_H_usage","Res_H_0", "Res_H_usage", "delta_FBD", "delta_Res","rel_FBD", "rel_Res", "ParmValue", "Parm")
    print(output)
    
    plotnames <- c(bquote("fracimp"~Parameter), bquote("propresimp"~Parameter), bquote(psi~Parameter),
                   bquote(beta["AA"]~Parameter), bquote(beta["HI"]~Parameter), bquote(beta["HD"]~Parameter), bquote(beta["HH"]~Parameter),
                   bquote(phi~Parameter), bquote(kappa~Parameter), bquote(alpha~Parameter), bquote(zeta~Parameter), bquote(r["H"]~Parameter), bquote(r["A"]~Parameter), 
                   bquote(mu["H"]~Parameter), bquote(mu["A"]~Parameter))[[j]]
    
    p1 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(delta_FBD))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, max(na.omit(output$delta_FBD)*1.1)), expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    p2 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(delta_Res))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    p3 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(rel_FBD))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, max(na.omit(output$rel_FBD)*1.1)), expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    p4 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(rel_Res))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous( expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    return(list(p1,p2,p3,p4, output))
  })
}

#deltaFBD
pdelta_FBD <- plot_grid(plot_grid(suppplotlist[[1]][[1]], suppplotlist[[2]][[1]], suppplotlist[[3]][[1]],suppplotlist[[4]][[1]], suppplotlist[[5]][[1]], 
                                  suppplotlist[[6]][[1]], suppplotlist[[7]][[1]], suppplotlist[[8]][[1]], suppplotlist[[9]][[1]], suppplotlist[[10]][[1]], suppplotlist[[11]][[1]],
                                  suppplotlist[[12]][[1]], suppplotlist[[13]][[1]], suppplotlist[[14]][[1]],  suppplotlist[[15]][[1]], nrow = 5, ncol =3), scale=0.95) + 
  draw_label("delta_FBD", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

pdelta_res <- plot_grid(plot_grid(suppplotlist[[1]][[2]], suppplotlist[[2]][[2]], suppplotlist[[3]][[2]],suppplotlist[[4]][[2]], suppplotlist[[5]][[2]], 
                                  suppplotlist[[6]][[2]], suppplotlist[[7]][[2]], suppplotlist[[8]][[2]], suppplotlist[[9]][[2]], suppplotlist[[10]][[2]], suppplotlist[[11]][[2]],
                                  suppplotlist[[12]][[2]], suppplotlist[[13]][[2]], suppplotlist[[14]][[2]],  suppplotlist[[15]][[2]], nrow = 5, ncol =3), scale=0.95) + 
  draw_label("delta_Res", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

prel_FBD <- plot_grid(plot_grid(suppplotlist[[1]][[3]], suppplotlist[[2]][[3]], suppplotlist[[3]][[3]],suppplotlist[[4]][[3]], suppplotlist[[5]][[3]], 
                                  suppplotlist[[6]][[3]], suppplotlist[[7]][[3]], suppplotlist[[8]][[3]], suppplotlist[[9]][[3]], suppplotlist[[10]][[3]], suppplotlist[[11]][[3]],
                                  suppplotlist[[12]][[3]], suppplotlist[[13]][[3]], suppplotlist[[14]][[3]], suppplotlist[[15]][[3]], nrow = 5, ncol =3), scale=0.95) + 
  draw_label("rel_FBD", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

prel_res <- plot_grid(plot_grid(suppplotlist[[1]][[4]], suppplotlist[[2]][[4]], suppplotlist[[3]][[4]],suppplotlist[[4]][[4]], suppplotlist[[5]][[4]], 
                                  suppplotlist[[6]][[4]], suppplotlist[[7]][[4]], suppplotlist[[8]][[4]], suppplotlist[[9]][[4]], suppplotlist[[10]][[4]], suppplotlist[[11]][[4]],
                                  suppplotlist[[12]][[4]], suppplotlist[[13]][[4]], suppplotlist[[14]][[4]],  suppplotlist[[14]][[4]],  nrow = 5, ncol =3), scale=0.95) + 
  draw_label("rel_Res", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

ggsave(pdelta_FBD, filename = "delta_FBD_parm_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")
ggsave(pdelta_res, filename = "delta_Res_parm_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")

ggsave(prel_FBD, filename = "rel_FBD_parm_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")
ggsave(prel_res, filename = "rel_Res_parm_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")

# LHS - Run Full Parameters ---------------------------------------------------------------------

#First thing is to select how many times we want to sample (h) - this is based on the number of "partitions" in our distribution 
parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
          betaAA = MAP_amp["betaAA", 2],
          betaHD = MAP_amp["betaHD", 2], betaHI = MAP_amp["betaHI", 2], phi = MAP_amp["phi", 2], 
          kappa = MAP_amp["kappa", 2], alpha = MAP_amp["alpha", 2], zeta = MAP_amp["zeta", 2], 
          psi = 0.656, fracimp = EU_cont, propres_imp = EU_res, eta = 0.0534)

#This is without tau

h <- 500
lhs <- maximinLHS(h, length(parms))

#this generates a scaling factor (0,1) - using a uniform distribution for every parameter -random values are sampled from each subsection (of h sections - going vertically)
#Uniform distribution can be transformed into any distribution using q... function (e.g qnorm) - different columns can have different distributions 
#the qnorm function - you give it a probability (from the LHS matrix) - put in the parameters of the distribution - and then it gives you the  

parmdetails <- data.frame("parms" = c("ra", "rh" ,"ua", "uh",
                                      "betaAA" ,"betaHH" ,"betaHD" ,"betaHI",
                                      "phi", "kappa", "alpha", "zeta", 
                                      "psi", "fracimp" , "propres_imp" ),
                          "lbound" = c(600^-1, 55^-1, 2400^-1, 288350^-1, 
                                       parms[["betaAA"]]/10, parms[["betaHH"]]/10, parms[["betaHD"]]/10,  parms[["betaHI"]]/10, 
                                       parms[["phi"]]/10, parms[["kappa"]]/10, 0, parms[["zeta"]]/10,
                                       0, 0, 0),
                          "ubound" = c(6^-1, 0.55^-1, 24^-1, 2883.5^-1, 
                                        parms[["betaAA"]]*10, parms[["betaHH"]]*10, parms[["betaHD"]]*10,  parms[["betaHI"]]*10, 
                                        parms[["phi"]]*10, parms[["kappa"]]*10, 1, parms[["zeta"]]*10,
                                        1, 1, 1))

#I want to multiply each column with the difference between the lower and upperbound for the particular parameter of interest 

lhsscaled <- data.frame(matrix(nrow = nrow(lhs), ncol = ncol(lhs)))
colnames(lhsscaled) <- unique(parmdetails[,1])

for(i in 1:length(parmdetails[,1])) {
  lhsscaled[,i] <- lhs[,i]*(parmdetails[i,3] - parmdetails[i,2]) + parmdetails[i,2]
}

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
modelrunlhs <- data.frame(matrix(nrow = h, ncol = nrow(parmdetails) + 4)) #+4 because I also want to keep track of 4 extra outcome measures
colnames(modelrunlhs) <- c(unique(parmdetails[,1]), "deltaFBD", "deltaRes", "relFBD", "relRes")

for(i in 1:h) {
  temptau <- matrix(nrow = 2, ncol = 2)
  parmslhs <- as.list(lhsscaled[i,])
  
  for(j in 1:2) { 
    parmslhstau <- c(parmslhs, "tau" =  c(0, UK_amp)[j])
    
    out <- runsteady(y = init, func = amrimp, times = c(0, Inf), parms = parmslhstau)
    
    temptau[j,] <- c((rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000,
              rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))

    print(temptau)
  }
  modelrunlhs[i,] <- c(parmslhs, 
                       temptau[1,1] - temptau[2,1], 
                       temptau[1,2] - temptau[2,2],
                       temptau[1,1]/temptau[2,1],
                       1 - temptau[1,2]/temptau[2,2])  
  
  print(paste0(round((i/h)*100, digits = 2), "%"))
}

prcc_fbd <- pcc(modelrunlhs[,1:15], modelrunlhs[,16], nboot = 100, rank=TRUE)
prcc_res <- pcc(modelrunlhs[,1:15], modelrunlhs[,17], nboot = 100, rank=TRUE)

prcc_fbd_rel <- pcc(modelrunlhs[,1:15], modelrunlhs[,18], nboot = 100, rank=TRUE)
prcc_res_rel <- pcc(modelrunlhs[,1:15], modelrunlhs[,19], nboot = 100, rank=TRUE)


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
  labs(title = bquote(bold("Sensitivity Analysis of " ~ Delta["FBD"])), x ="Model Parameters", y = "PRCC")

#Plotting Delta_RES PRCC
plotdf_res <- data.frame("parm" = rownames(prcc_res[[7]]), as.data.frame(prcc_res[[7]]))
#plotdf_res$parm <- c("ImpInf","ResPropInf", "psi", plotdf_res$parm[4:15])
plotdf_res$parm <- factor(plotdf_res$parm, levels = plotdf_res$parm) 

p_res <- ggplot(plotdf_res, aes(x = parm, y = original)) + geom_hline(yintercept = 0, col ="red", size = 1.05) + geom_point(stat = "identity", size = 3) + 
  theme_bw() + geom_errorbar(aes(ymin=min..c.i., ymax=max..c.i.), width=.1) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
                                                                  axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15)) +
  labs(title = bquote(bold("Sensitivity Analysis of " ~ Delta["RES"])), x ="Model Parameters", y = "PRCC")

#Plotting Deltas
PRCC_plot <- ggarrange(p_fbd, p_res, nrow = 2, ncol = 1,
                       align = "v", labels = c("A","B"), font.label = c(size = 20)) 

ggsave(PRCC_plot, filename = "LHS_PRCC.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")

#Plotting rel_FBD PRCC
plotrel_fbd <- data.frame("parm" = rownames(prcc_fbd_rel[[7]]), as.data.frame(prcc_fbd_rel[[7]]))
#plotrel_fbd$parm <- c("ImpInf","ResPropInf", "psi", plotrel_fbd$parm[4:15])
plotrel_fbd$parm <- factor(plotrel_fbd$parm, levels = plotrel_fbd$parm) 

p_relfbd <- ggplot(plotrel_fbd, aes(x = parm, y = original)) + geom_hline(yintercept = 0, col ="red", size = 1.05) + geom_point(stat = "identity", size = 3) + 
  theme_bw() + geom_errorbar(aes(ymin=min..c.i., ymax=max..c.i.), width=.1) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
                                                                  axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15)) +
  labs(title = bquote(bold("Sensitivity Analysis of relFBD")), x ="Model Parameters", y = "PRCC")

#Plotting rel_RES PRCC
plotrel_res <- data.frame("parm" = rownames(prcc_res_rel[[7]]), as.data.frame(prcc_res_rel[[7]]))
#plotrel_res$parm <- c("ImpInf","ResPropInf", "psi", plotrel_res$parm[4:15])
plotrel_res$parm <- factor(plotrel_res$parm, levels = plotrel_res$parm) 

p_relres <- ggplot(plotrel_res, aes(x = parm, y = original)) + geom_hline(yintercept = 0, col ="red", size = 1.05) + geom_point(stat = "identity", size = 3) + 
  theme_bw() + geom_errorbar(aes(ymin=min..c.i., ymax=max..c.i.), width=.1) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
                                                                  axis.text=element_text(size=12), axis.title =element_text(size=12), title = element_text(size=12)) +
  labs(title = bquote(bold("Sensitivity Analysis of relRes")), x ="Model Parameters", y = "PRCC")

#Plotting rel's
PRCC_plot <- ggarrange(p_relfbd, p_relres, nrow = 2, ncol = 1,
                       align = "v", labels = c("A","B"), font.label = c(size = 20)) 

ggsave(PRCC_plot, filename = "LHS_PRCC_rel.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")

# AMR Wrapper Function - for eFAST ----------------------------------------------------

ode_function_relfbd <- function(x) {
  times <- seq(0,3000, by = 10) 
  init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  tau_range <- c(0, UK_amp)
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms = c(ra = x[z,1], rh = x[z,2], ua = x[z,3], uh = x[z,4], 
                       betaAA = x[z,5], betaHH = x[z,6],
                       betaHD = x[z,7], betaHI = x[z,8], phi = x[z,9], 
                       kappa = x[z,10], alpha = x[z,11], zeta = x[z,12], 
                       psi = x[z,13], fracimp = x[z,14], propres_imp = x[z,15])
    print(paste0(round(z/nrow(x), digits  = 3)*100,"%"))
    temp1 <- matrix(nrow = 2, ncol = 2)
    
    for (i in 1:length(tau_range)) {
      temp <- vector()
      parmstemp <- c(parms, tau = tau_range[i])
      out <- ode(y = init, func = amrimp, times = times, parms = parmstemp)
      temp[1] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
      temp[2] <- rounding(out[nrow(out),7])/(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))
      temp1[i,] <- temp
    }
    return_vec[z] <- as.numeric(temp1[1,1]/temp1[2,1])
  }
  return(return_vec)
  #c(temp1[1,1]/temp1[2,1], temp1[1,2]/temp1[2,2], temp1[1,1] - temp1[2,1], temp1[1,2] - temp1[2,2])
  #relFBD, relRES, deltaFBD, deltaRES
}

ode_function_relres <- function(x) {
  times <- seq(0,3000, by = 10) 
  init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  tau_range <- c(0, UK_amp)
  print(x)
  
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms = c(ra = x[z,1], rh = x[z,2], ua = x[z,3], uh = x[z,4], 
                       betaAA = x[z,5], betaHH = x[z,6],
                       betaHD = x[z,7], betaHI = x[z,8], phi = x[z,9], 
                       kappa = x[z,10], alpha = x[z,11], zeta = x[z,12], 
                       psi = x[z,13], fracimp = x[z,14], propres_imp = x[z,15])
    print(paste0(round(z/nrow(x), digits  = 3)*100,"%"))
    print(paste0(round(z/nrow(x), digits  = 3)*100,"%"))
    
    temp1 <- matrix(nrow = 2, ncol = 2)
    
    for (i in 1:length(tau_range)) {
      temp <- vector()
      parmstemp <- c(parms, tau = tau_range[i])
      out <- ode(y = init, func = amrimp, times = times, parms = parmstemp)
      temp[1] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
      temp[2] <- rounding(out[nrow(out),7])/(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))
      temp1[i,] <- temp
    }
    return_vec[z] <- 1 - (as.numeric(temp1[1,2]/temp1[2,2]))
    
    if(is.nan(temp1[1,2]/temp1[2,2])) {
      return_vec[z] <-  0
    }
    #print(paste0(temp1[1,2], " , ", temp1[2,2]))
    #print(return_vec[z])
  }
  return(return_vec)
  #c(temp1[1,1]/temp1[2,1], temp1[1,2]/temp1[2,2], temp1[1,1] - temp1[2,1], temp1[1,2] - temp1[2,2])
  #relFBD, relRES, deltaFBD, deltaRES
}

# eFAST - Run the Model -------------------------------------------------------------------

factors <- c("ra", "rh", "ua", "uh", 
             "betaAA", "betaHH",
             "betaHD", "betaHI", "phi", 
             "kappa", "alpha", "zeta",
             "psi", "fracimp", "propres_imp")

start_time <- Sys.time()

testfbd <- fast99(model = ode_function_relfbd, factors = factors, n = 1000, 
                  q.arg = list(list(min=600^-1, max=6^-1), 
                               list(min=55^-1, max=0.55^-1), 
                               list(min=2400^-1, max=24^-1),
                               list(min=288350^-1, max=2883.5^-1),
                               list(min=parms[["betaAA"]]/10, max=parms[["betaAA"]]*10),
                               list(min=parms[["betaHH"]]/10, max=parms[["betaHH"]]*10),
                               list(min=parms[["betaHD"]]/10, max=parms[["betaHD"]]*10),
                               list(min=parms[["betaHI"]]/10, max=parms[["betaHI"]]*10),
                               list(min=parms[["phi"]]/10, max=parms[["phi"]]*10),
                               list(min=parms[["kappa"]]/10, max=parms[["kappa"]]*10),
                               list(min=0.001, max=1),
                               list(min=parms[["zeta"]]/10, max=parms[["zeta"]]*10),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1)))

testres <- fast99(model = ode_function_relres, factors = factors, n = 1000, 
                  q.arg = list(list(min=600^-1, max=6^-1), 
                               list(min=55^-1, max=0.55^-1), 
                               list(min=2400^-1, max=24^-1),
                               list(min=288350^-1, max=2883.5^-1),
                               list(min=parms[["betaAA"]]/10, max=parms[["betaAA"]]*10),
                               list(min=parms[["betaHH"]]/10, max=parms[["betaHH"]]*10),
                               list(min=parms[["betaHD"]]/10, max=parms[["betaHD"]]*10),
                               list(min=parms[["betaHI"]]/10, max=parms[["betaHI"]]*10),
                               list(min=parms[["phi"]]/10, max=parms[["phi"]]*10),
                               list(min=parms[["kappa"]]/10, max=parms[["kappa"]]*10),
                               list(min=0.001, max=1),
                               list(min=parms[["zeta"]]/10, max=parms[["zeta"]]*10),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1),
                               list(min=0.001, max=1)))

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
  geom_col(color = "black",position= "stack") + scale_x_discrete(expand = c(0, 0.5)) + scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(legend.position=c(0.25, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5)) + 
  scale_fill_manual(labels = c("Total Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = bquote(bold("eFAST total/first order sensitivity indices for RelFBD")),
       x ="", y = "Sensitivity Index") 

p_efast_res <- ggplot(plotdata_res, aes(fill = variable, x =  X, y = value)) + theme_bw() + 
  geom_col(color = "black",position= "stack") + scale_x_discrete(expand = c(0, 0.7)) + scale_y_continuous(limits = c(0,0.6), expand = c(0, 0)) +
  theme(legend.position=c(0.8, 0.85), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5)) + 
  scale_fill_manual(labels = c("Total Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = bquote(bold("eFAST total/first order sensitivity indices for RelRes")),
       x ="", y = "Sensitivity Index") 

#Plotting rel's

PRCC_plot <- ggarrange(p_efast_fbd, p_efast_res, nrow = 2, ncol = 1,
                       align = "v", labels = c("A","B"), font.label = c(size = 20), common.legend = TRUE,
                       legend = "bottom") 

ggsave(PRCC_plot, filename = "eFAST_relative.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in")

#### Combined plot for the RelRes Measure ####


p_relres <- p_relres + labs(title = bquote(bold("Sensitivity Analysis of relRes")), x ="", y = "PRCC")
p_efast_res <- p_efast_res + labs(title = bquote(bold("eFAST total/first order sensitivity indices for RelRes")),
                                  x ="Model Parameters", y = "Sensitivity Index") 
comb_relres <- ggarrange(p_relres, p_efast_res, nrow = 2, ncol = 1)
  
ggsave(comb_relres, filename = "sens_PRCC_eFAST_relres.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")

