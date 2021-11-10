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

#### Approximate Bayesian Computation - SMC ####

prior.non.zero<-function(par, lm.low, lm.upp){
  prod(sapply(1:7, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

#Return the sum of squares between resistance and the model output
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
           abs(tauoutput$ResPropHum[tauoutput$tau == UK_amp_usage] - 0.185),
           abs(tauoutput$ICombA[tauoutput$tau == UK_amp_usage] - UK_cont),
           abs(tauoutput$ResPropAnim[tauoutput$tau == UK_amp_usage] - UK_amp_res)))
}

# Single Particle of the Model Fit ---------------------------------------------

singlerun <- function(x, G, init.state, distanceABC, fitmodel, thetaparm, epsilon, 
                      tau_range, data, lm.low, lm.upp, w.old, sigma, res.old, N) {
  i <- 0
  m <- 0
  w.new <- 0

  while(i <= 1) {
    
    m <- m + 1
    
    if(G == 1) {
      d_betaAA <- runif(1, min = 0, max = 0.035)
      d_phi <- runif(1, min = 0, max = 0.5)
      d_kappa <- runif(1, min = 0, max = 100)
      d_alpha <- rbeta(1, 1.5, 8.5)
      d_zeta <- runif(1, 0, 0.005)
      d_betaHD <- runif(1, 0, 0.05)
      d_betaHI <- runif(1, 0, 0.01)
    } else { 
      p <- sample(seq(1,N),1,prob = w.old) # check w.old here
      par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
      d_betaAA<-par[1]
      d_phi<-par[2]
      d_kappa<-par[3]
      d_alpha<-par[4]
      d_zeta <- par[5]
      d_betaHD <- par[6]
      d_betaHI <- par[7]
    }
    
    new.parms = c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHI)
    
    if(prior.non.zero(new.parms, lm.low, lm.upp)) {
      
      thetaparm[c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHD", "betaHI")] <- new.parms
      
      dist_mod <- computeDistanceABC_ALEX(distanceABC, fitmodel, tau_range, thetaparm, init.state, data)
      
      if((dist_mod[1] <= epsilon[["dist"]][G]) && (dist_mod[2] <= epsilon[["foodH"]][G]) && (dist_mod[3] <= epsilon[["AMRH"]][G]) && 
         (dist_mod[4] <= epsilon[["foodA"]][G]) && (dist_mod[5] <= epsilon[["AMRA"]][G]) && (!is.na(dist_mod))) {
        
        if(G==1){
          w.new <- 1
        } else {
          w1 <- prod(c(sapply(c(1:3,5:7), function(b) dunif(new.parms[b], min=lm.low[b], max=lm.upp[b])),
                       dbeta(new.parms[4], 1.5, 8.5))) 
          w2 <- sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(new.parms, mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
          w.new <- w1/w2
        }
        i <- i + 1
        return(list(dist_mod, m, new.parms, w.new))
      }
    }
  }
}

# ABC-SMC Function --------------------------------------------------------

ABC_algorithm <- function(N, G, distanceABC, fitmodel, tau_range, init.state, data, epsilon, lm.low, lm.upp, thetaparm)  {
  out <- list()
  
  for(g in 1:G) {
    
    print(paste0("Generation ", g, " | Time: ", Sys.time()))
    
    if(g == 1) {
      sigma <- 0
      res.old <- 0
      w.old <- 0
    }
    
    clusterExport(cl, varlist = c("amrimp", "computeDistanceABC_ALEX", "prior.non.zero", "sum_square_diff_dist",
                                  "melt_amp_pigs", "UK_amp_res", "UK_amp_usage", "UK_cont"))
    
    particles <- parLapply(cl, 
                           1:N, 
                           singlerun, 
                           G = g, 
                           init.state = init.state,
                           distanceABC = sum_square_diff_dist,
                           fitmodel = amrimp, 
                           thetaparm = thetaparm, 
                           epsilon = epsilon,
                           tau_range = melt_amp_pigs$Usage,
                           data = melt_amp_pigs,
                           lm.low = lm.low,
                           lm.upp = lm.upp,
                           w.old = w.old, 
                           sigma = sigma, 
                           res.old = res.old,
                           N = N)
    
    dat_dist <- as.matrix(do.call(rbind, lapply(particles, "[[", 1)))
    dat_nruns <- do.call(sum, lapply(particles, "[[", 2))
    res.new <- as.matrix(do.call(rbind, lapply(particles, "[[", 3)))
    w.new <- as.matrix(do.call(rbind, lapply(particles, "[[", 4)))
    
    sigma <- cov(res.new) 
    res.old <- res.new
    w.old <- w.new/sum(w.new)
    
    out[[g]] <- list(dat_nruns, dat_dist, res.old, w.old)
    colnames(res.old) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHD", "betaHI")
    write.csv(res.old, file = paste("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/ABC_post_amppigs_",g,".csv",sep=""), row.names=FALSE)
    }
  return(out)
}

# Running the Model Fit ---------------------------------------------------
detectCores()

cl <- makeCluster(14, type="SOCK")

clusterEvalQ(cl, {c(library("rootSolve"), library("tmvtnorm"))})

test <- ABC_algorithm(N = 250,
                      G = 6,
                      distanceABC = sum_square_diff_dist, 
                      fitmodel = amrimp, 
                      tau_range = melt_amp_pigs$Usage, 
                      init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
                      data = melt_amp_pigs, 
                      epsilon = list("dist" =  c(5, 5, 5, 5, 5, 5),
                                     "foodH" = c(0.593, 0.593*0.75, 0.593*0.5, 0.593*0.3, 0.593*0.2, 0.593*0.1),
                                     "AMRH" =  c(0.185, 0.185*0.75, 0.185*0.5, 0.185*0.3, 0.185*0.2, 0.185*0.1),
                                     "foodA" = c(UK_cont, UK_cont*0.75, UK_cont*0.5, UK_cont*0.3, UK_cont*0.2, UK_cont*0.1),
                                     "AMRA" =  c(UK_amp_res, UK_amp_res*0.75, UK_amp_res*0.5, UK_amp_res*0.3, UK_amp_res*0.2, UK_amp_res*0.1)), 
                      lm.low = c(0, 0, 0, 0, 0, 0, 0), 
                      lm.upp = c(0.035, 0.5, 100, 1, 0.005, 0.05, 0.01), 
                      thetaparm = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = UK_food_usage,
                                    fracimp = EU_cont, propres_imp = EU_res, eta = 0.0534))

stopCluster(cl)

# Showing Posterior Dists -------------------------------------------------

post_dist_names <- grep("ABC_post_amppigs_",
                        list.files("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data"), value = TRUE)

setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")

post_dist <- lapply(post_dist_names, read.csv)

post_dist <- mapply(cbind, post_dist, "gen" = sapply(1:length(post_dist), function(x) paste0("gen_", x)), 
                    SIMPLIFY=F)
post_dist <- do.call("rbind", post_dist)

maps_est <- map_estimate(post_dist[post_dist$gen == tail(unique(post_dist$gen),1),][,1:8])

p_list <- list()

for(i in 1:(length(post_dist)-1)) {
  p_list[[i]] <- local ({
    name_exp <- post_dist[,c(i,8)]
    p <- ggplot(name_exp, aes(x= name_exp[,1], fill=gen)) + geom_density(alpha=.5) + 
      #geom_vline(xintercept = maps_est[i,2], size = 1.2, col = "red") +
        scale_x_continuous(expand = c(0, 0), name = colnames(post_dist)[-8][i]) + 
      scale_y_continuous(expand = c(0, 0)) +
      theme(legend.text=element_text(size=14),axis.text=element_text(size=14),
            axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
    return(p)
  })
}

