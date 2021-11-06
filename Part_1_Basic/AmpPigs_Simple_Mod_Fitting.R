library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("rootSolve"); library("parallel")

rm(list=ls())
#setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Data/Salmonella_Pigs")
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data")


# Single Model ------------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      0.5*zeta*Sa*(1-alpha) - 0.5*zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + 0.5*zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + 0.5*zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - 
      psi*(betaHD*Isa*Sh) - 
      psi*(1-alpha)*(betaHD*Ira*Sh) - 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh = betaHH*Ish*Sh + psi*betaHD*Isa*Sh + 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*psi*(betaHD*Ira*Sh) + 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    CumS = betaHH*Ish*Sh + psi*betaHD*Isa*Sh + (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) 
    CumR = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*psi*(betaHD*Ira*Sh) + (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh)
    
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

country_data_imp <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/FullData_2021_v1_trim.csv") #This is data for pigs 
country_data_imp$Corrected_Usage_18 <- country_data_imp$Corrected_Usage_18/100
country_data_imp$Foodborne_Carriage_2019 <- country_data_imp$Foodborne_Carriage_2019/100
country_data_imp[,12:13] <- country_data_imp[,12:13]/1000

UK_amp_res <- rowMeans(dataamp_pigs_raw[dataamp_pigs_raw$Country == "United Kingdom",12:15], na.rm = T)
UK_amp_usage <- rowMeans(dataamp_pigs_raw[dataamp_pigs_raw$Country == "United Kingdom",17:19], na.rm = T)/1000
UK_cont <- country_data_imp$Foodborne_Carriage_2019[country_data_imp$Country_of_Origin == "UK Origin"]
UK_food_usage <- country_data_imp$Corrected_Usage_18[country_data_imp$Country_of_Origin == "UK Origin"]

#Use the mean for the EU as the parameters (minus the UK) - only the main importers 

EU_cont <- mean(country_data_imp$Foodborne_Carriage_2019[2:10])
EU_res <- mean(country_data_imp$Prop_Amp_Res[2:10])

#### Approximate Bayesian Computation - SMC ####

prior.non.zero<-function(par, lm.low, lm.upp){
  prod(sapply(1:8, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
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
           abs(tauoutput$ICombH[tauoutput$tau == UK_amp_usage] - 0.593),
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
      d_phi <- runif(1, min = 0, max = 0.1)
      d_kappa <- runif(1, min = 0, max = 10)
      d_alpha <- rbeta(1, 1.5, 8.5)
      d_zeta <- runif(1, 0, 0.005)
      d_betaHD <- runif(1, 0, 0.001)
      d_betaHH <- runif(1, 0, 0.01)
      d_betaHI <- runif(1, 0, 0.0002)
    } else { 
      p <- sample(seq(1,N),1,prob = w.old) # check w.old here
      par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
      d_betaAA<-par[1]
      d_phi<-par[2]
      d_kappa<-par[3]
      d_alpha<-par[4]
      d_zeta <- par[5]
      d_betaHD <- par[6]
      d_betaHH <- par[7]
      d_betaHI <- par[8]
    }
    
    new.parms = c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHH, d_betaHI)
    
  
    if(prior.non.zero(new.parms, lm.low, lm.upp)) {
      
      thetaparm[c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHD", "betaHH", "betaHI")] <- new.parms
      
      dist <- computeDistanceABC_ALEX(distanceABC, fitmodel, tau_range, thetaparm, init.state, data)
      
      if((dist[1] <= epsilon[["dist"]][G]) && (dist[2] <= epsilon[["foodH"]][G]) && (dist[3] <= epsilon[["AMRH"]][G]) && 
         (dist[4] <= epsilon[["foodA"]][G]) && (dist[5] <= epsilon[["AMRA"]][G]) && (!is.na(dist))) {
        
        if(G==1){
          w.new <- 1
        } else {
          w1 <- prod(c(sapply(c(1:3,5:8), function(b) dunif(new.parms[b], min=lm.low[b], max=lm.upp[b])),
                       dbeta(new.parms[4], 1.5, 8.5))) 
          w2 <- sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(new.parms, mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
          w.new <- w1/w2
        }
        i <- i + 1
        return(list(dist, m, new.parms, w.new))
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
    
    print("test")
    
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
    colnames(res.old) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHD", "betaHH", "betaHI")
    write.csv(res.old, file = paste("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data/Test_051121/ABC_post_amppigs_",g,".csv",sep=""), row.names=FALSE)
  
    }
  return(out)
}

# Running the Model Fit ---------------------------------------------------
  
cl <- makeCluster(7, type="SOCK")

clusterEvalQ(cl, {c(library("rootSolve"), library("tmvtnorm"))})

test <- ABC_algorithm(N = 100,
                      G = 1,
                      distanceABC = sum_square_diff_dist, 
                      fitmodel = amrimp, 
                      tau_range = melt_amp_pigs$Usage, 
                      init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
                      data = melt_amp_pigs, 
                      epsilon = list("dist" <-  c(2, 1.75, 1.5, 1.25, 1, 0.9),
                                     "foodH" <- c(0.593, 0.593*0.75, 0.593*0.6, 0.593*0.5, 0.593*0.4, 0.593*0.2),
                                     "AMRH" <-  c(0.185, 0.185*0.75, 0.185*0.6, 0.185*0.5, 0.185*0.4, 0.185*0.2),
                                     "foodA" <- c(UK_cont, UK_cont*0.75, UK_cont*0.6, UK_cont*0.5, UK_cont*0.4, UK_cont*0.2),
                                     "AMRA" <-  c(UK_amp_res, UK_amp_res*0.75, UK_amp_res*0.6, UK_amp_res*0.5, UK_amp_res*0.4, UK_amp_res*0.2)), 
                      lm.low = c(0, 0, 0, 0, 0, 0, 0, 0), 
                      lm.upp = c(0.035, 0.1, 10, 1, 0.005, 0.001, 0.01, 0.0002), 
                      thetaparm = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = UK_food_usage,
                                    fracimp = EU_cont, propres_imp = EU_res))

stopCluster(cl)


















#Run the fit - This is where I will build the ABC-SMC Approach

ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data) {
  N_ITER_list <- list()
  for(g in 1:G) {
    i <- 1
    dist_data <- data.frame(matrix(nrow = 1000, ncol = 6))
    N_ITER <- 1
    
    while(i <= N) {
      
      N_ITER <- N_ITER + 1
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.035)
        d_phi <- runif(1, min = 0, max = 0.1)
        d_kappa <- runif(1, min = 0, max = 10)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 0.005)
        
        d_betaHD <- runif(1, 0, 0.001)
        d_betaHH <- runif(1, 0, 0.01)
        d_betaHI <- runif(1, 0, 0.0002)
        
      } else{ 
        p <- sample(seq(1,N),1,prob = w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp, algorithm = "gibbs")
        d_betaAA<-par[1]
        d_phi<-par[2]
        d_kappa<-par[3]
        d_alpha<-par[4]
        d_zeta <- par[5]
        
        d_betaHD <- par[6]
        d_betaHH <- par[7]
        d_betaHI <- par[8]
      }
      
      if(prior.non.zero(c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHH, d_betaHI))) {
        m <- 0
        
        thetaparm <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = 0, psi = 0.656,
                       
                       fracimp = EU_cont,
                       
                       propres_imp = EU_res,
                       
                       betaAA = d_betaAA, betaHH = d_betaHH, betaHD = d_betaHD,
                       betaHI = d_betaHI,
                       
                       phi = d_phi, kappa = d_kappa, alpha = d_alpha, zeta = d_zeta)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data)
        print(dist)

        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_foodH[g]) && (dist[3] <= epsilon_AMRH[g]) && 
           (dist[4] <= epsilon_foodA[g]) && (dist[5] <= epsilon_AMRA[g]) && (!is.na(dist))) {
          
          # Store results
          res.new[i,]<-c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHH, d_betaHI) 
          dist_data[i,] <- dist
          print(res.new[i,])

          
          # Calculate weights
          if(g==1){
            
            w.new[i] <- 1
            
          } else {
            w1<-prod(c(sapply(c(1:3,5:8), function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])),
                       dbeta(res.new[i,4], 1.5, 8.5))) 
            w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
            w.new[i] <- w1/w2
          }
          # Update counter
          print(paste0('Generation: ', g, ", particle: ", i,", weights: ", w.new[i]))
          
          i <- i+1
        }
      }
    }#
    N_ITER_list[[g]] <- list(N_ITER, dist_data)
    sigma <- cov(res.new) 
    res.old <- res.new
    print(res.old)
    w.old <- w.new/sum(w.new)
    colnames(res.new) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHD", "betaHH","betaHI")
    write.csv(res.new, file = paste("PART1_Simple_FIT_AMP_",g,".csv",sep=""), row.names=FALSE)
    ####
  }
  return(N_ITER_list)
}

N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 0, 0, 0)
lm.upp <- c(0.035, 0.1, 10, 1, 0.005, 0.001, 0.01, 0.0002)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=8,nrow=N)
res.new<-matrix(ncol=8,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <-  c(2, 1.75, 1.5, 1.25, 1, 0.9)
epsilon_foodH <- c(3.26, 3.26*0.75, 3.26*0.6, 3.26*0.5, 3.26*0.4, 3.26*0.2)
epsilon_AMRH <-  c(0.185, 0.185*0.75, 0.185*0.6, 0.185*0.5, 0.185*0.4, 0.185*0.2)
epsilon_foodA <- c(0.017173052, 0.017173052*0.75, 0.017173052*0.6, 0.017173052*0.5, 0.017173052*0.4, 0.017173052*0.2)
epsilon_AMRA <-  c(0.3333333, 0.3333333*0.75, 0.3333333*0.6, 0.3333333*0.5, 0.3333333*0.4, 0.3333333*0.2)

dist_save <- ABC_algorithm(N = 1000, 
                           G = 6,
                           sum.stats = summarystatprev, 
                           distanceABC = sum_square_diff_dist, 
                           fitmodel = amrimp, 
                           tau_range = country_data_gen$scaled_sales_amp, 
                           init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, 
                                          Sh = 1,
                                          Ish = 0, Irh = 0), 
                           times = seq(0, 2000, by = 50), 
                           data = country_data_gen)

end_time <- Sys.time(); end_time - start_time

saveRDS(dist_save, file = "dist_simple_amp.rds")

# Check Posterior ---------------------------------------------------------

amp_post <- do.call(rbind,
                    lapply(list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data", pattern = "Simple_FIT")[1:6], read.csv))
amp_post$gen <- as.vector(sapply(1:6, 
                                 function(x) rep(paste0("gen_",x), 1000)))

l_amp_post <- lapply(1:length(colnames(amp_post)[-1]), function(x) melt(amp_post, id.vars = "gen", measure.vars = colnames(amp_post)[x])[,c(1,3)])

names = c(expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")),
          expression(paste("Rate of Antibiotic-Resistant to Antibiotic-Sensitive Reversion (", phi, ")")),
          expression(paste("Efficacy of Antibiotic-Mediated Animal Recovery (", kappa, ")")),
          expression(paste("Transmission-related Antibiotic Resistant Fitness Cost (", alpha, ")")),
          expression(paste("Background Infection Rate (", zeta, ")")),
          expression(paste("Rate of Domestic Animal-to-Human Transmission (", beta[HD], ")")),
          expression(paste("Rate of Domestic Human-to-Human Transmission (", beta[HH], ")")),
          expression(paste("Rate of EU Animal-to-Human Transmission (", beta[HI], ")")))

amp_p_list <- list()

for(i in 1:length(names)){ 
  amp_p_list[[i]] <- ggplot(l_amp_post[[i]], aes(x=value, fill= gen)) + geom_density(alpha=.5) + 
    scale_x_continuous(expand = c(0, 0), name = names[i]) + 
    scale_y_continuous(expand = c(0, 0), name = "") +
    labs(fill = NULL) + scale_fill_discrete(labels = sapply(1:length(unique(amp_post$gen)), function(x) paste0("Generation ", x)))+
    theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
}


p_comb <- ggarrange(amp_p_list[[1]],amp_p_list[[2]],amp_p_list[[3]],amp_p_list[[4]],
          amp_p_list[[5]],amp_p_list[[6]],amp_p_list[[7]],amp_p_list[[8]], ncol = 2, nrow = 4)

ggsave(p_comb, filename = "post_simple.png",dpi = 300, type = "cairo", width = 12, height = 12, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Figures")
