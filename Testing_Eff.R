library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate"); library("rootSolve")
library("profvis")

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

country_data_imp <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model Fit Data/FullData_2021_v1_trim.csv") #This is data for pigs 
country_data_imp$Foodborne_Carriage_2019 <- country_data_imp$Foodborne_Carriage_2019/100
country_data_imp$Corrected_Usage_18 <- country_data_imp$Corrected_Usage_18/100
country_data_imp[,12:13] <- country_data_imp[,12:13]/1000

country_data_gen <- read.csv("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model Fit Data/res_sales_generalfit.csv") #This is data for pigs 
country_data_gen[,13:14] <- country_data_gen[,13:14]/1000

UK_tet <- country_data_gen$scaled_sales_tet[country_data_gen$Country == "United Kingdom"]
UK_amp <- country_data_gen$scaled_sales_amp[country_data_gen$Country == "United Kingdom"]

country_data_gen <- country_data_gen[country_data_gen$num_test_amp >= 10,]

plot(country_data_gen$scaled_sales_tet, country_data_gen$propres_tet, ylim = c(0,1))
plot(country_data_gen$scaled_sales_amp, country_data_gen$propres_amp, ylim = c(0,1))

# ABC Functions -----------------------------------------------------------

#### Approximate Bayesian Computation - SMC ####


summarystatprev <- function(prev) {
  return(prev$propres_amp)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, data) {
  tau_range <- c(0, tau_range, UK_amp)
  tauoutput <- matrix(nrow = length(tau_range), ncol = 5)
  
  for (i in 1:length(tau_range)) {
    temp <- numeric(5)
    
    parms2 = thetaparm
    parms2["tau"] = tau_range[i]
    
    out <- runsteady(y = init.state, func = fitmodel, parms = parms2, times = c(0, Inf))
    
    temp[1] <- tau_range[i]
    temp[2] <- out[[1]][["Isa"]] + out[[1]][["Ira"]]
    temp[3] <- sum(out[[1]][c(5:28)])*100000
    temp[4] <- out[[1]][["Ira"]]/ (out[[1]][["Isa"]] + out[[1]][["Ira"]])
    temp[5] <- sum(out[[1]][seq(6,28, by = 2)]) /  sum(out[[1]][c(5:28)])
    tauoutput[i,] <- temp
  }
  tauoutput <- data.frame(tauoutput)
  
  colnames(tauoutput) <- c("tau", "ICombA", "ICombH","propres_amp", "ResPropHum") 
  
  return(c(distanceABC(list(sum.stats), data, 
                       tauoutput[(!tauoutput$tau == UK_amp & !tauoutput$tau == 0),]),
           abs(tauoutput$ICombH[tauoutput$tau == UK_amp] - 3.26),
           abs(tauoutput$ResPropHum[tauoutput$tau == UK_amp] - 0.185),
           abs(tauoutput$ICombA[tauoutput$tau == UK_amp] - 0.017173052),
           abs(tauoutput$ICombA[tauoutput$tau == 0]),
           abs(tauoutput$propres_amp[tauoutput$tau == UK_amp] - 0.1111111)))
}

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

#Where G is the number of generations
prior.non.zero <- function(par){
  prod(sapply(1:10, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}


ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, data, data_match) {
  N_ITER_list <- list()
  
  fit_parms <- c("betaAA", "betaHH", "betaHD", "betaHI_EU", "imp_nEU", "propres_impnEU", "phi", "kappa", "alpha", "zeta")
  thetaparm <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = 0.656,
                 
                 fracimp1 = as.numeric(data_match[2,"Normalised_Usage_2018"]), fracimp2 = as.numeric(data_match[3,"Normalised_Usage_2018"]), fracimp3 = as.numeric(data_match[4,"Normalised_Usage_2018"]), 
                 fracimp4 = as.numeric(data_match[5,"Normalised_Usage_2018"]), fracimp5 = as.numeric(data_match[6,"Normalised_Usage_2018"]), fracimp6 = as.numeric(data_match[7,"Normalised_Usage_2018"]), 
                 fracimp7 = as.numeric(data_match[8,"Normalised_Usage_2018"]), fracimp8 = as.numeric(data_match[9,"Normalised_Usage_2018"]), 
                 fracimp9 = as.numeric(data_match[10,"Normalised_Usage_2018"]), fracimp_nEU = 1 - sum(as.numeric(country_data_imp[2:10,"Normalised_Usage_2018"])),
                 
                 imp1 = data_match[2,"Foodborne_Carriage_2019"], imp2 = data_match[3,"Foodborne_Carriage_2019"], imp3 = data_match[4,"Foodborne_Carriage_2019"], imp4 = data_match[5,"Foodborne_Carriage_2019"], 
                 imp5 = data_match[6,"Foodborne_Carriage_2019"], imp6 = data_match[7,"Foodborne_Carriage_2019"], imp7 = data_match[8,"Foodborne_Carriage_2019"], imp8 = data_match[9,"Foodborne_Carriage_2019"],
                 imp8 = data_match[9,"Foodborne_Carriage_2019"], imp9 = data_match[10,"Foodborne_Carriage_2019"],
                 
                 propres_imp1 = data_match[2,"Prop_Amp_Res"], propres_imp2 = data_match[3,"Prop_Amp_Res"], propres_imp3 = data_match[4,"Prop_Amp_Res"], propres_imp4 = data_match[5,"Prop_Amp_Res"], 
                 propres_imp5 = data_match[6,"Prop_Amp_Res"], propres_imp6 = data_match[7,"Prop_Amp_Res"], propres_imp7 = data_match[8,"Prop_Amp_Res"], propres_imp8 = data_match[9,"Prop_Amp_Res"], 
                 propres_imp9 = data_match[10,"Prop_Amp_Res"])
  
  for(g in 1:G) {
    i <- 1
    dist_data <- data.frame(matrix(nrow = G, ncol = 6))
    N_ITER <- 1
    
    while(i <= N) {
      
      N_ITER <- N_ITER + 1
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.05)
        d_phi <- runif(1, min = 0, max = 0.1)
        d_kappa <- runif(1, min = 0, max = 10)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 0.001)
        
        d_betaHD <- runif(1, 0, 0.001)
        d_betaHH <- runif(1, 0, 0.025)
        d_betaHI_EU <- runif(1, 0, 0.0005)
        d_imp_nEU <- runif(1, 0, 1)
        d_propres_impnEU <- runif(1, 0, 1)
        
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
        d_betaHI_EU <- par[8]
        d_imp_nEU <- par[9]
        d_propres_impnEU <- par[10]
      }
      
      if(prior.non.zero(c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHH, d_betaHI_EU, d_imp_nEU, d_propres_impnEU))) {
        m <- 0
        
        thetaparm[fit_parms] <- c(d_betaAA, d_betaHH, d_betaHD, d_betaHI_EU, d_imp_nEU, d_propres_impnEU, d_phi, d_kappa, d_alpha, d_zeta)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, data)
        print(dist)
        
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_foodH[g]) && (dist[3] <= epsilon_AMRH[g]) && 
           (dist[4] <= epsilon_foodA[g]) && (dist[5] < 0.95) && (dist[6] <= epsilon_AMRA[g]) && (!is.na(dist))) {
          
          # Store results
          res.new[i,]<-c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHH, d_betaHI_EU, d_imp_nEU, d_propres_impnEU) 
          dist_data[i,] <- dist
          print(res.new[i,])
          
          # Calculate weights
          if(g==1){
            
            w.new[i] <- 1
            
          } else {
            w1<-prod(c(sapply(c(1:3,5:10), function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])),
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
    colnames(res.new) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHD", "betaHH","betaHI_EU", "imp_nEU", "propres_impnEU")
    ####
  }
  return(N_ITER_list)
}

N <- 1 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
lm.upp <- c(0.05, 0.1, 10, 1, 0.001, 0.001, 0.025, 0.0005, 1, 1)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=10,nrow=N)
res.new<-matrix(ncol=10,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <-  c(2)
epsilon_foodH <- c(3.26)
epsilon_AMRH <-  c(0.185)
epsilon_foodA <- c(0.017173052)
epsilon_AMRA <-  c(0.3333333)

dist_save <- ABC_algorithm(N = 1, 
                           G = 1,
                           sum.stats = summarystatprev, 
                           distanceABC = sum_square_diff_dist, 
                           fitmodel = amrimp, 
                           tau_range = country_data_gen$scaled_sales_amp, 
                           init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, 
                                          Sh = 1,
                                          IshDA = 0,IrhDA = 0,
                                          IshA1 = 0,IrhA1 = 0,
                                          IshA2 = 0,IrhA2 = 0,
                                          IshA3 = 0,IrhA3 = 0,
                                          IshA4 = 0,IrhA4 = 0,
                                          IshA5 = 0,IrhA5 = 0,
                                          IshA6 = 0,IrhA6 = 0,
                                          IshA7 = 0,IrhA7 = 0,
                                          IshA8 = 0,IrhA8 = 0,
                                          IshA9 = 0,IrhA9 = 0,
                                          IshAnEU = 0,IrhAnEU = 0,
                                          IshH = 0, IrhH = 0),  
                           data = country_data_gen,
                           data_match = country_data_imp)

