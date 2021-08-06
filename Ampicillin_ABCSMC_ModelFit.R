library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs"); library("Surrogate")

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
    
    dSa = ua + ra*(Isa + Ira) + kappa*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    
    dIsa = betaAA*Isa*Sa + phi*Ira - kappa*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(1-Sh) - uh*Sh - 
      
      betaHH*Sh*(IshDA+IshA1+IshA2+IshA3+IshH) - 
      (1-alpha)*betaHH*Sh*(IrhDA+IrhA1+IrhA2+IrhA3+IrhH) - 
      
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
      
      (1-psi)*(imp_nEU)*(betaHI_nEU*fracimp_nEU*(1-propres_impnEU)*Sh) - 
      (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_nEU*fracimp_nEU*propres_impnEU*Sh)
    
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
    
    dIshAnEU = (1-psi)*(imp_nEU)*(betaHI_nEU*fracimp_nEU*(1-propres_impnEU)*Sh) - rh*IshAnEU - uh*IshAnEU 
    dIrhAnEU = (1-psi)*(imp_nEU)*(1-alpha)*(betaHI_nEU*fracimp_nEU*propres_impnEU*Sh) - rh*IrhAnEU - uh*IrhAnEU  
    
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

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data) {
  tauoutput <- matrix(nrow = 0, ncol = 5)
  tau_range <- c(tau_range, UK_amp)
  
  for (i in 1:length(tau_range)) {
    temp <- matrix(NA, nrow = 1, ncol = 5)
    
    parms2 = c(ra = thetaparm[["ra"]], rh = thetaparm[["rh"]], ua = thetaparm[["ua"]], uh = thetaparm[["uh"]], tau = tau_range[i], psi = thetaparm[["psi"]],
               
               fracimp1 = thetaparm[["fracimp1"]], fracimp2 = thetaparm[["fracimp2"]], fracimp3 = thetaparm[["fracimp3"]], fracimp4 = thetaparm[["fracimp4"]], 
               fracimp5 = thetaparm[["fracimp5"]], fracimp6 = thetaparm[["fracimp6"]], fracimp7 = thetaparm[["fracimp7"]], fracimp8 = thetaparm[["fracimp8"]], 
               fracimp9 = thetaparm[["fracimp9"]], fracimp_nEU = thetaparm[["fracimp_nEU"]],
               
               imp1 = thetaparm[["imp1"]], imp2 = thetaparm[["imp2"]], imp3 = thetaparm[["imp3"]], imp4 =thetaparm[["imp4"]], 
               imp5 = thetaparm[["imp5"]], imp6 = thetaparm[["imp6"]], imp7 = thetaparm[["imp7"]], imp8 = thetaparm[["imp8"]],
               imp8 = thetaparm[["imp8"]], imp9 = thetaparm[["imp9"]],
               
               propres_imp1 = thetaparm[["propres_imp1"]], propres_imp2 = thetaparm[["propres_imp2"]], propres_imp3 = thetaparm[["propres_imp3"]], propres_imp4 = thetaparm[["propres_imp4"]], 
               propres_imp5 = thetaparm[["propres_imp5"]], propres_imp6 = thetaparm[["propres_imp6"]], propres_imp7 = thetaparm[["propres_imp7"]], propres_imp8 = thetaparm[["propres_imp8"]], 
               propres_imp9 = thetaparm[["propres_imp9"]],
               
               betaAA = thetaparm[["betaAA"]], betaHH = thetaparm[["betaHH"]], betaHD = thetaparm[["betaHD"]],
               betaHI_EU = thetaparm[["betaHI_EU"]], betaHI_nEU = thetaparm[["betaHI_nEU"]], 
               
               imp_nEU = thetaparm[["imp_nEU"]], propres_impnEU = thetaparm[["propres_impnEU"]],
               phi = thetaparm[["phi"]], kappa = thetaparm[["kappa"]], alpha = thetaparm[["alpha"]], zeta = thetaparm[["zeta"]])
    
    out <- ode(y = init.state, func = fitmodel, times = times, parms = parms2)
    
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (out[nrow(out),"Isa"] + out[nrow(out),"Ira"])
    temp[1,3] <- (sum(out[nrow(out),seq(6, 29)]))*100000
    temp[1,4] <- out[nrow(out),"Ira"]/ (out[nrow(out),"Isa"] + out[nrow(out),"Ira"])
    temp[1,5] <- sum(out[nrow(out),seq(7, 29, by =2)]) /  sum(out[nrow(out),seq(6, 29)])
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau",  "ICombA", "ICombH",  "propres_amp", "ResPropHum") 
  
  
  return(c(distanceABC(list(sum.stats), data, 
                       tauoutput[!tauoutput$tau == UK_amp,]),
           abs(tauoutput$ICombH[tauoutput$tau == UK_amp] - 3.26),
           abs(tauoutput$ResPropHum[tauoutput$tau == UK_amp] - 0.185),
           abs(tauoutput$ICombA[tauoutput$tau == UK_amp] - 0.017173052),
           abs(tauoutput$propres_amp[tauoutput$tau == UK_amp] - 0.1111111)))
}

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

#Where G is the number of generations
prior.non.zero<-function(par){
  prod(sapply(1:11, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}


ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data, data_match) {
  N_ITER_list <- list()
  for(g in 1:G) {
    i <- 1
    dist_data <- data.frame(matrix(nrow = 1000, ncol = 5))
    N_ITER <- 1
    
    while(i <= N) {
      
      N_ITER <- N_ITER + 1
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.01)
        d_phi <- runif(1, min = 0, max = 0.2)
        d_kappa <- runif(1, min = 0, max = 50)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 0.001)
        
        d_betaHD <- runif(1, 0, 0.002)
        d_betaHH <- runif(1, 0, 0.01)
        d_betaHI_EU <- runif(1, 0, 0.002)
        d_betaHI_nEU <- runif(1, 0, 0.0025)
        d_imp_nEU <- runif(1, 0, 1)
        d_propres_impnEU <- runif(1, 0, 1)
        
      } else{ 
        p <- sample(seq(1,N),1,prob = w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
        print(par)
        d_betaAA<-par[1]
        d_phi<-par[2]
        d_kappa<-par[3]
        d_alpha<-par[4]
        d_zeta <- par[5]
        
        d_betaHD <- par[6]
        d_betaHH <- par[7]
        d_betaHI_EU <- par[8]
        d_betaHI_nEU <- par[9]
        d_imp_nEU <- par[10]
        d_propres_impnEU <- par[11]
      }
      
      
      if(prior.non.zero(c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHH, d_betaHI_EU, d_betaHI_nEU, d_imp_nEU, d_propres_impnEU))) {
        m <- 0
        
        thetaparm <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = tau_range[i], psi = 0.656,
                       
                       fracimp1 = data_match[2,"Corrected_Usage_18"], fracimp2 = data_match[3,"Corrected_Usage_18"], fracimp3 = data_match[4,"Corrected_Usage_18"], fracimp4 = data_match[5,"Corrected_Usage_18"], 
                       fracimp5 = data_match[6,"Corrected_Usage_18"], fracimp6 = data_match[7,"Corrected_Usage_18"], fracimp7 = data_match[8,"Corrected_Usage_18"], fracimp8 = data_match[9,"Corrected_Usage_18"], 
                       fracimp9 = data_match[10,"Corrected_Usage_18"], fracimp_nEU = 1 - sum(country_data_imp[1:10,"Corrected_Usage_18"]),
                       
                       imp1 = data_match[2,"Foodborne_Carriage_2019"], imp2 = data_match[3,"Foodborne_Carriage_2019"], imp3 = data_match[4,"Foodborne_Carriage_2019"], imp4 = data_match[5,"Foodborne_Carriage_2019"], 
                       imp5 = data_match[6,"Foodborne_Carriage_2019"], imp6 = data_match[7,"Foodborne_Carriage_2019"], imp7 = data_match[8,"Foodborne_Carriage_2019"], imp8 = data_match[9,"Foodborne_Carriage_2019"],
                       imp8 = data_match[9,"Foodborne_Carriage_2019"], imp9 = data_match[10,"Foodborne_Carriage_2019"],
                       
                       propres_imp1 = data_match[2,"Prop_Amp_Res"], propres_imp2 = data_match[3,"Prop_Amp_Res"], propres_imp3 = data_match[4,"Prop_Amp_Res"], propres_imp4 = data_match[5,"Prop_Amp_Res"], 
                       propres_imp5 = data_match[6,"Prop_Amp_Res"], propres_imp6 = data_match[7,"Prop_Amp_Res"], propres_imp7 = data_match[8,"Prop_Amp_Res"], propres_imp8 = data_match[9,"Prop_Amp_Res"], 
                       propres_imp9 = data_match[10,"Prop_Amp_Res"],
                       
                       betaAA = d_betaAA, betaHH = d_betaHH, betaHD = d_betaHD,
                       betaHI_EU = d_betaHI_EU, betaHI_nEU = d_betaHI_nEU, 
                       
                       imp_nEU = d_imp_nEU, propres_impnEU = d_propres_impnEU,
                       phi = d_phi, kappa = d_kappa, alpha = d_alpha, zeta = d_zeta)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data)
        print(dist)
        
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_foodH[g]) && (dist[3] <= epsilon_AMRH[g]) && 
           (dist[4] <= epsilon_foodA[g]) && (dist[5] <= epsilon_AMRA[g]) && (!is.na(dist))) {
          
          # Store results
          res.new[i,]<-c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHH, d_betaHI_EU, d_betaHI_nEU, d_imp_nEU, d_propres_impnEU) 
          dist_data[i,] <- dist
          print(res.new[i,])
          
          # Calculate weights
          if(g==1){
            
            w.new[i] <- 1
            
          } else {
            w1<-prod(c(sapply(c(1:3,5:11), function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])),
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
    print(N_ITER_list)
    sigma <- cov(res.new) 
    res.old <- res.new
    print(res.old)
    w.old <- w.new/sum(w.new)
    colnames(res.new) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "betaHD", "betaHH","betaHI_EU", "betaHI_nEU", "imp_nEU", "propres_impnEU")
    write.csv(res.new, file = paste("complexmodel_ABC_SMC_gen_amp_",g,".csv",sep=""), row.names=FALSE)
    ####
  }
  return(N_ITER_list)
}

N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
lm.upp <- c(0.01, 0.2, 50, 1, 0.001, 0.002, 0.01, 0.002, 0.0025, 1, 1)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=11,nrow=N)
res.new<-matrix(ncol=11,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <-  c(2, 1.75, 1.5, 1.25, 1, 0.9, 0.8)
epsilon_foodH <- c(3.26, 3.26*0.75, 3.26*0.6, 3.26*0.5, 3.26*0.25, 3.26*0.2, 3.26*0.15, 3.26*0.1)
epsilon_AMRH <-  c(0.185, 0.185*0.75, 0.185*0.6, 0.185*0.5, 0.185*0.25, 0.185*0.2, 0.185*0.15, 0.185*0.1)
epsilon_foodA <- c(0.017173052, 0.017173052*0.75, 0.017173052*0.6, 0.017173052*0.5, 0.017173052*0.25, 0.017173052*0.2, 0.017173052*0.15, 0.017173052*0.1)
epsilon_AMRA <-  c(0.3333333, 0.3333333*0.75, 0.3333333*0.6, 0.3333333*0.5, 0.3333333*0.25, 0.3333333*0.2, 0.3333333*0.15, 0.3333333*0.1)

dist_save <- ABC_algorithm(N = 1000, 
              G = 7,
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
              times = seq(0, 2000, by = 50), 
              data = country_data_gen,
              data_match = country_data_imp)

end_time <- Sys.time(); end_time - start_time

saveRDS(dist_save, file = "dist_amp_list.rds")

# Looking at Intermediate Posterior ---------------------------------------

post1_amp <- read.csv(tail(list.files(path = "C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data", pattern = "^complexmodel_ABC_SMC_gen_amp.*?\\.csv"), 1))

par(mfrow = c(3,4))
plot(density(post1_amp$betaAA))
plot(density(post1_amp$phi))
plot(density(post1_amp$kappa))
plot(density(post1_amp$alpha))
plot(density(post1_amp$zeta))
plot(density(post1_amp$betaHD))
plot(density(post1_amp$betaHH))
plot(density(post1_amp$betaHI_EU))
plot(density(post1_amp$betaHI_nEU))
plot(density(post1_amp$imp_nEU))
plot(density(post1_amp$propres_impnEU))

# Plotting the Resulting Model ----------------------------------------------

MAP_parms <- data.frame("parms" = colnames(post1_amp), "mean" = colMeans(post1_amp))

parmtau <- seq(0, 0.014, by = 0.001)

init = c(Sa=0.98, Isa=0.01, Ira=0.01, 
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
         IshH = 0, IrhH = 0)

output1 <- data.frame(matrix(ncol = 28, nrow = 0))
times <- seq(0, 10000, by = 1)

plot_analysis <- list()

for (j in 1:2) {
  
  parmtau <- list(country_data_gen$scaled_sales_amp,
                  seq(0, 0.014, by = 0.001))[[j]]
  
  plot_analysis[[j]] <- local({
    
    for (i in 1:length(parmtau)) {
      
      temp <- data.frame(matrix(NA, nrow = 1, ncol=28))
      
      parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = parmtau[i], psi = 0.656,
                
                fracimp1 = country_data_imp[2,"Corrected_Usage_18"], fracimp2 = country_data_imp[3,"Corrected_Usage_18"], fracimp3 = country_data_imp[4,"Corrected_Usage_18"], fracimp4 = country_data_imp[5,"Corrected_Usage_18"], 
                fracimp5 = country_data_imp[6,"Corrected_Usage_18"], fracimp6 = country_data_imp[7,"Corrected_Usage_18"], fracimp7 = country_data_imp[8,"Corrected_Usage_18"], fracimp8 = country_data_imp[9,"Corrected_Usage_18"], 
                fracimp9 = country_data_imp[10,"Corrected_Usage_18"], fracimp_nEU = 1 - sum(country_data_imp[1:10,"Corrected_Usage_18"]),
                
                imp1 = country_data_imp[2,"Foodborne_Carriage_2019"], imp2 = country_data_imp[3,"Foodborne_Carriage_2019"], imp3 = country_data_imp[4,"Foodborne_Carriage_2019"], imp4 = country_data_imp[5,"Foodborne_Carriage_2019"], 
                imp5 = country_data_imp[6,"Foodborne_Carriage_2019"], imp6 = country_data_imp[7,"Foodborne_Carriage_2019"], imp7 = country_data_imp[8,"Foodborne_Carriage_2019"], imp8 = country_data_imp[9,"Foodborne_Carriage_2019"],
                imp8 = country_data_imp[9,"Foodborne_Carriage_2019"], imp9 = country_data_imp[10,"Foodborne_Carriage_2019"],
                
                propres_imp1 = country_data_imp[2,"Prop_Amp_Res"], propres_imp2 = country_data_imp[3,"Prop_Amp_Res"], propres_imp3 = country_data_imp[4,"Prop_Amp_Res"], propres_imp4 = country_data_imp[5,"Prop_Amp_Res"], 
                propres_imp5 = country_data_imp[6,"Prop_Amp_Res"], propres_imp6 = country_data_imp[7,"Prop_Amp_Res"], propres_imp7 = country_data_imp[8,"Prop_Amp_Res"], propres_imp8 = country_data_imp[9,"Prop_Amp_Res"], 
                propres_imp9 = country_data_imp[10,"Prop_Amp_Res"],
                
                betaAA = MAP_parms[which(MAP_parms == "betaAA"),2], betaHH = MAP_parms[which(MAP_parms == "betaHH"),2], betaHD =  MAP_parms[which(MAP_parms == "betaHD"),2],
                betaHI_EU =  MAP_parms[which(MAP_parms == "betaHI_EU"),2], betaHI_nEU = MAP_parms[which(MAP_parms == "betaHI_nEU"),2], 
                
                imp_nEU =  MAP_parms[which(MAP_parms == "imp_nEU"),2], propres_impnEU =  MAP_parms[which(MAP_parms == "propres_impnEU"),2],
                phi =  MAP_parms[which(MAP_parms == "phi"),2], kappa =  MAP_parms[which(MAP_parms == "kappa"),2], alpha =  MAP_parms[which(MAP_parms == "alpha"),2], 
                zeta =  MAP_parms[which(MAP_parms == "zeta"),2])
      
      out <<- ode(y = init, func = amrimp, times = times, parms = parms)
      
      temp[1,1] <- parmtau[i]
      temp[1,2:13] <- out[nrow(out),seq(7, 29, by = 2)]/sum(out[nrow(out),seq(7, 29, by = 2)])
      temp[1,14:25] <- sapply(seq(6,28, by = 2), function(x) out[nrow(out),x]) + sapply(seq(7,29, by = 2), function(x) out[nrow(out),x]) 
      temp[1,26] <- sum(out[nrow(out),seq(7, 29, by = 2)])/sum(out[nrow(out),seq(6, 29)]) 
      temp[1,27] <- out[nrow(out), 4] / sum(out[nrow(out),3:4]) 
      temp[1,28] <- sum(out[nrow(out),3:4]) 
      
      print(paste0("Run ", j, ": ",temp[1,2]))
      output1 <- rbind.data.frame(output1, temp)
    }
    
    colnames(output1)[1:28] <- c("tau", 
                                 
                                 "PropRes_Domestic","PropRes_Netherlands","PropRes_Ireland","PropRes_Germany","PropRes_France","PropRes_Spain","PropRes_Italy",
                                 "PropRes_Belgium","PropRes_Poland","PropRes_Denmark","PropRes_NonEU", "PropRes_Human",
                                 
                                 "TotInf_Domestic","TotInf_Netherlands","TotInf_Ireland","TotInf_Germany","TotInf_France","TotInf_Spain","TotInf_Italy",
                                 "TotInf_Belgium","TotInf_Poland","TotInf_Denmark","TotInf_NonEU", "TotInf_Human", 
                                 
                                 "Res_PropHum", "Res_PropAnim", "Anim_Inf")
    return(output1)
  })
}

res_plotdata <- melt(plot_analysis[[2]], id.vars = c("tau"), measure.vars = colnames(plot_analysis[[1]])[2:13]) 
trim_names <- as.factor(sapply(strsplit(as.character(res_plotdata$variable), split = "_", fixed = TRUE), function(x) x[2]))
res_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

inf_plotdata <- melt(plot_analysis[[2]], id.vars = c("tau"), measure.vars = colnames(plot_analysis[[1]])[14:25]) 
inf_plotdata$variable <- factor(trim_names, levels = unique(trim_names))

#Model Output - Attribution Plots 

res_comb <- ggplot(res_plotdata, aes(fill = variable, x = tau, y = value)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'))  +
  labs(x ="Domestic Ampicillin Usage in Fattening Pig (g/PCU)", y = "Proportion of Attributable Human Resistance", fill = "Resistance Source")  

inf_comb <- ggplot(inf_plotdata, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() +  
  geom_col(color = "black",position= "stack", width  = 0.001) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 5)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Total Human Salmenollosis (per 100,000)", fill = "Resistance Source") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_vline(xintercept = UK_amp, size = 1.2, col = "red", lty = 2)

# Fit to Data -------------------------------------------------------------

plot_check <- data.frame("tau" = plot_analysis[[1]]$tau, "country" = country_data_gen$Country, 
                         "model_estim" = plot_analysis[[1]]$Res_PropAnim, "observ" = country_data_gen$propres_tet, 
                         "model_estim_inf" = plot_analysis[[1]]$Anim_Inf)

country_data_gen <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Data/res_sales_generalfit.csv") #This is data for pigs 
country_data_gen[,13:14] <- country_data_gen[,13:14]/1000

#Animal Resistance vs Animal Livestock Antibiotic Usage
anim_plot <- ggplot(plot_check, mapping = aes(x = tau)) + geom_point(aes(y = observ), size = 2) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Domestic Livestock Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp, y = country_data_gen$propres_amp[country_data_gen$Country == "United Kingdom"]) , col = "red", size = 5)

#Animal Infection vs Animal Livestock Antibiotic Usage

anim_plot_inf <- ggplot(plot_check, mapping = aes(x = tau)) + geom_line(aes(y = model_estim_inf), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.035)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Domestic Livestock Contaminated") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp, y = country_data_imp$Foodborne_Carriage_2019[country_data_imp$Country_of_Origin == "UK Origin"]) , col = "red", size = 5)

#Animal Resistance vs Animal Livestock Antibiotic Usage
plot_check_hum <- plot_check; plot_check_hum$model_estim <- plot_analysis[[1]]$Res_PropHum

hum_plot <- ggplot(plot_check_hum, mapping = aes(x = tau)) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Human Infections Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  geom_point(aes(x = UK_amp, y = 0.3108) , col = "red", size = 5)

ggplot(plot_check_hum, mapping = aes(x = tau)) + geom_line(aes(y = model_estim), size = 1.2) + theme_bw() +  
  scale_x_continuous(expand = c(0, 0), limits = c(0,0.015)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  labs(x ="Ampicillin Sales in Fattening Pig (g/PCU)", y = "Proportion of Human Infections Salmonella Resistant") +
  theme(legend.position= "bottom", legend.text=element_text(size=12), legend.title =element_text(size=12), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + geom_vline(xintercept = UK_amp, col = "red", size = 1.2, lty = 2)
