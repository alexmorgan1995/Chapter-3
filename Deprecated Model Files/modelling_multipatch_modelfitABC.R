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
    
    dSh = uh + rh*(1-Sh) - 
      
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

country_data_gen <- read.csv("C:/Users/amorg/Documents/PhD/Chapter_3/Data/res_sales_generalfit.csv") #This is data for pigs 
country_data_gen <- country_data_gen[country_data_gen$num_test_amp >= 10,]

plot(country_data_gen$scaled_sales_tet, country_data_gen$propres_tet)
plot(country_data_gen$scaled_sales_amp, country_data_gen$propres_amp)

# ABC Functions -----------------------------------------------------------


#### Approximate Bayesian Computation - Rejection Algorithm ####

summarystatprev <- function(prev) {
  return(prev$propres_tet)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}
#

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data) {
  tauoutput <- matrix(nrow = 0, ncol = 5)
  tau_range <- c(tau_range, 0.01287957)
  
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
  colnames(tauoutput) <- c("tau",  "ICombA", "ICombH",  "propres_tet", "ResPropHum") 
  
  
  return(c(distanceABC(list(sum.stats), data, 
                       tauoutput[!tauoutput$tau == 0.01287957,]),
           abs(tauoutput$ICombH[tauoutput$tau == 0.01287957] - 3.26),
           abs(tauoutput$ResPropHum[tauoutput$tau == 0.01287957] - 0.35),
           abs(tauoutput$ICombA[tauoutput$tau == 0.01287957] - 0.017173052),
           abs(tauoutput$propres_tet[tauoutput$tau == 0.01287957] - 0.3333333)))
}

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

#Where G is the number of generations
prior.non.zero<-function(par){
  prod(sapply(1:10, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data, data_match) {
  for(g in 1:G) {
    i <- 1
    while(i <= N) {
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.002)
        d_phi <- runif(1, min = 0, max = 0.1)
        d_kappa <- runif(1, min = 0, max = 7.5)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 0.002)
        
        d_betaHD <- runif(1, 0, 0.002)
        d_betaHI_EU <- runif(1, 0, 0.002)
        d_betaHI_nEU <- runif(1, 0, 0.002)
        d_imp_nEU <- runif(1, 0, 1)
        d_propres_impnEU <- runif(1, 0, 1)
        
      } else{ 
        p <- sample(seq(1,N),1,prob= w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
        d_betaAA<-par[1]
        d_phi<-par[2]
        d_kappa<-par[3]
        d_alpha<-par[4]
        d_zeta <- par[5]
        
        d_betaHD <- par[6]
        d_betaHI_EU <- par[7]
        d_betaHI_nEU <- par[8]
        d_imp_nEU <- par[9]
        d_propres_impnEU <- par[10]
      }
      if(prior.non.zero(c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHI_EU, d_betaHI_nEU, d_imp_nEU, d_propres_impnEU))) {
        m <- 0
        
        thetaparm <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, tau = tau_range[i], psi = 0.656,
                       
                       fracimp1 = data_match[2,"Corrected_Usage_18"], fracimp2 = data_match[3,"Corrected_Usage_18"], fracimp3 = data_match[4,"Corrected_Usage_18"], fracimp4 = data_match[5,"Corrected_Usage_18"], 
                       fracimp5 = data_match[6,"Corrected_Usage_18"], fracimp6 = data_match[7,"Corrected_Usage_18"], fracimp7 = data_match[8,"Corrected_Usage_18"], fracimp8 = data_match[9,"Corrected_Usage_18"], 
                       fracimp9 = data_match[10,"Corrected_Usage_18"], fracimp_nEU = 1 - sum(country_data_imp[1:10,"Corrected_Usage_18"]),
                       
                       imp1 = data_match[2,"Foodborne_Carriage_2019"], imp2 = data_match[3,"Foodborne_Carriage_2019"], imp3 = data_match[4,"Foodborne_Carriage_2019"], imp4 = data_match[5,"Foodborne_Carriage_2019"], 
                       imp5 = data_match[6,"Foodborne_Carriage_2019"], imp6 = data_match[7,"Foodborne_Carriage_2019"], imp7 = data_match[8,"Foodborne_Carriage_2019"], imp8 = data_match[9,"Foodborne_Carriage_2019"],
                       imp8 = data_match[9,"Foodborne_Carriage_2019"], imp9 = data_match[10,"Foodborne_Carriage_2019"],
                       
                       propres_imp1 = data_match[2,"Prop_Tet_Res"], propres_imp2 = data_match[3,"Prop_Tet_Res"], propres_imp3 = data_match[4,"Prop_Tet_Res"], propres_imp4 = data_match[5,"Prop_Tet_Res"], 
                       propres_imp5 = data_match[6,"Prop_Tet_Res"], propres_imp6 = data_match[7,"Prop_Tet_Res"], propres_imp7 = data_match[8,"Prop_Tet_Res"], propres_imp8 = data_match[9,"Prop_Tet_Res"], 
                       propres_imp9 = data_match[10,"Prop_Tet_Res"],
                       
                       betaAA = d_betaAA, betaHH = 0.00001, betaHD = d_betaHD,
                       betaHI_EU = d_betaHI_EU, betaHI_nEU = d_betaHI_nEU, 
                       
                       imp_nEU = d_imp_nEU, propres_impnEU = d_propres_impnEU,
                       phi = d_phi, kappa = d_kappa, alpha = d_alpha, zeta = d_zeta)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data)
        print(dist)
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_foodH[g]) && (dist[3] <= epsilon_AMRH[g]) && 
           (dist[4] <= epsilon_foodA[g]) && (dist[5] <= epsilon_AMRA[g]) && (!is.na(dist))) {
          
          # Store results
          res.new[i,]<-c(d_betaAA, d_phi, d_kappa, d_alpha, d_zeta, d_betaHD, d_betaHI_EU, d_betaHI_nEU, d_imp_nEU, d_propres_impnEU) 
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
    sigma <- cov(res.new) 
    res.old <- res.new
    print(res.old)
    w.old <- w.new/sum(w.new)
    colnames(res.new) <- c("betaAA", "phi", "kappa", "alpha", "zeta", "d_betaHD", "d_betaHI_EU", "d_betaHI_nEU", "d_imp_nEU", "d_propres_impnEU")
    write.csv(res.new, file = paste("complexmodel_ABC_SMC_gen_tet_",g,".csv",sep=""), row.names=FALSE)
    ####
  }
}

N <- 1000 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
lm.upp <- c(0.002, 0.1, 7.5, 1, 0.002, 0.002, 0.002, 0.002, 1, 1)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=10,nrow=N)
res.new<-matrix(ncol=10,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <- c(4.5, 4.25, 4.1, 4, 3.9)
epsilon_foodH <- c(3.26*0.75, 3.26*0.6, 3.26*0.5, 3.26*0.25, 3.26*0.1)
epsilon_AMRH <- c(0.35*0.75, 0.35*0.6, 0.35*0.5, 0.35*0.25, 0.35*0.1)
epsilon_foodA <- c(0.017173052*0.75, 0.017173052*0.6, 0.017173052*0.5, 0.017173052*0.25, 0.017173052*0.1)
epsilon_AMRA <- c(0.3333333*0.75, 0.3333333*0.6, 0.3333333*0.5, 0.3333333*0.25, 0.3333333*0.1)

ABC_algorithm(N = 1000, 
              G = 5,
              sum.stats = summarystatprev, 
              distanceABC = sum_square_diff_dist, 
              fitmodel = amrimp, 
              tau_range = country_data_gen$scaled_sales_tet, 
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
              times = seq(0, 2000, by = 100), 
              data = country_data_gen,
              data_match = country_data_imp)

end_time <- Sys.time(); end_time - start_time

# Looking at Intermediate Posterior ---------------------------------------

post1_tet <- read.csv("complexmodel_ABC_SMC_gen_tet_3.csv") #This is data for pigs 

plot(density(post1_tet$betaAA))
plot(density(post1_tet$phi))
plot(density(post1_tet$kappa))
plot(density(post1_tet$alpha))
plot(density(post1_tet$zeta))
plot(density(post1_tet$d_betaHD))
plot(density(post1_tet$d_betaHI_EU))
plot(density(post1_tet$d_betaHI_nEU))
plot(density(post1_tet$d_imp_nEU))
plot(density(post1_tet$d_propres_impnEU))


# Plotting the Resulting Model ----------------------------------------------


times <- seq(0,30000, by = 1) 

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, 
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
          IshA10 = 0,IrhA10 = 0,
          IshH = 0, IrhH = 0)

parms = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.0031, betaHH = 0.00001, tau = 0.0123,
          betaHD = (0.00001),  
          
          betaHI1 = (0.00001)/10, betaHI2 = (0.00001)/10, betaHI3 = (0.00001)/10,  betaHI4 = (0.00001)/10, betaHI5 = (0.00001)/10, betaHI6 = (0.00001)/10,
          betaHI7 = (0.00001)/10, betaHI8 = (0.00001)/10, betaHI9 = (0.00001)/10, betaHI10 = (0.00001)/10,
          
          phi = 0.03, theta = 1.1, alpha = 0.3, zeta = 0.0001, psi = 0.656, 
          
          fracimp1 = 0.076719577, fracimp2 = 0.046997389, fracimp3 = 0.022156129, fracimp4 = 0.116216976, fracimp5 = 0.083127572, fracimp6 = 0.057626057,
          fracimp7 = 0.04588015, fracimp8 = 0.010497067, fracimp9 = 0.05, fracimp10 = 0.2, 
          
          imp1 = 0.1642031, imp2 = 0.1356062, imp3 = 0.1355515, imp4 = 0.134563, imp5 = 0.104059, imp6 = 0.08818593,
          imp7 = 0.07600404, imp8 = 0.05488995, imp9 = 0.03716985, imp10 = 0.06976744,
          
          propres_imp1 = 0.347826087, propres_imp2 = 0.258064516, propres_imp3 = 0.266666667, propres_imp4 = 0.588235294, propres_imp5 = 0.674698795, 
          propres_imp6 = 0.54822335, propres_imp7 = 0.436893204, propres_imp8 = 0.4, propres_imp9 = 0.456, propres_imp10 = 0.5)


out <- ode(y = init, func = amrimp, times = times, parms = parms)

tail(out)

colnames(output1) <- c("tau", "ICombH","IResRat")
output1$IResRat[output1$tau == 0]

sum(out[nrow(out), seq(7,29,by = 2)])/ sum(out[nrow(out),6:29])


