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

#### Data Import ####
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

#### Approximate Bayesian Computation - Rejection Algorithm ####

summarystatprev <- function(prev) {
  return(prev$ResPropAnim)
}

sum_square_diff_dist <- function(sum.stats, data.obs, model.obs) {
  sumsquare <- sapply(sum.stats, function(x) {
    sumsquare <- abs((x(data.obs) - x(model.obs))^2)
  })
  return(sum(sumsquare))
}

computeDistanceABC_ALEX <- function(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data) {
  tauoutput <- matrix(nrow = 0, ncol = 4)
  tau_range <- append(tau_range, 0.0122887)
  for (i in 1:length(tau_range)) {
    temp <- matrix(NA, nrow = 1, ncol = 4)
    parms2 = c(fracimp = thetaparm[["fracimp"]], propres_imp =  thetaparm[["propres_imp"]], psi = thetaparm[["psi"]],
               ra = thetaparm[["ra"]], rh =  thetaparm[["rh"]], ua = thetaparm[["ua"]], uh = thetaparm[["uh"]], 
               betaAA = thetaparm[["betaAA"]], betaHH = thetaparm[["betaHH"]], 
               betaHI = thetaparm[["betaHI"]], betaHD = thetaparm[["betaHD"]], 
               phi = thetaparm[["phi"]], tau = tau_range[i], theta = thetaparm[["theta"]], 
               alpha = thetaparm[["alpha"]], zeta = thetaparm[["zeta"]])
    out <- ode(y = init.state, func = amrimp, times = times, parms = parms2)
    temp[1,1] <- tau_range[i]
    temp[1,2] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
    temp[1,3] <- (rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4])))
    temp[1,4] <- (rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    tauoutput <- rbind(tauoutput, temp)
  }
  tauoutput <- data.frame(tauoutput)
  colnames(tauoutput) <- c("tau", "ICombH", "ResPropAnim", "ResPropHum")  
  return(c(distanceABC(list(sum.stats), data, tauoutput[!tauoutput$tau == 0.0122887,]),
           abs(tauoutput$ICombH[tauoutput$tau == 0.0122887] - 3.26),
           abs(tauoutput$ResPropHum[tauoutput$tau == 0.0122887] - 0.35)))
}

#Run the fit - This is where I will build the ABC-SMC Approach

start_time <- Sys.time()

#Where G is the number of generations
prior.non.zero<-function(par){
  prod(sapply(1:9, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}

ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, times, data) {
  for(g in 1:G) {
    i <- 1
    while(i <= N) {
      if(g==1) {
        d_betaAA <- runif(1, min = 0, max = 0.2)
        d_betaHI <- runif(1, min = 0, max = 0.0001)
        d_betaHD <- runif(1, min = 0, max = 0.0001)
        d_phi <- runif(1, min = 0, max = 0.04)
        d_theta <- runif(1, min = 0, max = 2)
        d_alpha <- rbeta(1, 1.5, 8.5)
        d_zeta <- runif(1, 0, 0.15)
        d_fracimp <- runif(1, 0, 1)
        d_propres_imp <- runif(1, 0, 1)
      } else{ 
        p <- sample(seq(1,N),1,prob= w.old) # check w.old here
        par <- rtmvnorm(1,mean=res.old[p,], sigma=sigma, lower=lm.low, upper=lm.upp)
        d_betaAA<-par[1]
        d_betaHI <- par[2]
        d_betaHD <- par[3]
        d_phi<-par[4]
        d_theta<-par[5]
        d_alpha<-par[6]
        d_zeta <- par[7]
        d_fracimp <- par[8]
        d_propres_imp <- par[9]
      }
      if(prior.non.zero(c(d_betaAA, d_betaHI, d_betaHD, d_phi, d_theta, d_alpha, d_zeta, d_fracimp, d_propres_imp))) {
        m <- 0
        thetaparm <- c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = d_betaAA, betaHH = 0.00001, 
                       betaHI = d_betaHI, betaHD = d_betaHD, phi = d_phi, theta = d_theta, alpha = d_alpha, zeta = d_zeta,
                       fracimp = d_fracimp, propres_imp = d_propres_imp, psi = 0.7)
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, times, data)
        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_food[g]) && (dist[3] <= epsilon_AMR[g]) && (!is.na(dist))) {
          # Store results
          res.new[i,]<- c(d_betaAA, d_betaHI, d_betaHD, d_phi, d_theta, d_alpha, d_zeta, d_fracimp, d_propres_imp)  
          # Calculate weights
          w1<-prod(c(sapply(c(1:5,7:9), function(b) dunif(res.new[i,b], min=lm.low[b], max=lm.upp[b])),
                     dbeta(res.new[i,6], 1.5, 8.5))) 
          if(g==1){
            
            w2<-1
            
          } else {
            w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=lm.low, upper=lm.upp)))
          }
          w.new[i] <- w1/w2
          # Update counter
          print(paste0('Generation: ', g, ", particle: ", i,", weights: ", w.new[i]))
          print(dist)
          i <- i+1
        }
      }
    }#
    sigma <- cov(res.new) 
    res.old <- res.new
    print(res.old)
    w.old <- w.new/sum(w.new)
    colnames(res.new) <- c("betaAA","betaHI","betaHD", "phi", "theta", "alpha", "zeta", "fracimp", "propres_imp")
    write.csv(res.new, file = paste("results_ABC_SMC_gen_tet_",g,".csv",sep=""), row.names=FALSE)
    ####
  }
}

N <- 50 #(ACCEPTED PARTICLES PER GENERATION)

lm.low <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
lm.upp <- c(0.2, 0.0001, 0.0001, 0.04, 2, 1, 0.15, 1, 1)

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=9,nrow=N)
res.new<-matrix(ncol=9,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)

epsilon_dist <- c(2, 1.5, 1, 0.75, 0.55)
epsilon_food <- c(3.26*0.2, 3.26*0.15, 3.26*0.1, 3.26*0.05, 3.26*0.01)
epsilon_AMR <- c(0.35*0.2, 0.35*0.15, 0.35*0.1, 0.35*0.05,  0.35*0.01)

ABC_algorithm(N = 50, 
              G = 5,
              sum.stats = summarystatprev, 
              distanceABC = sum_square_diff_dist, 
              fitmodel = amr, 
              tau_range = datatetraanim$pig_tetra_sales, 
              init.state = c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0), 
              times = seq(0, 2000, by = 100), 
              data = datatetraanim)

end_time <- Sys.time(); end_time - start_time

#### Test Data ####
data1 <- cbind(read.csv("results_ABC_SMC_gen_tet_1.csv", header = TRUE), "group" = "data1")
data2 <- cbind(read.csv("results_ABC_SMC_gen_tet_2.csv", header = TRUE), "group" = "data2")
data3 <- cbind(read.csv("results_ABC_SMC_gen_tet_3.csv", header = TRUE), "group" = "data3")
data4 <- cbind(read.csv("results_ABC_SMC_gen_tet_4.csv", header = TRUE), "group" = "data4") 
data5 <- cbind(read.csv("results_ABC_SMC_gen_tet_5.csv", header = TRUE), "group" = "data5") 

map_phi <- mean(data4[,"phi"]) 
map_theta <- mean(data4[,"theta"]) 
map_betaAA <- mean(data4[,"betaAA"]) 
map_betaHD <- mean(data4[,"betaHD"]) 
map_betaHI <- mean(data4[,"betaHI"]) 
map_alpha <- mean(data4[,"alpha"]) 
map_zeta <- mean(data4[,"zeta"]) 
map_fracimp <- mean(data4[,"fracimp"]) 
map_propres_imp <- mean(data4[,"propres_imp"]) 


plot(density(data4[,"betaHI"]))
