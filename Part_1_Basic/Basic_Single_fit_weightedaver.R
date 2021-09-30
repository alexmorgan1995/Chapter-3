library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs")
library("rootSolve")

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
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - 
      psi*(betaHD*Isa*Sh) - 
      psi*(1-alpha)*(betaHD*Ira*Sh) - 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh = betaHH*Ish*Sh + psi*betaHD*Isa*Sh + 
      (1-psi)*(betaHI*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*psi*(betaHD*Ira*Sh) + 
      (1-psi)*(1-alpha)*(betaHI*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
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

#Use the mean for the EU as the parameters (minus the UK) - only the main importers 

#EU_cont <- mean(country_data_imp$Foodborne_Carriage_2019[2:10])
#EU_cont <- sum(country_data_imp[-1,"Foodborne_Carriage_2019"] * as.numeric(country_data_imp[-1,"Normalised_Usage_2018"]))
EU_cont <- sum(country_data_imp[-1,"Foodborne_Carriage_2019"] * 
                 as.numeric(country_data_imp[-1,"Normalised_Usage_2018"]))/sum(as.numeric(country_data_imp[-1,"Normalised_Usage_2018"]))

#EU_res <- mean(country_data_imp$Prop_Amp_Res [2:10])
EU_res <- sum(country_data_imp[-1,"Prop_Amp_Res"] * 
                as.numeric(country_data_imp[-1,"Normalised_Usage_2018"]))/sum(as.numeric(country_data_imp[-1,"Normalised_Usage_2018"]))

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
  tauoutput <- data.frame(matrix(nrow = length(tau_range), ncol = 5))
  
  for (i in 1:length(tau_range)) {
    
    parms2 = thetaparm
    parms2["tau"] = tau_range[i]
    
    out <- runsteady(y = init.state, func = fitmodel, parms = parms2, times = c(0, Inf))
    
    tauoutput[i,] <- c(tau_range[i], 
                       out[[1]][["Isa"]] + out[[1]][["Ira"]],
                       (out[[1]][["Ish"]] + out[[1]][["Irh"]])*100000,
                       out[[1]][["Ira"]] / (out[[1]][["Isa"]] + out[[1]][["Ira"]]),
                       out[[1]][["Irh"]] /  (out[[1]][["Ish"]] + out[[1]][["Irh"]]))
  }

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
  prod(sapply(1:8, function(a) as.numeric((par[a]-lm.low[a]) > 0) * as.numeric((lm.upp[a]-par[a]) > 0)))
}


ABC_algorithm <- function(N, G, sum.stats, distanceABC, fitmodel, tau_range, init.state, data) {
  N_ITER_list <- list()
  
  fit_parms <- c("betaAA", "betaHH", "betaHD", "betaHI", 
                 "phi", "kappa", "alpha", "zeta")
  thetaparm <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, psi = 0.656,
                 fracimp = EU_cont, propres_imp = EU_res)
  
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
        
        thetaparm[fit_parms] <- c(c(d_betaAA, d_betaHH, d_betaHD, d_betaHI, d_phi, d_kappa, d_alpha, d_zeta))
        
        dist <- computeDistanceABC_ALEX(sum.stats, distanceABC, fitmodel, tau_range, thetaparm, init.state, data)
        print(dist)

        if((dist[1] <= epsilon_dist[g]) && (dist[2] <= epsilon_foodH[g]) && (dist[3] <= epsilon_AMRH[g]) && 
           (dist[4] <= epsilon_foodA[g]) && (dist[5] < 0.95) && (dist[6] <= epsilon_AMRA[g]) && (!is.na(dist))) {
          
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
    write.csv(res.new, file = paste("TEST_WEIGHT_PART1_Simple_FIT_AMP_",g,".csv",sep=""), row.names=FALSE)
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
                           data = country_data_gen)

end_time <- Sys.time(); end_time - start_time

saveRDS(dist_save, file = "test_weightdist_simple_amp.rds")

# Check Posterior ---------------------------------------------------------

amp_post <- do.call(rbind,
                    lapply(list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data", pattern = "TEST_WEIGHT_PART1"), read.csv))
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


# Attribution Plot --------------------------------------------------------

post_amp <- read.csv(tail(list.files(path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Data/fit_data", pattern = "TEST_WEIGHT_PART1")[1:6], 1))
MAP_amp <- map_estimate(post_amp)

MAP_amp <- data.frame("Parameters" = colnames(post_amp), "MAP_Estimate" = colMeans(post_amp))



#Baseline 

parmtau <- seq(0,0.01, by = 0.0005)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
icombhdata <- data.frame(matrix(ncol = 8, nrow = 0))
times <- seq(0, 200000, by = 100)

for(j in 1:3) {
  output1 <- data.frame(matrix(ncol = 8, nrow = length(parmtau)))
  for (i in 1:length(parmtau)) {
    temp <- data.frame(matrix(nrow = 1, ncol =8))
    
    parms2 = c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, 
               
               betaAA = MAP_amp[MAP_amp$Parameter == "betaAA", 2], betaHH = MAP_amp[MAP_amp$Parameter == "betaHH", 2]*0.4, tau = parmtau[i],
               betaHD = MAP_amp[MAP_amp$Parameter == "betaHD", 2], betaHI = MAP_amp[MAP_amp$Parameter == "betaHI", 2], phi = MAP_amp[MAP_amp$Parameter == "phi", 2], 
               kappa = MAP_amp[MAP_amp$Parameter == "kappa", 2], alpha = MAP_amp[MAP_amp$Parameter == "alpha", 2], zeta = MAP_amp[MAP_amp$Parameter == "zeta", 2], 
               psi = c(0.656, 1, 0.1)[j], fracimp = 0.5, propres_imp = EU_res)
    
    out <- ode(y = init, func = amrimp, times = times, parms = parms2)
    temp[,1] <- parmtau[i]
    temp[,2] <- rounding(out[nrow(out),5]) 
    temp[,3] <- rounding(out[nrow(out),6]) 
    temp[,4] <- rounding(out[nrow(out),7])
    temp[,5] <- temp[3] + temp[4]
    temp[,6] <- signif(as.numeric(temp[4]/temp[5]), digits = 3)
    temp[,7] <- rounding(out[nrow(out),4]) / (rounding(out[nrow(out),3]) + rounding(out[nrow(out),4]))
    temp[,8] <- c("baseline","import_none", "import_90")[j]
    output1[i,] <- temp
  }
  icombhdata <- rbind(icombhdata, output1)
}

colnames(icombhdata) <- c("tau", "SuscHumans","InfHumans","ResInfHumans","ICombH","IResRat","IResRatA", "group")

plotdata <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[1],],
                 id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_base <- ggplot(plotdata, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() + 
  geom_vline(xintercept = UK_amp, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.0005) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,10), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[1]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF")) +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)", title = "Food from Domestic Sources = 65.6% (Baseline)") 


plotdata_imp_none <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[2],],
                          id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_imp_none <- ggplot(plotdata_imp_none, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() + 
  geom_vline(xintercept = UK_amp, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,10), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[2]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF"))  +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)", title = "Food from Domestic Sources = 90%") 

plotdata_imp_90 <- melt(icombhdata[icombhdata$group == unique(icombhdata$group)[3],],
                        id.vars = c("tau"), measure.vars = c("ResInfHumans","InfHumans")) 

p_imp_90 <- ggplot(plotdata_imp_90, aes(fill = variable, x = tau, y = value*100000)) + theme_bw() + 
  geom_vline(xintercept = UK_amp, alpha = 0.3, size = 2) + 
  geom_col(color = "black",position= "stack", width  = 0.001) + scale_x_continuous(expand = c(0, 0.0005)) + 
  scale_y_continuous(limits = c(0,10), expand = c(0, 0))  + 
  geom_text(label= c(round(icombhdata$IResRat[icombhdata$group == unique(icombhdata$group)[3]],digits = 2),rep("",length(parmtau))),vjust=-0.5, hjust = 0.05,
            position = "stack", angle = 45) +
  theme(legend.position=c(0.75, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm')) + 
  scale_fill_manual(labels = c("Antibiotic-Resistant Infection", "Antibiotic-Sensitive Infection"), values = c("#F8766D", "#619CFF"))  +
  labs(x ="Generic Antibiotic Usage (g/PCU)", y = "Infected Humans (per 100,000)", title = "Food from Domestic Sources = 10%") 

ggarrange(p_base, p_imp_none, p_imp_90, nrow  = 3, ncol = 1, common.legend = TRUE, legend = "bottom")
