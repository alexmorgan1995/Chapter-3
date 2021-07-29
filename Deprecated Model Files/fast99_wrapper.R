rm(list=ls())
library("deSolve"); library("ggplot2"); library("sensitivity"); library("reshape2"); library("ggpubr")
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Figures")


# AMR Import Function + Rounding --------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - 
      usage_dom*(betaHA*Isa*Sh) - usage_dom*(1-alpha)*(betaHA*Ira*Sh) - 
      (1-usage_dom)*(betaHA*fracimp*(1-propres_imp)*Sh) - (1-usage_dom)*(1-alpha)*(betaHA*fracimp*propres_imp*Sh) - uh*Sh 
    
    dIsh = betaHH*Ish*Sh + usage_dom*betaHA*Isa*Sh + (1-usage_dom)*(betaHA*fracimp*(1-propres_imp)*Sh) - rh*Ish - uh*Ish 
    
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*usage_dom*(betaHA*Ira*Sh) + (1-usage_dom)*(1-alpha)*(betaHA*fracimp*propres_imp*Sh) - rh*Irh - uh*Irh 
    
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

# AMR Wrapper Function ----------------------------------------------------

ode_function_relfbd <- function(x) {
  times <- seq(0,3000, by = 1) 
  init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  tau_range <- c(0, 0.0123)
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms <- c(ra = x[z,1], rh = x[z,2], ua = x[z,3], uh = x[z,4], betaAA = x[z,5], betaHH = x[z,6], 
               betaHA = x[z,7], phi = x[z,8], theta = x[z,9], alpha = x[z,10], zeta = x[z,11], usage_dom = x[z,12], fracimp = x[z,13], propres_imp = x[z,14])
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
  times <- seq(0,3000, by = 1) 
  init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
  tau_range <- c(0, 0.0123)
  print(x)
  
  return_vec <- vector()
  
  for(z in 1:nrow(x)) {
    
    parms <- c(ra = x[z,1], rh = x[z,2], ua = x[z,3], uh = x[z,4], betaAA = x[z,5], betaHH = x[z,6], 
               betaHA = x[z,7], phi = x[z,8], theta = x[z,9], alpha = x[z,10], zeta = x[z,11], usage_dom = x[z,12], fracimp = x[z,13], propres_imp = x[z,14])
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

    return_vec[z] <- as.numeric(temp1[1,2]/temp1[2,2])
    
    if(is.nan(temp1[1,2]/temp1[2,2])) {
      return_vec[z] <-  1
    }
    
    #print(paste0(temp1[1,2], " , ", temp1[2,2]))
    #print(return_vec[z])
  }
  
  
  return(return_vec)
  #c(temp1[1,1]/temp1[2,1], temp1[1,2]/temp1[2,2], temp1[1,1] - temp1[2,1], temp1[1,2] - temp1[2,2])
  #relFBD, relRES, deltaFBD, deltaRES
}

# eFAST -------------------------------------------------------------------

factors <- c("ra","rh", "ua", "uh", "betaAA", "betaHH", "betaHA", "phi", "theta", "alpha","zeta" ,"usage_dom" , "fracimp","propres_imp")

start_time <- Sys.time()


testfbd <- fast99(model = ode_function_relfbd, factors = factors, n = 1000, 
               q.arg = list(list(min=600^-1, max=6^-1), 
                            list(min=55^-1, max=0.55^-1), 
                            list(min=2400^-1, max=24^-1),
                            list(min=288350^-1, max=2883.5^-1),
                            list(min=0.0029, max=0.29),
                            list(min=0.000001, max=0.0001),
                            list(min=0.000001, max=0.0001),
                            list(min=0.00131, max=0.131),
                            list(min=0.113, max=11.3),
                            list(min=0.001, max=1),
                            list(min=0.00497, max=0.497),
                            list(min=0.001, max=1),
                            list(min=0.001, max=1),
                            list(min=0.001, max=1)))


testres <- fast99(model = ode_function_relres, factors = factors, n = 1000, 
               q.arg = list(list(min=600^-1, max=6^-1), 
                            list(min=55^-1, max=0.55^-1), 
                            list(min=2400^-1, max=24^-1),
                            list(min=288350^-1, max=2883.5^-1),
                            list(min=0.0029, max=0.29),
                            list(min=0.000001, max=0.0001),
                            list(min=0.000001, max=0.0001),
                            list(min=0.00131, max=0.131),
                            list(min=0.113, max=11.3),
                            list(min=0.001, max=1),
                            list(min=0.00497, max=0.497),
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
  geom_col(color = "black",position= "stack") + scale_x_discrete(expand = c(0, 0.5)) + scale_y_continuous(limits = c(0,0.6), expand = c(0, 0)) +
  theme(legend.position=c(0.25, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5)) + 
  scale_fill_manual(labels = c("Total Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = "eFAST total/first order sensitivity indices for RelFBD",
       x ="", y = "Sensitivity Index") 

p_efast_res <- ggplot(plotdata_res, aes(fill = variable, x =  X, y = value)) + theme_bw() + 
  geom_col(color = "black",position= "stack") + scale_x_discrete(expand = c(0, 0.5)) + scale_y_continuous(limits = c(0,0.6), expand = c(0, 0)) +
  theme(legend.position=c(0.25, 0.875), legend.text=element_text(size=12), legend.title = element_blank(), axis.text=element_text(size=12), 
        axis.title.y=element_text(size=12), axis.title.x= element_text(size=12), plot.margin = unit(c(0.35,1,0.35,1), "cm"),
        legend.spacing.x = unit(0.3, 'cm'), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5)) + 
  scale_fill_manual(labels = c("Total Order", "First Order"), values = c("lightgrey", "darkgrey")) +
  labs(title = "eFAST total/first order sensitivity indices for RelRes",
       x ="", y = "Sensitivity Index") 


#Plotting rel's

PRCC_plot <- ggarrange(p_efast_fbd, p_efast_res, nrow = 2, ncol = 1,
                       align = "v", labels = c("A","B"), font.label = c(size = 20), common.legend = TRUE,
                       legend = "bottom") 

ggsave(PRCC_plot, filename = "eFAST_relative.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in")
