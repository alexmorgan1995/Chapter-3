library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_2/Chapter2_Fit_Data/FinalData")

#Function to remove negative prevalence values and round large DP numbers
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}

# Single Model ------------------------------------------------------------

amrimp <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (betaAH*Ish*Sa) - (1-alpha)*(betaAH*Irh*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + betaAH*Ish*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAH*Irh*Sa + (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - (1-usage_dom)*(betaimp*fracimp*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh + (betaimp*(fracimp*(1-propres_imp))*Sh) - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) + (betaimp*(fracimp*propres_imp)*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIsa,dIra,dSh,dIsh,dIrh)))
  })
}


# Running the Basic Model -------------------------------------------------

#I want to do a baseline run and a model run where I test sampling the uniform distribution of %Domestic and Uniform proportion of FBD in importing patches 

#Baseline 

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)
tauseq <- c(seq(0,0.02, by = 0.001), 0.0123)

output1 <- data.frame(matrix(nrow = length(tauseq), ncol = 5)); colnames(output1) <- c("tau","SensH", "ResH","ICombH", "ResPropH")

for (i in 1:length(tauseq)) {
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, tau = tauseq[i], zeta = 0.0497,
             usage_dom = 0.6,
             betaimp = 0.00001,
             fracimp = 0.5,
             propres_imp = 0.5)
  
  out <- ode(y = init, func = amrimp, times = times, parms = parms2)
  temp <- c(parms2[["tau"]],
            rounding(out[nrow(out),6]),
            rounding(out[nrow(out),7]),
            rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]),
            rounding(out[nrow(out),7]) / (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
  output1[i,] <- temp
}

fin_base <- c(abs(output1[1,3] - output1[output1$tau == 0.0123,3]),
              abs(output1[1,4] - output1[output1$tau == 0.0123,4]))


meltoutpu1 <- melt(output1, id.vars = "tau", measure.vars = c("ResH","SensH"))

ggplot(meltoutpu1, aes(x = tau, y = value, fill = variable)) + geom_bar(position="stack", stat="identity")

#Testing Distributions

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)
output_fin <- data.frame(matrix(nrow = 0, ncol = 24))
plotlist <- list()

parmtau <- seq(0, 0.035, 0.001)

n = 1000


for(z in 1:3) {
  
  parms2 = c(ra = 60^-1, rh =  (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaAH = 0.00001, betaHH = 0.00001, 
             betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, tau = tauseq[i], zeta = 0.0497,
             betaimp = 0.00001,
             fracimp = 0,
             propres_imp = 0)
  
  for(i in 1:n) {
    
    parms[grep("tau_imp", names(parms))] <- rbeta(10, 1, 4)
    parms[grep("psi_imp", names(parms))] <- rbeta(10, 1, 1)
    
    parms2 <- parms
    
    for(k in 1:length(parms[grep("tau_imp", names(parms))])) {
      parms2[grep("tau_imp", names(parms))][k] <- 
        (parms[grep("tau_imp", names(parms))][k]/  sum(parms[grep("tau_imp", names(parms))])) * 
        ((parms[["tau_dom"]] * length(parms[grep("tau_imp", names(parms))])) - parms[["tau_dom"]])
      
      parms2[grep("psi_imp", names(parms))][k] <- 
        (parms[grep("psi_imp", names(parms))][k]/  sum(parms[grep("psi_imp", names(parms))])) * (1-parms[["psi_dom"]])
    }
    
    output1 <- data.frame(matrix(nrow = 0, ncol = 2))
    
    for (j in 1:2) {
      parms2[["tau_dom"]] <- c(0, parms[["tau_dom"]])[j]
      out <- ode(y = init, func = amrmulti, times = times, parms = parms2)
      temp <- c(rounding(out[nrow(out),36]) + rounding(out[nrow(out),37])*100000,
                rounding(out[nrow(out),37]) / (rounding(out[nrow(out),36]) + rounding(out[nrow(out),37])))
      output1[j,] <- temp
    }
    
    fin <- c(parms2[grep("tau", names(parms))], parms2[grep("psi", names(parms))],
             abs(output1[1,1] - output1[2,1]),
             abs(output1[1,2] - output1[2,2]))
    
    output_fin[i,] <- fin
    
    print(paste0((i/n)*100, " - psi_dom = ", parms2[["psi_dom"]]))
  }
  
  colnames(output_fin) <- c(names(parms2[grep("tau", names(parms))]), names(parms2[grep("psi", names(parms))]),
                            "diff_food", "diff_resis")
  
  p1 <- ggplot(data = output_fin, aes(x = diff_food, y = diff_resis)) + geom_point() + theme_bw() + 
    geom_point(x = fin_base[1]*100000, y = fin_base[2], size = 2, col = "red") + scale_x_continuous(limits = c(0,1.5)) + scale_y_continuous(limits = c(0,1)) +
    theme(legend.text=element_text(size=14),  axis.text=element_text(size=14),
          axis.title.y=element_text(size=14),axis.title.x= element_text(size=14), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
  
  plotlist[[z]] <- list(p1, output_fin)
}

