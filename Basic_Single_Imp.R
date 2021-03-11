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
    dSa = ua + ra*(Isa + Ira) + theta*tau*Isa - (betaAA*Isa*Sa) - (1-alpha)*(betaAA*Ira*Sa) - ua*Sa -
      zeta*Sa*(1-alpha) - zeta*Sa 
    dIsa = betaAA*Isa*Sa + phi*Ira - theta*tau*Isa - tau*Isa - ra*Isa - ua*Isa + zeta*Sa
    dIra = (1-alpha)*betaAA*Ira*Sa + tau*Isa - phi*Ira - ra*Ira - ua*Ira + zeta*Sa*(1-alpha)
    
    dSh = uh + rh*(Ish+Irh) - (betaHH*Ish*Sh) - (1-alpha)*(betaHH*Irh*Sh) - (betaHA*Isa*Sh) - (1-alpha)*(betaHA*Ira*Sh) - (1-usage_dom)*(betaHA*fracimp*Sh) - uh*Sh 
    dIsh = betaHH*Ish*Sh + betaHA*Isa*Sh + (betaHA*(betaHA*(1-propres_imp))*Sh) - rh*Ish - uh*Ish 
    dIrh = (1-alpha)*(betaHH*Irh*Sh) + (1-alpha)*(betaHA*Ira*Sh) + (betaHA*(betaHA*propres_imp)*Sh) - rh*Irh - uh*Irh 
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


# Testing for Monotonicity  -----------------------------------------------

#The aim of this section is to look at the relationship of changing each variable on delta_FBD and delta_res


parmdetails <- rbind(data.frame("Parameter" = "fracimp", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "propres_imp", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "usage_dom", "Value" = seq(0, 1, by = 1/100)),
                     
                     data.frame("Parameter" = "betaAA", "Value" = seq(0, 0.29, by = 0.29/100)),
                     data.frame("Parameter" = "betaHA", "Value" = seq(0, 0.0001, by = 0.0001/100)),
                     data.frame("Parameter" = "betaHH", "Value" = seq(0, 0.0001, by = 0.0001/100)),
                     data.frame("Parameter" = "phi", "Value" = seq(0, 0.131, by = 0.131/100)),
                     data.frame("Parameter" = "theta", "Value" = seq(0, 11.3, by = 11.3/100)),
                     data.frame("Parameter" = "alpha", "Value" = seq(0, 1, by = 1/100)),
                     data.frame("Parameter" = "zeta", "Value" = seq(0, 0.497, by = 0.497/100)),
                     data.frame("Parameter" = "rh", "Value" = seq(0.01, 0.55^-1, by = 0.55^-1/100)),
                     data.frame("Parameter" = "ra", "Value" = seq(0, 6^-1, by = 6^-1/100)),
                     data.frame("Parameter" = "uh", "Value" = seq(0, 2883.5^-1, by = 2883.5^-1/100)),
                     data.frame("Parameter" = "ua", "Value" = seq(0, 24^-1, by = 24^-1/100)))

times <- seq(0,30000, by = 100) 
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
tau_range <- c(0, 0.0123)


parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaAH = 0.00001, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, usage_dom = 0.6, fracimp = 0.5, propres_imp = 0.5)

suppplotlist <- list()

for (j in 1:length(unique(parmdetails[,1]))) { 
  
  suppplotlist[[j]] <- local ({ 
    output <- data.frame()
    
    for (x in 1:length(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]),2])) { #for the individual parameter values in the sequence
      temp1 <- data.frame()
      
      for (i in 1:length(tau_range)) {
        temp <- data.frame(matrix(nrow = 0, ncol=3))
        parmstemp <- c(parms, tau = tau_range[i])
        parmstemp[as.character(unique(parmdetails[,1])[j])] <- parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]
        out <- ode(y = init, func = amrimp, times = times, parms = parmstemp)
        temp[1,1] <- (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000
        temp[1,2] <- rounding(out[nrow(out),7])/(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))
        temp[1,3] <- as.character(parmdetails[parmdetails == as.character(unique(parmdetails[,1])[j]), 2][x]) #what is the parameter value used
        temp[1,4] <- as.character(unique(parmdetails[,1])[j]) # what is the parameter explored 
        temp1 <- rbind.data.frame(temp1, temp)
      }
      
      output <- rbind(output, stringsAsFactors = FALSE,
                      c(as.numeric(temp1[1,1]), #FBD at tau = 0
                        as.numeric(temp1[2,1]), #FBD at tau = 0.0123
                        as.numeric(temp1[1,2]), #ResH at tau = 0
                        as.numeric(temp1[2,2]), #ResH at tau = 0.0123 
                        as.numeric(temp1[1,1] - temp1[2,1]), #delta_FBD
                        as.numeric(temp1[1,2] - temp1[2,2]), #delta_Res
                        as.numeric(temp1[i,3]), #Parm Value
                        as.factor(temp1[i,4]))) #Parm Name 
      
      print(paste0("Parameter ",unique(parmdetails[,1])[j], " | ", round(x/101, digits = 2)*100,"%" ))
    }
    
    colnames(output)[1:7] <- c("FBD_H_0", "FBD_H_123","Res_H_0", "Res_H_123", "delta_FBD", "delta_Res", "ParmValue", "Parm")
    
    output <- output[!is.nan(output$PercInc) & !is.nan(output$RelInc) & !is.infinite(output$PercInc),] # Remove 0s and NAs - models which don't lift off 
    
    plotnames <- c(bquote(Imp["FBD"]~Parameter), bquote(Imp["Res"]~Parameter), bquote(psi["Dom"]~Parameter),
                   bquote(beta["AA"]~Parameter), bquote(beta["HA"]~Parameter), bquote(beta["HH"]~Parameter), bquote(beta["AH"]~Parameter), 
                   bquote(phi~Parameter), bquote(theta~Parameter), bquote(alpha~Parameter), bquote(zeta~Parameter), bquote(r["H"]~Parameter), bquote(r["A"]~Parameter), 
                   bquote(mu["H"]~Parameter), bquote(mu["A"]~Parameter))[[j]]
    
    p1 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(delta_FBD))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, max(output$delta_FBD)*1.1), expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    p2 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(delta_Res))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, max(output$delta_Res)*1.1), expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    #if(unique(parmdetails[,1])[j] == "betaHA"){
    #  p2 <- p2 + geom_vline(xintercept = 8.0e-06, col = "red", size  = 0.7, lty = 3)}
    #if(unique(parmdetails[,1])[j] == "zeta"){
    #  p2 <- p2 + geom_vline(xintercept = 0.011084193, col = "red", size  = 0.7, lty = 3)}
    #if(unique(parmdetails[,1])[j] == "ra"){
    #  p2 <- p2 + geom_vline(xintercept = 0.038333333, col = "red", size  = 0.7, lty = 3)}
    #if(unique(parmdetails[,1])[j] == "ua"){
    #  p2 <- p2 + geom_vline(xintercept = 0.0254166667, col = "red", size  = 0.7, lty = 3)}
    #if(unique(parmdetails[,1])[j] == "rh"){
     # p2 <- p2 + geom_vline(xintercept = 0.22818182, col = "red", size  = 0.7, lty = 3)}
    
    return(list(p1,p2))
  })
}

#Absolute Diff
pabdiff <- plot_grid(plot_grid(suppplotlist[[1]][[1]], suppplotlist[[2]][[1]], suppplotlist[[3]][[1]],suppplotlist[[4]][[1]], suppplotlist[[5]][[1]], 
                               suppplotlist[[6]][[1]], suppplotlist[[7]][[1]], suppplotlist[[8]][[1]], suppplotlist[[9]][[1]], suppplotlist[[10]][[1]], suppplotlist[[11]][[1]],
                               suppplotlist[[12]][[1]], nrow = 4, ncol =3), scale=0.95) + 
  draw_label("% Change in ICombH Relative to Baseline Usage", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

ggsave(pabdiff, filename = "Sensitivity_RelInc.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_2/Figures/Redraft_v1")


