library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Figures")

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
    
    colnames(output)[1:8] <- c("FBD_H_0", "FBD_H_123","Res_H_0", "Res_H_123", "delta_FBD", "delta_Res", "ParmValue", "Parm")
    print(output)
    
    plotnames <- c(bquote(Imp["FBD"]~Parameter), bquote(Imp["Res"]~Parameter), bquote(psi["Dom"]~Parameter),
                   bquote(beta["AA"]~Parameter), bquote(beta["HA"]~Parameter), bquote(beta["HH"]~Parameter),
                   bquote(phi~Parameter), bquote(theta~Parameter), bquote(alpha~Parameter), bquote(zeta~Parameter), bquote(r["H"]~Parameter), bquote(r["A"]~Parameter), 
                   bquote(mu["H"]~Parameter), bquote(mu["A"]~Parameter))[[j]]
    
    p1 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(delta_FBD))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(0, max(output$delta_FBD)*1.1), expand = c(0, 0)) +
      labs(x = plotnames) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.title.y=element_blank())
    
    p2 <- ggplot(output, aes(x = as.numeric(ParmValue), y = as.numeric(delta_Res))) + theme_bw() + geom_line(lwd = 1.02, col ="darkblue") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous( expand = c(0, 0)) +
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
    
    return(list(p1,p2, output))
  })
}

#deltaFBD
pdelta_FBD <- plot_grid(plot_grid(suppplotlist[[1]][[1]], suppplotlist[[2]][[1]], suppplotlist[[3]][[1]],suppplotlist[[4]][[1]], suppplotlist[[5]][[1]], 
                               suppplotlist[[6]][[1]], suppplotlist[[7]][[1]], suppplotlist[[8]][[1]], suppplotlist[[9]][[1]], suppplotlist[[10]][[1]], suppplotlist[[11]][[1]],
                               suppplotlist[[12]][[1]], suppplotlist[[13]][[1]], suppplotlist[[14]][[1]], nrow = 5, ncol =3), scale=0.95) + 
  draw_label("delta_FBD", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)

pdelta_res <- plot_grid(plot_grid(suppplotlist[[1]][[2]], suppplotlist[[2]][[2]], suppplotlist[[3]][[2]],suppplotlist[[4]][[2]], suppplotlist[[5]][[2]], 
                                  suppplotlist[[6]][[2]], suppplotlist[[7]][[2]], suppplotlist[[8]][[2]], suppplotlist[[9]][[2]], suppplotlist[[10]][[2]], suppplotlist[[11]][[2]],
                                  suppplotlist[[12]][[2]], suppplotlist[[13]][[2]], suppplotlist[[14]][[2]], nrow = 5, ncol =3), scale=0.95) + 
  draw_label("delta_Res", x=  0, y=0.5, vjust= 1.5, angle=90, size = 12)


ggsave(pdelta_FBD, filename = "delta_FBD_parm_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")
ggsave(pdelta_res, filename = "delta_Res_parm_mono.png", dpi = 300, type = "cairo", width = 5, height = 7, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures")


# LHS ---------------------------------------------------------------------

#First thing is to select how many times we want to sample (h) - this is based on the number of "partitions" in our distribution 

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, usage_dom = 0.6, fracimp = 0.5, propres_imp = 0.5)
#This is without tau

h <- 500
lhs <- maximinLHS(h, length(parms))

#this generates a scaling factor (0,1) - using a uniform distribution for every parameter -random values are sampled from each subsection (of h sections - going vertically)
#Uniform distribution can be transformed into any distribution using q... function (e.g qnorm) - different columns can have different distributions 
#the qnorm function - you give it a probability (from the LHS matrix) - put in the parameters of the distribution - and then it gives you the  

parmdetails <- data.frame("parms" = c("fracimp", "propres_imp" ,"usage_dom","betaAA" ,"betaHA" ,"betaHH" ,"phi" ,"theta",
                                 "alpha", "zeta", "rh", "ra" ,"uh" , "ua" ),
                          "lbound" = c(0, 0, 0, 0.0029, 0.000001, 0.000001,  0.00131, 0.113,
                                       0, 0.00497, 55^-1, 600^-1, 288350^-1, 2400^-1),
                          "ubound" = c(1, 1, 1, 0.29, 0.0001, 0.0001,  0.131, 11.3,
                                       1, 0.497, 0.55^-1, 6^-1, 2883.5^-1, 24^-1))

#I want to multiply each column with the difference between the lower and upperbound for the particular parameter of interest 

lhsscaled <- data.frame(matrix(nrow = nrow(lhs), ncol = ncol(lhs)))
colnames(lhsscaled) <- unique(parmdetails[,1])

for(i in 1:length(parmdetails[,1])) {
  lhsscaled[,i] <- lhs[,i]*(parmdetails[i,3] - parmdetails[i,2]) + parmdetails[i,2]
}

init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)
modelrunlhs <- data.frame(matrix(nrow = h, ncol = nrow(parmdetails) + 2)) #+2 because I also want to keep track of 2 extra outcome measures
colnames(modelrunlhs) <- c(unique(parmdetails[,1]), "deltaFBD", "deltaRes")

for(i in 1:h) {
  temptau <- matrix(nrow = 2, ncol = 2)
  parmslhs <- as.list(lhsscaled[i,])

  for(j in 1:2) { 
    parmslhstau <- c(parmslhs, "tau" =  c(0, 0.0123)[j])
    out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
    temp <- c((rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]))*100000,
              rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])))
    temptau[j,] <-temp
  }
  modelrunlhs[i,] <- c(parmslhs, temptau[2,1] - temptau[1,1], temptau[1,2] - temptau[2,2])  
  
  print(paste0(round((i/h)*100, digits = 2),"%" ))
}


prcc_fbd <- pcc(modelrunlhs[,1:14], modelrunlhs[,15], nboot = 100, rank=TRUE)
prcc_res <- pcc(modelrunlhs[,1:14], modelrunlhs[,16], nboot = 100, rank=TRUE)

#Plotting Delta_FBD PRCC
plotdf_fbd <- data.frame("parm" = rownames(prcc_fbd[[7]]), as.data.frame(prcc_fbd[[7]]))
plotdf_fbd$parm <- factor(plotdf_fbd$parm, levels = plotdf_fbd$parm) 

p_fbd <- ggplot(plotdf_fbd, aes(x = parm, y = original)) + geom_hline(yintercept = 0, col ="red", size = 1.05) + geom_point(stat = "identity", size = 3) + 
  theme_bw() + geom_errorbar(aes(ymin=min..c.i., ymax=max..c.i.), width=.1) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
                                                                  axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15)) +
  labs(title = bquote(bold("Sensitivity Analysis of " ~ Delta["FBD"])), x ="Model Parameters", y = "PRCC")

#Plotting Delta_RES PRCC
plotdf_res <- data.frame("parm" = rownames(prcc_res[[7]]), as.data.frame(prcc_res[[7]]))
plotdf_res$parm <- factor(plotdf_res$parm, levels = plotdf_res$parm) 

p_res <- ggplot(plotdf_res, aes(x = parm, y = original)) + geom_hline(yintercept = 0, col ="red", size = 1.05) + geom_point(stat = "identity", size = 3) + 
  theme_bw() + geom_errorbar(aes(ymin=min..c.i., ymax=max..c.i.), width=.1) +
  scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
                                                                  axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15)) +
  labs(title = bquote(bold("Sensitivity Analysis of " ~ Delta["RES"])), x ="Model Parameters", y = "PRCC")

PRCC_plot <- ggarrange(p_fbd, p_res, nrow = 2, ncol = 1,
                      align = "v", labels = c("A","B"), font.label = c(size = 20)) 

ggsave(PRCC_plot, filename = "LHS_PRCC.png", dpi = 300, type = "cairo", width = 10, height = 10, units = "in")

#Scatter Plot for LHS and Outcome Measure
plot(modelrunlhs$alpha,modelrunlhs$deltaRes)
t2 <- data.frame("rankparm" = rank(modelrunlhs$betaHA),
                 "delta_fbd" = modelrunlhs$deltaFBD)

plot(t2$rankparm,t2$delta_fbd)

# Uncertainty Analysis  ---------------------------------------------------
#This section will involve me overlaying the relationship between antibiotic usage and human resistance and foodborne disease
#With the every single model run obtained from the LHS sample 

#Run model for baseline parameters
#Run LHS parameter set and make sure I track fro both human FDB and Res
#Overlay 

parms <- c(ra = 60^-1, rh = (5.5^-1), ua = 240^-1, uh = 28835^-1, betaAA = 0.029, betaHH = 0.00001, 
           betaHA = (0.00001), phi = 0.0131, theta = 1.13, alpha = 0.43, zeta = 0.0497, usage_dom = 0.6, fracimp = 0.5, propres_imp = 0.5)

tau <- seq(0,0.02, by = 0.0005)
init <- c(Sa=0.98, Isa=0.01, Ira=0.01, Sh=1, Ish=0, Irh=0)
times <- seq(0, 200, by = 1)

#Baseline Run 
base_anal <- data.frame(matrix(nrow = length(tau), ncol = 3))
colnames(base_anal) <- c("FBD", "Res", "tau")

for(j in 1:length(tau)) {
  parmslhstau <- c(parms, "tau" =  tau[j])
  out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
  temp <- c(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]),
            rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])),
            tau[j])
  base_anal[j,] <- temp
  print(paste0(round((j/length(tau))*100, digits = 2),"%" ))
}

#Uncertainty Analysis


parmdetails <- data.frame("parms" = c("fracimp", "propres_imp" ,"usage_dom","betaAA" ,"betaHA" ,"betaHH" ,"phi" ,"theta",
                                      "alpha", "zeta", "rh", "ra" ,"uh" , "ua" ),
                          "lbound" = c(0, 0, 0, 0.0029, 0.000001, 0.000001,  0.00131, 0.113,
                                       0, 0.00497, 55^-1, 600^-1, 288350^-1, 2400^-1),
                          "ubound" = c(1, 1, 1, 0.29, 0.0001, 0.0001,  0.131, 11.3,
                                       1, 0.497, 0.55^-1, 6^-1, 2883.5^-1, 24^-1))
h <- 500; lhs <- maximinLHS(h, length(parms))
lhsscaled <- data.frame(matrix(nrow = nrow(lhs), ncol = ncol(lhs)))
colnames(lhsscaled) <- unique(parmdetails[,1])
for(i in 1:length(parmdetails[,1])) {
  lhsscaled[,i] <- lhs[,i]*(parmdetails[i,3] - parmdetails[i,2]) + parmdetails[i,2]
}

uncert_anal <- data.frame(matrix(nrow = h*length(tau), ncol = 4)) #+2 because I also want to keep track of 2 extra outcome measures
colnames(uncert_anal) <- c("FBD", "Res", "tau","run_n")
t <- 0

for(i in 1:h) {
  parmslhs <- as.list(lhsscaled[i,])
  
  for(j in 1:length(tau)) { 
    t <- t + 1 
    parmslhstau <- c(parmslhs, "tau" =  tau[j])
    out <- ode(y = init, func = amrimp, times = times, parms = parmslhstau)
    temp <- c(rounding(out[nrow(out),6]) + rounding(out[nrow(out),7]),
              rounding(out[nrow(out),7])/ (rounding(out[nrow(out),6]) + rounding(out[nrow(out),7])),
              tau[j],
              i)
    uncert_anal[t,] <- temp
  }
  print(paste0(round((i/h)*100, digits = 2),"%" ))
}

#Combining the two 
uncert_anal$group <- "LHS Runs"; base_anal$run_n <- uncert_anal$run_n[nrow(uncert_anal)] +1; base_anal$group <- "Baseline"
combdata <- rbind(uncert_anal, base_anal)

#Plotting

p1_FBD <- ggplot(combdata, aes(x = tau, y = FBD, group = run_n, col = group, size = group, alpha = group)) + geom_line() +
  scale_y_continuous(limits = c(0, 0.003), expand = c(0, 0)) +  scale_x_continuous(expand = c(0, 0)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
        axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15),
        legend.position = "bottom", legend.text = element_text(size=14), legend.title = element_blank()) +
  labs(title = "Uncertainty Analysis (FBD) with LHS runs", x ="Livestock Antibiotic Usage", y = "Human Foodborne Disease") +
  scale_color_manual(values = c("red", "black")) + scale_size_manual(values = c(1.2, 1)) + scale_alpha_manual(values = c(1, 0.2))

p2_FBD <- ggplot(combdata, aes(x = tau, y = FBD, group = run_n, col = group, size = group, alpha = group)) + geom_line() +
  scale_y_continuous(limits = c(0.00003, 0.00005), expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), 
        axis.text=element_text(size=10), axis.title = element_blank(), title = element_blank(),legend.position = "none") +
  labs(title = "Uncertainty Analysis (FBD) with LHS runs", x ="Livestock Antibiotic Usage", y = "Human Foodborne Disease") +
  scale_color_manual(values = c("red", "black")) + scale_size_manual(values = c(1.2, 1)) + scale_alpha_manual(values = c(1, 0.2))

p_FBD_comb <- p1_FBD + annotation_custom(ggplotGrob(p2_FBD), xmin = 0.0075, xmax = 0.02, 
                                   ymin = 0.002, ymax = 0.003)

p_Res <- ggplot(combdata, aes(x = tau, y = Res, group = run_n, col = group, size = group, alpha = group)) + geom_line() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +  scale_x_continuous(expand = c(0, 0)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm"), axis.text.x = element_text(angle = 45, vjust = 0.65, hjust=0.5),
        axis.text=element_text(size=14), axis.title =element_text(size=14), title = element_text(size=15),
        legend.position = "bottom", legend.text = element_text(size=14), legend.title = element_blank()) +
  labs(title = "Uncertainty Analysis (Res) with LHS runs", x ="Livestock Antibiotic Usage", y = "Proportion of Resistance") +
  scale_color_manual(values = c("red", "black")) + scale_size_manual(values = c(1.2, 1)) + scale_alpha_manual(values = c(1, 0.2))


unc_plot <- ggarrange(p_FBD_comb, p_Res, nrow = 2, ncol = 1,
                       align = "v", labels = c("A","B"), font.label = c(size = 20)) 

ggsave(unc_plot, filename = "uncertainty_anal.png", dpi = 300, type = "cairo", width = 10, height = 12, units = "in")

