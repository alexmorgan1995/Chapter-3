library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2"); library("sensitivity")
library("bayestestR"); library("tmvtnorm"); library("ggpubr"); library("cowplot"); library("lhs")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data")

# Import in last generation -----------------------------------------------

amp <- read.csv(tail(list.files("C:/Users/amorg/Documents/PhD/Chapter_3/Models/fit_data", pattern = "ICOMBHTEST"), 1))



plot_list <- list()

for(i in 1:length(colnames(amp))) {
  plot_list[[i]] <- ggplot(amp, aes_string(x = colnames(amp)[i])) + geom_density(alpha = 0.5, fill = "red") +
    geom_vline(xintercept = c(1.887899e-02, 2.882880e-02, 4.492843, 2.862765e-01, 9.247680e-05, 9.660318e-05, 2.882880e-02, 9.473187e-05,
                              6.665432e-01, 2.032954e-02)[i], col = "red", size = 1.2, lty = 2) + scale_x_continuous(name = colnames(amp)[i])
}

do.call(ggarrange,plot_list)
