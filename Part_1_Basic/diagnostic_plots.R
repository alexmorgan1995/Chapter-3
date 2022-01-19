library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1")

# Pairs Plot --------------------------------------------------------------

setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1")

post_amppigs <- read.csv(tail(grep("amppigs_gen", list.files("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1"), value = TRUE), 1))

plot_lower <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + scale_x_continuous(expand = c(0,0))  + scale_y_continuous(expand = c(0,0)) + 
    stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
    scale_fill_gradientn(colours=viridis::viridis(100))
  return(p)
}

plot_diag <- function(data, mapping){
  p <- ggplot(data = data, mapping = mapping) + scale_x_continuous(expand = c(0,0))  + 
    geom_density(fill = "grey", alpha = 0.3, size = 1.2) + theme_bw()
  return(p)
}

plot_amppig <- GGally::ggpairs(post_amppigs, lower=list(continuous=plot_lower), diag = list(continuous = plot_diag)) + theme_bw()

ggsave(plot_amppig, filename = "pairs_plot_amppigs.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures/New_Figures")

# Diagnostic Plots -------------------------------------------------------

setwd("C:/Users/amorg/Documents/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part1")

tetRDS <- readRDS(list.files(pattern = ".rds"))

#Create a list of the RDSs

p_diag_list <- list()

pre_plist <- list()
#Distances

case_RDS <- tetRDS

dist_dat <- data.frame(dist = sapply(1:8, function(x) case_RDS[[x]][[1]]), "gen" = sapply(1:8, function(x) paste0("gen", x)))
dist_dat$accept_perc <- (1000/dist_dat$dist)*100

#Summary Statistics 

sum_diag <- as.data.frame(cbind(t(sapply(1:8, function(x) colMeans(case_RDS[[x]][[2]]))/c(1, 0.593,0.2075134, 0.02865329, 0.4166667)),
                                t(sapply(1:8, function(x) apply(case_RDS[[x]][[2]], 2, min))/c(1, 0.593,0.2075134, 0.02865329, 0.4166667)),
                                t(sapply(1:8, function(x) apply(case_RDS[[x]][[2]], 2, max))/c(1, 0.593,0.2075134, 0.02865329, 0.4166667))))

colnames(sum_diag) <- c("mean_dist","mean_ICombH","mean_ResPropHum", "mean_ICombA", "mean_ResPropAnim",
                        "low_dist","low_ICombH","low_ResPropHum", "low_ICombA","low_ResPropAnim",
                        "high_dist","high_ICombH","high_ResPropHum","high_ICombA","high_ResPropAnim")

sum_diag$gen <- seq(1,8)

pre_plist[[1]] <- ggplot(sum_diag, aes(y = mean_dist, x = gen)) + 
  geom_ribbon(aes(ymin = low_dist, ymax = high_dist), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2)+
  scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), name = "Average Sum of Squared Distances") 

pre_plist[[2]] <-  ggplot(sum_diag, aes(y = mean_ICombH, x = gen)) + 
  geom_ribbon(aes(ymin = low_ICombH, ymax = high_ICombH), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2)+
  scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), limits = c(0,1), name = "Distance from Target ICombH Value")

pre_plist[[3]] <- ggplot(sum_diag, aes(y = mean_ResPropHum, x = gen)) + 
  geom_ribbon(aes(ymin = low_ResPropHum, ymax = high_ResPropHum), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2) +
  scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), name = "Distance from Target ResPropHum Value") + coord_cartesian(ylim  = c(0,1))

pre_plist[[4]] <-  ggplot(sum_diag, aes(y = mean_ICombA, x = gen)) + 
  geom_ribbon(aes(ymin = low_ICombA, ymax = high_ICombA), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2)+
  scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), limits = c(0,1), name = "Distance from Target ICombA Value")

pre_plist[[5]] <- ggplot(sum_diag, aes(y = mean_ResPropAnim, x = gen)) + 
  geom_ribbon(aes(ymin = low_ResPropAnim, ymax = high_ResPropAnim), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2) +
  scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), name = "Distance from Target ResPropAnim Value") + coord_cartesian(ylim  = c(0,1))


diag_plots <- ggarrange(pre_plist[[1]], pre_plist[[2]], pre_plist[[3]], 
                        pre_plist[[4]], pre_plist[[5]], ncol = 3, nrow = 2)

diag_plots <- annotate_figure(diag_plots, top = text_grob("Ampicillin Usage in Fattening Pigs", 
                                      color = "black", size = 16))


ggsave(diag_plots, filename = "diag_plots.png", dpi = 300, type = "cairo", width = 10, height = 8, units = "in",
       path = "C:/Users/amorg/Documents/PhD/Chapter_3/Figures/New_Figures")
