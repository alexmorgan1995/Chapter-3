library("deSolve"); library("ggplot2"); library("plotly"); library("reshape2")
library("bayestestR"); library("tmvtnorm"); library("ggpubr")

rm(list=ls())
setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha/inch")

# Posterior Distributions -------------------------------------------------

post_dist_names <- grep("complex",
                        list.files("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Model_Fit_Data/Part2/betaha/inch"), value = TRUE)

#final_amp_post <- read.csv(tail(grep(list.files(), pattern =  "amppigs_gen", value = TRUE),1))

post_dist <- lapply(post_dist_names, read.csv)
post_dist <- mapply(cbind, post_dist, "gen" = sapply(1:length(post_dist), function(x) paste0("gen_", x)), 
                    SIMPLIFY=F)
post_dist <- do.call("rbind", post_dist)

maps_est <- data.frame("Parameters" = colnames(post_dist[post_dist$gen == tail(unique(post_dist$gen),1),][1:8]), 
                       "MAP_Estimate" = colMeans(post_dist[post_dist$gen == tail(unique(post_dist$gen),1),][1:8]))

p_list <- list()

for(i in 1:(length(post_dist)-1)) {
  p_list[[i]] <- local ({
    name_exp <- post_dist[,c(i,9)]
    
    dens <- c()
    for (j in 1:8){
      dens[j] <- max(density(name_exp[name_exp$gen == unique(name_exp$gen)[j],1])[[2]])
    }
    
    
    max <- max(density(post_dist[,c(i,9)][,1])[[2]])
    
    p <- ggplot(name_exp, aes(x= name_exp[,1], fill=gen)) + geom_density(alpha=.5) +  theme_bw()  +
      geom_vline(xintercept = maps_est[i,2], size = 1.2, col = "red") +
      scale_x_continuous(expand = c(0, 0), name = c(expression(paste("Rate of Animal-to-Animal Transmission (", beta[AA], ")")),
                                                    expression(paste("Rate of Resistance Reversion (", phi, ")")),
                                                    expression(paste("Efficacy of Antibiotic-Mediated Recovery (", kappa, ")")), 
                                                    expression(paste("Antibiotic-Resistant Fitness Cost (", alpha, ")")), 
                                                    expression(paste("Background Infection Rate (", zeta, ")")),
                                                    expression(paste("Rate of Animal-to-Human Transmission (", beta[HA], ")")),
                                                    expression(paste("Proportion of Contaminated Imports (", Frac[Imp], ")")),
                                                    expression(paste("Proportion of Ampicillin-Resistant Cont Imports (", PropRes[Imp], ")")))[i]) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0,max(dens)*1.2), name = " ") +  
      theme(legend.text=element_text(size=10), axis.text.x=element_text(size=10),axis.ticks.y=element_blank(), axis.text.y=element_blank(),
            axis.title.y=element_text(size=10), axis.title.x= element_text(size=10), plot.margin = unit(c(0.25,0.4,0.15,0.55), "cm"),
            plot.title = element_text(size = 12, vjust = 3, hjust = 0.5, face = "bold")) + 
      labs(fill = "Generation")
    
    return(p)
  })
}

abc <- ggarrange(p_list[[1]], p_list[[2]], 
                 p_list[[3]], p_list[[4]], 
                 p_list[[5]], p_list[[6]], 
                 p_list[[7]], p_list[[8]], 
                 nrow = 4, ncol =2, 
                 font.label = c(size = 20), common.legend = TRUE, legend = "bottom",
                 align = "hv", vjust = 1.05)

ggsave(abc, filename = "ABC_SMC_Post_het.png", dpi = 300, width = 9, height = 8, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

# Pairs Plot --------------------------------------------------------------

final_amp_post <- read.csv(tail(grep(list.files(), pattern =  "complex", value = TRUE),1))

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

plot_amppig <- GGally::ggpairs(final_amp_post, lower=list(continuous=plot_lower), diag = list(continuous = plot_diag)) + theme_bw()

ggsave(plot_amppig, filename = "pairs_plot_amppig_hetero.png", dpi = 300, type = "cairo", width = 8, height = 8, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")

# Diagnostic Plots -------------------------------------------------------

ampRDS <- lapply(list.files(pattern = ".RDS"), readRDS)[[1]]

#Create a list of the RDSs

p_diag_list <- local({
  pre_plist <- list()
  #Distances
  
  dist_dat <- data.frame("dist" = sapply(1:8, function(x) ampRDS[[x]][[1]]), "gen" = sapply(1:8, function(x) paste0("gen", x)))
  dist_dat$accept_perc <- (1000/dist_dat$dist)*100
  
  #Summary Statistics 
  
  sum_diag <- as.data.frame(cbind(t(sapply(1:8, function(x) colMeans(ampRDS[[x]][[2]]))/c(1, 0.593, 0.2075134, 0.02865329, 0.4166667)),
                                  t(sapply(1:8, function(x) apply(ampRDS[[x]][[2]], 2, min))/c(1, 0.593, 0.2075134, 0.02865329, 0.4166667)),
                                  t(sapply(1:8, function(x) apply(ampRDS[[x]][[2]], 2, max))/c(1, 0.593, 0.2075134, 0.02865329, 0.4166667))))
  
  colnames(sum_diag) <- c("mean_dist","mean_IncH","mean_ResPropHum", "mean_LiveCont", "mean_ResPropAnim",
                          "low_dist","low_IncH","low_ResPropHum", "low_LiveCont", "low_ResPropAnim",
                          "high_dist","high_IncH","high_ResPropHum", "high_LiveCont", "high_ResPropAnim")
  
  sum_diag$gen <- seq(1,8)
  
  pre_plist[[1]] <- ggplot(sum_diag, aes(y = mean_dist, x = gen)) + 
    geom_ribbon(aes(ymin = low_dist, ymax = high_dist), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2)+
    scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), name = "Average Sum of Squared Distances") 
  
  pre_plist[[2]] <-  ggplot(sum_diag, aes(y = mean_IncH, x = gen)) + 
    geom_ribbon(aes(ymin = low_IncH, ymax = high_IncH), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2)+
    scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), limits = c(0,1), name = "Distance from Target Incidence Value")
  
  pre_plist[[3]] <- ggplot(sum_diag, aes(y = mean_ResPropHum, x = gen)) + 
    geom_ribbon(aes(ymin = low_ResPropHum, ymax = high_ResPropHum), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2) +
    scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), name = "Distance from Target ResPropHum Value") + coord_cartesian(ylim  = c(0,1))
  
  pre_plist[[4]] <-  ggplot(sum_diag, aes(y = mean_LiveCont, x = gen)) + 
    geom_ribbon(aes(ymin = low_LiveCont, ymax = high_LiveCont), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2)+
    scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), limits = c(0,1), name = "Distance from Target Livestock Cont Value")
  
  pre_plist[[5]] <- ggplot(sum_diag, aes(y = mean_ResPropAnim, x = gen)) + 
    geom_ribbon(aes(ymin = low_ResPropAnim, ymax = high_ResPropAnim), inherit.aes = TRUE , alpha = 0.5) + geom_line(size = 1.2) +
    scale_x_continuous(expand = c(0,0), name = "Generation") + scale_y_continuous(expand = c(0,0), name = "Distance from Target ResPropAnim Value") + coord_cartesian(ylim  = c(0,1))

  return(pre_plist)
})


diag_plots <- ggarrange(p_diag_list[[1]], p_diag_list[[2]], p_diag_list[[3]], p_diag_list[[4]],
                        p_diag_list[[5]], NULL,
                        ncol = 2, nrow = 3)

ggsave(diag_plots, filename = "diag_plots_heterofit.png", dpi = 300, width = 8, height = 11, units = "in",
       path = "//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Chapter_3/Models/Chapter-3/Figures")
