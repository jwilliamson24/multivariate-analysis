#
# Multivariate final project code
# fall 2024
# a little messy, trying different things
#

rm(list=ls())
setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multivariate-analysis")

library(vegan)
library(rpart)
library(rpart.plot)
library(party)
library(randomForest)
library(ggplot2)
library(viridis)
source("biostats.R")

#site-level data
#dat <- readRDS("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
dat <- readRDS("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
row.names(dat) <- dat[,1]

#dat2 <- subset(dat, year=="2024")
sals <- dat[26:27]
env <- dat[1:25]

drop <- c("lat","long","stand","tree_farm","landowner","site_id","year","weather")
env <- env[,!(colnames(env) %in% drop)]

env_cont <- env[,-1]

drop <- c("jul_date","veg_cov","fwd_cov","dwd_count","size_cl","decay_cl","char_cl","length_cl" )
env_subset <- env_cont[,!(colnames(env_cont) %in% drop)]

env_std <- decostand(env_cont, "standardize") #Z-scores the data in each column

### classification tree -------------------------------------------------------------------------

    # OSS presence/absence
    oss_PA <- ifelse(sals$oss > 0, "Present", "Absent")
    
    oss.tree <- rpart(oss_PA ~ ., data=env_subset, minsplit=2, xval=5)
    rpart.plot(oss.tree)
    plotcp(oss.tree)
    
    oss.tree.prune <- prune(oss.tree, 0.044)
    rpart.plot(oss.tree.prune)
    
    # using all env continuous variables
    oss.tree2 <- rpart(oss_PA ~ ., data=env_cont, minsplit=2, xval=5)
    rpart.plot(oss.tree2)
    plotcp(oss.tree2)
    
    oss.tree.prune2 <- prune(oss.tree2, 0.055)
    rpart.plot(oss.tree.prune2, cex=1.3)


### random forest -------------------------------------------------------------------------

    oss.forest <- randomForest(as.factor(oss_PA) ~ ., data=env_cont, ntree = 5000, mtry = 5, 
                               importance=TRUE, keep.forest=FALSE, na.action=na.omit)
    oss.forest  
    importance(oss.forest)
    varImpPlot(oss.forest)
    
    ForestData <- as.data.frame(importance(oss.forest))
    ForestData <- ForestData[order(ForestData[,1]),]
    ForestData$Var.Names <- row.names(ForestData)
    colnames(ForestData) <- c("Absent","Present","MeanDec","IncNodePurity","Var.Names")
    ForestData
    
    ggplot(ForestData, aes(x=Var.Names, y=MeanDec)) +
      geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=MeanDec), color="skyblue") +
      geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
      theme_light() +
      coord_flip() +
      theme(
        legend.position = c(1, 1), 
        legend.justification = c(1, 1),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 16)
      ) +
    labs(
      title = "Variable Importance from Random Forest Model",
      x = "Environmental Variables",
      y = "Mean Decrease in Accuracy",
      size = "Node Purity"
    )

# plot with importance size scaled by color
    ggplot(ForestData, aes(x = Var.Names, y = MeanDec)) +
      geom_segment(aes(x = Var.Names, xend = Var.Names, y = 0, yend = MeanDec, 
                       color = ifelse(MeanDec > 18, "Above 18", 
                                      ifelse(MeanDec >= 10, "10-18", "Below 10"))),
                   show.legend = FALSE) +
      geom_point(aes(size = IncNodePurity, 
                     color = ifelse(MeanDec > 18, "Above 18", 
                                    ifelse(MeanDec >= 10, "10-18", "Below 10"))), alpha = 0.6,
                 show.legend = FALSE) +
      theme_light() +
      coord_flip() +
      scale_color_manual(values = c("Above 18" = "blue", "10-18" = "#4CAF50", "Below 10" = "#FFC107")) +
      labs(color = "Range") +
    theme(
      legend.position = c(1, 1), 
      legend.justification = c(1, 1),
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      text = element_text(size = 14)
    ) +
      labs(
        title = "Variable Importance from Random Forest Model",
        x = "Environmental Variables",
        y = "Mean Decrease in Accuracy",
        size = "Node Purity"
      )
    
    # above plot with legend
    ggplot(ForestData, aes(x = Var.Names, y = MeanDec)) +
      geom_segment(aes(x = Var.Names, xend = Var.Names, y = 0, yend = MeanDec, 
                       color = ifelse(MeanDec > 18, "Above 18", 
                                      ifelse(MeanDec >= 10, "10-18", "Below 10"))), 
                   show.legend = FALSE) +
      geom_point(aes(size = IncNodePurity, 
                     color = ifelse(MeanDec > 18, "Above 18", 
                                    ifelse(MeanDec >= 10, "10-18", "Below 10"))), 
                 alpha = 0.6) +
      theme_light() +
      coord_flip() +
      scale_color_manual(values = c("Above 18" = "blue", "10-18" = "#66C2A5", "Below 10" = "#FFD700")) +
      labs(color = "Range", size = "Node Purity") +
      theme(
        text = element_text(size = 16)
      ) +
      labs(
        title = "Variable Importance from Random Forest Model",
        x = "Environmental Variables",
        y = "Mean Decrease in Accuracy",
        size = "Node Purity"
      )




### treatment classification efforts -------------------------------------------------------------------------

### tbRDA
    # these might be messy, original code in homework 6 rmd
    
    #site.sc <- scores(tbRDA, scaling=2, display="sites")
    #spe.sc <- scores(tbRDA, scaling=2, display="species")
    #env.sc <- scores(tbRDA, scaling=2, display="bp")
    
    #env_data_hellinger <- sqrt(scale(env_subset, center=TRUE, scale=TRUE))
    
    tbRDA <- rda(sal.hel ~ ., env_std)
    summary(tbRDA)
    
    plot(tbRDA, 
         choices=c(1,2), 
         type="n", 
         scaling=2, 
         main="tbRDA of Hellinger-transformed Sal Count",
         xlab="RDA 1", 
         ylab="RDA 2")	
    for (i in 1:length(groups)){
      dim_choice <- site.sc[rownames(site.sc) %in% rownames(env[env$trt == groups[i], ]), ]
      #dim_choice <- site.sc$sites[env$trt==groups[i],]
      points(dim_choice[,1], dim_choice[,2], 
             pch=19, 
             cex=1.4, 
             col=pt_col[i])
    }
    text(env.sc$biplot*2, 
         row.names(env.sc$biplot), 
         col="darkgray")
    arrows(0, 0, env.sc$biplot[,1]*1.7, env.sc$biplot[,2]*1.7, 
           lwd=2, 
           length=0.1, 
           col="darkgray")
    text(spe.sc$species*2, 
         row.names(spe.sc$species))
    arrows(0, 0, spe.sc$species[,1]*1.5, spe.sc$species[,2]*1.5, 
           lwd=2, 
           length=0.1)
    legend(x="bottomright", 
           legend=groups, 
           col=pt_col, 
           pch=19)


# code from chat gpt trying to fix plotting errors 
    # Plot RDA sites
    plot(tbRDA, 
         choices=c(1,2), 
         type="n", 
         scaling=1, 
         main="tbRDA of Hellinger-transformed Sal Count",
         xlab="RDA 1", 
         ylab="RDA 2") 
    
    # Add points for each group
    for (i in 1:length(groups)) {
      dim_choice <- site.sc[rownames(site.sc) %in% rownames(env[env$trt == groups[i], ]), ]
      points(dim_choice[, 1], dim_choice[, 2], 
             pch=19, 
             cex=1.4, 
             col=pt_col[i])
    }
    
    # Add text and arrows for environmental variables
    text(env.sc[ , 1] * 2, 
         env.sc[ , 2] * 2, 
         row.names(env.sc), 
         col="darkgray")
    arrows(0, 0, env.sc[ , 1] * 1.7, env.sc[ , 2] * 1.7, 
           lwd=2, length=0.1, col="darkgray")
    
    # Add text and arrows for species
    # Assuming spe.sc is a matrix with species coordinates in columns 1 and 2
    text(spe.sc[, 1] * 2, 
         spe.sc[, 2] * 2, 
         row.names(spe.sc))
    arrows(0, 0, spe.sc[, 1] * 1.5, spe.sc[, 2] * 1.5, 
           lwd=2, length=0.1)
    
    # Add legend for groups
    legend(x="bottomright", 
           legend=groups, 
           col=pt_col, 
           pch=19)


### PCoA -------------------------------------------------------------------------
# original code in lab 7 copy 

    env_std_subset <- env_std[,c("temp","dwd_cov","soil_moist","stumps","logs","decay_cl","canopy_cov")]	
    env_euc <- vegdist(env_std_subset, method="euclidean")
    
    env.pcoa <- cmdscale(env_euc, 
                         k=5, 
                         eig=TRUE, 
                         add=T)

    spe.sc <- wascores(env.pcoa$points[,1:2], env_subset)
    
    vec.sp <- envfit(as.data.frame(env.pcoa$points), env_subset, perm=1000)

#look for eigenvalues info
    100*env.pcoa$eig/sum(env.pcoa$eig)
    
    # Extract eigenvalues
    eigenvalues <- env.pcoa$eig
    
    # Calculate the percentage of variance explained by each axis
    percent_explained <- (eigenvalues / sum(eigenvalues)) * 100
    
    # View the percentages for the first few axes
    percent_explained[1:5]
    
    
# plot  
    groups <- levels(factor(dat$trt))
    pt_col <- viridis(length(groups))
    site.sc <- scores(env.pcoa, choices=c(1,2))
    
    plot(site.sc[, 1:2],  # First two dimensions
         main = "Env PCoA", 
         xlab = "PCoA 1", 
         ylab = "PCoA 2", 
         pch = 19)
    
    for (i in 1:length(groups))
    {
      dim_choice <- site.sc[dat$trt==groups[i],]
      points(dim_choice[,1], dim_choice[,2], 
             pch=19, 
             cex=1.4, 
             col=pt_col[i])
      
      # Calculate and add convex hull for the group
      chull_points <- dim_choice[chull(dim_choice[, 1], dim_choice[, 2]), ]
      
      polygon(chull_points[, 1], chull_points[, 2], 
              border=pt_col[i], 
              col=adjustcolor(pt_col[i], alpha.f = 0.3),  # Semi-transparent fill
              lwd=2)
      }
    
    text(spe.sc*1.5, row.names(spe.sc),
         cex = 1.5)
    arrows(0, 0, spe.sc[,1]*1.4, spe.sc[,2]*1.4, 
           lwd=2, 
           length=0.1,
           cex= 1.5)
    legend(x="bottomleft", 
           legend=levels(factor(dat$trt)), 
           col=pt_col[1:6], 
           pch=19, 
           cex=1.2)






