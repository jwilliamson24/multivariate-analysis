# -------------------------------------------------------------------------------------------------------
##
## randomforest.R 
##
## Random forest code, continued effort from multivariate class
##
## Jasmine Williamson
## Date Created: 01-15-2025
##
## 
## settings -----------------------------------------------------------------------------------------------
    
    rm(list=ls())
    setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multivariate-analysis")

    library(randomForest)
    library(ggplot2)
 
## load data--------------------------------------------------------------------------------------------------

# site-level data
    dat <- readRDS("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
    row.names(dat) <- dat[,1]
    
    sals <- dat[26:27]
    env <- dat[1:25]
    
    drop <- c("lat","long","stand","tree_farm","landowner","site_id","year","weather")
    env <- env[,!(colnames(env) %in% drop)]
    
    env_cont <- env[,-1]
    
    #took out elevation, jul date (we already know theyre elevationally/temporally specific) 
    drop <- c("jul_date","elev")
    env_subset <- env_cont[,!(colnames(env_cont) %in% drop)]

# sal presences absence
    oss_PA <- ifelse(sals$oss > 0, "Present", "Absent")
    enes_PA <- ifelse(sals$enes > 0, "Present", "Absent")
    
## oss random forest --------------------------------------------------------------------------------------------------

#all variables    
    oss.forest <- randomForest(as.factor(oss_PA) ~ ., data=env_cont, ntree = 5000, mtry = 5, 
                               importance=TRUE, keep.forest=FALSE, na.action=na.omit)
    oss.forest
    table(oss_PA)
    importance(oss.forest)
    # The error rate is pretty bad (29.92%), and most of that is due to the absences (62%). I tried tuning the mtry 
    # and ntree to improve the matrix but I havent figured out a way to make it better yet. 
    
#subset of variables
    oss.forest.sub <- randomForest(as.factor(oss_PA) ~ ., data=env_subset, ntree = 5000, mtry = 5, 
                               importance=TRUE, keep.forest=FALSE, na.action=na.omit)
    
## full variables plot -----------------------------------------------------------------------------------------------------

#with full variables
    varImpPlot(oss.forest)
    
    ossForestData <- as.data.frame(importance(oss.forest))
    ossForestData <- ossForestData[order(ossForestData[,1]),]
    ossForestData$Var.Names <- row.names(ossForestData)
    colnames(ossForestData) <- c("Absent","Present","MeanDec","IncNodePurity","Var.Names")
    ossForestData

    #ggplot     
    p1 <- ggplot(ossForestData, aes(x = Var.Names, y = MeanDec)) +
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
        text = element_text(size = 14)
      ) +
      labs(
        title = "Variable Importance from Random Forest Model",
        x = "Environmental Variables",
        y = "Mean Decrease in Accuracy",
        size = "Node Purity"
      )
 
ggsave(filename = "oss_varimp.png", plot = p, device = "png",
       path = )


## subset variables plot --------------------------------------------------------------------------------------------------------
    
    ossForestData.sub <- as.data.frame(importance(oss.forest.sub))
    ossForestData.sub <- ossForestData.sub[order(ossForestData.sub[,1]),]
    ossForestData.sub$Var.Names <- row.names(ossForestData.sub)
    colnames(ossForestData.sub) <- c("Absent","Present","MeanDec","IncNodePurity","Var.Names")
    
    #ggplot     
    p2 <- ggplot(ossForestData.sub, aes(x = Var.Names, y = MeanDec)) +
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
        text = element_text(size = 14)
      ) +
      labs(
        title = "Variable Importance from Random Forest Model",
        x = "Environmental Variables",
        y = "Mean Decrease in Accuracy",
        size = "Node Purity"
      )
    
## enes random forest --------------------------------------------------------------------------------------------------
    
    enes.forest <- randomForest(as.factor(enes_PA) ~ ., data=env_cont, ntree = 5000, mtry = 5, 
                               importance=TRUE, keep.forest=FALSE, na.action=na.omit)
    enes.forest
    table(enes_PA)
    importance(enes.forest)  
    
    # The error rate is  
 
       
#subset of variables
    enes.forest.sub <- randomForest(as.factor(oss_PA) ~ ., data=env_subset, ntree = 5000, mtry = 5, 
                                   importance=TRUE, keep.forest=FALSE, na.action=na.omit)
    
## full variables plot ------------------------------------------------------------------------------------------------
    
    varImpPlot(enes.forest)
    
    enesForestData <- as.data.frame(importance(enes.forest))
    enesForestData <- enesForestData[order(enesForestData[,1]),]
    enesForestData$Var.Names <- row.names(enesForestData)
    colnames(enesForestData) <- c("Absent","Present","MeanDec","IncNodePurity","Var.Names")
    enesForestData

#variable importance ggplot    
    p3 <- ggplot(enesForestData, aes(x = Var.Names, y = MeanDec)) +
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
        text = element_text(size = 14)
      ) +
      labs(
        title = "Variable Importance from Random Forest Model",
        x = "Environmental Variables",
        y = "Mean Decrease in Accuracy",
        size = "Node Purity"
      )

    
## subset variables plot -----------------------------------------------------------------------------------------------------

    enesForestData.sub <- as.data.frame(importance(enes.forest.sub))
    enesForestData.sub <- enesForestData.sub[order(enesForestData.sub[,1]),]
    enesForestData.sub$Var.Names <- row.names(enesForestData.sub)
    colnames(enesForestData.sub) <- c("Absent","Present","MeanDec","IncNodePurity","Var.Names")
    
    #ggplot     
    p4 <- ggplot(enesForestData.sub, aes(x = Var.Names, y = MeanDec)) +
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
        text = element_text(size = 14)
      ) +
      labs(
        title = "Variable Importance from Random Forest Model",
        x = "Environmental Variables",
        y = "Mean Decrease in Accuracy",
        size = "Node Purity"
      )
    
        