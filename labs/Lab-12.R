###############################################################
###Lab-12.R									###
###Classification and Regression Trees				###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: November 12, 2024						###
###############################################################

## Lab Objectives

		#Classification and Regression Trees (CART) are powerful tools in multivariate statistics, allowing us to model relationships 
			#between a dependent variable and one or more independent variables. In this lab, we will explore various tree-based 
			#methods, including classification trees, multi-state classification trees, regression trees, and random forests, using 
			#our class dataset of fish species abundance and environmental variables.

## Setting up the R Workspace and Preparing the Data

	#Load the necessary packages and prepare the dataset as we have done previously.

		library(vegan)
		library(rpart)
		library(rpart.plot)
		library(party)
		library(randomForest)
		library(ggplot2)

    #setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
    setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multivariate-analysis")
    source("Biostats.R")

    dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)

	#The CART framework is flexible enough that we do not need to worry about transforming or standardizing our data; however, 
		#you may choose to do so if you feel that it will improve the interpretation of your results. For instance, for our class 
		#dataset, we might choose to standardize fish abundances by sampling effort as we have done previously. For the sake of 
		#simplicity, I will use raw abundances in the examples below.

  		spp_N <- colSums(dat[,16:ncol(dat)])
  		spp_0 <- subset(spp_N, spp_N == 0)
  		omit <- names(spp_0)

  		dat2 <- dat[,!(colnames(dat) %in% omit)]
	
  		dat3 <- dat2[rowSums(dat2[,16:ncol(dat2)]) >0, ]
  
  		dat3$Herbaceous[is.na(dat3$Herbaceous)] <- 0 
  		dat3$Ann_Herb[is.na(dat3$Ann_Herb)] <- 0
	
  		dat3 <- dat3[complete.cases(dat3$SiteLength),]
  
  		dat_final <- dat3
	
  		fish <- dat_final[,16:ncol(dat_final)]
  		env <- dat_final[,1:15]
	
  		fish_red <- drop.var(fish, min.fo=1)
  
  		fish_dens <- fish_red
  		for(i in 1:nrow(fish_red)){	
    			fish_dens[i,] <- fish_red[i,]/dat_final$SiteLength[i]
  			}
  
  		fish_dens_log <- log(fish_dens + 1)
  
  		drop <- c("Latitude","Longitude","SiteLength","SiteWidth","SurfaceArea")
		env <- env[,!(colnames(env) %in% drop)]
		env_cont <- env[,!(colnames(env) %in% c("SMU","Pop","NLCD_Cat"))]
	
		env <- env[,!(colnames(env) %in% c("Ave_Max_D","Ann_Herb"))]
		env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]

## Classification Tree for Species Presence/Absence

	#A classification tree is a decision tree that is used to classify observations into predefined categories based on a set of predictor 
		#variables. The process involves recursively splitting the dataset into subsets based on the values of the predictor variables, 
		#which leads to a tree-like model structure. Each internal node represents a decision based on a predictor variable and each leaf 
		#node represents a classification outcome. The goal is to create a model that accurately predicts the target variable while 
		#minimizing misclassification.

	#Let's use redband trout presence/absence as an example:

		Trout_PA <- ifelse(fish$TROUT_RB > 0, "Present", "Absent")

		Trout.tree <- rpart(Trout_PA ~ ., data=env_cont, minsplit=2, xval=5)

		rpart.plot(Trout.tree)

	#The output of the classification tree shows how the data is partitioned based on the predictor variables. Each split indicates a 
		#decision rule, leading to final classifications of "Present" or "Absent." The tree visually represents the paths taken to 
		#reach these classifications, and the complexity of the tree indicates how many splits were necessary.

	#We can already see that this tree may be overfit (at the very least, it's challenging to interpret!), meaning it captures noise in 
		#the training data rather than the underlying distribution. This can lead to poor performance when applied to new data in a 
		#predictive capacity, so it's crucial to prune the tree to simplify the model while retaining predictive power.

	#The complexity parameter (cp) provides information about how we can best prune the tree. A smaller cp allows the tree to grow more 
		#complex, while a larger cp results in a simpler tree that may generalize better to unseen data. A good choice of cp for pruning 
		#is often the leftmost value for which the mean lies significantly below the horizontal line representing the 1 SE of the minimum.
		
		# 0.044 is the value to use according to this logic (also according to melanie, i asked her)
		
			plotcp(Trout.tree)
			

			Trout.tree.prune <- prune(Trout.tree, 0.076)
			rpart.plot(Trout.tree.prune)
			
			#trying with 0.044.... looks the same
			Trout.tree.prune <- prune(Trout.tree, 0.076)
			rpart.plot(Trout.tree.prune)

	#Now we can easily interpret how the branches in the tree lead to different classifications! Redband trout are more likely to be 
		#present in high canopy (\>43%), higher gradient (\>1.3%), and high elevation (\>1614 m) habitats.

### Multi-state Categorical Response

	#Classification trees can also handle multiple states or categories. In our case, we can extend the classification to include an 
		#"Abundant" category for sites with \>100 individuals.

		Trout_MS <- ifelse(fish$TROUT_RB > 0, "Present", "Absent")
		Trout_MS[fish$TROUT_RB > 100] = "Abundant"

		Trout.tree.ms <- rpart(Trout_MS ~ ., data=env_cont, minsplit=2, xval=5)
			rpart.plot(Trout.tree.ms)

	#Just like before, we will need to assess the complexity of this new tree and prune it accordingly to improve its predictive performance.

  		plotcp(Trout.tree.ms)

		Trout.tree.ms.prune <- prune(Trout.tree.ms, 0.03)
		rpart.plot(Trout.tree.ms.prune)

	#Notice that in the pruned version of the tree, the "Abundant" category does not appear prominently. This suggests that the splitting 
		#criteria may not have sufficiently distinct attributes to warrant a separate branch for this category, potentially indicating a 
		#lack of information gain from the predictors concerning abundance levels.

## Regression Trees

	#Regression trees, unlike classification trees, are used to predict a continuous outcome variable. They operate on similar principles 
		#as classification trees, partitioning the dataset into regions where the average response value is estimated based on the 
		#predictor variables. The tree structure helps to understand how predictors influence the response variable, and allows for 
		#complex, non-linear relationships.

			Trout.tree2 <- rpart(fish$TROUT_RB ~ ., data=env_cont, minsplit=2, xval=5)
			
			rpart.plot(Trout.tree2)

			plotcp(Trout.tree2)

			Trout.tree2.prune <- prune(Trout.tree2, 0.06)
			rpart.plot(Trout.tree2.prune)

	#The final, pruned output shows the splits based on our environmental variables, with terminal nodes indicating predicted values for 
		#redband trout abundance. By interpreting the pruned tree, we can discern which variables are most influential for determining 
		#trout abundance and how they interact. This tree essentially validates what we learned from the classification tree: low gradient 
		#(\<0.75%) and low canopy (\<2%) sites are less likely to support redband trout.

### Alternative using the `party` package

	#The `party` package offers an alternative approach to tree modeling through conditional inference trees. This method uses statistical 
		#tests to decide splits, ensuring that the selected predictors provide significant information gain at each node. This can be 
		#particularly useful for avoiding overfitting.

		Fish.rt <- ctree(fish$TROUT_RB ~ ., data=env_cont)
			plot(Fish.rt)

	#This time herbaceous cover (\> 21% \<) plays a stronger role in determining trout abundances, as does gradient (\>1.8%). I am not 
		#as satisfied with the performance of this tree; however, the box plots are a nice addition. If you wanted to, you could produce 
		#these box plots for a traditional regression tree using your own, custom R code.

## Random Forest Classification

	#Random forests are an ensemble learning method that builds multiple decision trees and merges them together to obtain a more accurate 
		#and stable prediction. They are particularly useful for handling overfitting issues that can arise with single trees, especially 
		#in high-dimensional spaces.

	#Random forests use bootstrapping to create subsets of the data and build multiple trees on these subsets, aggregating their 
		#predictions for final output. This approach enhances robustness and accuracy.

		Trout.forest <- randomForest(as.factor(Trout_PA) ~ ., data=env_cont, ntree = 5000, mtry = 5, importance=TRUE, keep.forest=FALSE, na.action=na.omit)
		Trout.forest

	#The output of the random forest model provides an overall accuracy measure and can be interpreted similarly to a single tree, but 
		#with the added benefit of reduced variance. When conducted using categorical data, the random forest function produces a 
		#confusion matrix of the data. The error rate of 17.29% we see for our data is quite high. An optimal error rate would be 10% 
		#or lower. We can see that most of the errors are occurring for "Absence" sites.

  		varImpPlot(Trout.forest)

	#The variable importance plot shows the relative importance of each environmental variable in the random forest model. This helps 
		#identify which features contribute most to the predictions. We already know that elevation, gradient, and canopy play a 
		#singificant role in determining redband trout abundances.

	#**Mean Decrease in Accuracy** shows how much the model's accuracy decreases when the variable is permuted (shuffled). Higher values 
		#indicate that the variable is more important for maintaining accuracy.

	#**Mean Decrease in Gini** measures how much each variable reduces the impurity (Gini impurity or entropy) on average across all trees. 
		#Larger values suggest greater importance.

	#To better visualize variable importance, we can create a more aesthetically pleasing plot that highlights the contribution of each 
		#predictor to the model's performance.

  		ForestData <- as.data.frame(importance(Trout.forest))
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
    				legend.position="bottom",
    				panel.grid.major.y = element_blank(),
    				panel.border = element_blank(),
    				axis.ticks.y = element_blank()
  				)

	#For regression tasks, random forests can similarly enhance predictive performance by averaging predictions from multiple trees.

		Trout.forest2 <- randomForest(fish$TROUT_RB ~ ., data=env_cont, ntree = 5000, mtry = 5, importance=TRUE, keep.forest=FALSE, na.action=na.omit)
		Trout.forest2
			
		ForestData <- as.data.frame(importance(Trout.forest2))
		ForestData <- ForestData[order(ForestData[,1]),]
		ForestData$Var.Names <- row.names(ForestData)
		colnames(ForestData) <- c("MeanDec","IncNodePurity","Var.Names")
		ForestData

		ggplot(ForestData, aes(x=Var.Names, y=MeanDec)) +
  			geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=MeanDec), color="skyblue") +
  			geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  			theme_light() +
  			coord_flip() +
  			theme(
    				legend.position="bottom",
    				panel.grid.major.y = element_blank(),
    				panel.border = element_blank(),
    				axis.ticks.y = element_blank()
  				)	

	#One limitation of the random forests procedure is that, because it is a tree-averaging procedure, it does not produce a tidy 
		#decision tree the way CART does. It can, however, be used in a predictive capacity to classify new data points. This 
		#makes CART/RF a great tool if you're looking to assign a classification system to new objects (e.g., assign habitat types to 
		#a given area based on environmental features, identify areas where certain species are/are not likely to occur).

