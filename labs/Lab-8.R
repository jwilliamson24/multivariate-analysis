###############################################################
###Lab-8.R									###
###Data Visualization and Interpretation				###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 24, 2024						###
###############################################################

## Lab Objectives

	#It is essential to effectively interpret and communicate ordination results. Understanding the relationships 
		#among species, sites, and environmental variables can be tricky, so our goal here is to develop creative 
		#strategies for representing these complex patterns. Effective data visualization allows us to reveal 
		#structure, identify outliers, and communicate insights clearly. In this lab, we will learn to:

		#-   Creatively present ordination output to highlight ecological trends
		#-   Use color, scaling, and other graphical techniques to add interpretive depth
		#-   Combine clustering and ordination to reveal hidden patterns in the data
		#-   Integrate external predictor variables to enhance interpretation of multivariate results

## Setting up the R Workspace and Preparing the Data

	#Load the necessary packages and prepare the dataset as we have done previously.

  		library(vegan)
  		library(viridis)
  		library(cluster)

      setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
      source("Biostats.R")
      source("coldiss.R")
      
      dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
  		sub_dat <- subset(dat, SMU=="Malheur")

	#As you may have guessed, we're going to subset the data to include only the Malheur sites.

  		spp_N <- colSums(sub_dat[,16:ncol(sub_dat)])
  		spp_0 <- subset(spp_N, spp_N == 0)
  		omit <- names(spp_0)

  		dat2 <- sub_dat[,!(colnames(sub_dat) %in% omit)]
	
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

## Variance Partitioning

	#As we already know, Principal Component Analysis (PCA) reduces dimensionality while preserving as much 
		#of the data’s variance as possible. In this step, we will examine how much variance in the dataset 
		#is explained by each principal component and/or each descriptor.

	#Each axis represents a certain percentage of the total variance. By focusing on the first few axes, we 
		#can simplify complex datasets for easier interpretation.

	#Run the PCA on the Hellinger-transformed fish density dataset:

		fish.hel <- decostand(fish_dens, "hellinger")
		fish.pca <- prcomp(fish.hel, scale=FALSE)	
			summary(fish.pca)

		pca.eigenvec(fish.pca, dim=6, digits=3, cutoff=0.1)

		site.sc <- scores(fish.pca)

	#Here's the basic biplot (we'll get to fancy visualizations in a second). What phenomenon do we see 
		#creeping through?

		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data", 
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		points(site.sc, 
			pch=20, 
			cex=1.4)
		text(site.sc, 
			rownames(dat_final), 
			pos=4, 
			cex=0.7)
		arrows(0, 0, fish.pca$rotation[,1]*0.8, fish.pca$rotation[,2]*0.8, 
			lwd=2, 
			length=0.1)
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))

	#Remember that the loadings are the square of the eigenvectors. Loadings allow us to see which descriptors 
		#(in this case, species) have the strongest influence on the formation of the principal components. 
		#We will focus on axes 1 and 2, as they explain the most variance.

		loadings <- fish.pca$rotation^2
		round(loadings, 2)

	#What do the loadings tell us about each species influence on the principal axes, especially axes 1 and 2? 
		#Interpretation of loadings helps identify which species are driving the patterns seen in the ordination 
		#plot. Species with high loadings on the first axis may be associated with environmental gradients, such 
		#as elevation or habitat type.

## Data Visualization

	#Now let's have some fun (as if we weren't having a great time in this class already)! Properly visualizing 
		#your data makes it much easier to conduct a thorough exploratory analysis and interpret your results. 
		#Here are some tips and tricks for making simple, visually pleasing, and easy to interpret biplots using 
		#the `ordiplot` function. This is intended for folks with a basic understanding of R. If you're an R whiz, 
		#feel free to get creative on your own! Don’t forget to share your code so everyone can replicate and learn 
		#from each other!

	#Option 1: Varying point size by descriptors

	#The `ordiplot` function is part of the `vegan` package and is used for plotting ordination results, particularly 
		#from PCA, NMDS, or CA. It offers flexibility in overlaying points, arrows, and text to visually represent 
		#your data in a way that emphasizes the relationships among variables.

	#In this example, we will vary the point size by redband trout abundance. To make the plot more interpretable, 
		#we first log-scale the data so that the point sizes reflect abundance intuitively. Before plotting, it’s 
		#a good idea to check the point size distribution with a histogram.

  		RBT_scale <- 1 + log(dat_final$TROUT_RB + 1)
    		hist(RBT_scale)

	#Next, we create a biplot with points sized according to trout abundance:

		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		points(site.sc, 
			pch=20, 
			cex=RBT_scale, 
			col="#50505050")
		text(site.sc, 
			rownames(dat_final), 
			pos=4, 
			cex=0.8)
		arrows(0, 0, fish.pca$rotation[,1]*0.8, fish.pca$rotation[,2]*0.8, 
			lwd=2, 
			length=0.1)
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))
		
		
		
		oss_scale <- 1 + log(dat2$oss + 1)
		hist(oss_scale)
		
		
		points(site.sc, 
		       pch=20, 
		       cex=RBT_scale, 
		       col="#50505050")

	#In this visualization, trout densities align well with the PC1 axis, demonstrating the strong influence of 
		#redband trout on that principal component.

	#Option 2: Color by Cluster

	#Let's timewarp back to a few weeks ago when we were learning about cluster analysis. Did you know ordination 
		#and clustering can be used in tandem (Legendre & Legendre 10.1 covers this nicely)? Combining ordination 
		#and clustering methods can reveal deeper insights into the structure of your data. Here, we combine cluster 
		#analysis with PCA to group our data by similarity and display it visually.

	#First, run the cluster analysis using Hellinger distances (which are the Euclidean distances among Hellinger-
		#transformed data points).
	
		fish.dist <- vegdist(fish.hel, "euclidean")
		fishcl.ward <- hclust(fish.dist, method = "ward.D2")

		plot(fishcl.ward, main="Ward's Minimum Variance Dendrogram", 
			xlab="Sites", 
			ylab="Hellinger Distance", 
			hang=-1)
			rect.hclust(fishcl.ward, k=5)

		fish_clus <- cutree(fishcl.ward, k=5)

	#Note that I only recommend this approach when you're truly going into your dataset blind. A lot of times, 
		#we already know something about the structure of our data, whether it be habitat type, population, 
		#treatment group, etc. It often makes more sense to interpret your data according to those known categorical 
		#predictors (which we'll do below!).

		groups <- levels(factor(fish_clus))
		pt_col <- viridis(length(groups))
		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[fish_clus==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
				}
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))

	#We can see that clusters divide nicely across the principal axes, indicating distinct groupings in the dataset.

	#For added visual clarity, dendrogram lines can be added:

		groups <- levels(factor(fish_clus))
		pt_col <- viridis(length(groups))
		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[fish_clus==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
				}
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))
		ordicluster(p, fishcl.ward, col="gray50")

	#Or even S.D. ellipses for each cluster (check out the `ordihull` function, too):

		groups <- levels(factor(fish_clus))
		pt_col <- viridis(length(groups))
		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[fish_clus==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
				}
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))
  		ordiellipse(fish.pca, groups=fish_clus, kind="sd", lwd=2, lty=2, col=pt_col)

	#Both these options (scaling by descriptor, clustering) rely exclusively on 'internal' variables 
		#(i.e., the descriptors that we used to run the ordination and/or clustering analysis in the first 
		#place). What about interpreting the structure of the descriptors in light of external variables?

	#Option 3: Color by Categorical Explanatory/Predictor Variable

	#Ordination plots can be enriched by relating them to external categorical factors, providing context and 
		#helping to interpret the underlying structure of the data. This technique allows us to map population 
		#or habitat data onto ordination results.

	#For example, for our class data, points can be colored by population.

		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))
		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[dat_final$Pop==groups[i],]
			points(dim_choice[,1], dim_choice[,2], pch=19, cex=1.4, col=pt_col[i])
			}
		ordiellipse(fish.pca, groups=dat_final$Pop, kind="sd", lwd=2, lty=2, col=pt_col)
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))
			legend(x="topright", legend=groups, col=pt_col, pch=19)

	#They can also be colored by habitat type (note that I've simplified the groups to avoid one-offs).

		NLCD_Cat2 <- dat_final$NLCD_Cat
			NLCD_Cat2[NLCD_Cat2 == "Wetlands-Woody"] = "Shrub-Scrub"
			NLCD_Cat2[NLCD_Cat2 == "Wetlands-Emergent"] = "Shrub-Scrub"
			NLCD_Cat2[NLCD_Cat2 == "Herbaceous-Grassland"] = "Shrub-Scrub"
			NLCD_Cat2[NLCD_Cat2 == "Pasture-Hay"] = "Shrub-Scrub"
		groups <- levels(factor(NLCD_Cat2))
		pt_col <- c("#3B528BFF","#5DC863FF")

		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[NLCD_Cat2==groups[i],]
			points(dim_choice[,1], dim_choice[,2], pch=19, cex=1.4, col=pt_col[i])
			}
		ordiellipse(fish.pca, groups=NLCD_Cat2, kind="sd", lwd=2, lty=2, col=pt_col)
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))
		legend(x="topright", legend=groups, col=pt_col, pch=19)

	#Option 4: Visualizing an Environmental Gradient

	#External numerical variables can also be visualized by scaling point sizes. For example, here we scale 
		#the points by elevation:

		El_scale <- 1 + (dat_final$Elev - min(dat_final$Elev))/120
			hist(El_scale)

		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))
		pt_col <- substr(pt_col,1,nchar(pt_col)-2)
		pt_col <- paste(pt_col, "80", sep = "")
		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		points(site.sc[,1], site.sc[,2], 
			pch=19, 
			cex=El_scale, 
			col="#50505050")
		arrows(0, 0, fish.pca$rotation[,1]*0.8, fish.pca$rotation[,2]*0.8, 
			lwd=2, 
			length=0.1)
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))

	#Now, let's step it up by using the `ordisurf` function to fit a smooth surface based on elevation. This 
		#will allow us to visualize the relationship between elevation and the PCA axes more clearly:

		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))
		pt_col <- substr(pt_col,1,nchar(pt_col)-2)
		pt_col <- paste(pt_col, "80", sep = "")
		p <- ordiplot(fish.pca, 
			type="n", 
			main="PCA of Hellinger-Transformed Fish Density Data",
			xlab="PC1 (50%)", 
			ylab="PC2 (29%)", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25,
			xlim=c(-0.75,1),
			ylim=c(-0.75,1))
		points(site.sc[,1], site.sc[,2], 
			pch=19, 
			cex=El_scale, 
			col="#50505050")
		arrows(0, 0, fish.pca$rotation[,1]*0.8, fish.pca$rotation[,2]*0.8, 
			lwd=2, 
			length=0.1)
		text(fish.pca$rotation[,1]*1, fish.pca$rotation[,2]*0.9, 
			row.names(fish.pca$rotation))
		ordisurf(fish.pca ~ Elev, 
			dat_final, 
			col="black", 
			add=TRUE, 
			select=FALSE)

	#The `ordisurf` function fits a smooth surface of elevation over the PCA plot, highlighting areas where 
		#elevation has a stronger influence on the ordination.

	#Note that for NMDS output, you can use the `MDSrotate` function to rotate the axes of an NMDS ordination 
		#plot to align them with an external environmental variable or gradient. This rotates the NMDS solution 
		#such that the variation in the external variable becomes more aligned with one of the ordination axes.

## Making Basic Inferences

	#We can further quantify this relationship by regressing the principal components against the external 
		#environmental factor (elevation). Let's see if elevation is significantly correlated with the first 
		#principal axis:

		summary(lm(site.sc[,1] ~ dat_final$Elev))

	#The summary shows there is a statistically significant relationship between elevation and the first principal 
		#component.

	#To visualize this relationship directly, let's plot elevation against the first principal axis and fit a 
		#regression line:

		plot(dat_final$Elev, site.sc[,1], 
	   		pch=20, 
	    	 	cex=2, 
	    	 	xlab="Elevation (m)", 
	    	 	ylab="PC1", 
	     		cex.axis=1.15, 
	     		cex.lab=1.25)
		abline(lm(site.sc[,1] ~ dat_final$Elev), 
		    	lty=2, 
		    	col="gray")

	#This plot provides a clear visual confirmation of the relationship between elevation and the first principal 
		#component, complementing the earlier `ordisurf` surface plot.


