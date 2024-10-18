###############################################################
###Lab-6.R									###
###Principal Component Analysis					###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 17, 2024						###
###############################################################

## Lab Objectives

	#The goal of principal component analysis (PCA) is to reduce the dimensionality of a dataset by 
		#transforming the original variables into a new set of uncorrelated variables called principal 
		#components. These components are ordered so that the first few retain most of the variation 
		#present in the original dataset. PCA helps to identify patterns, simplify data interpretation, 
		#and reduce noise by focusing on the most critical information. In this lab we will learn to:

		#-   Understand and apply PCA to environmental and species data.
		#-   Interpret PCA results.
		#-   Use the broken stick model to assess the significance of principal components.
		#-   Visualize PCA results using ordination biplots.
		#-   Compare the effects of PCA on raw and transformed species data.

## Setting up the R Workspace and Preparing the Data

	#We will load the necessary packages and prepare the dataset as we have done previously. You know what to do!

		library(vegan)
		library(viridis)
    setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
    dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
    source("Biostats.R")
    
  	#source("C:\\USGS_OCRU\\Teaching\\FW599_Multivariate_Statistics\\Data\\Biostats.R")
    #dat <- read.csv("C:\\USGS_OCRU\\Teaching\\FW599_Multivariate_Statistics\\Data\\Harney_Fishes_2007.csv", row.names = 1)

  	sub_dat <- subset(dat, SMU=="Malheur")

	#We will once again subset the data to include only the Malheur sites since the entire dataset is too 
		#large to effectively visualize.

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
    			fish_dens[i,] <- fish_red[i,]/env$SiteLength[i]
    			}
  
  		drop <- c("Latitude","Longitude","SiteLength","SiteWidth","SurfaceArea")
  		env <- env[,!(colnames(env) %in% drop)]
  		env_cont <- env[,!(colnames(env) %in% c("SMU","Pop","NLCD_Cat"))]
  
  		env <- env[,!(colnames(env) %in% c("Ave_Max_D","Ann_Herb"))]
  		env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]

## Principal Component Analysis: Environmental Data

	#The PCA is easily run using the `prcomp` function.This function computes the principal components by 
		#performing an eigen decomposition on the covariance matrix (or correlation matrix if the data 
		#are standardized). Setting `scale=TRUE` standardizes the data to have a mean of zero and a 
		#standard deviation of one, ensuring that differences in scale do not unduly influence the analysis. 
		#Note that we *do not* have to use the `decostand` function when `scale=TRUE`!

  		env.pca <- prcomp(env_cont, scale=TRUE)	
  		summary(env.pca)

	#The summary output shows the proportion of variance explained by each principal component, as well as 
		#the cumulative proportion. This helps to determine how many components are necessary to capture 
		#most of the variability in the data.

	#We can use the broken stick model to assess the importance of components. The broken stick model is a 
		#method for determining the significance of principal components by comparing the observed variance 
		#explained by each component to what would be expected by chance. Components explaining more variance 
		#than the broken stick model predicts are considered important.

  		screeplot(env.pca, bstick=TRUE, main="Broken Stick of PCs")

	#Note that just because the amount of variance explained by a principal component is less than what would 
		#be expected by random, doesn't mean the PC isn't important. Examine the total amount of variance 
		#explained by the first two PCs. If it's \>40 or 50%, you're probably on the right track. In this case, 
		#the total amount of variance explained by the first two PCs is about 63%, which is quite good.

	#We can examine the degree of influence of each environmental variable on the principal components using the 
		#`pca.eigenvec` function from the `Biostats` source code. It looks like herbaceous cover has stronger 
		#than expected loadings on PC2 and elevation, gradient, and canopy have strong loadings on PC1. This 
		#should be evident in the biplot.

  		pca.eigenvec(env.pca, dim=6, digits=3, cutoff=0.1)

	#The `biplot` function is the simplest way to look at the PCA output, but it's by far not the most elegant.

		biplot(env.pca, xlab="PC-1", ylab="PC-2")

	#The `ordiplot` function allows us to examine our PCA output with much more customization and detail. For example, 
		#we can scale site scores and descriptor vectors for easier visualization. We can also color points by 
		#"population" or other categorical variables. We will play around with this more next week.

  		site.sc <- scores(env.pca)
  		groups <- levels(factor(dat_final$Pop))
  		pt_col <- viridis(length(groups))

  		p <- ordiplot(env.pca, type="n", main="PCA of Malheur Environmental Data")
    		for (i in 1:length(groups)){
      		dim_choice <- site.sc[dat_final$Pop==groups[i],]
      		points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
			}
    			text(site.sc, rownames(env_cont), pos=4, cex=0.7)
    			abline(v=0, lty="dotted", col="gray")
    			abline(h=0, lty="dotted", col="gray")
    			arrows(0, 0, env.pca$rotation[,1]*3, env.pca$rotation[,2]*3, lwd=2, length=0.1)
    			text(env.pca$rotation[,1]*3.1, env.pca$rotation[,2]*3.1, row.names(env.pca$rotation))

###PCA of Hellinger-Transformed Species Data

	#Normally, using a PCA to analyze species abundance data results in a variety of issues that can muddy 
		#our interpretation of the data: the "double zero problem," heightened sensitivity to outliers, etc. 
		#An exception to this is transformation-based PCA, in which the data are transformed (usually with a 
		#Hellinger or Chord transformation) prior to analysis.

	#Let's see what the raw fish density PCA yields first. The procedure runs as above, except we are NOT 
		#standardizing our data using the `scale=TRUE` prompt.

		fish.pca <- prcomp(fish_dens, scale=FALSE)	
		summary(fish.pca)

		screeplot(fish.pca, bstick=TRUE, main="Broken Stick of PCs")

  		pca.eigenvec(fish.pca, dim=6, digits=3, cutoff=0.1)

		site.sc <- scores(fish.pca)
		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))

		p <- ordiplot(fish.pca, type="n", main="PCA of Malheur Fish Data (Untransformed)")
			for (i in 1:length(groups)){
				dim_choice <- site.sc[dat_final$Pop==groups[i],]
				points(dim_choice[,1], dim_choice[,2], 
					pch=19, 
					cex=1.4, 
					col=pt_col[i])
					}
			text(site.sc, rownames(fish), pos=4, cex=0.7)
			abline(v=0, lty="dotted", col="gray")
			abline(h=0, lty="dotted", col="gray")
			arrows(0, 0, fish.pca$rotation[,1]*0.8, fish.pca$rotation[,2]*0.8, lwd=2, length=0.1)
			text(fish.pca$rotation[,1]*0.9, fish.pca$rotation[,2]*0.9, row.names(fish.pca$rotation))

	#How do we feel about this result? The first principal component, which is strongly influenced by 
		#redband trout abundance, really seems to be explaining most of the variance.

	#Now let's try it with the Hellinger transformation:

  		fish_hel <- decostand(fish_dens, method="hellinger")

		hel.pca <- prcomp(fish_hel, scale=FALSE)	
		summary(hel.pca)

		screeplot(hel.pca, bstick=TRUE, main="Broken Stick of PCs")

  		pca.eigenvec(hel.pca, dim=6, digits=3, cutoff=0.1)

		site.sc <- scores(hel.pca)
		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))

		p <- ordiplot(hel.pca, type="n", main="PCA of Malheur Fish Data (Hellinger)")
			for (i in 1:length(groups)){
				dim_choice <- site.sc[dat_final$Pop==groups[i],]
				points(dim_choice[,1], dim_choice[,2], 
					pch=19, 
					cex=1.4, 
					col=pt_col[i])
					}
			text(site.sc, rownames(fish), pos=4, cex=0.7)
			abline(v=0, lty="dotted", col="gray")
			abline(h=0, lty="dotted", col="gray")
			arrows(0, 0, hel.pca$rotation[,1]*0.8, hel.pca$rotation[,2]*0.8, lwd=2, length=0.1)
			text(hel.pca$rotation[,1]*0.9, hel.pca$rotation[,2]*0.9, row.names(hel.pca$rotation))

	#We can see that the variance appears to be partitioned a little tidier this way; however, what phenomenon 
		#do you see starting to show through? Hint: you read about it in week 2!



