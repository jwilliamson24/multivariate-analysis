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

#1) Loading packages and source code
    
    library(vegan)
    library(viridis)
    #setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
    setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/multivariate-analysis")
    source("Biostats.R")
    
#2) Loading and subsetting site level data
    
    dat2 <- readRDS("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
    #dat2 <- readRDS("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
    row.names(dat2) <- dat2[,1]
    
    dat2 <- subset(dat2, year=="2024")
    
    sals <- dat2[,c("oss","enes")]
    
    drop <- c("lat","long","landowner","stand","tree_farm","year","weather","size_cl",
              "length_cl","site_id","oss","enes")
    env2 <- dat2[,!(colnames(dat2) %in% drop)]
    
    env_cont <- env2[, !colnames(env2) %in% "trt"]
    
    drop <- c("jul_date","veg_cov","fwd_cov","dwd_count","size_cl","decay_cl","char_cl","length_cl" )
    env_subset <- env_cont[,!(colnames(env_cont) %in% drop)]
    
    
#3) Standardizing salamanders by sampling area
    
    sal_dens <- sals
    for(i in 1:nrow(sals)){
      sal_dens[i,] <- sals[i,]/567
    }
    
    
#4) Transforming and standardizing the data as needed
    
    log_sal_cou <- log(sals + 1)
    # log_sal_dens <- log(sal_dens + 1)
    
    env_std <- decostand(env_cont, "standardize") #Z-scores the data in each column


## Principal Component Analysis: Environmental Data

	#The PCA is easily run using the `prcomp` function.This function computes the principal components by 
		#performing an eigen decomposition on the covariance matrix (or correlation matrix if the data 
		#are standardized). Setting `scale=TRUE` standardizes the data to have a mean of zero and a 
		#standard deviation of one, ensuring that differences in scale do not unduly influence the analysis. 
		#Note that we *do not* have to use the `decostand` function when `scale=TRUE`!

  		env.pca <- prcomp(env_subset, scale=TRUE)	
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
  		groups <- levels(factor(dat2$trt))
  		pt_col <- viridis(length(groups))

  		p <- ordiplot(env.pca, type="n", main="PCA of Malheur Environmental Data")
    		for (i in 1:length(groups)){
      		dim_choice <- site.sc[dat2$trt==groups[i],]
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

		sal.pca <- prcomp(sals, scale=FALSE)	
		summary(sal.pca)

		screeplot(sal.pca, bstick=TRUE, main="Broken Stick of PCs")

  		pca.eigenvec(sal.pca, dim=6, digits=3, cutoff=0.1)

		site.sc <- scores(sal.pca)
		groups <- levels(factor(dat2$trt))
		pt_col <- viridis(length(groups))

		p <- ordiplot(fish.pca, type="n", main="PCA of Sal Data (Untransformed)")
			for (i in 1:length(groups)){
				dim_choice <- site.sc[dat2$trt==groups[i],]
				points(dim_choice[,1], dim_choice[,2], 
					pch=19, 
					cex=1.4, 
					col=pt_col[i])
					}
			text(site.sc, rownames(sals), pos=4, cex=0.7)
			abline(v=0, lty="dotted", col="gray")
			abline(h=0, lty="dotted", col="gray")
			arrows(0, 0, sal.pca$rotation[,1]*0.8, sal.pca$rotation[,2]*0.8, lwd=2, length=0.1)
			text(sal.pca$rotation[,1]*0.9, sal.pca$rotation[,2]*0.9, row.names(sal.pca$rotation))

	#How do we feel about this result? The first principal component, which is strongly influenced by 
		#redband trout abundance, really seems to be explaining most of the variance.

	#Now let's try it with the Hellinger transformation:

  		sal_hel <- decostand(sal_dens, method="hellinger")

		hel.pca <- prcomp(sal_hel, scale=FALSE)	
		summary(hel.pca)

		screeplot(hel.pca, bstick=TRUE, main="Broken Stick of PCs")

  		pca.eigenvec(hel.pca, dim=6, digits=3, cutoff=0.1)

		site.sc <- scores(hel.pca)
		groups <- levels(factor(dat2$trt))
		pt_col <- viridis(length(groups))

		p <- ordiplot(hel.pca, type="n", main="PCA of Sal Data (Hellinger)")
			for (i in 1:length(groups)){
				dim_choice <- site.sc[dat2$trt==groups[i],]
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



