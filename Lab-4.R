###############################################################
###Lab-4.R									###
###Agglomerative Hierarchical Clustering				###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 10, 2024						###
###############################################################


## Lab Objectives

	#The goal of hierarchical cluster analysis is to classify objects such as species, habitats, or environmental 
		#variables, into clusters based on their similarities or dissimilarities. This allows ecologists to 
		#identify natural groupings and patterns within ecological data. Agglomerative hierarchical clustering 
		#starts with each object in its own cluster and iteratively merges clusters. 

	#In this lab, we will learn to:

		#Use an agglomerative hierarchical clustering approach to identify clusters in environmental and species data.
		#Apply different clustering methods to ecological data and identify the pros and cons of each one.
		#Use non-statistical and statistical methods to assess clustering results and determine the appropriate 
			#number of clusters.

## Setting up the R Workspace and Preparing the Data

	#We will load the necessary packages and prepare the dataset as we have done previously. This should be second 
		#nature by now!

	#1) Loading packages and source code

		library(vegan)
		library(cluster)

    setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")

	#2) Loading the data

      dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
      source("Biostats.R")
      
  		sub_dat <- subset(dat, SMU=="Malheur")

	#We will once again subset the data to include only the Malheur sites since the entire dataset is too large to 
		#effectively visualize.

	#3) Omitting species with zero observations and sites without fish

		spp_N <- colSums(sub_dat[,16:ncol(sub_dat)])
		spp_0 <- subset(spp_N, spp_N == 0)
		omit <- names(spp_0)

		dat2 <- sub_dat[,!(colnames(sub_dat) %in% omit)]
	
  		dat3 <- dat2[rowSums(dat2[,16:ncol(dat2)]) >0, ]

	#4) Dealing with missing data

		dat3$Herbaceous[is.na(dat3$Herbaceous)] <- 0 
		dat3$Ann_Herb[is.na(dat3$Ann_Herb)] <- 0
	
		dat3 <- dat3[complete.cases(dat3$SiteLength),]
		dat_final <- dat3

	#5) Splitting the data set into environmental variables and species abundances

		fish <- dat_final[,16:ncol(dat_final)]
		env <- dat_final[,1:15]

	#6) Dropping rare species

  		fish_red <- drop.var(fish, min.fo=1)

	#7) Standardizing species observed abundance by sampling effort

  		fish_dens <- fish_red
  		for(i in 1:nrow(fish_red)){	
    			fish_dens[i,] <- fish_red[i,]/env$SiteLength[i]
    			}

	#8) Selecting relevant environmental data without covarying factors

		drop <- c("Latitude","Longitude","SiteLength","SiteWidth","SurfaceArea")
		env <- env[,!(colnames(env) %in% drop)]
		env_cont <- env[,!(colnames(env) %in% c("SMU","Pop","NLCD_Cat"))]
  
		env <- env[,!(colnames(env) %in% c("Ave_Max_D","Ann_Herb"))]
		env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]

	#9) Checking for outliers (remember, we are ignoring them for now)

	#10) Transforming and standardizing the data as needed

  		log_fish_abu <- log(fish_red + 1)
  		log_fish_dens <- log(fish_dens + 1)

  		env_std <- decostand(env_cont, "max")

## Agglomerative Hierarchical Clustering Methods

	#The first step prior to conducting a cluster analysis is to produce an association matrix. We've already 
		#determined that, for the species abundance data, the Bray-Curtis matrix is a good choice.

  		fish.bra <- vegdist(log_fish_dens, "bray")
  		
  		
  
  		sals.bra <- vegdist(sals_red,"bray")

### Single Linkage

	#Single linkage clustering merges clusters based on the *minimum* distance between points. Although it is a 
		#simple and intuitive method, it often results in a "chaining effect," which can be observed in the results below:

  		fishcl.sin <- hclust(fish.bra, method = "single")

  		plot(fishcl.sin, main="Single Linkage Dendrogram",
			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)
  		
  		
  		sals.sin <- hclust(sals.bra, method = "single")
  		
  		plot(sals.sin, main="Single Linkage Dendrogram",
  		     xlab="Sites", 
  		     ylab="Bray-Curtis Dissimilarity", 
  		     hang=-1)

### Complete Linkage

	#Complete linkage clustering merges clusters based on the *maximum* distance between points. It is less 
		#prone to chaining and produces compact, spherical clusters.

  		fishcl.com <- hclust(fish.bra, method = "complete")

  		plot(fishcl.com, main="Complete Linkage Dendrogram",
			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)

### Average Linkage

	#Average linkage clustering (UPGMA) is a useful intermediate that merges clusters based on the average 
		#distance between points. The average value is weighted based on each cluster's size.

  		fishcl.ave <- hclust(fish.bra, method = "average")

  		plot(fishcl.ave, main="Average Linkage Dendrogram (UPGMA)",
    			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)

### Weighted Average Linkage

	#Weighted average linkage clustering (WPGMA) is similar to UPGMA except averages are weighted equally 
		#regardless of a cluster's size. This can be an advantage or disadvantage depending on whether it is 
		#appropriate to treat clusters equally.

  		fishcl.wpgma <- hclust(fish.bra, method = "mcquitty")

  		plot(fishcl.wpgma, main="Weighted Average Linkage Dendrogram (WPGMA)",
    			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)

### Centroid Linkage

	#Centroid linkage cluster analysis (UPGMC) merges clusters based on the distance between their centroids. 
		#It produces clusters based on the data's geometric center, which can be more representative of the 
		#cluster's overall location. This also means this method is less sensitive to outliers. The downside 
		#is that basing distances on the centroid can lead to reversals.

  		fishcl.cen <- hclust(fish.bra, method = "centroid")

  		plot(fishcl.cen, main="Centroid Linkage Dendrogram (UPGMC)",
    			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)

### Weighted Centroid Linkage

	#Weighted centroid linkage clustering (WPGMC) is similar to UPGMC except averages are weighted equally 
		#regardless of a cluster's size.

  		fishcl.wpgmc <- hclust(fish.bra, method = "median")

  		plot(fishcl.wpgmc, main="Weighted Centroid Linkage Dendrogram (WPGMC)",
    			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)

### Ward's Minimum Variance

	#Ward's minimum variance method merges clusters to achieve the smallest possible increase in the sum of 
		#squared distances within each cluster. It often results in balanced and easily interpretable clusters.

  		fishcl.ward <- hclust(fish.bra, method = "ward.D2")

  		plot(fishcl.ward, main="Ward's Minimum Variance Dendrogram",
    			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)

## Assessing Cluster Fit

### Agglomerative Coefficient

	#The agglomerative coefficient measures the strength of the clustering structure. It ranges from 0 to 1, 
		#where a value close to 1 indicates distinct, well-separated clusters, and a value close to 0 suggests 
		#a weak clustering structure.

  		coef.hclust(fishcl.sin)
  		coef.hclust(fishcl.com)
  		coef.hclust(fishcl.ave)
  		coef.hclust(fishcl.wpgma)
  		coef.hclust(fishcl.ward)

	#Note that the agglomerative coefficient produces an error for centroid linkage methods because there are reversals.

	#In this case Ward's method produces the highest agglomerative coefficient, so we can deduce that the 
		#resulting dendrogram is well-structured.

### Cophenetic Correlation Coefficient

	#The cophenetic correlation coefficient measures the correlation between the cophenetic distances 
		#obtained from the dendrogram and the original distances in the distance matrix. Values range from 0 – 1, 
		#where a value closer to 1 indicates a solution of high quality.

		hclus.cophenetic(fish.bra, fishcl.sin)
		hclus.cophenetic(fish.bra, fishcl.com)
		hclus.cophenetic(fish.bra, fishcl.ave)
		hclus.cophenetic(fish.bra, fishcl.ward)

	#Interestingly, Ward's method does not produce the best fit here. We should look at the dendrogram and 
		#determine for ourselves which makes the most sense. For this data set, Ward's method is perfectly adequate.

### How Many Clusters?

	#Commonly in ecology we want to know how many clusters reasonably exist in a dataset? This can be 
		#determined "non-statstically" by examining the dendrogram, or through various statistical methods.

	#The easiest is the "elbow method," by which the cluster sums of squares is plotted against the number of 
		#clusters and a natural break point is identified.

  		hclus.scree(fishcl.ward)

		plot(fishcl.ward, main="Ward's Minimum Variance Dendrogram", 
			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)
			rect.hclust(fishcl.ward, k=6)

		fishcl.class <- cutree(fishcl.ward, k=6)

		fish_dat_new <- cbind(fishcl.class, log_fish_dens)
		
		box.plots(fish_dat_new, by="fishcl.class")

