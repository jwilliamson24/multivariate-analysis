###############################################################
###Lab-4.R									###
###Agglomerative Hierarchical Clustering				###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 10, 2024						###
###############################################################


#i cant cluster my sites based on species assemblages; that's not what my data is meant to do
#i got a few to run just for the practice but theyre obviously quite ugly and not helpful


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

  #1) Loading packages and source code
  
      library(vegan)
      library(nomclust)
      library(BioStatR)
      library(RColorBrewer)
      library(gclus)
      
      setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
      source("Biostats.R")
      source("coldiss.R")
      
  
  #2) Loading and subsetting site level data
      
      dat2 <- readRDS("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
      row.names(dat2) <- dat2[,1]
      
      sals <- dat2[,c("oss","enes")]
      
      drop <- c("lat","long","landowner","stand","tree_farm","year","weather","size_cl",
                "length_cl","site_id","oss","enes")
      env <- dat2[,!(colnames(dat2) %in% drop)]
  
      env_cont <- env[, !colnames(env) %in% "trt"]


  #3) Standardizing salamanders by sampling area
      
      sal_dens <- sals 
      for(i in 1:nrow(sals)){	
          sal_dens[i,] <- sals[i,]/567
        }
      
 
  #4) Transforming and standardizing the data as needed
  		
  		log_sal_cou <- log(sals + 1)
  		log_sal_dens <- log(sal_dens + 1)
  		
  		env_std <- decostand(env_cont, "standardize") #Z-scores the data in each column
  		

## Agglomerative Hierarchical Clustering Methods

	#The first step prior to conducting a cluster analysis is to produce an association matrix. We've already 
		#determined that, for the species abundance data, the Bray-Curtis matrix is a good choice.

  		sals.bra <- vegdist(sals,"bray")
  		
  		sal.hel <- decostand(sals, "hellinger") 
  		
### Single Linkage

	#Single linkage clustering merges clusters based on the *minimum* distance between points. Although it is a 
		#simple and intuitive method, it often results in a "chaining effect," which can be observed in the results below:

  		fishcl.sin <- hclust(fish.bra, method = "single")

  		plot(fishcl.sin, main="Single Linkage Dendrogram",
			xlab="Sites", 
			ylab="Bray-Curtis Dissimilarity", 
			hang=-1)
  		
  		
  		sals.sin <- hclust(sal.hel, method = "single")
  		
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
  		
  		
  		salcl.com <- hclust(sal.hel, method = "complete")
  		
  		plot(salcl.com, main="Complete Linkage Dendrogram",
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

  		
  		salcl.ave <- hclust(sal.hel, method = "average")
  		
  		plot(salcl.ave, main="Average Linkage Dendrogram (UPGMA)",
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
  		
  		
  		salcl.wpgma <- hclust(sal.hel, method = "mcquitty")
  		
  		plot(salcl.wpgma, main="Weighted Average Linkage Dendrogram (WPGMA)",
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

