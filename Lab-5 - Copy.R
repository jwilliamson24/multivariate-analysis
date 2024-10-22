###############################################################
###Lab-5.R									###
###Divisive and Non-hierarchical Cluster Analysis		###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 15, 2024						###
###############################################################


#i cant cluster my sites based on species assemblages; that's not what my data is meant to do,
#so im not going to try to make any of the species clustering methods work. see the env clustering below


## Lab Objectives

	#The goal of cluster analysis is to classify objects such as species, habitats, or environmental 
		#variables, into clusters based on their similarities or dissimilarities. This allows ecologists 
		#to identify natural groupings and patterns within ecological data. Divisive hierarchical clustering 
		#starts with all objects in the same cluster and iteratively divides them into different clusters. 
		#Non-hierarchical clustering identifies clusters formed at the threshold of similarity without taking 
		#into account the hierarchical cluster structure. In this lab, we will learn to:

		#-   Use a divisive hierarchical clustering approach to identify clusters in environmental and 
			#species data.
		#-   Apply K-means partitioning to ecological data and recognize the advantages and disadvantages of 
			#this approach relative to hierarchical approaches.
		#-   Use indicator species analysis to identify species that are strongly associated with pre-determined 
			#site groups.

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
      
      dat2 <- subset(dat2, year=="2024")
      
      sals <- dat2[,c("oss","enes")]
      
      drop <- c("lat","long","landowner","stand","tree_farm","year","weather","size_cl",
                "length_cl","site_id","oss","enes")
      env2 <- dat2[,!(colnames(dat2) %in% drop)]
  
      env_cont <- env2[, !colnames(env2) %in% "trt"]


  #3) Standardizing salamanders by sampling area
      
      sal_dens <- sals 
      for(i in 1:nrow(sals)){	
          sal_dens[i,] <- sals[i,]/567
        }
      
 
  #4) Transforming and standardizing the data as needed
  		
  		log_sal_cou <- log(sals + 1)
  		log_sal_dens <- log(sal_dens + 1)
  		
  		env_std <- decostand(env_cont, "standardize") #Z-scores the data in each column
  		

# ## Setting up the R Workspace and Preparing the Data
# 
# 	#We will load the necessary packages and prepare the dataset as we have done previously. This time I'm 
# 		#hiding most of the code because you should know what to do by now.
# 
# 	#First load the necessary packages and data:
# 
# 		library(vegan)
# 		library(dplyr)
# 		library(ggplot2)
# 		library(cluster)
# 		library(indicspecies)
# 
# 		#source("C:\\USGS_OCRU\\Teaching\\FW599_Multivariate_Statistics\\Data\\Biostats.R")
#     setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
#     source("Biostats.R")
#     source("coldiss.R")
# 
#   		dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
# 	
#   		sub_dat <- subset(dat, SMU=="Malheur")
# 
# 	#We will once again subset the data to include only the Malheur sites since the entire dataset is too large 
# 		#to effectively visualize.
# 
#   		spp_N <- colSums(sub_dat[,16:ncol(sub_dat)])
#   		spp_0 <- subset(spp_N, spp_N == 0)
#   		omit <- names(spp_0)
# 
#   		dat2 <- sub_dat[,!(colnames(sub_dat) %in% omit)]
# 	
#   		dat3 <- dat2[rowSums(dat2[,16:ncol(dat2)]) >0, ]
#   
#   		dat3$Herbaceous[is.na(dat3$Herbaceous)] <- 0 
#   		dat3$Ann_Herb[is.na(dat3$Ann_Herb)] <- 0
# 	
#   		dat3 <- dat3[complete.cases(dat3$SiteLength),]
#   
#   		dat_final <- dat3
# 	
#   		fish <- dat_final[,16:ncol(dat_final)]
#   		env <- dat_final[,1:15]
# 	
#   		fish_red <- drop.var(fish, min.fo=1)
#   
#   		fish_dens <- fish_red
#   		for(i in 1:nrow(fish_red)){	
#     			fish_dens[i,] <- fish_red[i,]/env$SiteLength[i]
#     			}
#   
#   		drop <- c("Latitude","Longitude","SiteLength","SiteWidth","SurfaceArea")
#   		env <- env[,!(colnames(env) %in% drop)]
#   		env_cont <- env[,!(colnames(env) %in% c("SMU","Pop","NLCD_Cat"))]
#   
#   		env <- env[,!(colnames(env) %in% c("Ave_Max_D","Ann_Herb"))]
#   		env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]
# 	
#   		log_fish_abu <- log(fish_red + 1)
#   		log_fish_dens <- log(fish_dens + 1)
# 
#   		env_std <- decostand(env_cont, "max")	
# 
# ## Divisive Hierarchical Clustering (Divisive Analysis, or DIANA)
# 
# 	#As for agglomerative hierarchical clustering, the first step is to produce an association matrix. 
# 		#We've already determined that, for the species abundance data, the Bray-Curtis matrix is a good choice.
# 
#   		fish.bray <- vegdist(log_fish_dens, "bray")
# 
	#The clustering algorithm is carried out in R using the `diana` function.

  		diana_sals <- diana(as.dist(sal.hel))

	#Note that DIANA can also be carried out using the log-transformed abundance data with the code:
		#`diana_fish <- diana(log_fish_dens)`. This, may be advisable for very large datasets, where a
		#distance matrix will be substantially larger than the data matrix, thus requiring more computing power.

	#The cophenetic correlation coefficient can be calculated for a divisive analysis the same way it's
		#calculated
		#for an agglomerative analysis. Remember that values range from 0 – 1, where a value closer to 1
			#indicates a solution of high quality. How does this solution compare to our agglomerative solution(s)?

  		#hclus.cophenetic(fish.bray, diana_fish)

	#We can also use an elbow plot to examine the optimal number of clusters in the same way as for agglomerative
		#clustering.

  		hclus.scree(diana_sals)
  		
	#In this case, it looks like there's a bit of a natural break around five clusters, which makes sense when
		#we visually assess the dendrogram.

  		plot(diana_sals, which=2, main="DIANA Dendrogram", 
  		     xlab="Sites", 
  		     ylab="Hellinger", 
  		     hang=-1)
  		rect.hclust(diana_sals, k=5)
  		
  		
  		
# 
# 		fishcl.class <- cutree(diana_fish, k=5)
# 
# 		fish_dat_new <- cbind(fishcl.class, log_fish_dens)
# 			par(mfrow=c(3,3))
# 			box.plots(fish_dat_new, by="fishcl.class")
# 
#   		plot(log_fish_dens$TROUT_RB, log_fish_dens$DACE_UNID,
# 			col = fishcl.class,
# 			pch = 19,
#     			main = "DIANA Clustering with 5 Clusters",
# 			xlab="Redband trout",
# 			ylab="Dace")

## K-Means Partitioning

	#K-means partitioning or clustering is a non-hierarchical method used to partition data into k clusters by
		#minimizing within-cluster sum of squares error.

	#The desired number of clusters must be set a-priori by the user. Based on our previous analyses, we can
		#guess that somewhere between 4-6 clusters is optimal for this dataset. If we're still unsure, we can
		#examine an elbow plot and silhouette width plot te get a better sense for the optimal number of clusters.
      # dev.off()
  		# nhclus.scree(fish.bray, max.k = 20)

	#The Calinski-Harabasz (CH) Index, or the ratio of between-cluster separation to within-cluster dispersion
		#normalized by their degrees of freedom, is another way to determine the optimal number of clusters. A
		#higher CH value indicates better clustering.

# 		fishcl.kmeans.cas <- cascadeKM(fish.bray, inf.gr=2, sup.gr=10, iter=100)
# 		plot(fishcl.kmeans.cas, sortg=T)
# 		fishcl.kmeans.cas$results
# 
# 	#It looks like CH is maximized at four clusters.
# 
# 	#As for DIANA, K-means clustering can be performed either on the distance matrices or the raw data.
# 		#For species abundances, it's best practice to use the distance matrix.
# 
# 		fish.kmeans <- kmeans(fish.bray, centers=4, iter.max=10000, nstart=10)
# 		fish.kmeans
# 
#   		plot(log_fish_dens$TROUT_RB, log_fish_dens$DACE_UNID,
# 			col = fish.kmeans$cluster,
# 			pch = 19,
#   			main = "K-Means Clustering with 4 Clusters",
# 			xlab="Redband trout",
# 			ylab="Dace")

	#Because we'll have to provide pre-determined site groupings for indicator species analysis below,
		#let's examine clusters in the environmental data as well.

  	
  	env_std_subset <- env_std[,c("temp","dwd_cov","soil_moist","stumps","logs","decay_cl","canopy_cov")]	
  		
  	env.euc <- vegdist(env_std_subset, metric="euclidean")

		envcl.kmeans.cas <- cascadeKM(env.euc, inf.gr=2, sup.gr=10, iter=100)
		plot(envcl.kmeans.cas, sortg=T)
		envcl.kmeans.cas$results

		env.kmeans <- kmeans(env.euc, centers=4, iter.max=10000, nstart=10)
		env.kmeans
	
		env.class <- env.kmeans$cluster
		env_dat_new<-cbind(env.class, env_std_subset)
		par(mfrow=c(2,3))
		box.plots(env_dat_new, by="env.class")

	
		#In general, the sites split well by habitat type (shrub-scrub vs. forest), elevation, and gradient. 
		#If you remember, K-Means output is random depending on starting centroids. If we want consistent 
		#groupings, we can use agglomerative clustering to select sites.

## Indicator Species Analysis

	#Indicator species analysis is a very fun way to identify species that are strongly associated with 
		#specific groups of sites or environmental conditions. It is best used when you have *a priori* 
		#information about sites and are interested in identifying species associated with different systems, 
		#strata, treatment groups, etc.

	#Let's try applying indicator species analysis on treatments

 		site_num <- env2$trt
  		ind_sp <- multipatt(log_sal_dens, site_num, func = "IndVal.g", control = how(nperm=999))
  		summary(ind_sp)
	#Tells us that enes is associated with BU and UU; interesting because these are the sites that they
  		#wouldnt want to move out of, least disturbed
  		
  		
  	site_num <- dat2$landowner
  		ind_sp <- multipatt(log_sal_dens, site_num, func = "IndVal.g", control = how(nperm=999))
  		summary(ind_sp)
  #Enes assoc with BLM and ODF
  	
  #Enes is associated with the less harvested sites
  		
  		
  		
  #What if we use the analysis on sites that have been clustered into groups based on their full suite of environmental characteristics?
	#This takes a little bit more work, but we already know how to cluster. Lucky us!

		envcl.ward <- hclust(env.euc, method = "ward.D")

		plot(envcl.ward, main="Ward Dendrogram", 
			xlab="Sites", 
			ylab="Euclidean Distance", 
			hang=-1)
			rect.hclust(envcl.ward, k=5)

		ENV_CLUS <- cutree(envcl.ward, k=5)

		par(mfrow=c(2,2))

		plot(as.factor(ENV_CLUS), env_cont$canopy_cov, 
			ylim=c(0,4), 
			col = rgb(red = 0, green = 0, blue = 0, alpha = 0.35),
			ylab = "canopy", 
			xlab="Cluster")
			points(ENV_CLUS, env_cont$canopy_cov, 
				pch=19, 
				cex=1.5, 
				col = rgb(red = 0, green = 0, blue = 1, alpha = 0.4))

		plot(as.factor(ENV_CLUS), env_cont$dwd_cov, 
			ylim=c(0,4), 
			col = rgb(red = 0, green = 0, blue = 0, alpha = 0.35),
			ylab = "dwdcov", 
			xlab="Cluster")
			points(ENV_CLUS, env_cont$dwd_cov, 
				pch=19, 
				cex=1.5, 
				col = rgb(red = 0, green = 0, blue = 1, alpha = 0.4))

		plot(as.factor(ENV_CLUS), env_cont$decay_cl, 
			ylim=c(0,4), 
			col = rgb(red = 0, green = 0, blue = 0, alpha = 0.35),
			ylab = "decay", 
			xlab="Cluster")
			points(ENV_CLUS, env_cont$elev, 
				pch=19, 
				cex=1.5, 
				col = rgb(red = 0, green = 0, blue = 1, alpha = 0.4))

		plot(as.factor(ENV_CLUS), env_cont$temp, 
			ylim=c(0,100), 
			col = rgb(red = 0, green = 0, blue = 0, alpha = 0.35),
			ylab = "temp", 
			xlab="Cluster")
			points(ENV_CLUS, env_cont$temp, 
				pch=19, 
				cex=1.5, 
				col = rgb(red = 0, green = 0, blue = 1, alpha = 0.4))

  		site_num <- ENV_CLUS
  		ind_sp <- multipatt(log_sal_dens, site_num, func = "IndVal.g", control = how(nperm=999))
			summary(ind_sp)

	#oss assoc with clusters 2, 3, 4
			#2: no canopy, high dwd, med decay
			#3: high canopy, med dwd,high decay
			#4: no canopy, med decay, med dwdcov


