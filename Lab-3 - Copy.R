###############################################################
###Lab-3.R									###
###Association Matrices							###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 8, 2024						###
###############################################################


## Lab Objectives

	#Association matrices form the basis for most ordination and clustering techniques. Thus, selection of an 
		#appropriate measure of association is critical for the effective interpretation of multivariate ecological 
		#data, and using the wrong measure of association can have disastrous consequences. 

	#In this lab, we will learn to:

	#Compute, examine, and compare (dis)similarity coefficients (Q-mode analysis) for:
		#Species presence/absence data
		#Species abundance data
		#Multi-scale environmental data
	#Compute dependence matrices (R-mode analysis) among:
		#Species
		#Environmental varibales

## Setting up the R Workspace and Preparing the Data

	#We will load the necessary packages and prepare the dataset as we have done previously, including:

	#1) Loading packages and source code

		library(vegan)
		library(nomclust)
		library(BioStatR)
		library(RColorBrewer)
    library(gclus)

    setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
    source("Biostats.R")
    source("coldiss.R")
    
		# source("C:\\USGS_OCRU\\Teaching\\FW599_Multivariate_Statistics\\Data\\Biostats.R")
		# source("C:\\USGS_OCRU\\Teaching\\FW599_Multivariate_Statistics\\Data\\coldiss.R")

	#2) Loading site level data
    dat2 <- readRDS("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
    row.names(dat2) <- dat2[,1]
    
    sals <- dat2[,c("oss","enes")]
    
    drop <- c("lat","long","landowner","stand","tree_farm","year","weather","size_cl",
              "length_cl","site_id","oss","enes")
    env <- dat2[,!(colnames(dat2) %in% drop)]
    
    #env_cont <- env[,c("elev","temp","hum","soil_moist","canopy_cov","veg_cov",
    #                   "dwd_cov","fwd_cov","jul_date","dwd_dens","decay_cl","char_cl")]
    env_cont <- env[, !colnames(env) %in% "trt"]
    

	#10) Transforming and standardizing the data as needed
    
    
    log_sals <- log(sals + 1)
    env_std <- decostand(env_cont, "standardize") #Z-scores the data in each column
      

## Q-Mode: (Dis)similarity and Distance Matrices

### (Dis)similarity Coefficients for Species Presence/Absence (Binary) Data

	#Analyses can be performed on binary (0-1) data, including for species presence/absence data when binary 
		#values are the only data available or when abundances are irrelevant.

	#Since our class dataset is quantitative, we'll first convert it to a presence/absence format using the 
		#following transformation:

  		sal_occ <- data.trans(sals, method="power", exp=0, plot=F)

	#The simple matching coefficient is a symmetric measure of similarity that is computed as the number of 
		#double 0s and double 1s among a pair of sites divided by the total number of descriptors. Because it 
		#is a symmetric measure (i.e., double 0s are counted as perfectly similar), it is inappropriate for 
		#species data.

		sal.sim <- sm(sal_occ)

	#In contrast, Jaccard's similarity is an asymmetric measure of similarity that is computed as the number 
		#of double 1's divided by the number of descriptors excluding double 0s. It is generally robust for 
		#use with species presence/absence data. Note that `vegan` calculates these coefficients as a 
		#dissimilarity coefficients.

  		sal.jac <- vegdist(sal_occ, method="jaccard") 

	#By including double 0s, the simple matching coefficient tends to underestimate the degree of dissimilarity 
		#between two sites when compared to Jaccard's coefficient.

		plot(sal.jac, sal.sim, 
			xlab="Jaccard's coefficient", 
			ylab="Simple Matching coefficient",
			pch=21,
			col="black")
			abline(0, 1, col="darkgray")
	
	#Sorensen's similarity is calculated similarly to Jaccard's coefficient, but gives double weight to the number 
		#of double 1s. As such, it trends toward greater similarity values. Note that it is equivalent to the 
		#percentage difference or Bray Curtis method for quantitative (or abundance) data.

		sal.sor <- vegdist(sal_occ, method="bray") 

		plot(sal.jac, sal.sor, 
			xlab="Jaccard's coefficient", 
			ylab="Sorensen's coefficient",
			pch=21,
			col="black")
			abline(0, 1, col="darkgray")

	#The `coldiss` function is helpful for visualizing ordered similarities among sites (or objects), and is a 
		#precursor for identifying clusters.

  		coldiss(sal.jac, nc=5, byrank=FALSE, diag=TRUE)

	#Another helpful coefficient for presence/absence data is **Ochiai's similarity**, which is an asymmetric 
		#coefficient that is closely related to the chord, Hellinger, and log-chord transformations.

### (Dis)similarity Coefficients for Species Abundance Data

	#Species abundances generally require analysis using asymmetrical dissimilarity metrics. The four most 
		#frequently used asymmetrical coefficients are the percentage difference (i.e., Bray-Curtis), chord, 
		#Hellinger, and chi-square distances. Each one has its own advantages and disadvantages.

	#The percentage difference or Bray-Curtis dissimilarity metric is usually computed from log-transformed data, 
		#although it can be computed from raw data. This metric gives the same importance to absolute differences 
		#in abundance, irrespective of their order of magnitude. That is, it weighs abundant and rare species similarly.

		sal.bray <- vegdist(log_sals, method="bray")

	#The chord distance is a Euclidean distance computed on site (object) vectors that have been normalized to 
		#length 1 (also referred to as the chord transformation). The log-chord distance is simply the chord 
		#distance applied to log-transformed abundance data.

		sal.cho <- decostand(sals, "normalize") 
		sal.cho <- dist(sal.cho)

	#Similarly, the Hellinger distance is a Euclidean distance between site (object) vectors where abundance 
		#values are first divided by the site total abundance and the result is square-root transformed 
		#(also referred to as the Hellinger transformation).

		sal.hel <- decostand(sals, "hellinger") 
		sal.hel <- dist(sal.hel)

	#Note that both the chord and Hellinger transformations can be used prior to ordination using a Euclidean 
		#distance matrix (as in PCA or RDA) in a procedure referred to as transformation-based PCA, RDA, or 
		#K-means clustering. This has a similar outcome as using a chord or Hellinger transformation prior to 
		#use in PCoA or NMDS.

	#The chord and Hellinger distances are closely related to one-another.

		plot(sal.cho, sal.hel, 
	  		xlab="Chord distance", 
	  		ylab="Hellinger distance",
	  		pch=21,
	  		col="black")
	  		abline(0, 1, col="darkgray")

	#Bray-Curtis tends to demonstrate greater dissimilarity for some sites than the Hellinger distance. 
		#Why do you think this is?

		plot(sal.bray, sal.hel, 
	  		xlab="Bray-Curtis distance", 
	  		ylab="Hellinger distance",
	  		pch=21,
	  		col="black")
	  		abline(0, 1, col="darkgray")
	  
  		coldiss(sal.bray, nc=5, byrank=FALSE, diag=TRUE)
  		coldiss(sal.hel, nc=5, byrank=FALSE, diag=TRUE)

	#The chi-square distance is weighed by the relative frequency of occurrence. This means that it gives 
		#higher weights to rare than to common species. It is only recommended when rare species are considered 
		#to be good indicators of special ecological conditions.

	  sal.chi <- vegdist(log_sals, method="chi")

### Similarity Coefficients for Mixed Data Types and Environmental Data

	#For variables with a clear interpretation of double zeros, the symmetrical Euclidean distance matrix is the 
		#most straightforward coefficient. It has no upper limit and is thus strongly influenced by the scale 
		#of each descriptor. Therefore, it is recommended that all environmental data be standardized prior to analysis.

		env.euc <- vegdist(env_std, metric="euclidean")

	#Gower's coefficient is useful when the dataset contains categorical (or nominal) data. Gower's similarity has 
		#been designed to handle data of varying mathematical types such that each variable is treated separately. 
		#This is also a symmetrical coefficient with a structure similar to that of the simple matching coefficient 
		#when variables are binary or categorical.

		env.gower <- daisy(env_cont, metric="gower")

  		plot(env.gower, env.euc, 
			xlab="Gower dissimilarity", 
			ylab="Euclidean distance of standardized data",
			pch=21,
			col="black")
			abline(0, 1, col="darkgray")

## R-Mode: Identifying Associations Among Species and Environmental Data

### Species Associations

	#For binary species data, the Jaccard and Sorensen coefficients can also be used in R-mode to discern 
		#species associations.

  		sal_pa.t <- t(sal_occ)
	
  		sal.t.jac <- vegdist(sal_pa.t, "jaccard")

  		coldiss(sal.t.jac, diag=TRUE)

	#Parametric and non-parametric correlation coefficients are used to compare species distributions and 
		#associations in space and time. This includes the parametric Pearson's correlation coefficient and 
		#the non-parametric (i.e., rank-based) Spearman's r or Kendall's tau. The chi-square distance 
		#can also be used in R-mode for species abundance data.

  		sals.t <- t(log_sals)

  		sal.pears <- cor(sals)
  
  		sal.ken <- cor(sals, method="kendall")
  
  		sal.t.chi <- decostand(sals.t, "chi.square")
  		sal.t.chi2 <- dist(sal.t.chi) #Chi-square

  		coldiss(sal.pears, diag=TRUE)
  		coldiss(sal.ken, diag=TRUE)
  		coldiss(sal.t.chi2, diag=TRUE)

### Environmental Data

	#Quantitative variables such as environmental data must be dimensionally homogeneous (i.e., standardized) 
		#prior to R-mode analysis. For variables with largely linear relationships, Pearson's correlation coefficient 
		#is sufficient.
  	drop <- c("jul_date","canopy_cov","veg_cov","fwd_cov","stumps","logs","decay_cl","char_cl" )
  	env_std_subset <- env_std[,!(colnames(env_std) %in% drop)]
  		
		env.pearson <- cor(env_std_subset)
		env.o <- order.single(env.pearson)

		pairs(env_std_subset[, env.o],
			lower.panel = panel.smooth,
			upper.panel = panel.cor,
			diag.panel = panel.hist,
			main = "Pearson Correlation Matrix")

	#A multi-panel display of pairwise relationships among environmental variables can be a great addition to our 
		#exploratory analysis in Lab #1.

	#For quantitative variables with non-linear relationships, the non-parametric rank correlation coefficients 
		#such as Spearman's r or Kendall's tau may perform better.

		env.ken <- cor(env_std_subset, method="kendall")
		env.o <- order.single(env.ken)

		pairs(env_std_subset[, env.o],
			lower.panel = panel.smooth,
			upper.panel = panel.cor, 
			method = "kendall",
			diag.panel = panel.hist,
			main = "Kendall Correlation Matrix")

