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

	#2) Loading the data

  		dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
    
  		sub_dat <- subset(dat, SMU=="Malheur")

	#In this case, we will subset the data to include only the Malheur sites since the entire dataset is too 
		#large to effectively visualize.

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

## Q-Mode: (Dis)similarity and Distance Matrices

### (Dis)similarity Coefficients for Species Presence/Absence (Binary) Data

	#Analyses can be performed on binary (0-1) data, including for species presence/absence data when binary 
		#values are the only data available or when abundances are irrelevant.

	#Since our class dataset is quantitative, we'll first convert it to a presence/absence format using the 
		#following transformation:

  		fish_occ <- data.trans(fish_red, method="power", exp=0, plot=F)

	#The simple matching coefficient is a symmetric measure of similarity that is computed as the number of 
		#double 0s and double 1s among a pair of sites divided by the total number of descriptors. Because it 
		#is a symmetric measure (i.e., double 0s are counted as perfectly similar), it is inappropriate for 
		#species data.

		fish.sim <- sm(fish_occ)

	#In contrast, Jaccard's similarity is an asymmetric measure of similarity that is computed as the number 
		#of double 1's divided by the number of descriptors excluding double 0s. It is generally robust for 
		#use with species presence/absence data. Note that `vegan` calculates these coefficients as a 
		#dissimilarity coefficients.

  		fish.jac <- vegdist(fish_occ, method="jaccard") 

	#By including double 0s, the simple matching coefficient tends to underestimate the degree of dissimilarity 
		#between two sites when compared to Jaccard's coefficient.

		plot(fish.jac, fish.sim, 
			xlab="Jaccard's coefficient", 
			ylab="Simple Matching coefficient",
			pch=21,
			col="black")
			abline(0, 1, col="darkgray")
	
	#Sorensen's similarity is calculated similarly to Jaccard's coefficient, but gives double weight to the number 
		#of double 1s. As such, it trends toward greater similarity values. Note that it is equivalent to the 
		#percentage difference or Bray Curtis method for quantitative (or abundance) data.

		fish.sor <- vegdist(fish_occ, method="bray") 

		plot(fish.jac, fish.sor, 
			xlab="Jaccard's coefficient", 
			ylab="Sorensen's coefficient",
			pch=21,
			col="black")
			abline(0, 1, col="darkgray")

	#The `coldiss` function is helpful for visualizing ordered similarities among sites (or objects), and is a 
		#precursor for identifying clusters.

  		coldiss(fish.jac, nc=5, byrank=FALSE, diag=TRUE)

	#Another helpful coefficient for presence/absence data is **Ochiai's similarity**, which is an asymmetric 
		#coefficient that is closely related to the chord, Hellinger, and log-chord transformations.

### (Dis)similarity Coefficients for Species Abundance Data

	#Species abundances generally require analysis using asymmetrical dissimilarity metrics. The four most 
		#frequently used asymmetrical coefficients are the percentage difference (i.e., Bray-Curtis), chord, 
		#Hellinger, and chi-square distances. Each one has its own advantages and disadvantages.

	#The percentage difference or Bray-Curtis dissimilarity metric is usually computed from log-transformed data, 
		#although it can be computed from raw data. This metric gives the same importance to absolute differences 
		#in abundance, irrespective of their order of magnitude. That is, it weighs abundant and rare species similarly.

		fish.bray <- vegdist(log_fish_dens, method="bray")

	#The chord distance is a Euclidean distance computed on site (object) vectors that have been normalized to 
		#length 1 (also referred to as the chord transformation). The log-chord distance is simply the chord 
		#distance applied to log-transformed abundance data.

		fish.cho <- decostand(fish_dens, "normalize") 
		fish.cho <- dist(fish.cho)

	#Similarly, the Hellinger distance is a Euclidean distance between site (object) vectors where abundance 
		#values are first divided by the site total abundance and the result is square-root transformed 
		#(also referred to as the Hellinger transformation).

		fish.hel <- decostand(fish_dens, "hellinger") 
		fish.hel <- dist(fish.hel)

	#Note that both the chord and Hellinger transformations can be used prior to ordination using a Euclidean 
		#distance matrix (as in PCA or RDA) in a procedure referred to as transformation-based PCA, RDA, or 
		#K-means clustering. This has a similar outcome as using a chord or Hellinger transformation prior to 
		#use in PCoA or NMDS.

	#The chord and Hellinger distances are closely related to one-another.

		plot(fish.cho, fish.hel, 
	  		xlab="Chord distance", 
	  		ylab="Hellinger distance",
	  		pch=21,
	  		col="black")
	  		abline(0, 1, col="darkgray")

	#Bray-Curtis tends to demonstrate greater dissimilarity for some sites than the Hellinger distance. 
		#Why do you think this is?

		plot(fish.bray, fish.hel, 
	  		xlab="Bray-Curtis distance", 
	  		ylab="Hellinger distance",
	  		pch=21,
	  		col="black")
	  		abline(0, 1, col="darkgray")
	  
  		coldiss(fish.bray, nc=5, byrank=FALSE, diag=TRUE)
  		coldiss(fish.hel, nc=5, byrank=FALSE, diag=TRUE)

	#The chi-square distance is weighed by the relative frequency of occurrence. This means that it gives 
		#higher weights to rare than to common species. It is only recommended when rare species are considered 
		#to be good indicators of special ecological conditions.

		fish.chi <- vegdist(log_fish_abu, method="chi")

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

  		fish_pa.t <- t(fish_occ)
	
  		fish.t.jac <- vegdist(fish_pa.t, "jaccard")

  		coldiss(fish.t.jac, diag=TRUE)

	#Parametric and non-parametric correlation coefficients are used to compare species distributions and 
		#associations in space and time. This includes the parametric Pearson's correlation coefficient and 
		#the non-parametric (i.e., rank-based) Spearman's r or Kendall's tau. The chi-square distance 
		#can also be used in R-mode for species abundance data.

  		fish_dens.t <- t(log_fish_dens)

  		fish.pears <- cor(fish_dens)
  
  		fish.ken <- cor(fish_dens, method="kendall")
  
  		fish.t.chi <- decostand(fish_dens.t, "chi.square")
  		fish.t.chi2 <- dist(fish.t.chi) #Chi-square

  		coldiss(fish.pears, diag=TRUE)
  		coldiss(fish.ken, diag=TRUE)
  		coldiss(fish.t.chi2, diag=TRUE)

### Environmental Data

	#Quantitative variables such as environmental data must be dimensionally homogeneous (i.e., standardized) 
		#prior to R-mode analysis. For variables with largely linear relationships, Pearson's correlation coefficient 
		#is sufficient.

		env.pearson <- cor(env_std)
		env.o <- order.single(env.pearson)

		pairs(env_std[, env.o],
			lower.panel = panel.smooth,
			upper.panel = panel.cor,
			diag.panel = panel.hist,
			main = "Pearson Correlation Matrix")

	#A multi-panel display of pairwise relationships among environmental variables can be a great addition to our 
		#exploratory analysis in Lab #1.

	#For quantitative variables with non-linear relationships, the non-parametric rank correlation coefficients 
		#such as Spearman's r or Kendall's tau may perform better.

		env.ken <- cor(env_std, method="kendall")
		env.o <- order.single(env.ken)

		pairs(env_std[, env.o],
			lower.panel = panel.smooth,
			upper.panel = panel.cor, 
			method = "kendall",
			diag.panel = panel.hist,
			main = "Kendall Correlation Matrix")

