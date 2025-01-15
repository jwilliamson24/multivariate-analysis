###############################################################
###Lab-2.R									###
###Data Transformation and Standardization			###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 3, 2024						###
###############################################################

### Lab Objectives

	#Many frequentist or *least squares* statistical methods are dependent upon data meeting the 
		#assumptions of normality and homoscedasticity. The same is true for multivariate methods 
		#that rely on linear relationships among variables (e.g., Euclidean distances). In ecology, 
		#variables are rarely symmetrically distributed around their means. Furthermore, environmental 
		#variables may be multi-scalar, requiring standardization to facilitate robust comparisons. 
		
	#In this lab, we will learn to:

	#-   Examine the distribution and structure of the data
	#-   Evaluate the need for data transformation and standardization

### Setting up the R Workspace and Preparing the Data

	#First, we will load the necessary packages and prepare the dataset as in Lab #1.

	#Loading packages and source code:

    		library(vegan)
    		library(pastecs)
    		library(corrplot)
    		library(ggplot2)
    		library(ggpubr)

	#Loading data:

		#setwd("U:\\FW599_Multivariate\\Data\\")
		dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
		source("Biostats.R")

	#Remember, rows are named by sampling site and row names cannot be duplicates! Each site or "object" 
		#should have a unique ID.

	#Omitting species with zero observations:

		spp_N <- colSums(dat[,16:48])
		spp_0 <- subset(spp_N, spp_N == 0)
		omit <- names(spp_0)

		dat2 <- dat[,!(colnames(dat) %in% omit)]

	#Omitting sites without fish:

		dat3 <- dat2[rowSums(dat2[,16:36]) >0, ]

	#Treat missing environmental data:

		dat3$Herbaceous[is.na(dat3$Herbaceous)] <- 0 
		dat3$Ann_Herb[is.na(dat3$Ann_Herb)] <- 0

		dat3 <- dat3[complete.cases(dat3$SiteLength),]
		dat_final <- dat3

	#Split the datasets:

		fish <- dat_final[,16:36]
		env <- dat_final[,1:15]

	#Removing very rare species:

  		fish_red <- drop.var(fish, min.fo=2)

	#Standardizing by unit effort (segment length) to ind/m:

  		fish_dens <- fish_red
  		for(i in 1:nrow(fish_red)){	
    		fish_dens[i,] <- fish_red[i,]/env$SiteLength[i]     }

	#Trimming the environmental data set:

    		drop <- c("Latitude","Longitude","SiteLength","SiteWidth","SurfaceArea") 
    		env <- env[,!(colnames(env) %in% drop)] 
    		head(env) 
    		stat.desc(env)

  		env_cont <- env[,!(colnames(env) %in% c("SMU","Pop","NLCD_Cat"))]
  		head(env_cont)

	#Checking for multivariate outliers:

  		mv.outliers(fish_dens, method = "euclidean", sd.limit=2)

	#Remember, unusually high abundances of speckled dace seem to be driving species trends at a 
		#few sites, but we're leaving them in for now.

  		mv.outliers(env_cont, method = "euclidean", sd.limit=2)

	#None of the environmental outlier sites overlap with the species outlier sites.

	#Remove co-varying environmental factors:

  		env <- env[,!(colnames(env) %in% c("Ave_Max_D","Ann_Herb"))]
  		env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]

### Distribution of Data: Fish Abundance

	#Let's first re-examine how fish abundances are distributed.

	#A barplot of the distribution:

  		ac <- table(unlist(fish))
  
		barplot(ac, 
			las = 1, 
			xlab = "Abundance class", 
			ylab = "Frequency", 
			col = gray(length(ac): 0/length(ac)),
			ylim=c(0,100)
			)

	#The data are right skewed and also what we would call "zero-skewed." The overabundance of 
		#20-values is driven by speckled dace and brook trout (it's an error in the dataset that 
		#I couldn't completely rectify...we can ignore it for now).

	#We can see that the skewness of the dataset indicates it will need to be transformed!

### Transformation and Standardization: Fish

	#Examine the distributional properties of the data using the `uv.plots()` function. 
		#They will definitely need to be transformed.

		uv.plots(fish)
		
	#You can do this item-by-item. For example, for redband trout:

		par(mfrow = c(2,2))
		hist(fish_dens$TROUT_RB, 
			xlab="Fish (ind/m)", 
			ylab="Frequency", 
			main="Raw")
		hist(sqrt(fish_dens$TROUT_RB), 
			xlab="Fish (ind/m)", 
			ylab="Frequency", 
			main="Square Root")
		hist(log(fish_dens$TROUT_RB+1), 
			xlab="Fish (ind/m)", 
			ylab="Frequency", 
			main="Log + 1")
		hist(log(fish_dens$TROUT_RB+0.1), 
			xlab="Fish (ind/m)", 
			ylab="Frequency", 
			main="Log + 0.1")

	#We also can use the data.trans() function to play with transformations:

		data.trans(fish_dens, method="log", plot=F) #For heavily skewed data
		data.trans(fish_dens, method="power", exp=0.5, plot=F) #For slightly skewed data
		data.trans(fish_dens, method="asin", plot=F) #For proportional data (scaled 0-1, not relevant here)

	#We can use the `decostand()` function in the "vegan" package to standardize the data. 
		#See `?decostand` for more information.

	#The coefficient of variation value (cv) in the `stat.desc()` function will provide clues as to 
		#whether standardization is necessary. If cv \< 50, standardization won't make a difference.
		#If cv \> 100 consider standardizing. Since cv values are very low for the species data set, 
		#let's skip standardizing for now.

		dev.off()
		
		boxplot(fish_dens$TROUT_RB,
			sqrt(fish_dens$TROUT_RB),
			log1p(fish_dens$TROUT_RB),
			las = 1,
			main = "Simple Transformations",
			names = c("raw data","sqrt","log"),
			col = "#FF5050"
			)

	#For future analyses, we will log-transform the data as this appears to be the best option for correcting skew 
		#and outliers. If your data has zeros, you will need to "log+1" transform it. Another option (shown above)
		#is to base the value you use on the scale of the data. In this case, we've tried 0.1. Although this 
		#adjustment does a better job of correcting skew, it results in negative values, which can be challenging 
		#to interpret for some analyses. For now, let's plan on just using the log+1 transformation, even though we 
		#know it isn't perfect.

		fish_log <- log(fish_dens + 1)

### Transformation and Standardization: Environmental Data

	#As above, first let's examine the distributional properties of the data using the `uv.plots()` 
		#function or with an item-by-item histogram.

		par(mfrow = c(3,2))
		hist(env$Max_Depth, 
			xlab="Max Depth (m)", 
			main=NA)
		hist(env$Gradient, 
			xlab="Gradient (%)", 
			main=NA)
		hist(env$Elev, 
			xlab="Elevation (m)", 
			main=NA)
		hist(env$Canopy, 
			xlab="Canopy Cover (%)", 
			main=NA)
		hist(env$Herbaceous, 
			xlab="Herbaceous Cover (%)", 
			main=NA)

	#It looks like some variables are normally distributed and some aren't.

		env.log <- data.trans(env_cont, method="log", plot=F)
		env.power <- data.trans(env_cont, method="power", exp=0.5, plot=F)
		env.asin <- data.trans(env_cont, method="asin", plot=F)

	#What about standardization?

  		stat.desc(env_cont)

	#It looks like cv values are still low; however, we know that these variables are mixed-scale!

		env.scal <- decostand(env_cont, "max") #Standardization by max value of each column

		env.relsp <- decostand(env_cont, "total", MARGIN=2) #Standardization by column totals

		env.zscore <- decostand(env_cont, "standardize") #Z-scores the data in each column

		env.rel <- decostand(env_cont, "total") #Standardization by max value of each site/row

		env.norm <- decostand(env_cont, "normalize") #Give a length of 1 to each row vector (chord 
			#transformation), very useful for Euclidean distances (PCA, RDA) and can be used on 
			#log-transformed data

		env.hel <- decostand(env_cont, "hellinger") #Square root of relative values per site, obtained 
			#by applying chord transformation to square root transformed data

		env.chi <- decostand(env_cont, "chi.square") #Double standardization by columns and rows

		env.wis <- wisconsin(env_cont) #Wisconsin standardization, values are ranged by column maxima 
			#and then by site totals

	#Let's see what each of these transformations looks like for the "Max Depth (m)" variable:

		par(mfrow = c(2,2))
		boxplot(env_cont$Max_Depth,
			env.power$Max_Depth,
			env.log$Max_Depth,
			las = 1,
			main = "Simple Transformations",
			names = c("raw data","sqrt","log"),
			col = "#FF5050"
			)
		boxplot(env.scal$Max_Depth,
			env.relsp$Max_Depth,
   			env.zscore$Max_Depth,
			las = 1,
			main = "Standardizations by Variable",
			names = c("max","total","Z-score"),
			col = "#9BE1AF"
			)
		boxplot(env.hel$Max_Depth,
			env.rel$Max_Depth,
			env.norm$Max_Depth,
			las = 1,
			main = "Standardizations by Sites",
			names = c("Hellinger","total","norm"),
			col = "#46AFAA"
			)
		boxplot(env.chi$Max_Depth,
			env.wis$Max_Depth,
			las = 1,
			main = "Double Standardizations",
			names = c("Chi-square","Wisconsin"),
			col = "#FAD223"
			)

	#Since the goal is to standardize across environmental variables with different units, consider standardizing 
		#by column. Z-scoring is one of the most commonly used options, but what works best for your data 
		#may be different depending on the circumstances.

