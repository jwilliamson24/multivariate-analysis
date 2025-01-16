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
        setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
        source("Biostats.R")

#Loading subplot level data:
    # dat1 <- read.csv("habitat.occu.complete.csv", row.names = 1)
    # dat <- readRDS("habitat.occu.complete.rds")
    # row.names(dat) <- dat[,1]
     
    
#Loading site level data
    dat2 <- readRDS("site_level_df.rds")
    row.names(dat2) <- dat2[,1]
    sals <- dat2[19:20]
    env <- dat2[1:18]
    
    drop <- c("lat","long")
    env <- env[,!(colnames(env) %in% drop)]
    
    env_cont <- env[,c("elev","temp","hum","soil_moist","canopy_cov","veg_cov","dwd_cov","fwd_cov","jul_date")]

	
# #Subplot-level data
#     sals <- dat[,24:29]
#     env <- dat[,1:22]
# 
#     sals_red <- drop.var(sals, min.fo=5) #Removing very rare species:
# 
# 	#Standardizing by unit effort (segment length) to ind/m:
#   # 		fish_dens <- fish_red
#   # 		for(i in 1:nrow(fish_red)){	
#   #   		fish_dens[i,] <- fish_red[i,]/env$SiteLength[i]     }
# 
# 	#Trimming the environmental data set:
# 
#   		drop <- c("date","lat","long")
#   		env <- env[,!(colnames(env) %in% drop)]
#   		head(env)
#   		stat.desc(env)
#   		
#   		#We can also parse out continuous variables first
#   		
#   		env_cont <- env[,c("elev","temp","hum","soil_moist_avg","canopy_cov","veg_cov","dwd_cov","fwd_cov")]
#   		head(env_cont)

#Start Exploring Data
	#Checking for multivariate outliers:

  		saloutlier <- mv.outliers(sals, method = "euclidean", sd.limit=2)
  		envoutlier <- mv.outliers(env_cont, method = "euclidean", sd.limit=2)
  		intersect(rownames(saloutlier),rownames(envoutlier))

	#There is one site that repeats in both outlier dfs. "12213 _ 1 _ 2023"

	#Remove co-varying environmental factors:

  		# env <- env[,!(colnames(env) %in% c("Ave_Max_D"))]
  		# env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]

### Distribution of Data: Fish Abundance

	#Let's first re-examine how fish abundances are distributed.

	#A barplot of the distribution:

  	ac <- table(unlist(sals))
  
  	 dev.off()
		barplot(ac, 
			las = 1, 
			xlab = "Abundance class", 
			ylab = "Frequency", 
			col = gray(length(ac): 0/length(ac)),
			ylim=c(0,150)
			)

	#over 5000 zeros in the oss occu data
	#The data are right skewed and also what we would call "zero-skewed."

	#We can see that the skewness of the dataset indicates it will need to be transformed!

### Transformation and Standardization

	#Examine the distributional properties of the data using the `uv.plots()` function. 
		#They will definitely need to be transformed.

		uv.plots(sals)
		
	# oss and enes are the only spp with enough data
	#cloglog???????????????????????????????????????????????????????????	
		
	#You can do this item-by-item. 

		par(mfrow = c(3,2))
		hist(sals$oss, 
			xlab="count", 
			ylab="Frequency", 
			main="Raw")
		hist(sqrt(sals$oss), 
			xlab="count", 
			ylab="Frequency", 
			main="Square Root")
		hist(log(sals$oss), 
			xlab="count", 
			ylab="Frequency", 
			main="Log")
		hist(log(sals$oss+0.1), 
			xlab="count", 
			ylab="Frequency", 
			main="Log + 0.1")
		hist(log(sals$oss+1), 
		     xlab="count", 
		     ylab="Frequency", 
		     main="Log + 1")
		
	#so many zeros, makes it hard to figure this out. do i have enough data for this?

	#We also can use the data.trans() function to play with transformations:

		data.trans(sals, method="log", plot=F) #For heavily skewed data
		data.trans(sals, method="power", exp=0.5, plot=F) #For slightly skewed data
		data.trans(sals, method="asin", plot=F) #For proportional data (scaled 0-1, not relevant here)

	#We can use the `decostand()` function in the "vegan" package to standardize the data. 
		#See `?decostand` for more information.

	#The coefficient of variation value (cv) in the `stat.desc()` function will provide clues as to 
		#whether standardization is necessary. If cv \< 50, standardization won't make a difference.
		#If cv \> 100 consider standardizing. Since cv values are very low for the species data set, 
		#let's skip standardizing for now.

		stat.desc(sals)
		#coeff var is low so i guess standardizing these data doesnt make sense anyways
		
		dev.off()
		
		boxplot(sals$oss,
			sqrt(sals$oss),
			log1p(sals$oss),
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
		hist(env$elev, 
			xlab="elev", 
			main=NA)
		hist(env$temp, 
			xlab="temp", 
			main=NA)
		hist(env$soil_moist, 
			xlab="soil moist", 
			main=NA)
		hist(env$hum, 
			xlab="hum", 
			main=NA)
		hist(env$dwd_cov, 
			xlab="dwd", 
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
		boxplot(env_cont$temp,
			env.power$temp,
			env.log$temp,
			las = 1,
			main = "Simple Transformations",
			names = c("raw data","sqrt","log"),
			col = "#FF5050"
			)
		boxplot(env.scal$temp,
			env.relsp$temp,
   			env.zscore$temp,
			las = 1,
			main = "Standardizations by Variable",
			names = c("max","total","Z-score"),
			col = "#9BE1AF"
			)
		boxplot(env.hel$temp,
			env.rel$temp,
			env.norm$temp,
			las = 1,
			main = "Standardizations by Sites",
			names = c("Hellinger","total","norm"),
			col = "#46AFAA"
			)
		boxplot(env.chi$temp,
			env.wis$temp,
			las = 1,
			main = "Double Standardizations",
			names = c("Chi-square","Wisconsin"),
			col = "#FAD223"
			)

	#Since the goal is to standardize across environmental variables with different units, consider standardizing 
		#by column. Z-scoring is one of the most commonly used options, but what works best for your data 
		#may be different depending on the circumstances.

		par(mfrow = c(3,2))
		hist(env.zscore$elev, 
		     xlab="elev", 
		     main=NA)
		hist(env.zscore$temp, 
		     xlab="temp", 
		     main=NA)
		hist(env.zscore$soil_moist, 
		     xlab="soil moist", 
		     main=NA)
		hist(env.zscore$hum, 
		     xlab="hum", 
		     main=NA)
		hist(env.zscore$dwd_cov, 
		     xlab="dwd", 
		     main=NA)
		