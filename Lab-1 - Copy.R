###############################################################
###Lab-1.R									###
###Data Screening and Exploration				 	###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 1, 2024						###
###############################################################

### Lab Objectives

	#Data screening is an essential precursor to any data analysis. In this lab, we will learn to:

		#Load data and R packages
		#Generate summary statistics and screen data for errors
		#Evaluate and correct missing data
		#Check for outliers
		#Check for multicollinearity

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
    dat <- readRDS("site_level_matrix.rds")
    row.names(dat) <- dat[,1]
    sals <- dat[19:20]
    env <- dat[,-c(19:20)]
    
    drop <- c("lat","long","landowner","stand","tree_farm","year","weather")
    env <- env[,!(colnames(env) %in% drop)]
    
    env_cont <- env[,c("elev","temp","hum","soil_moist","canopy_cov","veg_cov","dwd_cov","fwd_cov","jul_date")]

### Check the Data Structure

	#Some helpful prompts:

		head(dat) 	# Displays first 6 lines of the data set
		tail(dat) 	#Displays last 6 lines of the data set
		dat[1:5, 1:10] 	#Displays i = 5 rows and j = 10 columns

		nrow(dat) #Number of rows (sites)
		ncol(dat) #Number of columns (habitat variables/species)

		str(dat) #Examine the structure of the data set including information about the types of columns 
				#and a preview of the data

		#stat.desc(dat) #From the 'pastecs' package, compute and display descriptive statistics for the data 
					#frame, including statistical summaries

### Splitting the Data

	#We can split the data into separate "fish" and "environmental" matrices. i.e., our "response" and "predictor" 
		#matrices.
# 
# 		sals <- dat[,24:29]
# 		env <- dat[,1:22]

	#Make sure to complete any of the aforementioned data processing/cleaning steps BEFORE you split the matrices! 
		#The exception to this is if you need to replace missing values differently for your response and predictor 
		#variables.

### Exploratory Data Analysis: Fish Abundance

	#Let's first examine how fish abundances are distributed.

		range(sals) #Min/max values

		apply(sals, 2, range) #Min/max values for each species

		ac <- table(unlist(sals2)) #Number of cases for each abundance class

	#A barplot of the distribution:

		barplot(ac, 
			las = 1, 
			xlab = "Abundance class", 
			ylab = "Frequency", 
			col = gray(length(ac): 0/length(ac)),
			ylim=c(0,100)
			)

	#The data are right skewed and also what we would call "zero-skewed." 
		#The skewness of the data indicates that they will likely need to be transformed!

	#At how many sites does each species occur?

		spe_pres <- apply(sals > 0, 2, sum) #Number of occurrences
		sort(spe_pres)
		# ANFE TAGR PLDU AMGR ENES  OSS 
		#   1    3    4    5   107  163 
		
		
		spe_relf <- 100*spe_pres/nrow(sals) #Relative frequency of occurrences
		round(sort(spe_relf), 1)
		# ANFE  TAGR  PLDU  AMGR  ENES  OSS 
		#  0.1   0.3  0.4   0.6   12.0  18.3 
		
		
		par(mfrow = c(1,2))
		hist(spe_pres,
			main = "Species Occurrences",
			right = FALSE,
			las = 1,
			xlab = "Number of occurrences",
			ylab = "Number of species",
			breaks = seq(0, max(spe_pres)+5, by=5),
			col = "#FF5050"
			)
		hist(spe_relf,
			main = "Species Relative Frequencies",
			right = FALSE,
			las = 1,
			xlab = "Frequency of occurrences (%)",
			ylab = "Number of species",
			breaks = seq(0, 100, by=5),
			col = "#FF5050"
			)

	#The `foa.plots()` function will produce a series of plots with this information as well

### When to drop "insufficient" or overabundant species

	#Depending on how your data were collected and your plans for analysis, you might consider dropping 
		#"rare" descriptors (species). This is especially important for analyses that ascribe a greater 
		#importance to rare species (e.g., Canonical Correspondence Analysis, which uses a Chi-square distance). 
		#A good rule of thumb for larger data sets is to omit species that have non-zero values in less than 5% 
		#of sites.

		testsals <- drop.var(sals, min.po=5)
  		ncol(sals)
  		ncol(testsals)

	#This actually omits quite a few species from our data set!

	#We can also drop species based on the number of non-zero values.

		testsals <- drop.var(sals, min.fo=10)
  		ncol(sals)
  		ncol(testsals)

	#You might also consider omitting over-abundant (generalist) species. If you know they occur pretty much 
		#everywhere and won't contribute anything to an analysis of community-level differences (I would 
		#caution folks against this, as it inches deeper into the realm of subjectivity).

		testsals <- drop.var(sals, max.po=90)
  		ncol(sals)
  		ncol(testsals)

	#In our case, we don't have any "generalist" species to remove.

	#Let's work with a very conservative cutoff for now.

		sals_red <- drop.var(sals, min.fo=2)


### Exploratory Data Analysis: Environmental Variables

	#Use the `summary()` function to produce a Summary of the environmental data.

	#There are many variables, but we're only interested in some of these as predictors

		drop <- c("date","lat","long")
		env <- env[,!(colnames(env) %in% drop)]
  			head(env)
		stat.desc(env)

	#We can also parse out continuous variables first

		env_cont <- env[,c("elev","temp","hum","soil_moist_avg","jul_date")]
		head(env_cont)

### Checking for Outliers

	#First, let's check for outliers in our environmental data. We can use the same multivariate 
		#procedure as for the species data above.

		mv.outliers(env_cont, method = "euclidean", sd.limit=2)

	#If any of these sites showed up as outliers for both the species and environmental data, I would 
		#consider omitting them, especially if the SD was \> 4. For now, we can keep them in and assess later.

### Checking for Covariance in Predictors

	#Next we'll want to determine if any of our predictors covary. This is important to know, even if the analysis 
		#we're using accounts for correlation among predictors.

		pairs(env_cont,
			panel = panel.smooth,
			main = "Bivariate Plots with Smooth Curves")

		P.corr <- cor(env_cont, method = "pearson", use = "complete.obs")
	  		round(P.corr, 2)
    
	  dev.off()
		corrplot(P.corr, 
			type = "upper", 
			order = "hclust", 
			tl.col = "black", 
			tl.srt = 45)

	#This indicates thattemp and hum covary substantially. julian date and
		# soil moisture also covary. We can remove co-varying factors if we want.


	#What about relationships between categorical and continuous variables? Let's look at among-basin differences.

		p1 <- ggplot(env) + geom_boxplot(aes(x = trt, y = soil_moist), fill = "#FF5050", alpha=0.8)
		p2 <- ggplot(env) + geom_boxplot(aes(x = trt, y = temp), fill = "#FF5050", alpha=0.8)
		p3 <- ggplot(env) + geom_boxplot(aes(x = trt, y = hum), fill = "#FF5050", alpha=0.8)
		p4 <- ggplot(env) + geom_boxplot(aes(x = trt, y = canopy_cov), fill = "#FF5050", alpha=0.8)
		p5 <- ggplot(env ) + geom_boxplot(aes(x = trt, y = dwd_cov), fill = "#FF5050", alpha=0.8)

  	ggarrange(p1, p2, p3, p4, p5, ncol=2, nrow=3)

