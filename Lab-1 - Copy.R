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

### Load the Necessary Packages

	#The `vegan` package is the quintessential community ecology package. It provides functions for ordination, 
		#diversity analysis, and other multivariate analyses commonly used in ecological studies.

	#The `pastecs` package includes tools for descriptive statistics, handling and analyzing time series data, 
		#and other useful statistical methods.

	#The `corrplot`, `ggplot2`, and `ggpubr` packages are used to eloquently visualize data.

    library(vegan)
		library(pastecs)
		library(corrplot)
		library(ggplot2)
		library(ggpubr)
    library(dplyr)

### Load the Data Set(s)

	#If you would like, you can set a working directory prior to loading your data. 

	#The first column in your csv should be the row names. If you have response and predictor variables, you 
		#can load them as two separate matrices; however, I find that data sets are easiest to clean when you 
		#load everything at once, as below.

		#setwd("U:\\FW599_Multivariate\\Data\\")
		#dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
    dat <- read.csv("sals.complete.csv")

	#In this case, rows are named by sampling site. As mentioned in the "notes" box above, row names cannot 
		#be duplicates! If this is an issue for you, you will need to assign a unique "Site ID" for each row.

	#This source code from Julian Olden has some useful functions for cleaning and processing multivariate data:

		source("biostats.R")

### Check the Data Structure

	#Some helpful prompts:

		head(dat) 	# Displays first 6 lines of the data set
		tail(dat) 	#Displays last 6 lines of the data set
		dat[1:5, 1:10] 	#Displays i = 5 rows and j = 10 columns

		nrow(dat) #Number of rows (sites)
		ncol(dat) #Number of columns (habitat variables/species)

		str(dat) #Examine the structure of the data set including information about the types of columns 
				#and a preview of the data

		stat.desc(dat) #From the 'pastecs' package, compute and display descriptive statistics for the data 
					#frame, including statistical summaries

# did the below steps already in oss-occu repo
### Data Cleaning and Preparation
### Omit Missing Data
### Splitting the Data

		sals <- dat
		site <- read.csv("site.complete.csv")
    subplot <- read.csv("subplot.complete.csv")

### Treatments for Missing Data
# no NA's, dealt with in oss-occu repo
    
    na.count <- sals %>%
      summarise(across(everything(), ~ sum(is.na(.))))
    print(na.count)
    
    na.count <- site %>%
      summarise(across(everything(), ~ sum(is.na(.))))
    print(na.count)
    
    na.count <- subplot %>%
      summarise(across(everything(), ~ sum(is.na(.))))
    print(na.count)

# not relevant or already did it
### Exploratory Data Analysis: Fish Abundance
### When to drop "insufficient" or overabundant species
### Scaling Data by Per-Unit Effort

### Checking for Outliers
# i only have two main species so i dont think this is something that will apply to me
 
      sals$spp <- as.factor(sals$spp)
      summary(sals$spp)
      # could remove other species here, counts too low to matter
      # AMGR  ANFE  ENES  OSS  PLDU  TAGR 
      #   5    1    138   257    5    3 


### Exploratory Data Analysis: Environmental Variables

      #lets delete the columns we arent interested in as predictors
      summary(site)
      rownames(site) <- site$site_id
      drop <- c("site_rep","date","elev","stand","year","landowner","site_id")
      site <- site[,!(colnames(site) %in% drop)]
      head(site)
      stat.desc(site)
      site_cont <- site[,c("temp","hum")]
      
      summary(subplot)
	    drop <- c("site_rep","date","lat","long","time","stand","year")
      subplot <- subplot[,!(colnames(subplot) %in% drop)]
      head(subplot)
      stat.desc(subplot)

### Checking for Outliers

	#First, let's check for outliers in our environmental data. We can use the same multivariate 
		#procedure as for the species data above.

		mv.outliers(site_cont, method = "euclidean", sd.limit=2)
		# avedist    sd
		# 12441 _ 1 _ 2023   37.297 2.039
		# 21159 _ 1 _ 2023   39.262 2.335
		# 12560 _ 1 _ 2024   39.493 2.370
		# 12776 _ 1 _ 2024   50.635 4.045
		# 33628 _ 1 _ 2024   50.272 3.990
		# 33651 _ 1 _ 2024   46.584 3.436
		# 408180 _ 1 _ 2024  38.583 2.233
		
	#not sure what to do with these

### Checking for Covariance in Predictors  ########################################################################
		
		# to do this i want to make a dataframe that has all site and subplot data in the same matrix

	#Next we'll want to determine if any of our predictors covary. This is important to know, even if the analysis 
		#we're using accounts for correlation among predictors.

		pairs(env_cont,
			panel = panel.smooth,
			main = "Bivariate Plots with Smooth Curves")

		P.corr <- cor(env_cont, method = "pearson", use = "complete.obs")
	  		round(P.corr, 2)

		corrplot(P.corr, 
			type = "upper", 
			order = "hclust", 
			tl.col = "black", 
			tl.srt = 45)

	#This indicates that Maximum Depth and Average Maximum Depth covary substantially. Percent Herbaceous and 
		#Annual Herbaceous plants also covary. We can remove co-varying factors if we want.

		env <- env[,!(colnames(env) %in% c("Ave_Max_D","Ann_Herb"))]
		env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]

	#What about relationships between categorical and continuous variables? Let's look at among-basin differences.

		p1 <- ggplot(env) + geom_boxplot(aes(x = SMU, y = Max_Depth), fill = "#FF5050", alpha=0.8)
		p2 <- ggplot(env) + geom_boxplot(aes(x = SMU, y = Gradient), fill = "#FF5050", alpha=0.8)
		p3 <- ggplot(env) + geom_boxplot(aes(x = SMU, y = Elev), fill = "#FF5050", alpha=0.8)
		p4 <- ggplot(env) + geom_boxplot(aes(x = SMU, y = Canopy), fill = "#FF5050", alpha=0.8)
		p5 <- ggplot(env) + geom_boxplot(aes(x = SMU, y = Herbaceous), fill = "#FF5050", alpha=0.8)

  	ggarrange(p1, p2, p3, p4, p5, ncol=2, nrow=3)

