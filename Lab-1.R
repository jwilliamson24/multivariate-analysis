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

### Load the Data Set(s)

	#If you would like, you can set a working directory prior to loading your data. 

	#The first column in your csv should be the row names. If you have response and predictor variables, you 
		#can load them as two separate matrices; however, I find that data sets are easiest to clean when you 
		#load everything at once, as below.

		#setwd("U:\\FW599_Multivariate\\Data\\")
		dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)

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

		#stat.desc(dat) #From the 'pastecs' package, compute and display descriptive statistics for the data 
					#frame, including statistical summaries

### Data Cleaning and Preparation

	#In the class data set, columns 16-48 are the response matrix: the fish community.

		head(dat[,16:48])

	#The response varibles, or environmental data include:

	#-   SMU (Categorical): The watershed or basin where each sampling site is located.
	#-   Pop (Categorical): The river system where each sampling site is located.
	#-   Latitude/Longitude: Each site's GPS coordinates.These values are mostly useful if we're assessing 
		#spatial autocorrelation.
	#-   SiteLength (m): The length of the sampled stream segment.
	#-   SiteWidth (m): The width of the sampled stream segment.
	#-   SurfaceArea (m\^2): The estimated total surface area of the sampled stream segment.
	#-   Max_Depth (m): The maximum observed depth of the sampled stream segment.
	#-   Ave_Max_D (m): The average maximum observed depth of the segment.
	#-   Gradient (%): Estimated stream slope.
	#-   Elevation (m): The site's elevation.
	#-   NLCD_Cat (Categorical): NLCD Land Cover category, 2008 <https://www.mrlc.gov/>
	#-   Canopy (%): Tree canopy cover, 2011 <https://www.mrlc.gov/>
	#-   Herbaceous (%): Herbaceous rangeland cover, 2016 <https://www.mrlc.gov/>
	#-   Ann_Herb (%): Annual herbaceous rangeland cover, 2016 <https://www.mrlc.gov/>

	#If you don't know how R has classified something, you can determine the class of a given object 
		#with `class()`. This is especially helpful when you're getting weird errors because R imported a 
		#numeric value as a character value, etc. You can do this for all variables at once using the `str()` 
		#and `stat.desc()` functions shown above.

		class(dat$SMU)

### Omit Missing Data

	#Let's start by finding and omitting/treating missing data points.

	#First, we can drop any species (columns) that have zero observations because they contribute nothing to the analysis.

		colSums(dat[,16:48])

	#We have quite a few fish species with no observations.

		spp_N <- colSums(dat[,16:48])
		spp_0 <- subset(spp_N, spp_N == 0)
		omit <- names(spp_0)

		dat2 <- dat[,!(colnames(dat) %in% omit)]

	#Fish columns are now 16-36

		length(omit)
		colSums(dat2[,16:36])

	#What about sites where no fish were observed at all? There are different ways to treat sites without 
		#observations, but for now, let's omit those sites.

		dat3 <- dat2[rowSums(dat2[,16:36]) >0, ]

	#An alternative to omitting fishless sites would be to add a "no fish" column to the response matrix. 
		#This can be problematic for some analyses (ordination), but works better for presence/absence data.

### Missing Environmental Data

	#We can treat missing environmental variables in a variety of ways. Some analyses will support missing 
		#values while others won't.

	#For Herbaceous and Annual Herbaceous cover, we can convert the NAs to 0s based on what we know about 
		#these data (values were only calculated where herbaceous plants were present).

		dat3$Herbaceous[is.na(dat3$Herbaceous)] <- 0 
		dat3$Ann_Herb[is.na(dat3$Ann_Herb)] <- 0

	#Because we'll need to scale the fish data (haha) by survey effort, let's omit any sites where SiteLength 
		#was not measured.

		dat3 <- dat3[complete.cases(dat3$SiteLength),]
		dat_final <- dat3

	#And there we have it, a cleaned data set!

### Splitting the Data

	#We can split the data into separate "fish" and "environmental" matrices. i.e., our "response" and "predictor" 
		#matrices.

		fish <- dat_final[,16:36]
		env <- dat_final[,1:15]

	#Make sure to complete any of the aforementioned data processing/cleaning steps BEFORE you split the matrices! 
		#The exception to this is if you need to replace missing values differently for your response and predictor 
		#variables.

### Treatments for Missing Data

	#If you have NAs/missing data in your species/response matrix, you can replace small amounts of missing data 
		#with the median or mean of the column.

	#For count data, using the median (default) option is best because it spits out an integer value

		testfish <- replace.missing(fish)

	#For the environmental/predictor matrix, a mean replacement might be best for a continuous scale

		testenv <- replace.missing(env, method="mean")	

	#It gives us an error because some of our variables are categorical. We can ignore this.

	#Now it's time for the fun stuff!

### Exploratory Data Analysis: Fish Abundance

	#Let's first examine how fish abundances are distributed.

		range(fish) #Min/max values

		apply(fish, 2, range) #Min/max values for each species

		ac <- table(unlist(fish)) #Number of cases for each abundance class

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

		spe_pres <- apply(fish > 0, 2, sum) #Number of occurrences
		sort(spe_pres)

		spe_relf <- 100*spe_pres/nrow(fish) #Relative frequency of occurrences
		round(sort(spe_relf), 1)

		par(mfrow = c(1,2))
		hist(spe_pres,
			main = "Species Occurrences",
			right = FALSE,
			las = 1,
			xlab = "Number of occurrences",
			ylab = "Number of species",
			breaks = seq(0, max(spe_pres), by=5),
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

		testfish <- drop.var(fish, min.po=5)
  		ncol(fish)
  		ncol(testfish)

	#This actually omits quite a few species from our data set!

	#We can also drop species based on the number of non-zero values.

		testfish <- drop.var(fish, min.fo=10)
  		ncol(fish)
  		ncol(testfish)

	#You might also consider omitting over-abundant (generalist) species. If you know they occur pretty much 
		#everywhere and won't contribute anything to an analysis of community-level differences (I would 
		#caution folks against this, as it inches deeper into the realm of subjectivity).

		testfish <- drop.var(fish, max.po=90)
  		ncol(fish)
  		ncol(testfish)

	#In our case, we don't have any "generalist" species to remove.

	#Let's work with a very conservative cutoff for now.

		fish_red <- drop.var(fish, min.fo=2)

	#We're just dropping yellow perch from the data set.

	#Remember to check that removing these species didn't create empty rows!

### Scaling Data by Per-Unit Effort

	#Our raw data are ordinal (counts) of species captured using electrofishing. In this case, ODFW 
		#did not sample the same length of stream at each site.

		hist(env$SiteLength, col="#FF5050")

	#We should consider standardizing by unit effort (segment length) to ind/m. This can be done using a 
		#for loop as follows:

		fish_dens <- fish_red   #Create new data frame

		for(i in 1:nrow(fish_red)){   #For each row (i)
			fish_dens[i,] <- fish_red[i,]/env$SiteLength[i]   #Divide the number of individuals by segment length
	  		}   #End Loop

### Checking for Outliers

	#There is no rule for what constitutes an "extreme value," so use your best judgment

	#Univariate outlier check (redband trout and speckled dace):

		uv.outliers(fish_dens, id="BASS_UNID:TROUT_RB", var="TROUT_RB", sd.limit=2)
		uv.outliers(fish_dens, id="BASS_UNID:TROUT_RB", var="DACE_SP", sd.limit=2)

	#Look at the last column in the output. This column tells us how many standard deviations observed 
		#abundance at a given site deviates from average. Anything \>2 standard deviations is considered 
		#an outlier; however, you may consider other criteria first. What about from a multivariate 
		#(multi-species) perspective?

		mv.outliers(fish_dens, method = "euclidean", sd.limit=2)

	#This procedure identified three sites with multivariate outliers. Note that these are the same sites 
		#that were identified for speckled dace. So it appears that unusally high abundances of speckled dace 
		#are driving species trends at these sites. If you were not interested in anomalies in dace abundance, 
		#you might consider removing these rows. Otherwise, you could wait until you delve deeper into your 
		#analysis to make that call. For now, let's note them and move on.

### Checking for covariance in abundant fish

	#You can use a correlation bi-plot to examine whether certain species tend to occur concurrently with each other. This may be visually helpful, but isn't as analytically crucial as checking for covariance in environmental (predictor) variables.

		abu_fish <- drop.var(fish_dens, min.po=5)

		P.corr <- cor(abu_fish, method = "pearson", use = "complete.obs")
		round(P.corr, 2)

		corrplot(P.corr, 
			type = "upper", 
			order = "hclust", 
			tl.col = "black", 
			tl.srt = 45)

### Exploratory Data Analysis: Environmental Variables

	#Use the `summary()` function to produce a Summary of the environmental data.

	#There are many variables, but we're only interested in some of these as predictors: 
		#SMU, Pop, Max_Depth, Ave_Max_D, Gradient, Elev, NLCD_Cat, Canopy, Herbaceous, and Ann_Herb

		drop <- c("Latitude","Longitude","SiteLength","SiteWidth","SurfaceArea")
		env <- env[,!(colnames(env) %in% drop)]
  			head(env)
		stat.desc(env)

	#We can also parse out continuous variables first

		env_cont <- env[,!(colnames(env) %in% c("SMU","Pop","NLCD_Cat"))]
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

