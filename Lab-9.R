###############################################################
###Lab-9.R									###
###Constrained Ordination Techniques				###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 29, 2024						###
###############################################################

## Lab Objectives

	#Constrained ordination is a multivariate statistical technique used to explore and quantify the 
		#relationships between response (usually species) variables and predictor (usually environmental) 
		#variables. These methods are crucial in ecological studies, as they help to uncover the driving 
		#environmental factors shaping species communities.

	#In this lab, we will learn to:

	#-   Use asymmetric constrained ordination techniques, such as **redundancy analysis (RDA)** and 
		#**canonical correspondence analysis (CCA)**, to explore species-environment relationships.
	#-   Explore **transformation-based RDA (tbRDA)** to improve ordination results when species abundances 
		#exhibit unimodal or skewed distributions.
	#-   Implement **distance-based redundancy analysis (dbRDA)** to handle nonlinear species-environment 
		#relationships using dissimilarity measures.
	#-   Perform model selection using stepwise approaches and permutation tests to identify significant 
		#environmental predictors.
	#-   Visualize the relationships between species and environmental variables using ordination plots 
		#(triplots).
	#-   Investigate **co-inertia analysis (CoIA)** for symmetric comparisons between species and 
		#environmental data matrices.

## Setting up the R Workspace and Preparing the Data

	#Load the necessary packages and prepare the dataset as we have done previously.

  		library(vegan)
  		library(viridis)
  		library(MASSExtra)
  		library(ade4)

      setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
      source("Biostats.R")
      
      dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)

	#This time let's use the full dataset to capture the entire suite of environmental conditions and species! 
		#Just a warning, it's going to look messy...

  		spp_N <- colSums(dat[,16:ncol(dat)])
  		spp_0 <- subset(spp_N, spp_N == 0)
  		omit <- names(spp_0)

  		dat2 <- dat[,!(colnames(dat) %in% omit)]
	
  		dat3 <- dat2[rowSums(dat2[,16:ncol(dat2)]) >0, ]
  
  		dat3$Herbaceous[is.na(dat3$Herbaceous)] <- 0 
  		dat3$Ann_Herb[is.na(dat3$Ann_Herb)] <- 0
	
  		dat3 <- dat3[complete.cases(dat3$SiteLength),]
  
  		dat_final <- dat3
	
  		fish <- dat_final[,16:ncol(dat_final)]
  		env <- dat_final[,1:15]
	
  		fish_red <- drop.var(fish, min.fo=1)
  
  		fish_dens <- fish_red
  		for(i in 1:nrow(fish_red)){	
    			fish_dens[i,] <- fish_red[i,]/dat_final$SiteLength[i]
  			}
  
  		fish_dens_log <- log(fish_dens + 1)
  
  		drop <- c("Latitude","Longitude","SiteLength","SiteWidth","SurfaceArea")
		env <- env[,!(colnames(env) %in% drop)]
		env_cont <- env[,!(colnames(env) %in% c("SMU","Pop","NLCD_Cat"))]
	
		env <- env[,!(colnames(env) %in% c("Ave_Max_D","Ann_Herb"))]
		env_cont <- env_cont[,!(colnames(env_cont) %in% c("Ave_Max_D","Ann_Herb"))]

## Redundancy Analysis (RDA)

	#Redundancy analysis (RDA) is a constrained ordination technique used to test directional hypotheses 
		#about the effect of a "predictor" (usually environmental) matrix **X** on a "response" (usually 
		#species abundance) matrix **Y**. In essence, RDA is a multivariate extension of linear regression 
		#that seeks to explain as much variation in species data as possible through linear combinations 
		#of environmental predictors.

	#Remember from lecture, we can do this manually through the following steps:

		#**Step 1)** Conduct a multivariate linear regression of **Y** on **X** to produce a matrix 
			#of fitted values.

		#**Step 2)** Run a principal component analysis (PCA) on the matrix of fitted values to produce 
			#canonical eigenvalues, eigenvectors, and axes.

	#Luckily, the `vegan` package streamlines this process using the `rda` function.

	#Let's start with a simple RDA of fish densities without any transformation.

  		RDA <- rda(fish_dens ~ ., env_cont)
    			summary(RDA)

	#At a first glance of the summary table, we can see that the simple RDA performs very poorly for 
		#the complete dataset. This is to be expected given that species abundances are more likely to 
		#have unimodal or otherwise nonlinear relationships with environmental factors than linear 
		#relationships. This is a known limitation of RDA when species distributions respond to 
		#environmental gradients in a nonlinear fashion (e.g., in hump-shaped patterns along gradients 
		#like temperature or elevation).

	#Looking at the "Partitioning of variance" table, the total percentage of constrained inertia (i.e., 
		#variance accounted for by the environmental predictors) is only 2.7%. That's pretty dismal! 
		#Of this tiny slice of constrained inertia, RDA axis 1 explains about 91% and axis 2 explains 
		#about 8%. This indicates that most of the explained variance is concentrated on the first axis, 
		#but the overall fit of the model is poor, suggesting that the species-environment relationships 
		#are not well captured by this linear model.

	#The `summary` function also spits out the species scores, site scores, and biplot scores (environmental 
		#scores) for the RDA triplot. Note that the default scaling for RDA in the `vegan` package is 
		#Scaling 2, which focuses on preserving the correlations between species and environmental 
		#variables. Scaling 1, on the other hand, emphasizes the distances between sampling sites. 
		#Scaling 2 is often preferred when you are more interested in understanding species-environment 
		#relationships, whereas Scaling 1 may be useful when interpreting site-level differences.

	#We have a few more tricks up our sleeve for evaluating this RDA:

  		RsquareAdj(RDA)
  		anova(RDA, permutations = how(nperm = 999))

	#The `RsquareAdj` function provides a "**raw**" and an "**adjusted**" **R-squared** value, which 
		#corrects the raw R-squared value for the number of predictors, giving a more realistic estimate 
		#of the variance explained by the model. This adjusted value is often lower than the raw R-squared 
		#but gives a better idea of the model's true fit.

	#The `anova` function performs a permutation test to assess the statistical significance of the RDA model.
		#It tests the null hypothesis that the environmental variables (or predictors) have no effect on 
		#species abundances (or response variables). If the permutation test produces a low p-value, it 
		#indicates that the environmental predictors are significantly related to species variation.

	#We can also check the significance of specific axes and environmental variables using the `anova` 
		#function as follows:

  		anova(RDA, by = "axis", permutations = how(nperm = 999))
  		anova(RDA, by = "terms", permutations = how(nperm = 999))

	#So, while none of the constrained ordination axes are significant, the environmental variable "gradient" 
		#does appear to have a marginally significant effect.

	#Two more functions, `goodness` and `vif.cca`, can be helpful for evaluating the model.

	#The `goodness` function calculates the **goodness of fit** for each species, showing how well the 
		#species' variation is explained by the environmental predictors. A higher value across axes 1 
		#and 2 indicates a stronger relationship.

	#The `vif.cca` function computes the **variance inflation factors (VIFs)** for each environmental 
		#predictor, identifying multicollinearity issues. High VIF values indicate that some predictors 
		#are strongly correlated with each other, which can inflate standard errors and make it harder to 
		#detect significant relationships.

  		goodness(RDA)
  		vif.cca(RDA)

	#In this case, the goodness of fit is very poor for all species.

	#We can use the basic plot function to evaluate the results of the RDA in a triplot.

		par(mfrow=c(1,2))
		plot(RDA,
			scaling = 1,
			display = c("sp","lc","cn"),
			main = "Triplot of RDA, Scaling = 1")
		plot(RDA,
			scaling = 2,
			display = c("sp","lc","cn"),
			main = "Triplot of RDA, Scaling = 2")

	#When visualizing RDA results in triplots, Scaling 1 may not perform well because certain species can 
		#exert an outsized influence, disproportionately affecting the ordination results. These species 
		#may dominate the variation, overshadowing patterns in other species, leading to skewed ordinations.

	#As always, we can modify plots to meet our specific needs. The same tricks and tips we learned in Lab 8 
		#apply for RDA/CCA triplots too!

	#For example, points can be colored by population:

  		groups <- levels(factor(env$SMU))
		pt_col <- viridis(length(groups))

		site.sc <- summary(RDA)[2]
		spe.sc <- summary(RDA)[1]
		env.sc <- summary(RDA)[4]

		plot(RDA, 
			choices=c(1,2), 
			type="n", 
			scaling=2, 
			main="RDA of Fish Density",
			xlab="RDA 1", 
			ylab="RDA 2")	
		for (i in 1:length(groups)){
			dim_choice <- site.sc$sites[env$SMU==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
			pch=19, 
			cex=1.4, 
			col=pt_col[i])
			}
		text(env.sc$biplot*4, 
			row.names(env.sc$biplot), 
			col="darkgray")
		arrows(0, 0, env.sc$biplot[,1]*3.7, env.sc$biplot[,2]*3.7, 
			lwd=2, 
			length=0.1, 
			col="darkgray")
		text(spe.sc$species*10, 
			row.names(spe.sc$species))
		arrows(0, 0, spe.sc$species[,1]*9, spe.sc$species[,2]*9, 
			lwd=2, 
			length=0.1)
		legend(x="bottomright", 
			legend=groups, 
			col=pt_col, 
			pch=19)

	#Or habitat type:

		NLCD_Cat2 <- dat_final$NLCD_Cat
		NLCD_Cat2[NLCD_Cat2 == "Crops-Cultivated"] = "Pasture-Hay"
		NLCD_Cat2[NLCD_Cat2 == "Water"] = "Wetlands-Emergent"
		
  		groups <- levels(factor(NLCD_Cat2))
		pt_col <- viridis(length(groups))

		site.sc <- summary(RDA)[2]
		spe.sc <- summary(RDA)[1]
		env.sc <- summary(RDA)[4]

		plot(RDA, 
			choices=c(1,2), 
			type="n", 
			scaling=2, 
			main="RDA of Fish Density",
			xlab="RDA 1", 
			ylab="RDA 2")	
		for (i in 1:length(groups)){
			dim_choice <- site.sc$sites[NLCD_Cat2==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
			pch=19, 
			cex=1.4, 
			col=pt_col[i])
			}
		text(env.sc$biplot*4, 
			row.names(env.sc$biplot), 
			col="darkgray")
		arrows(0, 0, env.sc$biplot[,1]*3.7, env.sc$biplot[,2]*3.7, 
			lwd=2, 
			length=0.1, 
			col="darkgray")
		text(spe.sc$species*10, 
			row.names(spe.sc$species))
		arrows(0, 0, spe.sc$species[,1]*9, spe.sc$species[,2]*9, 
			lwd=2, 
			length=0.1)
		legend(x="bottomright", 
			legend=groups, 
			col=pt_col, 
			pch=19)

### Transformation-based RDA

	#The simple RDA performed poorly due to unmet assumptions of linearity. Transformation-based RDA 
		#(tbRDA) applies a Hellinger (or chord) transformation to species data, which down-weights rare 
		#species and makes the data more suitable for linear methods like RDA. You'll probably remember 
		#from our tbPCA analysis that the Hellinger transformation is particularly useful for community 
		#data that contain many zeros or have skewed distributions.

		fish.hel <- decostand(fish_dens, "hellinger")

		tbRDA <- rda(fish.hel ~ ., env_cont)
			summary(tbRDA)

		anova(tbRDA, permutations = how(nperm = 999))
		anova(tbRDA, by = "axis", permutations = how(nperm = 999))
		anova(tbRDA, by = "terms", permutations = how(nperm = 999))
	
		RsquareAdj(tbRDA)
		goodness(tbRDA)
		vif.cca(tbRDA)

	#Ok, not bad. After transforming the data, the proportion of constrained inertia increased to 16%, 
		#and the "full" model became significant. Most predictor variables now show significant 
		#relationships with species variation, suggesting that the transformation better captures 
		#species-environment relationships.

  		groups <- levels(factor(env$SMU))
		pt_col <- viridis(length(groups))

		site.sc <- summary(tbRDA)[2]
		spe.sc <- summary(tbRDA)[1]
		env.sc <- summary(tbRDA)[4]

		plot(tbRDA, 
			choices=c(1,2), 
			type="n", 
			scaling=2, 
			main="tbRDA of Hellinger-transformed Fish Density",
			xlab="RDA 1", 
			ylab="RDA 2")	
		for (i in 1:length(groups)){
			dim_choice <- site.sc$sites[env$SMU==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
			pch=19, 
			cex=1.4, 
			col=pt_col[i])
			}
		text(env.sc$biplot*2, 
			row.names(env.sc$biplot), 
			col="darkgray")
		arrows(0, 0, env.sc$biplot[,1]*1.7, env.sc$biplot[,2]*1.7, 
			lwd=2, 
			length=0.1, 
			col="darkgray")
		text(spe.sc$species*2, 
			row.names(spe.sc$species))
		arrows(0, 0, spe.sc$species[,1]*1.5, spe.sc$species[,2]*1.5, 
			lwd=2, 
			length=0.1)
		legend(x="bottomright", 
			legend=groups, 
			col=pt_col, 
			pch=19)

	#Even though the model performed much better, we can see there is a pretty substantial "horseshoe" 
		#effect going on; however, some species-environment trends are becoming apparent. For instance, 
		#we can see that there's a strong relationship between redband trout and high-gradient habitats 
		#and between dace and highly herbaceous (i.e., low canopy cover) habitats. The Warner sucker is 
		#driving an elongated y-axis because it only occurs at certain (deeper) sites in the Warner Basin.

### Stepwise Model Selection

	#In constrained ordination methods like RDA or CCA, not all environmental variables may be significant. 
		#Stepwise model selection allows us to refine the model by either adding (forward selection) or 
		#removing (backward selection) predictors based on their statistical significance. This helps to 
		#identify the most important variables while reducing model complexity.

	#The `vegan` package allows us to conduct a stepwise forward or backward model selection procedure as 
		#follows:

  		mod0 <- rda(fish.hel ~ 1, env_cont)
  		step.forward <- ordistep(mod0,
			scope = formula(tbRDA),
			direction = "forward",
			permutations = how(nperm=999))
		RsquareAdj(step.forward)

		step.backward <- ordistep(tbRDA, permutations = how(nperm=999))
		RsquareAdj(step.backward)

		anova(step.backward, permutations = how(nperm = 999))
		anova(step.backward, by = "terms", permutations = how(nperm = 999))

	#If we have a variable we are particularly interested in (for exmple, gradient), we can force the model 
		#to keep it:

		mod0p <- rda(fish.hel ~ Gradient, env_cont)
		mod1p <- rda(fish.hel ~ . + Gradient, env_cont)
		step.forward.p <- ordistep(mod0p,
			scope = formula(mod1p),
			direction = "forward",
			permutations = how(nperm=999))
		RsquareAdj(step.forward.p)
	 
  		anova(step.forward.p, permutations = how(nperm = 999))
		anova(step.forward.p, by = "terms", permutations = how(nperm = 999))

	#Be cautious. The final model may contain different combinations of variables depending on the stepwise 
		#procedure you use! This is where a detailed knowledge of your study system may come into play.

### Distance-Based RDA

	#Distance-based RDA (dbRDA) extends traditional RDA by allowing for nonlinear species-environment 
		#relationships. It uses a dissimilarity matrix (e.g., Bray-Curtis) to perform constrained 
		#ordination on distance measures, making it more suitable for community data where species 
		#distributions do not follow linear patterns.

	#Let's perform a distance-based RDA using a Bray-Curtis distance function as we have in the past for 
		#principal coordinates analysis (PCoA) and nonmetric multidimensional scaling (NMDS).

  		dbRDA <- dbrda(fish_dens_log ~ ., env_cont, "bray")

		anova(dbRDA, permutations = how(nperm = 500))
		anova(dbRDA, by = "axis", permutations = how(nperm = 500))
		anova(dbRDA, by = "terms", permutations = how(nperm = 500))

		RsquareAdj(dbRDA)
		vif.cca(dbRDA)	

	#The `dbrda` function in `vegan` doesn't compute goodness-of-fit estimates for species because it 
		#operates on distances rather than raw data. As a result, dbRDA outputs can be interpreted through 
		#the significance of environmental variables, but species scores are not directly available for 
		#triplots. (**If anyone can find a way to force these or calculate them manually, please include 
		#in your homework this week and share with the class!**)

	#We can get around this limitation for now by sizing points according to certain species' densities, 
		#for example, redband trout.

  		groups <- levels(factor(env$SMU))
		pt_col <- viridis(length(groups))

		site.sc <- summary(dbRDA)[1]
		env.sc <- summary(dbRDA)[3]

		plot(dbRDA, 
			choices=c(1,2), 
			type="n", 
			scaling=2, 
			main="dbRDA of Log-Fish Density (Bray-Curtis)",
			xlab="RDA 1", 
			ylab="RDA 2")	
		for (i in 1:length(groups)){
			dim_choice <- site.sc$sites[env$SMU==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
			pch=19, 
			cex=(1 + 2*fish_dens_log$TROUT_RB), 
			col=pt_col[i])
			}
		text(env.sc$biplot*5, 
			row.names(env.sc$biplot), 
			col="darkgray")
		arrows(0, 0, env.sc$biplot[,1]*4.7, env.sc$biplot[,2]*4.7, 
			lwd=2, 
			length=0.1, 
			col="darkgray")
		legend(x="bottomright", 
			legend=groups, 
			col=pt_col, 
			pch=19)

## Canonical Correspondence Analysis (CCA)

	#Canonical Correspondence Analysis (CCA) is a constrained ordination technique that assumes unimodal 
		#species responses to environmental gradients. Unlike RDA, which assumes linear relationships, 
		#CCA captures nonlinear relationships, making it more appropriate when species exhibit bell-shaped 
		#distributions along environmental gradients.

	#The calls in `vegan` are nearly identical as for RDA. (note that I've hashtagged the `summary` function 
		#for the sake of document length)

		CCA <- cca(fish_dens_log ~ ., env_cont)
			summary(CCA)

		anova(CCA, permutations = how(nperm = 999))
		anova(CCA, by = "axis", permutations = how(nperm = 999))
		anova(CCA, by = "terms", permutations = how(nperm = 999))

		RsquareAdj(CCA)
		goodness(CCA)
		vif.cca(CCA)

	#When applied to our data, CCA performs fairly well, though the proportion of constrained inertia is 9%, 
		#indicating that species-environment relationships are still not fully captured. As with RDA, we 
		#can use triplots to visualize these relationships, and axes can reveal specific environmental 
		#variables driving species distributions.

  		groups <- levels(factor(env$SMU))
		pt_col <- viridis(length(groups))

		site.sc <- summary(CCA)[2]
		spe.sc <- summary(CCA)[1]
		env.sc <- summary(CCA)[4]

		plot(CCA, 
			choices=c(1,2), 
			type="n", 
			scaling=2, 
			main="CCA of Log-Fish Density",
			xlab="CCA 1", 
			ylab="CCA 2")	
		for (i in 1:length(groups)){
			dim_choice <- site.sc$sites[env$SMU==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
			pch=19, 
			cex=1.4, 
			col=pt_col[i])
			}
		text(env.sc$biplot*10, 
			row.names(env.sc$biplot), 
			col="darkgray")
		arrows(0, 0, env.sc$biplot[,1]*9.7, env.sc$biplot[,2]*9.7, 
			lwd=2, 
			length=0.1, 
			col="darkgray")
		text(spe.sc$species*5, 
			row.names(spe.sc$species))
		arrows(0, 0, spe.sc$species[,1]*4.5, spe.sc$species[,2]*4.5, 
			lwd=2, 
			length=0.1)
		legend(x="bottomright", 
			legend=groups, 
			col=pt_col, 
			pch=19)

	#What's going on with the Y-axis. Which species/environmental variable is driving this trend? What 
		#phenomenon or characteristic of CCA is this related to?

### RDA or CCA?

	#There is a useful trick for determining whether your response-predictor relationships are linear (RDA) 
		#or unimodal/non-linear (CCA/tbRDA/dbRDA). Remember how I told you "just don't" use detrended 
		#correspondence analysis (DCA)? Well, here is a situation where it comes in handy. DCA estimates 
		#the length of the ordination axis in terms of standard deviations of species scores. If the 
		#gradient length exceeds 4, a unimodal model like CCA is more appropriate. If it is less than 2, 
		#a linear model like RDA is preferable. Intermediate values between 2 and 4 suggest that either 
		#model could be used, but CCA is more likely.

  		decorana(fish_dens_log, ira=0)

	#For these data, the axis length is \>4, so we know a CCA, tbRDA, or dbRDA is going to outperform a 
		#simple RDA. Ideally, this check is performed *prior to* your analysis.

## Co-inertia Analysis (CoIA)

	#Co-inertia analysis (CoIA) is a symmetric ordination method used to explore the relationships between 
		#two datasets without assuming a directional relationship. CoIA maximizes the covariance between 
		#the two datasets, making it a powerful exploratory tool for identifying co-structures between 
		#species and environmental variables. Usually in ecology, we're testing directional hypotheses, 
		#so a symmetric method such as CoIA is less appropriate; however, using CoIA from an exploratory 
		#standpoint can reveal some pretty cool trends. Let's try it!

		dudi.fish <- dudi.pca(fish.hel, scale=FALSE, scannf=FALSE)
		dudi.env <- dudi.pca(env_cont, scale=TRUE, scannf = FALSE)

		coia <- coinertia(dudi.fish, dudi.env, scannf = FALSE, nf=2)
			summary(coia)
		randtest(coia, nrepet = 999)

  		plot(coia)

	#**Total Inertia:** A total inertia of 0.15 reflects the overall strength of the covariance or shared 
		#structure between the two datasets. This value indicates how much of the variance in both datasets 
		#is explained by their common structure. The value is relatively small, suggesting that while there 
		#is a relationship between fish these environmental variables, it's not very strong.

	#**RV Coefficient**: The RV coefficient of 0.19 is a measure of the global correlation between the 
		#two datasets. This value is pretty low (on a scale of 0 to 1), suggesting that, although there 
		#is some relationship between fish and environmental drivers, it is not particularly strong.

	#The biplot, which shows individual sample scores, demonstrates how the samples relate to one another 
		#in terms of their co-inertia. Samples that are closer together are more similar in their joint 
		#configuration across the two datasets. The arrows connecting the samples represent the difference 
		#between the fish and environmental data for each site. Short arrows indicate a good alignment 
		#between the two datasets, whereas longer arrows suggest discrepancies between the fish and 
		#environmental patterns for a given sample.

	#In the canonical weights plots, which display the contributions of variables to the shared structure, 
		#the environmental variable "Max Depth" has the highest weight on the first axis, suggesting it is 
		#the most significant environmental factor influencing the fish dataset. Other variables, such as 
		#elevation ("Elev") and herbaceous cover ("Herbaceous"), contribute less to the co-inertia. In the 
		#fish species dataset, Warner sucker has the most substantial influence, followed by redband trout 
		#and dace. This implies that Warner sucker are particularly important in explaining the co-inertia 
		#between the fish and environmental data.




