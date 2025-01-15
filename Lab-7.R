###############################################################
###Lab-7.R									###
###CA, PCA, and NMDS							###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 22, 2024						###
###############################################################

## Lab Objectives

	#As we learned in the previous lecture and in Lab 6, PCA is not well-suited to study species abundance data 
		#and other types of non-linearly distributed data without substantial transformation. Other types of 
		#"unconstrained" ordination techniques such as Correspondence Analysis (CA), Principal Coordinates 
		#Analysis (PCoA), and Nonmetric Multidimensional Scaling (NMDS) may improve our ability to explore 
		#ecological patterns in the data in reduced dimensional space. In this lab we will learn to:

		#-   Understand and apply these methods to species abundance data
		#-   Use the broken stick model to assess the significance of principal axes
		#-   Use permutation to assess the significance of species/descriptor loadings
		#-   Visualize output using ordination biplots

## Setting up the R Workspace and Preparing the Data

	#We will load the necessary packages and prepare the dataset as we have done previously. Because we are 
		#only examining the species abundance data, we do not need to mess with the environmental portion 
		#of the dataset.

		library(vegan)
		library(viridis)

    setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
    source("Biostats.R")
    source("coldiss.R")
    
    dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)
    sub_dat <- subset(dat, SMU=="Malheur")
    
  		# source("C:\\USGS_OCRU\\Teaching\\FW599_Multivariate_Statistics\\Data\\Biostats.R")
  		# 
  		# dat <- read.csv("C:\\USGS_OCRU\\Teaching\\FW599_Multivariate_Statistics\\Data\\Harney_Fishes_2007.csv", row.names = 1)

  	
	#We will once again subset the data to include only the Malheur sites because the entire dataset is too large 
		#to effectively visualize for learning purposes.

  		spp_N <- colSums(sub_dat[,16:ncol(sub_dat)])
  		spp_0 <- subset(spp_N, spp_N == 0)
  		omit <- names(spp_0)

  		dat2 <- sub_dat[,!(colnames(sub_dat) %in% omit)]
	
  		dat3 <- dat2[rowSums(dat2[,16:ncol(dat2)]) >0, ]
  
  		dat3$Herbaceous[is.na(dat3$Herbaceous)] <- 0 
  		dat3$Ann_Herb[is.na(dat3$Ann_Herb)] <- 0
	
  		dat3 <- dat3[complete.cases(dat3$SiteLength),]
  
  		dat_final <- dat3
	
  		fish <- dat_final[,16:ncol(dat_final)]
	
  		fish_red <- drop.var(fish, min.fo=1)
  
  		fish_dens <- fish_red
  		for(i in 1:nrow(fish_red)){	
    			fish_dens[i,] <- fish_red[i,]/dat_final$SiteLength[i]
    			}

	#For the sake of consistency, we will run our analyses using the log-transformed fish density dataset. 
		#For CA, it's technically better to use the raw count data; however, as you'll see in a moment, 
		#the method performs poorly for this dataset regardless.

  		fish_dens_log <- log(fish_dens + 1)

## Correspondence Analysis (CA)

	#Correspondence Analysis (CA) is an ordination technique specifically designed for count data, such as 
		#species abundance data, which often have a unimodal distribution. CA is a form of eigenanalysis, 
		#similar to PCA, but it uses a chi-squared distance metric instead of Euclidean distance. CA is 
		#often used to identify patterns in ecological data by examining how different sites or species 
		#are distributed across multiple dimensions.

	#In `vegan`, the CA is run using the `cca` function, which is the same function we'll use to run Canonical 
		#Correspondence Analysis (CCA), a constrained ordination technique that we'll learn about later.

  		fish.ca <- cca(fish_dens_log)
			summary(fish.ca)

	#The summary output shows the eigenvalues, species scores, and site scores for the first few principal axes. 
		#These scores represent the position of each species and site along the principal axes, which are the 
		#dimensions that explain the most variation in the data. The species and site scores can be interpreted 
		#similarly to loadings in PCA.

	#The proportion of variance accounted for by each axis, as determined from the eigenvalues, can be calcuated 
		#as follows:

  		100*eigenvals(fish.ca)/sum(eigenvals(fish.ca))

	#We can use the broken stick model to compare the observed eigenvalues against a random distribution of 
		#eigenvalues. This model helps to determine whether the variance explained by each axis is significant 
		#or if it could have occurred by chance.

  		screeplot(fish.ca, bstick=TRUE, main="Inertia, CA vs. Random")

	#Similar to PCA, there are two scaling methods that can be used to visualize the data. Scaling 1 preserves 
		#the chi-squared distance among sites, making it more suitable for interpreting relationships between 
		#sites. Scaling 2, on the other hand, preserves the chi-squared distance among species, making it more 
		#appropriate for interpreting relationships between species. We can compare the output as such:

		par(mfrow=c(1,2))
		plot(fish.ca, 
			scaling=1, 
			main="CA of Log-Fish Density, Scaling 1")
		plot(fish.ca, 
			main="CA of Log-Fish Density, Scaling 2")

	#In this case, any differences in the two scaling methods are overshadowed by the brook trout outlier, 
		#Silvies-12. Here is a slightly clearer visualization using the `ordiplot` function.

		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))
		site.sc <- scores(fish.ca, choices=c(1,2))

		dev.off()
		p <- ordiplot(fish.ca, 
			type="n", 
			main="CA of Log-Fish Density",
			xlab="Axis 1", 
			ylab="Axis 2", 
			xlim=c(-8,2))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc$sites[dat_final$Pop==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
				}
		text(site.sc$species*1.2, row.names(site.sc$species))
		arrows(0, 0, site.sc$species[,1]*1.1, site.sc$species[,2]*1.1, 
			lwd=2, 
			length=0.1)
		legend(x="bottomleft", 
			legend=levels(factor(dat_final$Pop)), 
			col=pt_col[1:6], 
			pch=19, 
			cex=1.2)

## Principal Coordinates Analysis (PCoA)

	#Principal Coordinates Analysis (PCoA), also known as metric multidimensional scaling (MDS), is a 
		#method used to explore similarities or dissimilarities in data. Unlike PCA, which directly 
		#reduces dimensionality based on Euclidean distances, PCoA can use any distance metric, making 
		#it more flexible for ecological data that may not adhere to Euclidean assumptions. PCoA aims to 
		#represent data in a reduced number of dimensions while preserving the original distances as 
		#accurately as possible.

	#PCoA is run using the `cmdscale` function. This function performs classical multidimensional scaling 
		#on a distance matrix, which must be computed beforehand. The `k` parameter specifies the number 
		#of dimensions to retain, and `eig=TRUE` allows us to examine the eigenvalues, which represent the 
		#variance explained by each principal coordinate.

		fish.bray <- vegdist(fish_dens_log, "bray")
		fish.pcoa <- cmdscale(fish.bray, 
			k=5, 
			eig=TRUE, 
			add=T)
		fish.pcoa

	#We can examine the proportion of explained variance the same way we did for CA:

  		100*fish.pcoa$eig/sum(fish.pcoa$eig)

	#Look at the broken stick plot to see the percentage of variation explained by the first few principal 
		#coordinates compared to what would be expected by chance (note that the `screeplot` function doesn't       ##### how to interpret this plot
		#work well with PCoA):

		plot(100*fish.pcoa$eig/sum(fish.pcoa$eig),
			type="b",
			lwd=2,
			col="blue",
			xlab="Principal Component from PCoA",
			ylab="Percent variation explained",
			main="Broken Stick Model")
		lines(bstick(length(fish.pcoa$eig))*100, 
			type="b", 
			lwd=2, 
			col="red")

	#A fun aspect of PCoA is that, even though the descriptors (species) are not implicit in the calculation 
		#of principal coordinates, we can calculate their loadings as follows *and* use permuatation to ascribe   ##### how to interpret this code
		#significance values to them!

		spe.sc <- wascores(fish.pcoa$points[,1:2], fish_dens_log)

		vec.sp <- envfit(as.data.frame(fish.pcoa$points), fish_dens_log, perm=1000)

	#As would be expected, redband trout and dace have "significant" loadings on the principal coordinates, and 
		#redside shiner is marginally significant.

		ordiplot(scores(fish.pcoa, choices=c(1,2)), 
			type="t", 
			main="PCoA of Log-Fish Density")
		abline(h=0, lty=3)
		abline(v=0, lty=3)
		points(scores(fish.pcoa, choices=c(1,2)), 
			pch=19)
		plot(vec.sp, 
			p.max=0.1, 
			col="blue")

	#And our fun, site-colored plot:

		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))
		site.sc <- scores(fish.pcoa, choices=c(1,2))

		p <- ordiplot(fish.pcoa, 
			type="n", 
			main="PCoA of Log-Fish Density",
			xlab="Axis 1", 
			ylab="Axis 2", 
			xlim=c(-0.8,1))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[dat_final$Pop==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
			}
		text(spe.sc*1.5, row.names(spe.sc))
		arrows(0, 0, spe.sc[,1]*1.4, spe.sc[,2]*1.4, 
			lwd=2, 
			length=0.1)
		legend(x="bottomleft", 
			legend=levels(factor(dat_final$Pop)), 
			col=pt_col[1:6], 
			pch=19, 
			cex=1.2)

	#What do we notice here? Still no trends by population, but sites do seem to partition out well among 
		#redband trout-dominant, dace-dominant, and shiner-dominant species profiles.

## Nonmentric Multidimensional Scaling (NMDS)

	#In this last section, we'll cover Nonmetric Multidimensional Scaling (NMDS). Like PCoA, NMDS aims 
		#to represent the dissimilarities between data points in a reduced number of dimensions, but it 
		#does so non-parametrically. NMDS ranks the dissimilarities between data points and attempts to find 
		#a configuration in reduced dimensional space that best preserves the rank order of these 
		#dissimilarities. NMDS is particularly useful for ecological data because it does not rely on 
		#assumptions about the distribution of the data.

	#NMDS is often applied using the Bray-Curtis dissimilarity metric, which is sensitive to species 
		#composition. Unlike PCA or PCoA, NMDS does not produce eigenvalues, so the solution must be 
		#evaluated based on the stress value, which indicates how well the ordination preserves the 
		#original dissimilarities.

		fish.nmds2 <- metaMDS(fish_dens_log, 
			distance="bray", 
			k=2, 
			autotransform=FALSE, 
			trymax=100)
		fish.nmds2

		fish.nmds3 <- metaMDS(fish_dens_log, 
			distance="bray", 
			k=3, 
			autotransform=FALSE, 
			trymax=100)
		fish.nmds3

	#The stressplot shows the relationship between the observed dissimilarities and the dissimilarities 
		#in the ordination space. A stress value below 0.2 indicates acceptable fit.

	#As you can see below, adding a third dimension substantially reduces stress. The stress for the 
		#two-dimensional solution is about 0.15, which is technically "acceptable" but not ideal. The 
		#three-dimensional solution has a stress value of about 0.08, which is quite good. We can examine 
		#how stress declines with the addition of more dimensions using the `nmds.scree` function.

		nmds.scree(fish_dens_log, 
			distance="bray",
			k=10, 
			autotransform=FALSE,
			trymax=50)

		stressplot(fish.nmds2, main="Shepard Plot")
	
		stressplot(fish.nmds3, main="Shepard Plot")

	#The Shepard plot is a diagnostic tool that shows how well the distances in the ordination space correspond 
		#to the original dissimilarities. Ideally, the points should lie close to a straight line, indicating 
		#a good fit.

	#For the sake of simplicity, we'll continue with the two-dimensional solution. Similar to PCoA, we can 
		#calculate species loadings and assess their significance using permutation tests.

		spe.sc <- wascores(as.data.frame(fish.nmds2$points[,1:2]), fish_dens_log)

		vec.sp <- envfit(as.data.frame(fish.nmds2$points), fish_dens_log, perm=1000)

	#The `envfit` function fits environmental vectors or factors onto an ordination and provides significance 
		#values for the loadings. Significant species loadings suggest that these species contribute 
		#substantially to the observed patterns. In this case, redband trout, dace, and redside shiner are 
		#significant, meaning they are important for structuring the community composition in the ordination space.

	#Next, we visualize the NMDS results. The ordination plot displays the sites in the reduced dimensional space, 
		#with species vectors indicating the direction and strength of species associations.

		plot(fish.nmds2, 
			choices=c(1,2),
			type = "text", 
			display = "sites",
			main="NMDS of Log-Fish Density")
		plot(vec.sp, 
			p.max=0.1, 
			col="blue")

	#In this plot, the arrows represent the direction and strength of species' influence on the ordination. Longer 
		#arrows indicate stronger associations with the NMDS axes.

	#We'll play around with visualization more in the next lab, but the ordination biplot can be scaled by 
		#goodness of fit, which is a measure of how well each site is represented in the reduced dimensional space:

		gof <- goodness(fish.nmds2)
		plot(fish.nmds2, 
			choices=c(1,2),
			type = "text", 
			display = "sites",
			main="NMDS of Log-Fish Density")
		points(fish.nmds2, 
			display="sites", 
			cex=gof*200)

	#In this visualization, larger points indicate better representation (lower stress) for those sites in 
		#the ordination space.

	#We can also highlight specific species abundances in the NMDS plot. For instance, let's focus on redband trout 
		#by scaling the points based on their log-transformed density.

		plot(fish.nmds2, 
			choices=c(1,2),
			type = "text", 
			display = "sites",
			main="NMDS of Log-Fish Density")	
		points(fish.nmds2, 
			cex=10*fish_dens_log$TROUT_RB)

	#This plot emphasizes the abundance of redband trout across sites, allowing us to see which sites are 
		#most strongly associated with this species.

	#How does the site-coded NMDS biplot compare to the PCoA plot from earlier?

		groups <- levels(factor(dat_final$Pop))
		pt_col <- viridis(length(groups))
		site.sc <- scores(fish.nmds2, choices=c(1,2))

		p <- ordiplot(fish.nmds2, 
			type="n", 
			main="NMDS of Log-Fish Density",
			xlab="Axis 1", 
			ylab="Axis 2")
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc$sites[dat_final$Pop==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
				}
		text(site.sc$species, row.names(site.sc$species))
		arrows(0, 0, site.sc$species[,1]*0.9, site.sc$species[,2]*0.9, 
			lwd=2, 
			length=0.1)
		legend(x="topright", 
			legend=levels(factor(dat_final$Pop)), 
			col=pt_col[1:6], 
			pch=19, 
			cex=1.2)	

