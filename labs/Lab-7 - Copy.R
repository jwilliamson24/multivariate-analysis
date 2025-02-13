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

#1) Loading packages and source code
    
    library(vegan)
    library(viridis)
    setwd("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data")
    source("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis/Biostats.R")
    
#2) Loading and subsetting site level data
    
    dat2 <- readRDS("~/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/oss-occu/data/site_level_matrix.rds")
    row.names(dat2) <- dat2[,1]
    
    dat2 <- subset(dat2, year=="2024")
    
    sals <- dat2[,c("oss","enes")]
    
    drop <- c("lat","long","landowner","stand","tree_farm","year","weather","size_cl",
              "length_cl","site_id","oss","enes")
    env2 <- dat2[,!(colnames(dat2) %in% drop)]
    
    env_cont <- env2[, !colnames(env2) %in% "trt"]
    
    drop <- c("jul_date","veg_cov","fwd_cov","dwd_count","size_cl","decay_cl","char_cl","length_cl" )
    env_subset <- env_cont[,!(colnames(env_cont) %in% drop)]
    
    
#3) Standardizing salamanders by sampling area
    
    sal_dens <- sals
    for(i in 1:nrow(sals)){
      sal_dens[i,] <- sals[i,]/567
    }
    
    
#4) Transforming and standardizing the data as needed
    
    log_sal_cou <- log(sals + 1)
    log_sal_dens <- log(sal_dens + 1)
    
    env_std <- decostand(env_cont, "standardize") #Z-scores the data in each column
    

## Correspondence Analysis (CA)

	#Correspondence Analysis (CA) is an ordination technique specifically designed for count data, such as 
		#species abundance data, which often have a unimodal distribution. CA is a form of eigenanalysis, 
		#similar to PCA, but it uses a chi-squared distance metric instead of Euclidean distance. CA is 
		#often used to identify patterns in ecological data by examining how different sites or species 
		#are distributed across multiple dimensions.

	#In `vegan`, the CA is run using the `cca` function, which is the same function we'll use to run Canonical 
		#Correspondence Analysis (CCA), a constrained ordination technique that we'll learn about later.

    
# doesnt work, Error in cca.default(sals) : 
#    all row sums must be >0 in the community data matrix 
  		sal.ca <- cca(sals)
  		
  		env.ca <- cca(env_cont)
			summary(env.ca)

	#The summary output shows the eigenvalues, species scores, and site scores for the first few principal axes. 
		#These scores represent the position of each species and site along the principal axes, which are the 
		#dimensions that explain the most variation in the data. The species and site scores can be interpreted 
		#similarly to loadings in PCA.

	#The proportion of variance accounted for by each axis, as determined from the eigenvalues, can be calcuated 
		#as follows:

  		100*eigenvals(env.ca)/sum(eigenvals(env.ca))

	#We can use the broken stick model to compare the observed eigenvalues against a random distribution of 
		#eigenvalues. This model helps to determine whether the variance explained by each axis is significant 
		#or if it could have occurred by chance.

  		screeplot(env.ca, bstick=TRUE, main="Inertia, CA vs. Random")

	#Similar to PCA, there are two scaling methods that can be used to visualize the data. Scaling 1 preserves 
		#the chi-squared distance among sites, making it more suitable for interpreting relationships between 
		#sites. Scaling 2, on the other hand, preserves the chi-squared distance among species, making it more 
		#appropriate for interpreting relationships between species. We can compare the output as such:

		par(mfrow=c(1,2))
		plot(env.ca, 
			scaling=1, 
			main="Env CA, Scaling 1")
		plot(env.ca, 
			main="Env CA, Scaling 2")

	#In this case, any differences in the two scaling methods are overshadowed by the brook trout outlier, 
		#Silvies-12. Here is a slightly clearer visualization using the `ordiplot` function.

		groups <- levels(factor(dat2$trt))
		pt_col <- viridis(length(groups))
		site.sc <- scores(env.ca, choices=c(1,2))

		dev.off()
		p <- ordiplot(env.ca, 
			type="n", 
			main="Env CA",
			xlab="Axis 1", 
			ylab="Axis 2", 
			xlim=c(-8,2))
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc$sites[dat2$trt==groups[i],]
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
			legend=levels(factor(dat2$trt)), 
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

		
#sal data doesnt work, Error in cmdscale(sal.bray, k = 5, eig = TRUE, add = T) : 
#NA values not allowed in 'd'
		sal.bray <- vegdist(log_sal_dens, "bray")
		sal.pcoa <- cmdscale(sal.bray, 
			k=5, 
			eig=TRUE, 
			add=T)
	
		
#env data		
		env_std_subset <- env_std[,c("temp","dwd_cov","soil_moist","stumps","logs","decay_cl","canopy_cov")]	
		env_euc <- vegdist(env_std_subset, method="euclidean")
		
		env.pcoa <- cmdscale(env_euc, 
		                     k=5, 
		                     eig=TRUE, 
		                     add=T)
		env.pcoa

	#We can examine the proportion of explained variance the same way we did for CA:

  		100*env.pcoa$eig/sum(env.pcoa$eig)

	#Look at the broken stick plot to see the percentage of variation explained by the first few principal 
		#coordinates compared to what would be expected by chance (note that the `screeplot` function doesn't 
		#work well with PCoA):

		plot(100*env.pcoa$eig/sum(env.pcoa$eig),
			type="b",
			lwd=2,
			col="blue",
			xlab="Principal Component from PCoA",
			ylab="Percent variation explained",
			main="Broken Stick Model")
		lines(bstick(length(env.pcoa$eig))*100, 
			type="b", 
			lwd=2, 
			col="red")

	#A fun aspect of PCoA is that, even though the descriptors (species) are not implicit in the calculation 
		#of principal coordinates, we can calculate their loadings as follows *and* use permuatation to ascribe 
		#significance values to them!

		spe.sc <- wascores(env.pcoa$points[,1:2], env_subset)

		vec.sp <- envfit(as.data.frame(env.pcoa$points), env_subset, perm=1000)


	#And our fun, site-colored plot:

		groups <- levels(factor(dat2$trt))
		pt_col <- viridis(length(groups))
		site.sc <- scores(env.pcoa, choices=c(1,2))

		plot(site.sc[, 1:2],  # First two dimensions
		     main = "Env PCoA", 
		     xlab = "PCoA 1", 
		     ylab = "PCoA 2", 
		     pch = 19)

		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[dat2$trt==groups[i],]
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
			legend=levels(factor(dat2$trt)), 
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

		nmds2 <- metaMDS(env_std_subset, 
			distance="euclidean", 
			k=2, 
			autotransform=FALSE, 
			trymax=100)
		nmds2

		nmds3 <- metaMDS(env_std_subset, 
			distance="euclidean", 
			k=3, 
			autotransform=FALSE, 
			trymax=100)
		nmds3

	#The stressplot shows the relationship between the observed dissimilarities and the dissimilarities 
		#in the ordination space. A stress value below 0.2 indicates acceptable fit.

	#As you can see below, adding a third dimension substantially reduces stress. The stress for the 
		#two-dimensional solution is about 0.15, which is technically "acceptable" but not ideal. The 
		#three-dimensional solution has a stress value of about 0.08, which is quite good. We can examine 
		#how stress declines with the addition of more dimensions using the `nmds.scree` function.

		nmds.scree(env_std_subset, 
			distance="euclidean",
			k=10, 
			autotransform=FALSE,
			trymax=50)

		stressplot(nmds2, main="Shepard Plot")
	
		stressplot(nmds3, main="Shepard Plot")

	#The Shepard plot is a diagnostic tool that shows how well the distances in the ordination space correspond 
		#to the original dissimilarities. Ideally, the points should lie close to a straight line, indicating 
		#a good fit.

	#For the sake of simplicity, we'll continue with the two-dimensional solution. Similar to PCoA, we can 
		#calculate species loadings and assess their significance using permutation tests.

		spe.sc <- wascores(as.data.frame(nmds2$points[,1:2]), env_std_subset)

		vec.sp <- envfit(as.data.frame(nmds2$points), env_std_subset, perm=1000)

	#The `envfit` function fits environmental vectors or factors onto an ordination and provides significance 
		#values for the loadings. Significant species loadings suggest that these species contribute 
		#substantially to the observed patterns. In this case, redband trout, dace, and redside shiner are 
		#significant, meaning they are important for structuring the community composition in the ordination space.

	#Next, we visualize the NMDS results. The ordination plot displays the sites in the reduced dimensional space, 
		#with species vectors indicating the direction and strength of species associations.

		plot(nmds3, 
			choices=c(1,2),
			type = "text", 
			display = "sites",
			main="Env NMDS")
		plot(vec.sp, 
			p.max=0.1, 
			col="blue")

	#In this plot, the arrows represent the direction and strength of species' influence on the ordination. Longer 
		#arrows indicate stronger associations with the NMDS axes.

	#We'll play around with visualization more in the next lab, but the ordination biplot can be scaled by 
		#goodness of fit, which is a measure of how well each site is represented in the reduced dimensional space:

		gof <- goodness(nmds2)
		plot(nmds2, 
			choices=c(1,2),
			type = "text", 
			display = "sites",
			main="NMDS2")
		points(nmds2, 
			display="sites", 
			cex=gof*200)

	#In this visualization, larger points indicate better representation (lower stress) for those sites in 
		#the ordination space.


	#How does the site-coded NMDS biplot compare to the PCoA plot from earlier?

		groups <- levels(factor(dat2$trt))
		pt_col <- viridis(length(groups))
		site.sc <- scores(nmds2, choices=c(1,2))

		env_fit <- envfit(nmds2, env_std_subset, permutations = 999)  # env_std_subset is your z-scored environmental data
		
		p <- ordiplot(nmds2, 
			type="n", 
			main="NMDS2",
			xlab="Axis 1", 
			ylab="Axis 2")
		for (i in 1:length(groups))
			{
			dim_choice <- site.sc[dat2$trt==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
				pch=19, 
				cex=1.4, 
				col=pt_col[i])
		}
		plot(env_fit, p.max = 0.05, col = "red", lwd = 1.5)  # Use p.max to only show significant arrows, if desired
		legend(x="topright", 
			legend=levels(factor(dat2$trt)), 
			col=pt_col[1:6], 
			pch=19, 
			cex=1.2)	

