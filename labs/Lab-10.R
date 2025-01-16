###############################################################
###Lab-10.R									###
###Statistical Inference - Part 1					###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: October 31, 2024						###
###############################################################

## Lab Objectives

	#In the last module, we learned about constrained ordination—a **direct comparison** (or **direct 
		#gradient analysis**) technique used to explore and quantify the relationships between a set 
		#of response variables and a set of predictor variables. In the next two modules, we'll learn 
		#some **indirect comparison** techniques (with the exception of linear discriminant analysis) 
		#for making these sorts of inferences. These methods are crucial in ecological studies, as they 
		#can be used to identify environmental factors influencing species distributions and community 
		#structures.

	#In this lab, we will learn about:

		#-   **Principal Component Regression (PCR):** Used to reduce data dimensionality, where predictor 
			#variables are highly collinear.
		#-   **Linear Discriminant Analysis (LDA):** A classification technique that identifies linear 
			#combinations of features to separate two or more classes, such as habitats or species 
			#assemblages.
		#-   **Multi-response Permutation Procedure (MRPP):** A non-parametric method that compares 
			#differences between predefined groups, assessing whether there is a significant difference 
			#in species assemblages among groups like clusters or habitat types.

	#These methods are essential for understanding response-predictor relationships for complex, multivariate 
			#data.

## Setting up the R Workspace and Preparing the Data

	#Load the necessary packages and prepare the dataset as we have done previously.

		library(vegan)
		library(viridis)
		library(MASSExtra)
		library(ade4)

    setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
    source("Biostats.R")

    dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)

	#Once again, we'll use the **full** dataset to capture the entire suite of environmental conditions 
		#and species.

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

## Principal Component Regression

	#**Principal Component Regression (PCR)** is used when predictor variables (environmental data in 
		#this case) are collinear, as it reduces the dimensionality by using Principal Component Analysis 
		#(PCA) prior to regression. If you remember, PCA creates new, *uncorrelated* variables (principal 
		#components or axes) based on the original data in reduced dimensional space.

	#The first step is to run a PCA of the data. In this example, the PCA is run on the standardized 
		#environmental variables and loadings are examined to identify key environmental drivers of the 
		#first principal component/axis, which accounts for the greatest proportion of variance in the 
		#data (about 45%).

		env.pca <- prcomp(env_cont, scale=TRUE)	
			summary(env.pca)
		
		site.sc <- scores(env.pca)
  		p <- ordiplot(env.pca, 
    			type="n", 
    			main="PCA of Standardized Environmental Data")
		points(site.sc[,1], site.sc[,2], 
		      pch=19, 
		      cex=1.4, 
		      col="gray")
		arrows(0, 0, env.pca$rotation[,1]*4, env.pca$rotation[,2]*4, 
		       lwd=2, 
		       length=0.1)
		text(env.pca$rotation[,1]*4.5, env.pca$rotation[,2]*4.5, row.names(env.pca$rotation))

	#Yes, you can use this method with PCoA, CA, or NMDS if those ordination approaches are more 
		#appropriate for your data!

	#As another reminder, the visual output shows sites plotted along the first two principal components, 
		#and environmental variables are represented as vectors indicating how strongly they load onto 
		#these axes. As before, we can see that elevation, gradient, and canopy all seem to have strong 
		#loadings on PC1.

	#We can verify this using the following bit of code, which performs a *permutational* **environmental 
		#vector fitting analysis**, showing which environmental variables are significantly correlated 
		#with the ordination axes.

		env_cont_std <- decostand(env_cont, method="standardize")

		vec.env <- envfit(as.data.frame(scores(env.pca)), env_cont_std, perm=1000)
			vec.env

	#The next step of PCR is to regress the response variable (in this case, log-fish density) against the 
		#PC(s) that account for the most variation.

	#Initial results indicate that trout densities are negatively correlated with PC1, meaning that as the 
		#values along the first principal axis increase, trout density decreases. We can also interpret 
		#this as redband trout being associated with high gradient, high elevation, high canopy sites.

		plot(site.sc[,1], fish_dens_log$TROUT_RB, 
			pch=20, 
			cex=2,
			xlab="PC1",
			ylab="Log-redband trout density")
			abline(lm(fish_dens_log$TROUT_RB ~ site.sc[,1]))
			summary(lm(fish_dens_log$TROUT_RB ~ site.sc[,1]))

		plot(site.sc[,2], fish_dens_log$TROUT_RB, 
			pch=20, 
			cex=2,
			xlab="PC1",
			ylab="Log-redband trout density")
			abline(lm(fish_dens_log$TROUT_RB ~ site.sc[,2]))
			summary(lm(fish_dens_log$TROUT_RB ~ site.sc[,2]))

	#The "ordisurf" plot backs this up.

		p <- ordiplot(env.pca, 
    			type="n", 
   			 main="PCA of Standardized Environmental Data")
		points(site.sc[,1], site.sc[,2], 
		      pch=19, 
		      cex=(1+fish_dens_log$TROUT_RB), 
		      col="black")
		arrows(0, 0, env.pca$rotation[,1]*4, env.pca$rotation[,2]*4, 
		       lwd=2, 
		       length=0.1)
		text(env.pca$rotation[,1]*4.5, env.pca$rotation[,2]*4.5, row.names(env.pca$rotation))
		ordisurf(env.pca ~ TROUT_RB, fish_dens_log, col="blue", add=TRUE, select=FALSE)

	#AND we can verify the significance of this relationship for all fish species using the same 
		#environmental vector fitting analysis as above!

		vec.sp <- envfit(as.data.frame(scores(env.pca)), fish_dens_log, perm=1000)
			vec.sp

	#The relationship for redband trout is significant at the P \< 0.0001 level.

## Linear Discriminant Analysis

	#**Linear Discriminant Analysis (LDA)** is used for classifying objects (in this case, sites) into 
		#groups based on descriptors (in this case, environmental variables). The process involves 
		#finding linear combinations of descriptors/variables that best separate predefined groups 
		#(e.g., habitat types or clusters). Because LDA assumes linear relationships between response 
		#and predictor variables, it's generally not ideal for species data.

	#LDA by cluster:

	#Hierarchical clustering (Ward's method) is used to group sites based on environmental similarity. 
		#The clusters will then be used as the categorical variables in the LDA.

		env.euc <- vegdist(env_cont_std, "euclidean")
	
		envcl.ward <- hclust(env.euc, method = "ward.D2")

		plot(envcl.ward, main="Ward's Minimum Variance Dendrogram", 
			xlab="Sites", 
			ylab="Euclidean Distance", 
			hang=-1)
			rect.hclust(envcl.ward, k=6)

		env_clus <- cutree(envcl.ward, k=6)

	#It is a good idea to check that the data satisfy the assumptions of **multivariate homogeneity of 
		#variance** (using the `betadisper` function) and **distinct group means** (using **Wilks' 
		#lambda**).

  	MHV <- betadisper(env.euc, as.factor(env_clus))
		permutest(MHV)
		
  	wgc <- manova(as.matrix(env_cont_std) ~ as.factor(env_clus))
		summary(wgc, test = "Wilks")

	#In this case, both of these assumptions are met. Note that LDA is somewhat resilient to the violation 
		#of these assumptions.

	#Next, Perform LDA to find linear discriminants that best separate the clusters. The effectiveness of 
		#LDA is assessed through a confusion matrix, showing the percentage of correct classifications 
		#(in this case 91% accuracy, which is pretty good).

		clus.lda <- lda(env_cont_std, env_clus)
			clus.lda

		clus.lda.pred <- predict(clus.lda)
			clus.scores <- clus.lda.pred$x
			clus.df <- as.data.frame(cbind(env_clus, clus.scores))

		conf.matrix <- table(env_clus, clus.lda.pred$class)
			conf.matrix
		sum(diag(conf.matrix))/sum(conf.matrix)	

	#An ordination biplot is used to visualize the clusters and the environmental variables that contribute 
		#most to this separation. What does the LDA output look like?

		groups <- levels(factor(env_clus))
		pt_col <- viridis(length(groups))
		site.sc <- clus.scores
		spe.sc <- clus.lda[4]

		p <- ordiplot(clus.lda, 
    			type="n", 
    			xlim=c(-6,9), 
    			ylim=c(-3,5), 
    			main="LDA of Environmental Data",
			xlab="LD1", 
			ylab="LD2", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25)
		for (i in 1:length(groups)){
			dim_choice <- clus.scores[env_clus==i,]
			points(dim_choice[,1], dim_choice[,2], 
		      pch=19, 
		      cex=1.4, 
		      col=pt_col[i])
			}
		text(spe.sc$scaling*4.5, row.names(spe.sc$scaling))
		arrows(0, 0, spe.sc$scaling[,1]*3.5, spe.sc$scaling[,2]*3.8, lwd=2, length=0.1)
		ordiellipse(clus.lda, 
    			groups=env_clus, 
    			kind="sd", 
    			lwd=2, 
    			lty=2, 
    			col=pt_col)

	#We can do the same analysis for pre-defined groups, like habitat type (note that I clumped some groups 
		#together into broader categories).

		NLCD_Cat2 <- env$NLCD_Cat
		NLCD_Cat2[NLCD_Cat2 == "Crops-Cultivated"] = "Agriculture"
		NLCD_Cat2[NLCD_Cat2 == "Pasture-Hay"] = "Agriculture"
		NLCD_Cat2[NLCD_Cat2 == "Water"] = "Wetlands"
		NLCD_Cat2[NLCD_Cat2 == "Wetlands-Emergent"] = "Wetlands"
		NLCD_Cat2[NLCD_Cat2 == "Wetlands-Woody"] = "Wetlands"
		table(NLCD_Cat2)

	#The habitat groupings also meet the assumptions of multivariate homogeneity of variance and distinct 
		#group means.

		MHV <- betadisper(env.euc, as.factor(NLCD_Cat2))
		permutest(MHV)

		wgc <- manova(as.matrix(env_cont_std) ~ as.factor(NLCD_Cat2))
		summary(wgc, test = "Wilks")

		hab.lda <- lda(env_cont_std, NLCD_Cat2)
		hab.lda

		hab.lda.pred <- predict(hab.lda)
		hab.scores <- hab.lda.pred$x
		hab.df <- as.data.frame(cbind(NLCD_Cat2, hab.scores))
		plot(factor(hab.df$NLCD_Cat2), as.numeric(hab.df$LD1))

		conf.matrix <- table(NLCD_Cat2, hab.lda.pred$class)
			conf.matrix
		sum(diag(conf.matrix))/sum(conf.matrix)

	#The classification accuracy of this model is 81%, which is not awful, but also not great.

		groups <- levels(factor(NLCD_Cat2))
		pt_col <- viridis(length(groups))
		site.sc <- hab.scores
		spe.sc <- hab.lda[4]

		p <- ordiplot(hab.lda, 
    			type="n", 
    			xlim=c(-5,5), 
    			ylim=c(-4,4), 
    			main="LDA of Environmental Data",
			xlab="LD1", 
			ylab="LD2", 
			cex.axis=1.15, 
			cex.lab=1.25, 
			cex.main=1.25)
		for (i in 1:length(groups)){
			dim_choice <- hab.scores[NLCD_Cat2==groups[i],]
			points(dim_choice[,1], dim_choice[,2], 
		      pch=19, 
		      cex=1.4, 
		      col=pt_col[i])
			}
		text(spe.sc$scaling*3, row.names(spe.sc$scaling))
		arrows(0, 0, spe.sc$scaling[,1]*2.5, spe.sc$scaling[,2]*2.8, lwd=2, length=0.1)
		ordiellipse(hab.lda, 
    			groups=NLCD_Cat2, 
    			kind="sd", 
    			lwd=2, 
    			lty=2, 
    			col=pt_col)

	#We can still see a lot of overlap among groups in the LDA ordination output.

## Multi-response Permutation Procedure

	#**Multi-response permutation procedure (MRPP)** is a nonparametric (permutational) method used to test 
		#whether there are significant differences among groups, in this case, environmental clusters and 
		#habitat types. The method is similar to an analysis of variance (ANOVA), but it does not assume 
		#normality. It can be run in R using just a few lines of code.

	#We already know from our previous analyses that environmental characteristics differ among clusters and 
		#habitats. The MRPP can explicitly test for significance (but does not produce ordination output 
		#in the same way LDA does). Here we can see that the significance for both effects is *P* \< 0.001, 
		#though the effect size (A; on a scale of 0-1) for the cluster analysis is modest (0.38) and for 
		#the habitat analysis is small (0.14).

		mrpp.clus <- mrpp(env_cont_std, factor(env_clus), distance="euclidean")
			mrpp.clus 

		mrpp.hab <- mrpp(env_cont_std, NLCD_Cat2, distance="euclidean")
			mrpp.hab

	#What about species assemblages? Do they differ significantly among environmental clusters?

		mrpp.fish <- mrpp(fish_dens_log, factor(env_clus), distance="bray")
			mrpp.fish

	#Yes! But the effects size (A = 0.07) is very small. We would expect nearly this amount of group 
		#heterogeneity by chance.

	#MRPP **does not** assume equal group sizes. It is a non-parametric test that compares differences 
		#between predefined groups based on their distances, and it can handle unequal group sizes 
		#effectively. The test focuses on the within-group distances and compares them to the between-group 
		#distances, without requiring groups to be the same size.

	#Next week we will learn about more (mostly non-parametric) methods for making ecological inferences 
		#and comparing multivariate data among groups!



