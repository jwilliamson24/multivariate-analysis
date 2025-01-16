###############################################################
###Lab-11.R									###
###Statistical Inference - Part 2					###
###FW 599 Multivariate Analysis of Ecological Data in R	###
###Author: Melanie Davis						###
###Date: November 5, 2024						###
###############################################################

## Lab Objectives

	#This lab module continues our exploration of **indirect gradient analysis** techniques, which allow us to 
		#make statistical inferences about the relationships between a set of response variables (such as 
		#species abundances or environmental measures) and a set of predictor variables (such as environmental 
		#gradients, experimental treatments, or spatial factors).

	#In this lab, we will learn about:

		#-   **Analysis of Similarities (ANOSIM):** A non-parametric statistical test used to determine whether 
			#there is a significant difference between two or more groups based on ranked dissimilarities or 
			#distances between samples.
		#-   **(Permutational) Multivariate Analysis of Variance (PERMANOVA):** A robust statistical technique 
			#used to compare multivariate datasets across different groups. It extends traditional ANOVA by 
			#assessing the significance of differences between group centroids in multivariate space, based 
			#on a distance matrix.
		#-   **Similarity Percentages (SIMPER):** A post-hoc method used to quantify which species or variables 
			#contribute most to the dissimilarities between groups.
		#-   **Mantel Test:** A non-parametric test used to assess the correlation between two distance matrices, 
			#such as comparing spatial distances and environmental distances or genetic differences and 
			#ecological distances.
		#-   **Procrustes Analysis:** A technique that evaluates the similarity between two datasets by aligning 
			#them in multidimensional space.

	#These methods are essential for understanding response-predictor relationships for complex, multivariate data.

## Setting up the R Workspace and Preparing the Data

	#Load the necessary packages and prepare the dataset as we have done previously.

		library(vegan)
		library(viridis)
		library(MASSExtra)
		library(ade4)
    
    setwd("C:/Users/jasmi/OneDrive/Documents/Academic/OSU/Git/multivariate-analysis")
    source("Biostats.R")
    
    dat <- read.csv("Harney_Fishes_2007.csv", row.names = 1)

	#Once again, we'll use the **full** dataset to capture the entire suite of environmental conditions and species.

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

## Analysis of Similarity

	#**Analysis of Similarity (ANOSIM)** is a non-parametric test used to evaluate whether there are significant 
		#differences in the composition of communities or multivariate datasets among predefined groups. It 
		#operates by comparing the ranked dissimilarities between and within groups, using a permutation-based 
		#approach to assess significance.

	#Because our data are already organized by Species Management Unit (SMU) or watershed, we can ask the question: 
		#*does fish community structure differ by SMU?*

	#The code to run an ANOSIM is very straightforward:

		SMU_anosim <- anosim(fish_dens_log, dat_final$SMU, permutations = 999, distance="bray")
			summary(SMU_anosim)

	#The main output of an ANOSIM analysis is the R statistic, which ranges from -1 to 1:

		#-   R ≈ 1 indicates that dissimilarities between groups are greater than within groups, suggesting 
			#strong differences between group compositions.
		#-   R ≈ 0 suggests no meaningful difference between groups.
		#-   R \< 0 indicates more variation within groups than between groups, which is rarely seen and 
			#usually points to an anomaly in the dataset or group structure.

	#The output here indicates that the fish community differs significantly among SMUs at the *P* = 0.002 level. 
		#This suggests that fish communities are structured differently depending on the SMU in which they are 
		#found. However, it is important to note that *P*-values can vary depending on the number of permutations 
		#used. Increasing the number of permutations improves the precision of the significance test, but it also 
		#increases computational time.

  		par(mfrow=c(1,2))
		plot(SMU_anosim)
	
	#This plot typically includes:
		#1.  Left Plot: A boxplot showing the distribution of rank dissimilarities within and between groups. 
			#If between-group ranks are higher than within-group ranks, this suggests group separation.
		#2.  Right Plot: A cumulative plot of the observed R values from the permutation tests, showing how the 
			#test statistic compares to the randomized data.

	#In this case, we would likely observe higher dissimilarities between SMUs than within SMUs, confirming the 
		#significant P-value. The cumulative plot helps visualize the separation by showing whether the observed 
		#R value stands out from the randomized distribution.

	#If we're more interested in habitat, we can ask a different question: *does fish community structure differ 
		#by habitat type?* This can be tested using pre-defined habitat categories such as the NLCD (National 
		#Land Cover Database) categories. To simplify, categories like "Crops-Cultivated" and "Pasture-Hay" were 
		#grouped into "Agriculture," and several water-related categories were grouped under "Wetlands":

		NLCD_Cat2 <- dat_final$NLCD_Cat
		NLCD_Cat2[NLCD_Cat2 == "Crops-Cultivated"] = "Agriculture"
		NLCD_Cat2[NLCD_Cat2 == "Pasture-Hay"] = "Agriculture"
		NLCD_Cat2[NLCD_Cat2 == "Water"] = "Wetlands"
		NLCD_Cat2[NLCD_Cat2 == "Wetlands-Emergent"] = "Wetlands"
		NLCD_Cat2[NLCD_Cat2 == "Wetlands-Woody"] = "Wetlands"
		table(NLCD_Cat2)

		NLCD_anosim <- anosim(fish_dens_log, NLCD_Cat2, permutations = 999, distance="bray")
			summary(NLCD_anosim)

	#We can see that fish community compositions are significantly shaped by habitat characteristics.

	#Another approach is to group habitats based on environmental characteristics using cluster analysis. This 
		#is useful when the habitats are not predefined and you wish to group sites based on similarities in 
		#environmental conditions:

  		env_cont_std <- decostand(env_cont, method="standardize")
		env.euc <- vegdist(env_cont_std, "euclidean")

		envcl.ward <- hclust(env.euc, method = "ward.D2")

		plot(envcl.ward, main="Ward's Minimum Variance Dendrogram", 
			xlab="Sites", 
			ylab="Euclidean Distance", 
			hang=-1)
			rect.hclust(envcl.ward, k=6)

		hclus.scree(envcl.ward, xlim=c(0,10))

		env_clus <- cutree(envcl.ward, k=6)

		env_dat_new <- cbind(env_clus, env_cont)

	#The scree plot helps identify an optimal number of clusters (in this case, 6), which can then be used to 
		#group sites and define habitat types. After grouping, these habitat types are assigned to the dataset 
		#for further analysis:

  		par(mfrow=c(2,3))
		box.plots(env_dat_new, by="env_clus")

		dat_final$HAB2 <- 0

		dat_final$HAB2[env_clus == "1"] = "High-Elev-Forest"

		dat_final$HAB2[env_clus == "4"] = "High-Gradient-Forest"

		dat_final$HAB2[env_clus == "6"] = "Low-Elev-Forest"

		dat_final$HAB2[env_clus == "2"] = "Low-Elev-Open"

		dat_final$HAB2[env_clus == "3"] = "Low-Elev-Pooled"

		dat_final$HAB2[env_clus == "5"] = "Low-Elev-Shrub-Scrub"

		hab_anosim <- anosim(fish_dens_log, dat_final$HAB2, permutations = 999, distance="bray")
			summary(hab_anosim)

	#Again, fish communities appear to be significantly different among the cluster-based habitat groups.

	#While ANOSIM is a useful tool, it has limitations:

		#-   **Assumption of ranked dissimilarities:** ANOSIM works with rank-order distances, which may not 
			#capture certain types of variation well, especially if data distributions are uneven.
		#-   **No interaction effects:** ANOSIM cannot model interaction effects, such as the combined effect of 
			#SMU and habitat type on community structure.
		#-   **Sensitivity to sampling effort:** ANOSIM's results can be sensitive to unequal sample sizes or 
			#sampling effort across groups, which can skew the interpretation of group differences.

	#Despite these limitations, ANOSIM provides a simple and effective way to test for differences in community 
		#structure across groups and is widely used in ecological and environmental studies.

## Multivariate Analysis of Variance

	#**Multivariate Analysis of Variance (MANOVA)** is an extension of ANOVA that allows for the analysis of 
		#multiple response variables simultaneously. MANOVA tests whether the mean vectors of several groups are 
			#significantly different from each other. It evaluates the effect of one or more independent variables
			#on multiple dependent variables, making it particularly useful in ecological studies where response 
			#variables (like species abundances) are often interdependent.

	#MANOVA assumes:
		#-   **Multivariate normality:** The response variables should be normally distributed within each group.
		#-   **Homogeneity of covariance matrices:** The variance-covariance matrices of the dependent variables 
			#should be the same across all groups.

	#These assumptions can be restrictive, particularly when working with ecological or environmental data that do 
		#not meet these conditions.

	#**Permutational Multivariate Analysis of Variance (PERMANOVA)**, on the other hand, is a more flexible and 
		#robust alternative that does *not* require assumptions of normality or homogeneity of variances. Instead, 
		#PERMANOVA tests for differences among group centroids based on a distance matrix (e.g., Bray-Curtis), 
		#using permutations to assess statistical significance.

	#Again, the code to run these functions is quite simple.

	#Before running a MANOVA, it's good practice to check the assumption of multivariate homogeneity of variances 
		#(MHV). This can be done using an MHV test, which can be visually assessed via plots:

  		fish.bray <- vegdist(fish_dens_log, "bray")

		MHV_smu <- betadisper(fish.bray, dat_final$SMU)
			permutest(MHV_smu)
			plot(MHV_smu, hull=FALSE, ellipse=TRUE)

	#The output from `permutest()` provides a test for multivariate dispersion, which is related to the assumption 
		#of homogeneity of variances. If the test is **non-significant**, it indicates that the assumption holds 
		#for the comparison of SMUs. Do our data violate this assumption?

	#We can also check for separation of means after running the MANOVA using the Wilks' Lambda statistic.

		smu_manova <- manova(as.matrix(fish_dens_log) ~ dat_final$SMU)
			summary(smu_manova)
			summary(smu_manova, test = "Wilks")

	#Because MANOVA relies on stricter assumptions, PERMANOVA is often the preferred method for ecological data. 
		#PERMANOVA can be implemented using the `adonis2()` function, which evaluates group differences based on 
		#a distance matrix. This method is more robust to violations of assumptions like normality and homogeneity 
		#of variances:

		smu_perm <- adonis2(fish.bray ~ dat_final$SMU, permutations = 999)
			smu_perm

	#The output provides an F statistic and a P-value, much like traditional ANOVA, but based on permutation 
		#testing. A significant result suggests that fish community structure varies significantly among SMUs.

	#The `adonis2` function can also handle interaction effects!

		int_perm <- adonis2(fish.bray ~ dat_final$SMU + dat_final$HAB2, permutations = 999)
			int_perm

		int2_perm <- adonis2(fish.bray ~ dat_final$SMU*dat_final$HAB2, permutations = 999)
			int2_perm

	#AND continuous predictor variables!

		elev_perm <- adonis2(fish.bray ~ dat_final$Elev, permutations = 999)
			elev_perm

		canopy_perm <- adonis2(fish.bray ~ dat_final$Canopy, permutations = 999)
			canopy_perm

	#This flexibility allows PERMANOVA to test for relationships between multivariate community structure and 
		#both categorical and continuous predictors.

## Similarity Percentages

	#**Similarity Percentages (SIMPER)** is a post-hoc analysis that helps identify which variables (here, species) 
		#contribute most to the observed differences among groups. It breaks down the overall dissimilarity 
		#between groups (based on a distance measure like Bray-Curtis) into contributions from individual 
		#species. This is particularly useful for ecological studies where researchers want to determine which 
		#species are driving differences in community composition.

	#For our class dataset, SIMPER can help answer the question: *Which species are driving among-SMU differences 
		#in community structure?*

		smu_simper <- simper(fish_dens_log, dat_final$SMU)
			summary(smu_simper)

	#The SIMPER output typically includes:
		#-   **Species contributions:** A list of species along with their contribution to the overall 
			#dissimilarity between group pairs. Species that contribute the most are responsible for driving 
			#the differences.
		#-   **Cumulative percentage:** The cumulative percentage contribution, showing how much of the 
			#dissimilarity is explained by the top contributing species. This helps determine which species 
			#are most influential in differentiating among groups.
		#-   **Average dissimilarity:** The average dissimilarity between pairs of groups based on species 
			#composition.

	#If a particular species (e.g., dace) consistently contributes more than others, it indicates that this 
		#species' abundance patterns differ significantly across management units.

## Mantel Test

	#The **Mantel Test** is a non-parametric test used to evaluate the correlation between two distance matrices, 
		#typically measuring ecological or genetic similarity and geographic or environmental distances. It is 
		#commonly used to assess whether community composition is influenced by spatial or environmental gradients.

	#The Mantel test operates by calculating the correlation (usually Pearson or Spearman) between corresponding 
		#entries in two distance matrices and then permutes the rows and columns of one matrix to generate a null 
		#distribution for the correlation coefficient. The test evaluates whether the observed correlation is 
		#greater than expected by chance.

	#In ecological studies, the Mantel test is often used to relate species composition (distance matrix of 
		#community dissimilarities) to environmental variables or geographic distances. For instance, it can 
		#address questions like: *Does geographic distance influence fish community structure?* or *Is there a 
		#correlation between environmental gradients and species composition?*

	#The code below shows how to apply the Mantel test to compare the geographic and environmental distances, and 
		#then geographic and species distances:

		dist.fish <- vegdist(fish_dens_log, method="bray")

		dist.env <- vegdist(env_cont, method="euclidean")

		geog <- cbind(dat_final$Latitude, dat_final$Longitude)
		dist.geog <- vegdist(geog, method="euclidean")

		genv_mantel <- mantel(dist.geog, dist.env, method="spearman", permutations=999, na.rm=TRUE)
			genv_mantel
	
		gfish_mantel <- mantel(dist.geog, dist.fish, method="spearman", permutations=999, na.rm=TRUE)
			gfish_mantel

	#Output includes:
		#-   *Correlation coefficient (r value):* This indicates the strength of the relationship between the 
			#two matrices (geographic and environmental, or geographic and fish community distances). A higher 
			#r value means a stronger correlation between the two matrices.
		#-   *P-value:* Derived from the permutation test. A low p-value (e.g., p \< 0.05) indicates that the 
			#observed correlation is significant and unlikely to have arisen by chance.

	#We can see that both environmental and species distances are associated with geographic distances in our class 
		#dataset.

## Procrustes Analysis

	#A **Procrustes Analysis** is used to compare two multivariate data sets by assessing their similarity in 
		#multivariate space. It minimizes the sum of squared differences between corresponding points in two 
		#datasets through translation, rotation, and scaling, effectively aligning the two configurations. In 
		#ecological applications, it is typically used to compare ordinations from different data sets, such as 
		#species composition and environmental variables, or different transformations of the same dataset.

	#The Procrustes analysis is often followed by a **Procrustes Randomization Test (PROTEST)**, which assesses the 
		#significance of the Procrustes fit using permutations.

	#In the example below, we first perform principal component analysis (PCA) on both environmental and species 
		#data, and then use Procrustes analysis to compare these ordination results.

		pca.env <- prcomp(env_cont, scale = TRUE)

		fish.hel <- decostand(fish_dens, method="hellinger")
		pca.hel <- prcomp(fish.hel, scale = FALSE)

		procr <- procrustes(pca.env, pca.hel, symmetric=FALSE)
			procr

		protest(pca.env, pca.hel, scores="sites", permutations = 999)

	#Procrustes analysis results include:
		#-   **Procrustes residuals:** These are the differences between corresponding points after translation, 
			#rotation, and scaling. Smaller residuals indicate better alignment between the two data sets.
		#-   **Protest output:** The `protest()` function tests the significance of the Procrustes fit. It 
			#provides a correlation-like statistic (m\^2), representing how well the two ordinations match, and 
			#a p-value from permutation tests. A low p-value suggests that the ordination structures are 
			#significantly similar.

	#Two types of plots can be used to interpret Procrustes analysis visually:

		plot(procr, kind = 1, type = "text")
		plot(procr, kind = 2)

	#-   **Plot 1:** This shows the Procrustes residuals as arrows. The length and direction of the arrows 
			#represent the degree of difference between corresponding points (sites) in the two ordination 
			#spaces. Shorter arrows indicate better alignment.
	#-   **Plot 2:** This shows the Procrustes residuals as bars.

	#This concludes our journey through the world of statistical inference in multivariate analysis. Throughout 
		#the lab, we explored techniques like ANOSIM, PERMANOVA, SIMPER, Mantel tests, and Procrustes analysis. 
		#These tools are essential for understanding complex relationships between environmental variables and 
		#species communities, enabling us to make informed ecological inferences.


