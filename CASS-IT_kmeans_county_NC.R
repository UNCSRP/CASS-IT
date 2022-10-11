#################################################################################################
#################################################################################################
#### CASS-IT: county level social and environmental exposures in NC
####
#### Calculating pairwise correlations between social and metals variables at county-level in NC
#### Conducting k-mean clustering analysis on social and metals variables 
#### 
#### INITIAL DRAFT:01/22/2021 
#### UPDATED: 01/27/2021, 02/09/2021, 02/15/2021, 02/16/2021, 04/27/2021, 09/21/2021, 10/03/2022
#### by Lauren Eaves                  
#################################################################################################
#################################################################################################

#################################################################################################
#################################################################################################
#### Step 1: Installing and activating appropriate packages in R 
#################################################################################################
#################################################################################################

sessionInfo()
rm(list=ls())

library(dplyr)
library(plyr)
library(tidyverse)
library(corrplot)
library(data.table)
library(Hmisc)
library(stats)
library(factoextra)
library(RColorBrewer)
library(NbClust)
library(tidycensus)
library(dplyr)
library(tidyverse)
library(corrplot)
library(data.table)
library(Hmisc)
library(sf)
library(ggplot2)
library(leaflet)
library(classInt)
library(biscale)
library(cowplot)
library(lsr)
library(missForest)

#################################################################################################
#################################################################################################
#### Step 2: Set working directory, output folder path, load data, clean IDs etc
#################################################################################################
#################################################################################################

# Set working directory
setwd("/Users/lauren_eaves/IEHS Dropbox/lauren Eaves/1_Projects/13_DMAC/2_Aim1/2_correlations_clustering/input")
getwd()
# Set output folder
Output_Folder <- ("/Users/lauren_eaves/IEHS Dropbox/lauren Eaves/1_Projects/13_DMAC/2_Aim1/2_correlations_clustering/output/vulnerabilityresourcesmetals")

# Create a current date variable to name outputfiles
cur_date <- paste0(str_replace_all(Sys.Date(),"-",""))

# Import county-level data
county <- read.csv(file="20220601_DMAC_NCcounties_metals_vulnerability_resources.csv")
dim(county)
colnames(county)
county[1]<-NULL


#################################################################################################
#################################################################################################
#### Step 3: Data cleaning, imputation and distributions 
#################################################################################################
#################################################################################################

#vulnerability variables:
#Social Vulnerability Index https://www.atsdr.cdc.gov/placeandhealth/svi/at-a-glance_svi.html :
#SVI_RPL_THEME1_SOCIECO = Percentile ranking for socioeconomic theme)
#SVI_RPL_THEME2_HH_DISB = Percentile ranking for household composition theme)
#SVI_RPL_THEME3_MINO = Percentile ranking for minority status/language theme)
#SVI_RPL_THEME4_HH_TRANS = Percentile ranking for housing/transportation theme)

vulnerabilityvars <-c("SVI_RPL_THEME1_SOCIECO","SVI_RPL_THEME2_HH_DISB", "SVI_RPL_THEME3_MINO","SVI_RPL_THEME4_HH_TRANS")

#Toxic metal variables: 
#NCWELL database: https://dataverse.unc.edu/dataset.xhtml?persistentId=doi:10.15139/S3/BDQG9O 
#mean levels of Arsenic, Lead, and Manganese
metalsvars <-c("Arsenic.Mean_avg", "Lead.Mean_avg","Manganese.Mean_avg")

#Resources variables 
#FEMA RAPT tool: social resources and health resources https://www.fema.gov/sites/default/files/documents/fema_rapt-user-guide-2022.pdf
#FCC county connections database: information resources https://www.fcc.gov/form-477-county-data-internet-access-services 
resourcesvars <- c("resources_social","resources_health","resources_info")

#There are 7 counties for which the resources_info = -9999 where the data has been masked for firm confidentiality
#Will impute these values using Random Forest Modelling (missForest package)

county <- county %>% 
  mutate(resources_info=ifelse(resources_info==-9999.0,NA,resources_info))%>% 
  column_to_rownames(var="county")

# Set a seed number for reproducibility
set.seed(123)

# Running the MissForect function, across all included variables
county_imp <- missForest(county, variablewise = FALSE, verbose=TRUE)
county_imp_df <- county_imp$ximp    # resulting imputed values included alongside covariates
county_imp$OOBerror   # model performance parameters

#Use imputed dataset from now on
county <- county_imp_df

#To make the resources variables more interpretable, we will negate them so they become higher = lack of resources
summary(county$resources_health)
summary(county$resources_info)
summary(county$resources_social)

county <- county %>% 
  mutate(resources_health = -resources_health) %>% 
  mutate(resources_info = -resources_info) %>% 
  mutate(resources_social = -resources_social) 

summary(county$resources_health)
summary(county$resources_info)
summary(county$resources_social)


#Calculate distributions for all the variables
min_max_etc <- list(
  min = ~min(.x, na.rm = TRUE), 
  max = ~max(.x, na.rm = TRUE),
  mean = ~mean(.x, na.rm=TRUE),
  median = ~median(.x, na.rm=TRUE),
  SD = ~sd(.x, na.rm=TRUE),
  IQR = ~IQR(.x, na.rm=TRUE)
)
minmaxetc <- county %>% 
  dplyr::summarise(across(.cols=everything(),  min_max_etc, .names = "{.col}.{.fn}"))
minmaxetc <- minmaxetc %>% 
  gather(key="Measure", value="Value") %>% 
  mutate(Variable=str_extract(Measure, "^([^.])+")) %>% 
  mutate(Stat=str_extract(Measure, "[^.]+$")) %>% 
  select(-Measure) %>% 
  spread(key="Stat",value="Value")

write.csv(minmaxetc, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_vulnerability_metal_resources_statedistributions.csv"), row.names= TRUE)

write.csv(county, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_input_data_to_Kmeans_pre_zscore_premetallog_transformation.csv"), row.names= TRUE)


#Log transform the metals variables 
county <- county %>% 
  mutate(Arsenic.Mean_avg=log2(Arsenic.Mean_avg)) %>% 
  mutate(Lead.Mean_avg=log2(Lead.Mean_avg)) %>% 
  mutate(Manganese.Mean_avg=log2(Manganese.Mean_avg))

write.csv(county, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_input_data_to_Kmeans_pre_zscore_transformation.csv"), row.names= TRUE)


#################################################################################################
#################################################################################################
#### Step 4: Calculate Spearman correlations pew-z score transformation 
#################################################################################################
#################################################################################################

county_forcorr<-as.matrix(county)

county_forcorr_spear <-rcorr(county_forcorr, type=c("spearman"))
#access correlation matrix 
county_forcorr_corrmat_spear <- as.matrix(county_forcorr_spear[[1]])
#access p values 
county_forcorr_pvalues_spear <-as.matrix(county_forcorr_spear[[3]])

county_forcorr_corrmat_spear <- as.data.frame(county_forcorr_corrmat_spear) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_corrmat_spear<-as.data.frame(t(county_forcorr_corrmat_spear)) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_corrmat_spear <-as.matrix(county_forcorr_corrmat_spear)


county_forcorr_pvalues_spear <- as.data.frame(county_forcorr_pvalues_spear) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_pvalues_spear<-as.data.frame(t(county_forcorr_pvalues_spear)) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_pvalues_spear <-as.matrix(county_forcorr_pvalues_spear)


corrplot(county_forcorr_corrmat_spear, method="color", tl.col="black",
         tl.cex = 0.75, tl.srt=45, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), p.mat=county_forcorr_pvalues_spear)

png(file = (paste0(Output_Folder,"/", cur_date, "_Spearmancorrelations_correlogram.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
corrplot(county_forcorr_corrmat_spear, method="color", tl.col="black",
         tl.cex = 0.75, tl.srt=45, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), p.mat=county_forcorr_pvalues_spear)

dev.off()

coeff <- as.data.frame(county_forcorr_corrmat_spear) %>% 
  rownames_to_column(var="Stressor 1") %>% 
  gather(key="Stressor 2", value="Spearman_coefficient",2:11) %>% 
  filter(Spearman_coefficient != 1)

pvalues <- as.data.frame(county_forcorr_pvalues_spear) %>% 
  rownames_to_column(var="Stressor 1") %>% 
  gather(key="Stressor 2", value="P-value", 2:11) 

joined <- inner_join(coeff, pvalues, by=c("Stressor 1", "Stressor 2"))

write.csv(joined, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_social_and_metal_Spearmancorrelations.csv"), row.names= TRUE)


#################################################################################################
#################################################################################################
#### Step 5: Calculate Pearson correlations on z-score transformed values 
#################################################################################################
#################################################################################################

county_forcorr<-as.matrix(scale(county))

county_forcorr_pear <-rcorr(county_forcorr, type=c("pearson"))
#access correlation matrix 
county_forcorr_corrmat_pear <- as.matrix(county_forcorr_pear[[1]])
#access p values 
county_forcorr_pvalues_pear <-as.matrix(county_forcorr_pear[[3]])


county_forcorr_corrmat_pear <- as.data.frame(county_forcorr_corrmat_pear) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_corrmat_pear<-as.data.frame(t(county_forcorr_corrmat_pear)) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_corrmat_pear <-as.matrix(county_forcorr_corrmat_pear)


county_forcorr_pvalues_pear <- as.data.frame(county_forcorr_pvalues_pear) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_pvalues_pear<-as.data.frame(t(county_forcorr_pvalues_pear)) %>% 
  dplyr::rename("Arsenic"="Arsenic.Mean_avg",
                "Lead"="Lead.Mean_avg",
                "Manganese"="Manganese.Mean_avg",
                "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                "Low Social resources" ="resources_social",
                "Low Health resources" = "resources_health",
                "Low Information resources" = "resources_info")
county_forcorr_pvalues_pear <-as.matrix(county_forcorr_pvalues_pear)

corrplot(county_forcorr_corrmat_pear, method="color", tl.col="black",
         tl.cex = 0.75, tl.srt=45, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), p.mat=county_forcorr_pvalues_pear)

png(file = (paste0(Output_Folder,"/", cur_date, "_Pearsoncorrelations_correlogram.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
corrplot(county_forcorr_corrmat_pear, method="color", tl.col="black",
         tl.cex = 0.75, tl.srt=45, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), p.mat=county_forcorr_pvalues_pear)
dev.off()


coeff <- as.data.frame(county_forcorr_corrmat_pear) %>% 
  rownames_to_column(var="Stressor 1") %>% 
  gather(key="Stressor 2", value="Pearson_coefficient",2:11) %>% 
  filter(Pearson_coefficient != 1)

pvalues <- as.data.frame(county_forcorr_pvalues_pear) %>% 
  rownames_to_column(var="Stressor 1") %>% 
  gather(key="Stressor 2", value="P-value", 2:11) 

joined <- inner_join(coeff, pvalues, by=c("Stressor 1", "Stressor 2"))

write.csv(joined, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_social_and_metal_Pearsoncorrelations.csv"), row.names= TRUE)

#################################################################################################
#################################################################################################
#### Step 6: K-means clustering algorithm
#################################################################################################
#################################################################################################
##Useful links for this analysis:

# https://www.datanovia.com/en/blog/types-of-clustering-methods-overview-and-quick-start-r-code/
# http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization#visualizing-dimension-reduction-analysis-outputs
# https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/

#For clustering, need to make sure that there are no missing values and need to scale variables  
#so that they are comparable 
sapply(county, function(x) sum(is.na(x)))
county <- as.data.frame(scale(county)) #scaled 
set.seed(123) #set seed for reproducibility 
write.csv(county, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_input_data_to_Kmeans_zscore_transformation.csv"), row.names= TRUE)


#K-means with K = 1
km1 <- kmeans(county, 1, nstart = 20) #20 iterations, 1 cluster
km1$size 
head(km1$centers) ## mean of each chemical in the cluster 
mean(county[,1])
km1$totss ## total sum of squares
km1$withinss ## within cluster sum of squares by cluster- within cluster variation for each cluster 
km1$tot.withinss ## total within cluster sum of squares -- this is what we want to minimize!
km1$betweenss ## between cluster sum of squares
Percent_bt_k1 = paste(round(100*(km1$betweenss/km1$totss),1), "%", sep = "")

#K-means with K = 2
km2 <- kmeans(county, 2, nstart = 20) 
km2$size 
head(km2$centers) #now seeing some separation of the means by the different clusters 
km2$totss 
km2$withinss
km2$tot.withinss #minimize this
km2$betweenss #now have a between cluster variation, maximize this 
Percent_bt_k2 = paste(round(100*(km2$betweenss/km2$totss),1), "%", sep = "")

#K-means with K = 3
km3 <- kmeans(county, 3, nstart = 20)
km3$size
head(km3$centers)
km3$totss 
km3$withinss 
km3$tot.withinss  
km3$betweenss 
Percent_bt_k3 = paste(round(100*(km3$betweenss/km3$totss),1), "%", sep = "")

#K-means with K = 4
km4 <- kmeans(county, 4, nstart = 20)
km4$size
head(km4$centers)
km4$totss 
km4$withinss 
km4$tot.withinss  
km4$betweenss 
Percent_bt_k4 = paste(round(100*(km4$betweenss/km4$totss),1), "%", sep = "")

#K-means with K = 5
km5 <- kmeans(county, 5, nstart = 20)
km5$size
head(km5$centers)
km5$totss 
km5$withinss 
km5$tot.withinss  
km5$betweenss 
Percent_bt_k5 = paste(round(100*(km5$betweenss/km5$totss),1), "%", sep = "")

#K-means with K = 6
km6 <- kmeans(county, 6, nstart = 20)
km6$size
head(km6$centers)
km6$totss 
km6$withinss 
km6$tot.withinss  
km6$betweenss 
Percent_bt_km6 = paste(round(100*(km6$betweenss/km6$totss),1), "%", sep = "")

#K-means with K = 7
km7 <- kmeans(county, 7, nstart = 20)
km7$size
head(km7$centers)
km7$totss 
km7$withinss 
km7$tot.withinss  
km7$betweenss 
Percent_bt_k7 = paste(round(100*(km7$betweenss/km7$totss),1), "%", sep = "")

#K-means with K = 8
km8 <- kmeans(county, 8, nstart = 20)
km8$size
head(km8$centers)
km8$totss 
km8$withinss 
km8$tot.withinss  
km8$betweenss 
Percent_bt_k8 = paste(round(100*(km8$betweenss/km8$totss),1), "%", sep = "")

#K-means with K = 9
km9 <- kmeans(county, 9, nstart = 20)
km9$size
head(km9$centers)
km9$totss 
km9$withinss 
km9$tot.withinss  
km9$betweenss 
Percent_bt_km9 = paste(round(100*(km9$betweenss/km9$totss),1), "%", sep = "")

#K-means with K = 10
km10 <- kmeans(county, 10, nstart = 20)
km10$size
head(km10$centers)
km10$totss 
km10$withinss 
km10$tot.withinss  
km10$betweenss 
Percent_bt_k10 = paste(round(100*(km10$betweenss/km10$totss),1), "%", sep = "")

# create an empty data frame
km.res <- data.frame(matrix(NA, 99, 4))

# use a loop to run k-means for 0 - 30 clusters and generate variance values
for (k in 1:99) {
  set.seed(123)
  kk <- kmeans(county, k, nstart = 20)
  kk$size <- toString(kk$size)
  km.res[k, ] <- cbind(k, kk$tot.withinss, kk$totss, kk$size) 
}


# name the variables created
names(km.res) <- c("clusters", "WithinSS", "TotSS", "ClusterSizes")

# create a variable for within cluster variance/total SS
km.res$WithinSS <- as.numeric(km.res$WithinSS)
km.res$TotSS <-as.numeric(km.res$TotSS)
km.res$clusters <-as.numeric(km.res$clusters)
km.res$PropWithin <- 100*(km.res$WithinSS/km.res$TotSS)

#output K-means results for all possible k 
write.csv(km.res, paste0(Output_Folder,"/", cur_date, "_Kmeans_allks_results.csv"), row.names= TRUE)



#################################################################################################
#################################################################################################
#### Step 7: K-means clustering algorithm: decide on optimal k through balancing the number of single
#### county clusters (don't want any) and the proportion of within cluster variance (want this to be minimal)
#################################################################################################
#################################################################################################

gapstatplot <- fviz_nbclust(county, kmeans, method = "gap_stat", nboot = 500, k.max = 5)+
  labs(subtitle = "Gap statistic method")
plot(gapstatplot)
silh <- fviz_nbclust(county, kmeans, method = "silhouette", nboot = 500, k.max =5)+
  labs(subtitle = "Silhouette method")
plot(silh)
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_gapstatisticplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(gapstatplot)
dev.off()
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_silhouetteplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(silh)
dev.off()
NbClust(data = county, diss = NULL, distance = "euclidean", min.nc = 3, max.nc = 5, 
        method = "kmeans")

# compare the within cluster variance/total SS to number of clusters
km.res %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_line() +
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
propwithinplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_point() + 
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=100,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))
plot(propwithinplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_propwithin_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(propwithinplot)
dev.off()

#compare number of single county clusters 
km.res <- km.res %>% 
  mutate(SinglecountyClusters = str_count(ClusterSizes, regex("\\s[1],\\s")))
km.res %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_line() +
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
singleclustersplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_point() + 
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=12,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))+
  scale_y_continuous(breaks=seq(0,12,2))
plot(singleclustersplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_singleclusters_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(singleclustersplot)
dev.off()

#################################################################################################
#################################################################################################
#### Step 8: K-means clustering algorithm: vizualizations of k=4 solution
#################################################################################################
#################################################################################################

#creating dataframe with the mean of each cluster for each metal in z-score transformed value
km4_centers <- as.data.frame.matrix(t(km4$centers)) #dataframe containing the cluster means for each metal
km4_centers$var <- row.names(km4_centers)
km4_centers <- km4_centers %>% 
  mutate(Group = ifelse(var %in% metalsvars,"Metal",
                        ifelse(var %in% vulnerabilityvars, "Vulnerability", "Resources"))) %>% 
  gather(key = "Cluster", value = "mean", -var, -Group) %>% 
  dplyr::rename(Mean_zscore = mean)

#bargraphs of each cluster with each component represented 
#install.packages("remotes")
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

plot_chem_means_zscore <- km4_centers %>% 
  mutate(Cluster = as.factor(Cluster),
         Cluster = fct_recode(Cluster, "Cluster 1" = "1",
                              "Cluster 2" = "2",
                              "Cluster 3" = "3",
                              "Cluster 4" = "4"),
         var = fct_recode(var, "Arsenic"="Arsenic.Mean_avg",
                          "Lead"="Lead.Mean_avg",
                          "Manganese"="Manganese.Mean_avg",
                          "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                          "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                          "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                          "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS",
                          "Low Social resources" ="resources_social",
                          "Low Health resources" = "resources_health",
                          "Low Information resources" = "resources_info"),
         Group = as.factor(Group)) %>%
  ggplot(aes(x = var, y = Mean_zscore, fill = Cluster)) + geom_col() +
  #geom_point(aes(y = Mean_zscore), size = 1) +
  facet_wrap(~ Cluster) + theme_bw() + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.text.y=element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.caption = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, size = 0.2) +
  labs(x = "Stressor",
       y = "Z-score Standardized Mean")+
  scale_fill_brewer(palette = "Set2", name = "Cluster") + 
  geom_bar_pattern(stat = "identity",
                   pattern_color = "black",
                   pattern_fill = "black",
                   pattern_density = 0.1,
                   pattern_spacing = 0.02,
                   aes(pattern = Group))+
  theme(legend.position = "none")

plot(plot_chem_means_zscore)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_clustersbargraph_zscoredvalues.png")), width = 8, height = 10, units = "in", pointsize = 12, res = 800)
plot(plot_chem_means_zscore)
dev.off()

#PCA based plot
pca_clusterplot <- fviz_cluster(km4, data = county,
                                ggtheme = theme_minimal(),
                                ellipse.type = "convex",
                                palette= "Set2",
                                main = "",
                                geom.var = c("point", "text"),
                                repel=TRUE, max.overlaps=Inf)
plot(pca_clusterplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_PCAclusterplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 800)
plot(pca_clusterplot)
dev.off()

#Map with coloring by clusters
#dataframe with counties and their clusters
countycluster <- cbind(county, cluster = km4$cluster) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="county") %>% 
  mutate(county=as.character(str_trim(county)))

#Get census/geographic info needed
#API key to access census data. Request your own API key at http://api.census.gov/data/key_signup.html 
#census_api_key("0661e062882e541f3caf807f0e7a8de5a2add7ae", install=TRUE, overwrite=TRUE) 
#readRenviron("~/.Renviron")
#v18 <- load_variables(2018, "acs5", cache = TRUE) #loading variables in ACS to search for 
census <-  tidycensus::get_estimates(state="NC", geography = "county",
                                     product = "population", year = 2018, geometry = TRUE) %>%
  dplyr::filter(variable=="POP") %>% 
  mutate(county=as.character(str_trim((str_remove(NAME,"County, North Carolina")))))%>% 
  as.data.frame(.) 

formap <- join(countycluster, census, by="county")

Cluster_map <-
  formap %>% 
  mutate(cluster=as.factor(cluster)) %>% 
  ggplot(aes(fill = cluster, geometry=geometry)) + 
  ggthemes::theme_map() + 
  geom_sf(color = "black", size=0.1) +
  scale_fill_brewer(palette = "Set2", name="cluster", direction = 1,
                    labels=c("1", "2","3","4","5")) +
  theme(legend.position="right") 
plot(Cluster_map)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_mapofNC_withclusters.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(Cluster_map)
dev.off()


#################################################################################################
#################################################################################################
#### Step 9: Run one way ANOVAs to compare differences in mean values of metals and SVIs between clusters 
#################################################################################################
#################################################################################################

#Null hypothesis: the means of the different groups are the same
#Alternative hypothesis: At least one sample mean is not equal to the others.

km4_anovas <- countycluster %>% select(-county) %>% 
  mutate(cluster=as.factor(cluster))

results.anova <- data.frame()
for (i in colnames(km4_anovas)[1:10]) {
  print(i)
  km4_anovas$var <- as.numeric(km4_anovas[[i]])
  anova <- aov(var ~ cluster, data = km4_anovas)
  output.anova <- summary(anova)
  results <- output.anova[[1]] %>% 
    dplyr::rename(p.ftest = "Pr(>F)") %>% 
    mutate(var=i,
           p.ftest.adjust=p.adjust(p.ftest,method="BH",n=10)) #p adjustment for conducting 10 different anovas
  results.anova <- rbind(results.anova, results) 
}
colnames(results.anova)

write.csv(results.anova, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_social_and_metal_ANOVAresults.csv"), row.names= TRUE)


#################################################################################################
#################################################################################################
#### Step 10: export datasheet with each county, raw and z-score  and cluster  
#################################################################################################
#################################################################################################

countycluster <- countycluster %>% select(county,cluster)
write.csv(countycluster, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_clusterassign.csv"), row.names= TRUE)

#################################################################################################
#################################################################################################
#### Step 11: FOR VULNERABILTY VARIABLES ONLY K-means clustering algorithm 
#################################################################################################
#################################################################################################

# Reset output folder to a folder containing only results based on vulnerability variables only analysis
Output_Folder <- ("/Users/lauren_eaves/IEHS Dropbox/lauren Eaves/1_Projects/13_DMAC/2_Aim1/2_correlations_clustering/output/vulnerability")

county_vulnerability <- county %>% 
  select(all_of(c(vulnerabilityvars)))


#K-means with K = 1
km1 <- kmeans(county_vulnerability, 1, nstart = 20) #20 iterations, 1 cluster
km1$size 
head(km1$centers) ## mean of each chemical in the cluster 
km1$totss ## total sum of squares
km1$withinss ## within cluster sum of squares by cluster- within cluster variation for each cluster 
km1$tot.withinss ## total within cluster sum of squares -- this is what we want to minimize!
km1$betweenss ## between cluster sum of squares
Percent_bt_k1 = paste(round(100*(km1$betweenss/km1$totss),1), "%", sep = "")

#K-means with K = 2
km2 <- kmeans(county_vulnerability, 2, nstart = 20) 
km2$size 
head(km2$centers) #now seeing some separation of the means by the different clusters 
km2$totss 
km2$withinss
km2$tot.withinss #minimize this
km2$betweenss #now have a between cluster variation, maximize this 
Percent_bt_k2 = paste(round(100*(km2$betweenss/km2$totss),1), "%", sep = "")

#K-means with K = 3
km3 <- kmeans(county_vulnerability, 3, nstart = 20)
km3$size
head(km3$centers)
km3$totss 
km3$withinss 
km3$tot.withinss  
km3$betweenss 
Percent_bt_k3 = paste(round(100*(km3$betweenss/km3$totss),1), "%", sep = "")

#K-means with K = 4
km4 <- kmeans(county_vulnerability, 4, nstart = 20)
km4$size
head(km4$centers)
km4$totss 
km4$withinss 
km4$tot.withinss  
km4$betweenss 
Percent_bt_k4 = paste(round(100*(km4$betweenss/km4$totss),1), "%", sep = "")

#K-means with K = 5
km5 <- kmeans(county_vulnerability, 5, nstart = 20)
km5$size
head(km5$centers)
km5$totss 
km5$withinss 
km5$tot.withinss  
km5$betweenss 
Percent_bt_k5 = paste(round(100*(km5$betweenss/km5$totss),1), "%", sep = "")

#K-means with K = 6
km6 <- kmeans(county_vulnerability, 6, nstart = 20)
km6$size
head(km6$centers)
km6$totss 
km6$withinss 
km6$tot.withinss  
km6$betweenss 
Percent_bt_km6 = paste(round(100*(km6$betweenss/km6$totss),1), "%", sep = "")

#K-means with K = 7
km7 <- kmeans(county_vulnerability, 7, nstart = 20)
km7$size
head(km7$centers)
km7$totss 
km7$withinss 
km7$tot.withinss  
km7$betweenss 
Percent_bt_k7 = paste(round(100*(km7$betweenss/km7$totss),1), "%", sep = "")

#K-means with K = 8
km8 <- kmeans(county_vulnerability, 8, nstart = 20)
km8$size
head(km8$centers)
km8$totss 
km8$withinss 
km8$tot.withinss  
km8$betweenss 
Percent_bt_k8 = paste(round(100*(km8$betweenss/km8$totss),1), "%", sep = "")

#K-means with K = 9
km9 <- kmeans(county_vulnerability, 9, nstart = 20)
km9$size
head(km9$centers)
km9$totss 
km9$withinss 
km9$tot.withinss  
km9$betweenss 
Percent_bt_km9 = paste(round(100*(km9$betweenss/km9$totss),1), "%", sep = "")

#K-means with K = 10
km10 <- kmeans(county_vulnerability, 10, nstart = 20)
km10$size
head(km10$centers)
km10$totss 
km10$withinss 
km10$tot.withinss  
km10$betweenss 
Percent_bt_k10 = paste(round(100*(km10$betweenss/km10$totss),1), "%", sep = "")

# create an empty data frame
km.res <- data.frame(matrix(NA, 99, 4))

# use a loop to run k-means for 0 - 30 clusters and generate variance values
for (k in 1:99) {
  set.seed(123)
  kk <- kmeans(county_vulnerability, k, nstart = 20)
  kk$size <- toString(kk$size)
  km.res[k, ] <- cbind(k, kk$tot.withinss, kk$totss, kk$size) 
}


# name the variables created
names(km.res) <- c("clusters", "WithinSS", "TotSS", "ClusterSizes")

# create a variable for within cluster variance/total SS
km.res$WithinSS <- as.numeric(km.res$WithinSS)
km.res$TotSS <-as.numeric(km.res$TotSS)
km.res$clusters <-as.numeric(km.res$clusters)
km.res$PropWithin <- 100*(km.res$WithinSS/km.res$TotSS)

#output K-means results for all possible k 
write.csv(km.res, paste0(Output_Folder,"/", cur_date, "_Kmeans_allks_results.csv"), row.names= TRUE)



#################################################################################################
#################################################################################################
#### Step 12: FOR VULNERABILTY VARIABLES ONLY K-means clustering algorithm: decide on optimal k through balancing the number of single
#### county clusters (don't want any) and the proportion of within cluster variance (want this to be minimal)
#################################################################################################
#################################################################################################

gapstatplot <- fviz_nbclust(county_vulnerability, kmeans, method = "gap_stat", nboot = 500, k.max = 5)+
  labs(subtitle = "Gap statistic method")
plot(gapstatplot)
silh <- fviz_nbclust(county_vulnerability, kmeans, method = "silhouette", nboot = 500, k.max =5)+
  labs(subtitle = "Silhouette method")
plot(silh)
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_gapstatisticplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(gapstatplot)
dev.off()
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_silhouetteplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(silh)
dev.off()
NbClust(data = county_vulnerability, diss = NULL, distance = "euclidean", min.nc = 3, max.nc = 5, 
        method = "kmeans")

# compare the within cluster variance/total SS to number of clusters
km.res %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_line() +
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
propwithinplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_point() + 
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=100,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))
plot(propwithinplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_propwithin_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(propwithinplot)
dev.off()

#compare number of single county clusters 
km.res <- km.res %>% 
  mutate(SinglecountyClusters = str_count(ClusterSizes, regex("\\s[1],\\s")))
km.res %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_line() +
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
singleclustersplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_point() + 
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=12,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))+
  scale_y_continuous(breaks=seq(0,12,2))
plot(singleclustersplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_singleclusters_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(singleclustersplot)
dev.off()

#################################################################################################
#################################################################################################
#### Step 13: FOR VULNERABILTY VARIABLES ONLY K-means clustering algorithm: vizualizations of k=4 solution
#################################################################################################
#################################################################################################

#creating dataframe with the mean of each cluster for each metal in z-score transformed value
km4_centers <- as.data.frame.matrix(t(km4$centers)) #dataframe containing the cluster means for each metal
km4_centers$var <- row.names(km4_centers)
km4_centers <- km4_centers %>% 
  gather(key = "Cluster", value = "mean", -var) %>% 
  dplyr::rename(Mean_zscore = mean)

plot_chem_means_zscore <- km4_centers %>% 
  mutate(Cluster = as.factor(Cluster),
         Cluster = fct_recode(Cluster, "Cluster 1" = "1",
                              "Cluster 2" = "2",
                              "Cluster 3" = "3",
                              "Cluster 4" = "4"),
         var = fct_recode(var, 
                          "Minority Status and Language"="SVI_RPL_THEME3_MINO",
                          "Household Composition and Disability"="SVI_RPL_THEME2_HH_DISB",
                          "Socioeconomic status "="SVI_RPL_THEME1_SOCIECO",
                          "Housing Type and Transport"="SVI_RPL_THEME4_HH_TRANS")) %>%
  ggplot(aes(x = var, y = Mean_zscore, fill = Cluster)) + geom_col() +
  #geom_point(aes(y = Mean_zscore), size = 1) +
  facet_wrap(~ Cluster) + theme_bw() + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.text.y=element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.caption = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, size = 0.2) +
  labs(x = "Stressor",
       y = "Z-score Standardized Mean")+
  scale_fill_brewer(palette = "Set2", name = "Cluster") + 
  theme(legend.position = "none")

plot(plot_chem_means_zscore)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_clustersbargraph_zscoredvalues.png")), width = 8, height = 10, units = "in", pointsize = 12, res = 800)
plot(plot_chem_means_zscore)
dev.off()

#PCA based plot
pca_clusterplot <- fviz_cluster(km4, data = county_vulnerability,
                                ggtheme = theme_minimal(),
                                ellipse.type = "convex",
                                palette= "Set2",
                                main = "",
                                geom.var = c("point", "text"),
                                repel=TRUE, max.overlaps=Inf)
plot(pca_clusterplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_PCAclusterplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 800)
plot(pca_clusterplot)
dev.off()

#Map with coloring by clusters
#dataframe with counties and their clusters
countycluster <- cbind(county_vulnerability, cluster = km4$cluster) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="county") %>% 
  mutate(county=as.character(str_trim(county)))

formap <- join(countycluster, census, by="county")

Cluster_map <-
  formap %>% 
  mutate(cluster=as.factor(cluster)) %>% 
  ggplot(aes(fill = cluster, geometry=geometry)) + 
  ggthemes::theme_map() + 
  geom_sf(color = "black", size=0.1) +
  scale_fill_brewer(palette = "Set2", name="cluster", direction = 1,
                    labels=c("1", "2","3","4","5")) +
  theme(legend.position="right") 
plot(Cluster_map)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_mapofNC_withclusters.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(Cluster_map)
dev.off()


#################################################################################################
#################################################################################################
#### Step 14: FOR VULNERABILTY VARIABLES ONLY Run one way ANOVAs to compare differences in mean values of metals and SVIs between clusters 
#################################################################################################
#################################################################################################

#Null hypothesis: the means of the different groups are the same
#Alternative hypothesis: At least one sample mean is not equal to the others.

km4_anovas <- countycluster %>% select(-county) %>% 
  mutate(cluster=as.factor(cluster))

results.anova <- data.frame()
for (i in colnames(km4_anovas)[1:4]) {
  print(i)
  km4_anovas$var <- as.numeric(km4_anovas[[i]])
  anova <- aov(var ~ cluster, data = km4_anovas)
  output.anova <- summary(anova)
  results <- output.anova[[1]] %>% 
    dplyr::rename(p.ftest = "Pr(>F)") %>% 
    mutate(var=i,
           p.ftest.adjust=p.adjust(p.ftest,method="BH",n=4)) #p adjustment for conducting 4 different anovas
  results.anova <- rbind(results.anova, results) 
}
colnames(results.anova)

write.csv(results.anova, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_ANOVAresults.csv"), row.names= TRUE)


#################################################################################################
#################################################################################################
#### Step 15: FOR VULNERABILTY VARIABLES ONLY export datasheet with each county, raw and z-score  and cluster  
#################################################################################################
#################################################################################################

countycluster <- countycluster %>% select(county,cluster)
write.csv(countycluster, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_clusterassign.csv"), row.names= TRUE)

#################################################################################################
#################################################################################################
#### Step 16: FOR RESOURCES VARIABLES ONLY K-means clustering algorithm 
#################################################################################################
#################################################################################################


# Reset output folder to a folder containing only results based on vulnerability variables only analysis
Output_Folder <- ("/Users/lauren_eaves/IEHS Dropbox/lauren Eaves/1_Projects/13_DMAC/2_Aim1/2_correlations_clustering/output/resources")

county_resources <- county %>% 
  select(all_of(c(resourcesvars)))

#K-means with K = 1
km1 <- kmeans(county_resources, 1, nstart = 20) #20 iterations, 1 cluster
km1$size 
head(km1$centers) ## mean of each chemical in the cluster 
km1$totss ## total sum of squares
km1$withinss ## within cluster sum of squares by cluster- within cluster variation for each cluster 
km1$tot.withinss ## total within cluster sum of squares -- this is what we want to minimize!
km1$betweenss ## between cluster sum of squares
Percent_bt_k1 = paste(round(100*(km1$betweenss/km1$totss),1), "%", sep = "")

#K-means with K = 2
km2 <- kmeans(county_resources, 2, nstart = 20) 
km2$size 
head(km2$centers) #now seeing some separation of the means by the different clusters 
km2$totss 
km2$withinss
km2$tot.withinss #minimize this
km2$betweenss #now have a between cluster variation, maximize this 
Percent_bt_k2 = paste(round(100*(km2$betweenss/km2$totss),1), "%", sep = "")

#K-means with K = 3
km3 <- kmeans(county_resources, 3, nstart = 20)
km3$size
head(km3$centers)
km3$totss 
km3$withinss 
km3$tot.withinss  
km3$betweenss 
Percent_bt_k3 = paste(round(100*(km3$betweenss/km3$totss),1), "%", sep = "")

#K-means with K = 4
km4 <- kmeans(county_resources, 4, nstart = 20)
km4$size
head(km4$centers)
km4$totss 
km4$withinss 
km4$tot.withinss  
km4$betweenss 
Percent_bt_k4 = paste(round(100*(km4$betweenss/km4$totss),1), "%", sep = "")

#K-means with K = 5
km5 <- kmeans(county_resources, 5, nstart = 20)
km5$size
head(km5$centers)
km5$totss 
km5$withinss 
km5$tot.withinss  
km5$betweenss 
Percent_bt_k5 = paste(round(100*(km5$betweenss/km5$totss),1), "%", sep = "")

#K-means with K = 6
km6 <- kmeans(county_resources, 6, nstart = 20)
km6$size
head(km6$centers)
km6$totss 
km6$withinss 
km6$tot.withinss  
km6$betweenss 
Percent_bt_km6 = paste(round(100*(km6$betweenss/km6$totss),1), "%", sep = "")

#K-means with K = 7
km7 <- kmeans(county_resources, 7, nstart = 20)
km7$size
head(km7$centers)
km7$totss 
km7$withinss 
km7$tot.withinss  
km7$betweenss 
Percent_bt_k7 = paste(round(100*(km7$betweenss/km7$totss),1), "%", sep = "")

#K-means with K = 8
km8 <- kmeans(county_resources, 8, nstart = 20)
km8$size
head(km8$centers)
km8$totss 
km8$withinss 
km8$tot.withinss  
km8$betweenss 
Percent_bt_k8 = paste(round(100*(km8$betweenss/km8$totss),1), "%", sep = "")

#K-means with K = 9
km9 <- kmeans(county_resources, 9, nstart = 20)
km9$size
head(km9$centers)
km9$totss 
km9$withinss 
km9$tot.withinss  
km9$betweenss 
Percent_bt_km9 = paste(round(100*(km9$betweenss/km9$totss),1), "%", sep = "")

#K-means with K = 10
km10 <- kmeans(county_resources, 10, nstart = 20)
km10$size
head(km10$centers)
km10$totss 
km10$withinss 
km10$tot.withinss  
km10$betweenss 
Percent_bt_k10 = paste(round(100*(km10$betweenss/km10$totss),1), "%", sep = "")
# create an empty data frame
km.res <- data.frame(matrix(NA, 99, 4))

# use a loop to run k-means for 0 - 30 clusters and generate variance values
for (k in 1:99) {
  set.seed(123)
  kk <- kmeans(county, k, nstart = 20)
  kk$size <- toString(kk$size)
  km.res[k, ] <- cbind(k, kk$tot.withinss, kk$totss, kk$size) 
}


# name the variables created
names(km.res) <- c("clusters", "WithinSS", "TotSS", "ClusterSizes")

# create a variable for within cluster variance/total SS
km.res$WithinSS <- as.numeric(km.res$WithinSS)
km.res$TotSS <-as.numeric(km.res$TotSS)
km.res$clusters <-as.numeric(km.res$clusters)
km.res$PropWithin <- 100*(km.res$WithinSS/km.res$TotSS)

#output K-means results for all possible k 
write.csv(km.res, paste0(Output_Folder,"/", cur_date, "_Kmeans_allks_results.csv"), row.names= TRUE)



#################################################################################################
#################################################################################################
#### Step 17: FOR RESOURCES VARIABLES ONLY K-means clustering algorithm: decide on optimal k through balancing the number of single
#### county clusters (don't want any) and the proportion of within cluster variance (want this to be minimal)
#################################################################################################
#################################################################################################

gapstatplot <- fviz_nbclust(county_resources, kmeans, method = "gap_stat", nboot = 500, k.max = 5)+
  labs(subtitle = "Gap statistic method")
plot(gapstatplot)
silh <- fviz_nbclust(county_resources, kmeans, method = "silhouette", nboot = 500, k.max =5)+
  labs(subtitle = "Silhouette method")
plot(silh)
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_gapstatisticplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(gapstatplot)
dev.off()
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_silhouetteplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(silh)
dev.off()
NbClust(data = county_resources, diss = NULL, distance = "euclidean", min.nc = 3, max.nc = 5, 
        method = "kmeans")

# compare the within cluster variance/total SS to number of clusters
km.res %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_line() +
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
propwithinplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_point() + 
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=100,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))
plot(propwithinplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_propwithin_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(propwithinplot)
dev.off()

#compare number of single county clusters 
km.res <- km.res %>% 
  mutate(SinglecountyClusters = str_count(ClusterSizes, regex("\\s[1],\\s")))
km.res %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_line() +
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
singleclustersplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_point() + 
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=12,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))+
  scale_y_continuous(breaks=seq(0,12,2))
plot(singleclustersplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_singleclusters_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(singleclustersplot)
dev.off()

#################################################################################################
#################################################################################################
#### Step 18: FOR RESOURCES VARIABLES ONLY K-means clustering algorithm: vizualizations of k=4 solution
#################################################################################################
#################################################################################################

#creating dataframe with the mean of each cluster for each metal in z-score transformed value
km4_centers <- as.data.frame.matrix(t(km4$centers)) #dataframe containing the cluster means for each metal
km4_centers$var <- row.names(km4_centers)
km4_centers <- km4_centers %>% 
  gather(key = "Cluster", value = "mean", -var) %>% 
  dplyr::rename(Mean_zscore = mean)

plot_chem_means_zscore <- km4_centers %>% 
  mutate(Cluster = as.factor(Cluster),
         Cluster = fct_recode(Cluster, "Cluster 1" = "1",
                              "Cluster 2" = "2",
                              "Cluster 3" = "3",
                              "Cluster 4" = "4"),
         var = fct_recode(var, 
                          "Low Social resources" ="resources_social",
                          "Low Health resources" = "resources_health",
                          "Low Information resources" = "resources_info")) %>%
  ggplot(aes(x = var, y = Mean_zscore, fill = Cluster)) + geom_col() +
  #geom_point(aes(y = Mean_zscore), size = 1) +
  facet_wrap(~ Cluster) + theme_bw() + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.text.y=element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.caption = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, size = 0.2) +
  labs(x = "Stressor",
       y = "Z-score Standardized Mean")+
  scale_fill_brewer(palette = "Set2", name = "Cluster") + 
  theme(legend.position = "none")

plot(plot_chem_means_zscore)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_clustersbargraph_zscoredvalues.png")), width = 8, height = 10, units = "in", pointsize = 12, res = 800)
plot(plot_chem_means_zscore)
dev.off()

#PCA based plot
pca_clusterplot <- fviz_cluster(km4, data = county_resources,
                                ggtheme = theme_minimal(),
                                ellipse.type = "convex",
                                palette= "Set2",
                                main = "",
                                geom.var = c("point", "text"),
                                repel=TRUE, max.overlaps=Inf)
plot(pca_clusterplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_PCAclusterplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 800)
plot(pca_clusterplot)
dev.off()

#Map with coloring by clusters
#dataframe with counties and their clusters
countycluster <- cbind(county_resources, cluster = km4$cluster) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="county") %>% 
  mutate(county=as.character(str_trim(county)))

formap <- join(countycluster, census, by="county")

Cluster_map <-
  formap %>% 
  mutate(cluster=as.factor(cluster)) %>% 
  ggplot(aes(fill = cluster, geometry=geometry)) + 
  ggthemes::theme_map() + 
  geom_sf(color = "black", size=0.1) +
  scale_fill_brewer(palette = "Set2", name="cluster", direction = 1,
                    labels=c("1", "2","3","4","5")) +
  theme(legend.position="right") 
plot(Cluster_map)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_mapofNC_withclusters.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(Cluster_map)
dev.off()


#################################################################################################
#################################################################################################
#### Step 19: FOR RESOURCES VARIABLES ONLY Run one way ANOVAs to compare differences in mean values of metals and SVIs between clusters 
#################################################################################################
#################################################################################################

#Null hypothesis: the means of the different groups are the same
#Alternative hypothesis: At least one sample mean is not equal to the others.

km4_anovas <- countycluster %>% select(-county) %>% 
  mutate(cluster=as.factor(cluster))

results.anova <- data.frame()
for (i in colnames(km4_anovas)[1:3]) {
  print(i)
  km4_anovas$var <- as.numeric(km4_anovas[[i]])
  anova <- aov(var ~ cluster, data = km4_anovas)
  output.anova <- summary(anova)
  results <- output.anova[[1]] %>% 
    dplyr::rename(p.ftest = "Pr(>F)") %>% 
    mutate(var=i,
           p.ftest.adjust=p.adjust(p.ftest,method="BH",n=3)) #p adjustment for conducting 3 different anovas
  results.anova <- rbind(results.anova, results) 
}
colnames(results.anova)

write.csv(results.anova, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_social_and_metal_ANOVAresults.csv"), row.names= TRUE)


#################################################################################################
#################################################################################################
#### Step 20: FOR RESOURCES VARIABLES ONLY export datasheet with each county, raw and z-score  and cluster  
#################################################################################################
#################################################################################################

countycluster <- countycluster %>% select(county,cluster)
write.csv(countycluster, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_clusterassign.csv"), row.names= TRUE)

#################################################################################################
#################################################################################################
#### Step 21: FOR METALS VARIABLES ONLY K-means clustering algorithm 
#################################################################################################
#################################################################################################


# Reset output folder to a folder containing only results based on vulnerability variables only analysis
Output_Folder <- ("/Users/lauren_eaves/IEHS Dropbox/lauren Eaves/1_Projects/13_DMAC/2_Aim1/2_correlations_clustering/output/metals")

county_metals <- county %>% 
  select(all_of(c(metalsvars)))


#K-means with K = 1
km1 <- kmeans(county_metals, 1, nstart = 20) #20 iterations, 1 cluster
km1$size 
head(km1$centers) ## mean of each chemical in the cluster 
km1$totss ## total sum of squares
km1$withinss ## within cluster sum of squares by cluster- within cluster variation for each cluster 
km1$tot.withinss ## total within cluster sum of squares -- this is what we want to minimize!
km1$betweenss ## between cluster sum of squares
Percent_bt_k1 = paste(round(100*(km1$betweenss/km1$totss),1), "%", sep = "")

#K-means with K = 2
km2 <- kmeans(county_metals, 2, nstart = 20) 
km2$size 
head(km2$centers) #now seeing some separation of the means by the different clusters 
km2$totss 
km2$withinss
km2$tot.withinss #minimize this
km2$betweenss #now have a between cluster variation, maximize this 
Percent_bt_k2 = paste(round(100*(km2$betweenss/km2$totss),1), "%", sep = "")

#K-means with K = 3
km3 <- kmeans(county_metals, 3, nstart = 20)
km3$size
head(km3$centers)
km3$totss 
km3$withinss 
km3$tot.withinss  
km3$betweenss 
Percent_bt_k3 = paste(round(100*(km3$betweenss/km3$totss),1), "%", sep = "")

#K-means with K = 4
km4 <- kmeans(county_metals, 4, nstart = 20)
km4$size
head(km4$centers)
km4$totss 
km4$withinss 
km4$tot.withinss  
km4$betweenss 
Percent_bt_k4 = paste(round(100*(km4$betweenss/km4$totss),1), "%", sep = "")

#K-means with K = 5
km5 <- kmeans(county_metals, 5, nstart = 20)
km5$size
head(km5$centers)
km5$totss 
km5$withinss 
km5$tot.withinss  
km5$betweenss 
Percent_bt_k5 = paste(round(100*(km5$betweenss/km5$totss),1), "%", sep = "")

#K-means with K = 6
km6 <- kmeans(county_metals, 6, nstart = 20)
km6$size
head(km6$centers)
km6$totss 
km6$withinss 
km6$tot.withinss  
km6$betweenss 
Percent_bt_km6 = paste(round(100*(km6$betweenss/km6$totss),1), "%", sep = "")

#K-means with K = 7
km7 <- kmeans(county_metals, 7, nstart = 20)
km7$size
head(km7$centers)
km7$totss 
km7$withinss 
km7$tot.withinss  
km7$betweenss 
Percent_bt_k7 = paste(round(100*(km7$betweenss/km7$totss),1), "%", sep = "")

#K-means with K = 8
km8 <- kmeans(county_metals, 8, nstart = 20)
km8$size
head(km8$centers)
km8$totss 
km8$withinss 
km8$tot.withinss  
km8$betweenss 
Percent_bt_k8 = paste(round(100*(km8$betweenss/km8$totss),1), "%", sep = "")

#K-means with K = 9
km9 <- kmeans(county_metals, 9, nstart = 20)
km9$size
head(km9$centers)
km9$totss 
km9$withinss 
km9$tot.withinss  
km9$betweenss 
Percent_bt_km9 = paste(round(100*(km9$betweenss/km9$totss),1), "%", sep = "")

#K-means with K = 10
km10 <- kmeans(county_metals, 10, nstart = 20)
km10$size
head(km10$centers)
km10$totss 
km10$withinss 
km10$tot.withinss  
km10$betweenss 
Percent_bt_k10 = paste(round(100*(km10$betweenss/km10$totss),1), "%", sep = "")
# create an empty data frame
km.res <- data.frame(matrix(NA, 99, 4))


# create an empty data frame
km.res <- data.frame(matrix(NA, 99, 4))

# use a loop to run k-means for 0 - 30 clusters and generate variance values
for (k in 1:99) {
  set.seed(123)
  kk <- kmeans(county_metals, k, nstart = 20)
  kk$size <- toString(kk$size)
  km.res[k, ] <- cbind(k, kk$tot.withinss, kk$totss, kk$size) 
}


# name the variables created
names(km.res) <- c("clusters", "WithinSS", "TotSS", "ClusterSizes")

# create a variable for within cluster variance/total SS
km.res$WithinSS <- as.numeric(km.res$WithinSS)
km.res$TotSS <-as.numeric(km.res$TotSS)
km.res$clusters <-as.numeric(km.res$clusters)
km.res$PropWithin <- 100*(km.res$WithinSS/km.res$TotSS)

#output K-means results for all possible k 
write.csv(km.res, paste0(Output_Folder,"/", cur_date, "_Kmeans_allks_results.csv"), row.names= TRUE)



#################################################################################################
#################################################################################################
#### Step 22: FOR METALS VARIABLES ONLY K-means clustering algorithm: decide on optimal k through balancing the number of single
#### county clusters (don't want any) and the proportion of within cluster variance (want this to be minimal)
#################################################################################################
#################################################################################################

gapstatplot <- fviz_nbclust(county_metals, kmeans, method = "gap_stat", nboot = 500, k.max = 5)+
  labs(subtitle = "Gap statistic method")
plot(gapstatplot)
silh <- fviz_nbclust(county_metals, kmeans, method = "silhouette", nboot = 500, k.max =5)+
  labs(subtitle = "Silhouette method")
plot(silh)
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_gapstatisticplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(gapstatplot)
dev.off()
png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_silhouetteplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(silh)
dev.off()
NbClust(data = county_metals, diss = NULL, distance = "euclidean", min.nc = 3, max.nc = 5, 
        method = "kmeans")

# compare the within cluster variance/total SS to number of clusters
km.res %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_line() +
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
propwithinplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_point() + 
  labs(y = "Proportion of Within over Total SS (%)", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=100,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))
plot(propwithinplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_propwithin_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(propwithinplot)
dev.off()

#compare number of single county clusters 
km.res <- km.res %>% 
  mutate(SinglecountyClusters = str_count(ClusterSizes, regex("\\s[1],\\s")))
km.res %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_line() +
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")
# look at a smaller subset of data
singleclustersplot <- km.res[1:30,]  %>% 
  ggplot(aes(x = clusters, y = SinglecountyClusters)) + geom_point() + 
  labs(y = "Number of Single county Clusters", 
       x = "Number of Clusters")+
  geom_vline(xintercept=4, colour="red")+
  annotate(x=4,y=12,label="K=4",vjust=2,geom="label")+
  scale_x_continuous(breaks=seq(0,30,5))+
  scale_y_continuous(breaks=seq(0,12,2))
plot(singleclustersplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_singleclusters_vs_numberk.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(singleclustersplot)
dev.off()

#################################################################################################
#################################################################################################
#### Step 23: FOR METALS VARIABLES ONLY K-means clustering algorithm: vizualizations of k=4 solution
#################################################################################################
#################################################################################################

#creating dataframe with the mean of each cluster for each metal in z-score transformed value
km4_centers <- as.data.frame.matrix(t(km4$centers)) #dataframe containing the cluster means for each metal
km4_centers$var <- row.names(km4_centers)
km4_centers <- km4_centers %>% 
  mutate(Group = ifelse(var %in% metalsvars,"Metal",
                        ifelse(var %in% socialvars, "Vulnerability", "Resources"))) %>% 
  gather(key = "Cluster", value = "mean", -var, -Group) %>% 
  dplyr::rename(Mean_zscore = mean)

plot_chem_means_zscore <- km4_centers %>% 
  mutate(Cluster = as.factor(Cluster),
         Cluster = fct_recode(Cluster, "Cluster 1" = "1",
                              "Cluster 2" = "2",
                              "Cluster 3" = "3",
                              "Cluster 4" = "4"),
         var = fct_recode(var, "Arsenic"="Arsenic.Mean_avg",
                          "Lead"="Lead.Mean_avg",
                          "Manganese"="Manganese.Mean_avg")) %>% 
  ggplot(aes(x = var, y = Mean_zscore, fill = Cluster)) + geom_col() +
  #geom_point(aes(y = Mean_zscore), size = 1) +
  facet_wrap(~ Cluster) + theme_bw() + 
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.text.y=element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.caption = element_text(size = 12, hjust = 0),
        strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, size = 0.2) +
  labs(x = "Stressor",
       y = "Z-score Standardized Mean")+
  scale_fill_brewer(palette = "Set2", name = "Cluster") + 
  theme(legend.position = "none")

plot(plot_chem_means_zscore)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_clustersbargraph_zscoredvalues.png")), width = 8, height = 10, units = "in", pointsize = 12, res = 800)
plot(plot_chem_means_zscore)
dev.off()

#PCA based plot
pca_clusterplot <- fviz_cluster(km4, data = county_metals,
                                ggtheme = theme_minimal(),
                                ellipse.type = "convex",
                                palette= "Set2",
                                main = "",
                                geom.var = c("point", "text"),
                                repel=TRUE, max.overlaps=Inf)
plot(pca_clusterplot)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_PCAclusterplot.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 800)
plot(pca_clusterplot)
dev.off()

#Map with coloring by clusters
#dataframe with counties and their clusters
countycluster <- cbind(county_metals, cluster = km4$cluster) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="county") %>% 
  mutate(county=as.character(str_trim(county)))

formap <- join(countycluster, census, by="county")

Cluster_map <-
  formap %>% 
  mutate(cluster=as.factor(cluster)) %>% 
  ggplot(aes(fill = cluster, geometry=geometry)) + 
  ggthemes::theme_map() + 
  geom_sf(color = "black", size=0.1) +
  scale_fill_brewer(palette = "Set2", name="cluster", direction = 1,
                    labels=c("1", "2","3","4","5")) +
  theme(legend.position="right") 
plot(Cluster_map)

png(file = (paste0(Output_Folder,"/", cur_date, "_Kmeans_mapofNC_withclusters.png")), width = 10, height = 10, units = "in", pointsize = 12, res = 300)
plot(Cluster_map)
dev.off()


#################################################################################################
#################################################################################################
#### Step 24: FOR METALS VARIABLES ONLY Run one way ANOVAs to compare differences in mean values of metals and SVIs between clusters 
#################################################################################################
#################################################################################################

#Null hypothesis: the means of the different groups are the same
#Alternative hypothesis: At least one sample mean is not equal to the others.

km4_anovas <- countycluster %>% select(-county) %>% 
  mutate(cluster=as.factor(cluster))

results.anova <- data.frame()
for (i in colnames(km4_anovas)[1:3]) {
  print(i)
  km4_anovas$var <- as.numeric(km4_anovas[[i]])
  anova <- aov(var ~ cluster, data = km4_anovas)
  output.anova <- summary(anova)
  results <- output.anova[[1]] %>% 
    dplyr::rename(p.ftest = "Pr(>F)") %>% 
    mutate(var=i,
           p.ftest.adjust=p.adjust(p.ftest,method="BH",n=3)) #p adjustment for conducting 3 different anovas
  results.anova <- rbind(results.anova, results) 
}
colnames(results.anova)

write.csv(results.anova, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_social_and_metal_ANOVAresults.csv"), row.names= TRUE)


#################################################################################################
#################################################################################################
#### Step 25: FOR METALS VARIABLES ONLY export datasheet with each county, raw and z-score  and cluster  
#################################################################################################
#################################################################################################

countycluster <- countycluster %>% select(county,cluster)
write.csv(countycluster, paste0(Output_Folder,"/", cur_date, "_DMAC_NCcounties_clusterassign.csv"), row.names= TRUE)

