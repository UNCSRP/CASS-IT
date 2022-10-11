# DMAC 
#######

getwd()

rm(list=ls())
library(dplyr)
library(readr)
library(ggcorrplot)
library(tidyverse)
library(tidyLPA)
library(janitor)


# LOAD DATASET #####
glimpse(dmac)

# Run and compare potential optimal model ======

dmac <- 
    dmac %>% 
    mutate(ars =  (Arsenic.Mean_avg - mean(Arsenic.Mean_avg))/sd(Arsenic.Mean_avg),
           lead = (Lead.Mean_avg - mean(Lead.Mean_avg))/sd(Lead.Mean_avg),
           mang = (Manganese.Mean_avg - mean(Manganese.Mean_avg))/sd(Manganese.Mean_avg),
           svi1 = (SVI_RPL_THEME1_SOCIECO - mean(SVI_RPL_THEME1_SOCIECO))/sd(SVI_RPL_THEME1_SOCIECO),
           svi2 = (SVI_RPL_THEME2_HH_DISB - mean(SVI_RPL_THEME2_HH_DISB))/sd(SVI_RPL_THEME2_HH_DISB),
           svi3 = (SVI_RPL_THEME3_MINO - mean(SVI_RPL_THEME3_MINO))/sd(SVI_RPL_THEME3_MINO),
           svi4 = (SVI_RPL_THEME4_HH_TRANS - mean(SVI_RPL_THEME4_HH_TRANS))/sd(SVI_RPL_THEME4_HH_TRANS),
           res1 = (resources_social - mean(resources_social))/sd(resources_social),
           res2 = (resources_health - mean(resources_health))/sd(resources_health),
           res3 = (resources_info - mean(resources_info))/sd(resources_info))
    

dmac_svi_metals <- dmac %>% 
    select(-county) %>% 
    estimate_profiles(1:6) 

## After running potential models, compare model fit indices
fit_indices <- get_fit(dmac_svi_metals)
compare_solutions(dmac_svi_metals, statistics = "BIC")


# Focus on 4-profile ====

m4_svi_metals <- dmac %>% 
    estimate_profiles(4, 
                      select_vars=c("ars", 
                                    "lead", 
                                    "mang", 
                                    "res1",
                                    "res2",
                                    "res3",
                                    "svi1", 
                                    "svi2",
                                    "svi3",
                                    "svi4")
                                       )  

## Plot visually for 4-profile solution
plot_profiles(m4_svi_metals,
               variables = NULL, 
               bw = FALSE,   alpha_range = c(0, 0.1)) 
 
## Plot bar graphs of means of each class 
get_estimates(m4_svi_metals)
estimates_m4_svi_metals <-  get_estimates(m4_svi_metals) 

### focus only means 
means <- estimates_m4_svi_metals %>% 
        filter(Category == "Means") 
    

### Plot 4-profile in ggplot2 
library(ggplot2)
p <- ggplot(data = means, mapping = aes(x = Parameter, y = Estimate)) 

p + geom_col() +
    facet_wrap(~Class) +
    theme_bw() + 
    labs(x = "Stressor",
         ) + 
    theme(
          axis.text.y = element_text(size = 4,  color = "black"),
          axis.text.x = element_text(size = 4,  color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          axis.title.y = element_text(size = 8, color = "black")) -> plot_final





