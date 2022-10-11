# DMAC 
## tract level

rm(list=ls())
library(dplyr)
library(readr)
library(tidyverse)
library(tidyLPA)
library(janitor)


# LOAD DATASET #####
glimpse(dmac)



# Run and compare potential optimal model ======
names(dmac)

# create z scores 
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
    select(-tract) %>% 
    estimate_profiles(1:6) 

# Focus on 3-profile ====

m3_svi_metals <- dmac %>% 
    estimate_profiles(3, 
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

## Plot visually for 3-profile solution
 plot_profiles(m3_svi_metals,
 variables = NULL, ci = 0.95, sd = TRUE,
 add_line = TRUE,
 rawdata = TRUE,  bw = FALSE,   alpha_range = c(0, 0.1))
 
## Plot bar graphs of means of each class 
get_estimates(m3_svi_metals)
estimates_m3_svi_metals <-  get_estimates(m3_svi_metals)

### focus only means 
means <- estimates_m3_svi_metals %>% 
        filter(Category == "Means") 
    
### Plot 3-profile in ggplot2 
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

getwd()
### Save as picture
ggsave("tract_level/lpa_plot_m3.png", plot = plot_final, 
       type = 'cairo', width = 4, height = 4, dpi = 300, 
       units = "in", bg = "#ffffff")    





