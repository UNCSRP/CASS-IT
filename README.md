# CASS-IT

This repository contains code and scripts correspond to the manuscript: "Generation of the Chemical and Social Stressors Integration Technique (CASS-IT) to identify areas of holistic public health concern: an application to North Carolina"

Abstract: Due to structural racism and income inequality, exposure to environmental chemicals is tightly linked to socioeconomic factors. In addition, exposure to psychosocial stressors, such as racial discrimination, as well as having limited resources, can increase susceptibility to environmentally induced disease. Yet, studies are often conducted separately in fields of social science and environmental science, reducing the potential for holistic risk estimates. To tackle this gap, we developed the Chemical and Social Stressors Integration Technique (CASS-IT) to integrate environmental chemical and social stressor datasets. The CASS-IT provides a framework to identify distinct geographic areas based on combinations of environmental chemical exposure, social vulnerability, and access to resources. It incorporates two data dimension reduction tools: k-means clustering and latent profile analysis. Here, the CASS-IT was applied to North Carolina (NC) as a case study. Environmental chemical data included toxic metals – arsenic, manganese, and lead – in private drinking well water. Social stressor data was captured by the CDC’s social vulnerability index’s four domains: socioeconomic status, household composition and disability, minority status and language, and housing type and transportation. Data on resources were derived from Federal Emergency Management Agency (FEMA’s) Resilience and Analysis Planning Tool, which generated measures of health resources, social resources, and information resources. The results highlighted 31 NC counties where exposure to both toxic metals and social stressors are elevated, and health resources are minimal; these are counties in which environmental justice is of utmost concern. A census-tract level analysis was also conducted to demonstrate the utility of CASS-IT at different geographical scales. The tract-level analysis highlighted specific tracts within counties of concern that are particularly high priority. In future research, the CASS-IT can be used to analyze United-States wide environmental datasets providing guidance for targeted public health interventions and reducing environmental disparities. 

# Scripts contained herein:
1. CASS-IT_kmeans_county_NC.R: k-means analysis at a county level 
2. CASS-IT_kmeans_tract_NC.R: k-means analysis at a census tract level 
3. CASS-IT_LPA_county_NC.R: LPA at a county level
4. CASS-IT_LPA_tract_NC.R: LPA at a census tract level 
# Data contained herein:
1. NCcounties_metals_vulnerability_resources: county level data used in scripts above
2. NCtracts_metals_vulnerability_resources: census tract level data used in scripts above



