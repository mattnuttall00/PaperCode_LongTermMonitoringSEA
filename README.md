# PaperCode_LongTermMonitoringSEA
Code for the manuscript submitted to Conservation Science and Practice (2021) - Long-term monitoring of wildlife populations for protected area management in Southeast Asia

## Overview
The R scripts here are the full scripts used for all the different analyses (details below). However, none of the scripts will run as we have not provided the raw data that are required for the analyses. This is because the data are sensitive - they contain the locations of rare and endangered species. The metadata are available on the Global Biodiversity Information Facility, and can be found here: https://doi.org/10.15468/37thhj. If you would like access to the raw data, please contact Olly Griffin from WCS Cambodia on ogriffin@wcs.org. Alternatively please get in touch with WCS Cambodia.

The analyses presented in the paper were conducted by Matt Nuttall, with support from Olly Griffin, Rachel Fewster, and Nils Bunnefeld. For any questions relating to the analyses, please contact Matt Nuttall (mattnuttall00@gmail.com).

## Analyses

### Annual population estimates
For the script used to produce the population point estimates, see the script "Conventional_DistSamp_abundance_estimates.R". This anlysis uses a conventional distance sampling framework. There is a second script in the folder titled "Analytical_approach_test_CDS.R". This script tested the incusion and exclusion of the geographical priority strata that are a feature in some of the years.

### Population trends
For the script that was used to produce the temporal population trends, see the script "Population_trends.R". This anlysis uses the detection function models from the above point estimate analysis, followed by Generalized Additive Models and bootstrapping to estimate long-term population trends. There is a second script in the folder titled "Rfunc.R". This script was written by Prof Rachel Fewster (co-author), and edited by Matt Nuttall. This script is required to run the Population_trends.R script. 

### Spatial models
For the script used to produce the relative abundance maps, see the script "Density_surface_models.R". 
