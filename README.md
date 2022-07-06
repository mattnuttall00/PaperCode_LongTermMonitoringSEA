# R code for the analyses within the paper titled "Long-term monitoring of wildlife populations for protected area management in Southeast Asia"

Code for the manuscript published in Conservation Science and Practice (2022) - Long-term monitoring of wildlife populations for protected area management in Southeast Asia. The publication can be found here: https://conbio.onlinelibrary.wiley.com/doi/full/10.1111/csp2.614

## Overview

Corresponding author: Matt Nuttall (mattnuttall00@gmail.com)

The R scripts here are the full scripts used for all the different analyses (details below). However, we have not provided the raw data that are required for the analyses in this repository. The metadata, plus raw data for primates, green peafowl, and red muntjac, are stored publicly on the Global Biodiversity Information Facility and can be found here: https://doi.org/10.15468/37thhj. Due to the potential risks of releasing precise location data for species that are locally rare, with limited distribution, that are possible targets for poachers, data for Gaur, Banteng, Eld's deer, and Sambar are not included in the public dataset, but are available on request from Olly Griffin / WCS Cambodia. 

These scripts are provided largely for the purposes of transparency and to hopefully be useful for researchers conducting similar analyses. If you would like access to any of the required input files, please contact either Matt Nuttall (mattnuttall00@gmail.com) or Olly Griffin (ogriffin@wcs.org) / Cain Agger (cagger@wcs.org). Up-to-date contact details for WCS Cambodia should be posted on the GBIF database linked above.


The analyses presented in the paper were conducted by Matt Nuttall, with support from Olly Griffin, Rachel Fewster, and Nils Bunnefeld. For any questions relating to the analyses, please contact Matt Nuttall (mattnuttall00@gmail.com).

## Analyses

### Annual population estimates
For the code used to produce the population point estimates, see the script "Conventional_DistSamp_abundance_estimates.R". This anlysis uses a conventional distance sampling framework. There is a second script in the folder titled "Analytical_approach_test_CDS.R". This script tested the incusion and exclusion of the geographical priority strata that are a feature in some of the years.

### Population trends
For the code used to produce the temporal population trends, see the script "Population_trends.R". This anlysis uses the detection function models from the above point estimate analysis, followed by Generalized Additive Models and bootstrapping to estimate long-term population trends. There is a second script in the folder titled "Rfunc.R". This script was written by Prof Rachel Fewster (co-author), and edited by Matt Nuttall. This script is required to run the Population_trends.R script. 

### Spatial models
For the code used to produce the relative abundance maps, see the script "Density_surface_models.R". 

## License

The contents of this repository are covered by the MIT license (Open Source Initiative compatible). Please see the "LICENSE.txt" file for details. 
