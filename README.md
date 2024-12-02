# eBird_elevational_migration_WG

This repository contains code for a manuscript exploring the elevational migration strategies of birds in the Western Ghats of India.

Below is a description of what each script (.R) in this repository contains:

- `01_filtering_ebird_data.R`:  
  This script is used to clean and filter the eBird Basic Dataset (EBD) downloaded from [eBird](https://science.ebird.org/en/use-ebird-data/download-ebird-data-products). This freely accessible dataset was used in all downstream analyses. Here, we apply spatial, temporal, and sampling effort filters to the EBD and extract elevation of each occurrence record.

- `02_bias_correction_weighted_bootstrap_estimates.R`:  
  This script was used to correct sampling effort biases in the filtered eBird data. We developed a novel approach that corrects these biases using fitted GAMs and weighted resampling to equate sampling effort across seasons and elevations.

- `03_method_validation.R`:  
  This script was used to validate our bias correction method against previously published methods used in Tsai et al. (2021). The data used in this script was not included in the "data" folder but can be found in the supporting information section of the paper mentioned above.

- `04_climatic_vars_prep.R`:  
  This script was used to calculate a species' tolerance to temperature, precipitation, and wind speed in the Western Ghats. We used mean monthly rasters from the CHELSA dataset and species occurrence records from our filtered eBird data to calculate the lower (2.5%), median (50%), and upper (97.5%) quantiles of the environmental variables of interest, quantifying a species' tolerance to these variables. The CHELSA data used in this script is not included in this repository due to its large size but can be downloaded from [CHELSA](https://chelsa-climate.org/downloads/).

- `05_data_prep_phyloglm.R`:  
  This script was used to prepare the climatic tolerance and species traits data along with median elevation estimates of species. These data were then used to fit seasonal Phylogenetic Logistic Regression models examining the effect of species traits on the elevational migration status (migrant/non-migrant) of a species.

- `06_pgls_models_distance.R`:  
  This script was used to fit seasonal Phylogenetic Generalized Least Squares models examining the effect of species traits on the elevational migration distance of a species.

- `07_figs.R`:  
  This script was used to create the figures accompanying the main text of the manuscript.
