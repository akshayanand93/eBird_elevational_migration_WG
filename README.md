# Avian Elevational Migration in the Western Ghats of India

This repository contains the code and associated files for a manuscript exploring the elevational migration strategies of birds in the Western Ghats of India. Below is a detailed breakdown of the repository's structure.

---

## Scripts

### `01_filtering_ebird_data.R`
- Cleans and filters the eBird Basic Dataset (EBD) downloaded from [eBird](https://science.ebird.org/en/use-ebird-data/download-ebird-data-products).
- Applies spatial, temporal, and sampling effort filters.
- Extracts elevation data for each occurrence record.

### `02_bias_correction_weighted_bootstrap_estimates.R`
- Corrects sampling effort biases in the filtered eBird data.
- Uses fitted GAMs and weighted resampling to equate sampling effort across seasons and elevations.

### `03_method_validation.R`
- Validates the bias correction method against published methods from [Tsai et al. (2021)](https://nsojournals.onlinelibrary.wiley.com/doi/10.1111/ecog.05196).
- Data for this script is not included in the repository but is available in the supporting information of the referenced paper.

### `04_climatic_vars_prep.R`
- Calculates species’ tolerance to temperature, precipitation, and wind speed using occurrence records and CHELSA climate data.
- Generates tolerance estimates (2.5%, 50%, 97.5% quantiles) for these variables.
- Note: CHELSA data is not included due to size constraints but can be downloaded from [CHELSA](https://chelsa-climate.org/downloads/).

### `05_data_prep_phyloglm.R`
- Prepares climatic tolerance and species trait data along with median elevation estimates.
- Fits seasonal Phylogenetic Logistic Regression (PhyloGLM) models to analyze the effects of species traits on elevational migration status.

### `06_pgls_models_distance.R`
- Fits seasonal Phylogenetic Generalized Least Squares (PGLS) models to examine the effect of species traits on elevational migration distances.

### `07_figs.R`
- Generates the figures accompanying the manuscript.

---

## Data

### Included Files
- `avonet_WG.csv`: Species trait data for the analyzed species pool.  
- `avonet-ebird.csv`: Full species trait dataset from [AVONET](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13898).  
- `ebird_taxonomy_v2022.csv`: Taxonomic classifications from eBird.  
- `India-Checklist_v8_2.csv`: Checklist of birds from the Indian subcontinent, used to identify Western Ghats endemics.  
- `resident_list.csv`: List of resident species used to build a consensus phylogenetic tree in `05_data_prep_phyloglm.R`.  
- `unique_chklst_locations_WG.csv`: Latitude and longitude of eBird checklists used in `01_filtering_ebird_data.R`.  

### Not Included (Requestable)
- `aster-gdem-wg.tiff`: ASTER GDEM raster file used to extract elevations in `01_filtering_ebird_data.R`.  
- `ebd_filtered_WG.txt`: Filtered eBird dataset created in `01_filtering_ebird_data.R`.  
- `ebd_IN_201303_202303_relAug-2023.txt`: Unfiltered eBird dataset.  
- `ebird_elev_residents_WG.csv`: Filtered dataset with elevation and seasonal occurrence data, created as output in `01_filtering_ebird_data.R`.  

---

## Output

### Data Files
- `cons_resident_tree_wg.tree`: Consensus phylogenetic tree of analyzed bird species, generated in `05_data_prep_phyloglm.R`.  
- `distance_models_results.csv`: Results of seasonal PGLS models for migration distance and traits, output of `06_pgls_models_distance.R`.  
- `eBird_climatic_vars_wg.csv`: Environmental tolerance estimates for analyzed bird species, output of `04_climatic_vars_prep.R`.  
- `elev_mig_climate_wg.csv`: Comprehensive dataset with elevational migration, traits, and environmental tolerances, created in `05_data_prep_phyloglm.R`.  
- `elev_mig_trait.csv`: Migration and species trait data, used to create `elev_mig_climate_wg.csv` and output from `05_data_prep_phyloglm.R`.  
- `elev_migration_values.csv`: Migration values for each species, output of `05_data_prep_phyloglm.R`.  
- `median_quantile_estimates.csv`: Seasonal median elevation estimates with confidence intervals, output of `02_bias_correction_weighted_bootstrap_estimates.R`.  
- `table_2_binom_mod_results.csv`: Results of PhyloGLM models for migration status, output of `05_data_prep_phyloglm.R`.  
- `table_mod_sel.csv`: Model selection table for seasonal PGLS models, output of `06_pgls_models_distance.R`.  

---

## Figures

All manuscript figures, including supplementary figures, are located in the `figs` folder. Figures are generated using `07_figs.R`.

---

## Shapefiles

- The `shp` folder contains the Western Ghats boundary shapefile used for spatial subsetting in `01_filtering_ebird_data.R`.

---

## Notes
- Data files not included in this repository due to size constraints can be requested from the repository owner.
