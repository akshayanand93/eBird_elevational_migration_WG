#clear environment and perform garbage collection
rm(list = ls())
gc()

#load required libraries
library(terra)
library(sf)
library(data.table)
library(tidyverse)

#list raster files from 2008 - 2018
ppt_file_list <- list.files("chelsa/pr", pattern = "200[8-9]|201[0-8]", full.names = TRUE)
temp_file_list <- list.files("chelsa/tas", pattern = "200[8-9]|201[0-8]", full.names = TRUE)
wind_file_list <- list.files("chelsa/sfcWind", pattern = "200[8-9]|201[0-8]", full.names = TRUE)

#create raster layer
ppt_raster <- rast(ppt_file_list)
temp_raster <-rast(temp_file_list)
wind_raster <- rast(wind_file_list)

#stack rasters
ppt_stack <- c(ppt_raster)
temp_stack <- c(temp_raster)
wind_stack <- c(wind_raster)

#read in ebird data
ebird_occs <- read.csv("ebird_elev_residents_WG.csv")

#select useful columns
ebird_occs <- ebird_occs[, c(6, 7, 14, 15)]

#convert to terra vector
ebird_occs_vect <- vect(ebird_occs, geom = c("longitude", "latitude"), 
                        crs = "epsg:4326")

#extract
ppt_wg <- terra::extract(ppt_stack, ebird_occs_vect)
temp_wg <- terra::extract(temp_stack, ebird_occs_vect)
wind_wg <- terra::extract(wind_stack, ebird_occs_vect)

#removing id column
ppt_wg <- ppt_wg[, -1]
temp_wg <- temp_wg[, -1]
wind_wg <- wind_wg[, -1]

#calculate mean, lower and upper quantiles for each species per envt variable
#precip
#join and subset
ebird_precip <- cbind(ebird_occs, ppt_wg) %>% 
  select(c(1, 5:136))

#convert to data table
ebird_precip <- as.data.table(ebird_precip)

#calculate mean, lower, median, and upper (hq) quantiles
ebird_precip_val <- ebird_precip[, .(
  precip_mean = mean(unlist(.SD), na.rm = TRUE),
  precip_lq = quantile(unlist(.SD), probs = 0.025, na.rm = TRUE),
  precip_med = quantile(unlist(.SD), probs = 0.50, na.rm = TRUE),
  precip_hq = quantile(unlist(.SD), probs = 0.975, na.rm = TRUE)
), by = scientific_name, .SDcols = 2:133]

#convert to mm
ebird_precip_val$precip_mean <- ebird_precip_val$precip_mean/100
ebird_precip_val$precip_lq <- ebird_precip_val$precip_lq/100
ebird_precip_val$precip_med <- ebird_precip_val$precip_med/100
ebird_precip_val$precip_hq <- ebird_precip_val$precip_hq/100

#calculate tolerance ranges (breadth)
ebird_precip_val[, precip_breadth := precip_hq - precip_lq]

#temp
#join and subset
ebird_temp <- cbind(ebird_occs, temp_wg) %>% 
  select(c(1, 5:136))

#convert to data table
ebird_temp <- as.data.table(ebird_temp)

#calculate mean, lower, median, and upper (hq) quantiles
ebird_temp_val <- ebird_temp[, .(
  temp_mean = mean(unlist(.SD), na.rm = TRUE), 
  temp_lq = quantile(unlist(.SD), probs = 0.025, na.rm = TRUE), 
  temp_med = quantile(unlist(.SD), probs = 0.50, na.rm = TRUE), 
  temp_hq = quantile(unlist(.SD), probs = 0.975, na.rm = TRUE)
), by = scientific_name, .SDcols = 2:133]

#convert to celcius
ebird_temp_val$temp_mean <- (ebird_temp_val$temp_mean/10) - 273.15
ebird_temp_val$temp_lq <- (ebird_temp_val$temp_lq/10) - 273.15
ebird_temp_val$temp_med <- (ebird_temp_val$temp_med/10) - 273.15
ebird_temp_val$temp_hq <- (ebird_temp_val$temp_hq/10) - 273.15

#calculate breadth
ebird_temp_val[, temp_breadth := temp_hq - temp_lq]

#wind
#join and subset
ebird_wind <- cbind(ebird_occs, wind_wg) %>% 
  select(c(1, 5:136))

#convert to data table
ebird_wind <- as.data.table(ebird_wind)

#calculate mean, lower, median and upper (hq) quantiles
ebird_wind_val <- ebird_wind[, .(
  wind_mean = mean(unlist(.SD), na.rm = TRUE), 
  wind_lq = quantile(unlist(.SD), probs = 0.025, na.rm = TRUE), 
  wind_med = quantile(unlist(.SD), probs = 0.50, na.rm = TRUE), 
  wind_hq = quantile(unlist(.SD), probs = 0.975, na.rm = TRUE)
), by = scientific_name, .SDcols = 2:133]

#calculate breadth
ebird_wind_val[, wind_breadth := wind_hq - wind_lq]

#join
climatic_vars_wg <- left_join(ebird_precip_val, ebird_temp_val, by = "scientific_name") %>% 
  left_join(., ebird_wind_val, by = "scientific_name")

#adding "_" to scientific_name for joining
climatic_vars_wg$scientific_name <- str_replace(climatic_vars_wg$scientific_name, " ", "_")

#rounding climatic vars
climatic_vars_wg[, (names(climatic_vars_wg)[-1]) := lapply(.SD, round, 2), .SDcols = -1]

#save data
fwrite(climatic_vars_wg, "output/eBird_climatic_vars_wg.csv", row.names = FALSE)