#load required libraries
library(auk)
library(lubridate)
library(sf)
library(tidyverse)
library(ggplot2)
library(terra)
library(raster)
library(data.table)
library(readr)
library(stringr)

#resolve namespace conflicts
select <- dplyr::select

#read in western ghats shapefile for spatial filtering
wg <- read_sf("shp/wg.shp")

#specify input and output files
f_in <- "data/ebd_IN_201303_202303_relAug-2023.txt" 
f_out <- "data/ebd_filtered_WG.txt"

#define filters
ebird_filters <- f_in %>%
  #reference file
  auk_ebd() %>% 
  #spatial filter
  auk_bbox(wg) %>%
  #filtering for 10 years of data with 3 full seasons per year
  auk_date(date = c("2013-03-01", "2023-02-28")) %>% 
  #protocol filter
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
  #effort distance filter
  auk_distance(distance = c(0, 1)) %>% 
  #duration filter
  auk_duration(duration = c(0, 300)) %>% 
  #include only complete checklists where observers reported all species identified
  auk_complete()

#check filters
ebird_filters 

#run filters
auk_filter(ebird_filters, file = f_out, overwrite = TRUE)

#read in filtered file
ebird_data <- read_ebd("data/ebd_filtered_WG.txt")

#sub-setting the ebird data to include only passerine birds
#read in the ebird taxonomy file
taxonomy <- read.csv("data/ebird_taxonomy_v2022.csv")

#clean capitalization and column names for merging
new_col_names <- colnames(taxonomy) %>%
  str_replace("SCI_NAME", "scientific_name") %>% 
  str_to_lower()
setnames(taxonomy, new_col_names)
 
#merge with ebird data
ebird_data <- merge(ebird_data, taxonomy, by = "scientific_name")

#subset to remove raptors and aquatic species
ebird_data <- filter(ebird_data, order1 %in% c("Passeriformes", "Bucerotiformes",
                               "Trogoniformes", "Piciformes", "Coraciiformes", 
                               "Cuculiformes", "Psittaciformes", "Columbiformes"))

#swallows and kingfishers are still in the data set, filter to remove these groups
ebird_data <- ebird_data %>%
  filter(!family %in% c("Alcedinidae (Kingfishers)", "Hirundinidae (Swallows)"))

#additional filtering
ebird_data <- ebird_data %>% 
  #keeping only approved observations
  filter(reviewed == 0 | approved == 1) %>% 
  #keep only single sampling id
  mutate(sampling_id = ifelse(is.na(group_identifier), 
                              sampling_event_identifier, 
                              group_identifier)) %>% 
  #include checklists with less than 10 observers
  filter(number_observers <= 10) %>%
  #removing repeats
  group_by(sampling_id,scientific_name) %>% 
  slice(1) %>% 
  ungroup %>%
  #add number of species column
  group_by(sampling_id) %>% 
  mutate(no_sp = n_distinct(scientific_name)) %>% 
  #adding year month day column
  mutate(date = as_date(observation_date),
         year = year(observation_date),
         month = month(observation_date),
         day = day(observation_date),
         DoY = yday(observation_date)) %>% 
  #removing observation count X = 1
  mutate(count = ifelse(observation_count == "X", "1",
                        observation_count)) %>% 
  ungroup()

#select useful columns
ebird_data <- ebird_data %>% 
  dplyr::select(scientific_name, common_name, sampling_id, date, 
         year, month, day, DoY, category.x, state, state_code,
         locality_id, latitude, longitude, time_observations_started,
         protocol_type, duration_minutes, effort_distance_km, effort_area_ha,
         no_sp, count)

#cleaning up variables
#changing effort distance to 0 for stationary checklists and converting to meters
ebird_data <- ebird_data %>% 
  mutate(effort_distance = if_else(protocol_type == "Stationary", 
                                    0, effort_distance_km))
ebird_data$effort_distance <- ebird_data$effort_distance*1000

#changing column name category.x to category
names(ebird_data)[names(ebird_data) == "category.x"] <- "category"

#adding elevation column
#read in ASTER GDEM tiff file
aster <- raster("data/aster-gdem-wg.tiff")

#create a spatial data frame from filtered ebird data
ebird_data_sf <- st_as_sf(ebird_data, coords = c("longitude", "latitude"), 
                       crs = 4326, remove = "F")

#final spatial filter to restrict checklists only within the Western Ghats bounds.
#the auk_bbox function creates a bounding box of the limits of the polygon, 
#including points outside the western ghats boundary
ebird_data_filt <- st_filter(ebird_data_sf, wg)

#extract elevation from ASTER GDEM 
elev_data <- raster::extract(aster, ebird_data_filt, buffer = 100, fun = mean)

#bind elevation with ebird data
ebird_elev <- cbind(ebird_data_filt, elev_data)

#rename elevation column
names(ebird_elev)[names(ebird_elev) == "dat_elev"] <- "elev"
ebird_elev$elev <- round(ebird_elev$elev, 2)

#adding a season column
seasons = function(x){
  if(x %in% 3:5) return("summer")
  if(x %in% 6:10) return("monsoon")
  if(x %in% c(11,12,1,2)) return("winter")
  
}

#apply function
ebird_elev$season <- sapply(ebird_elev$month, seasons)

#creating final data frame for analysis
dat <- ebird_elev %>% 
  dplyr::select(date, year, month, day, DoY, scientific_name, common_name, category, count, 
         sampling_id, state, state_code, locality_id, latitude, longitude, no_sp, 
         protocol_type, time_observations_started, duration_minutes, effort_distance, 
         elev, season) %>% 
  st_drop_geometry()

#create an object storing unique locations of all checklists in the filtered ebird data
locations <- dat %>%
  distinct(latitude, longitude, .keep_all = T) %>% 
  dplyr::select(locality_id, latitude, longitude)

#calculate sample size
sample_size_season <- dat %>% 
  group_by(scientific_name, season) %>% 
  summarise(n = n()) %>% 
  ungroup()
sample_size_season <- sample_size_season %>% 
  pivot_wider(names_from = season, values_from = n)

#create new species list based on 10 observations per season
sample_size_season <- sample_size_season %>% 
  filter(monsoon >= 10 & summer >= 10 & winter >=10)
all_species_list <- unique(sample_size_season$scientific_name)

#subset data to only include residents
#read in avonet data set
avonet <- read_csv("data/avonet-ebird.csv")

#the avonet data uses genus Streptopelia instead of Spilopelia
#changing the genus for filtering
avonet$Species2[avonet$Species2 == "Streptopelia chinensis"] <- "Spilopelia chinensis"
avonet$Species2[avonet$Species2 == "Streptopelia senegalensis"] <- "Spilopelia senegalensis"

#filtering avonet data for study species
avonet <- avonet %>% 
  filter(Species2 %in% all_species_list) %>% 
  dplyr::select(Species2, `Hand-Wing.Index`, Mass, Trophic.Niche, Migration)#Malabar Imperial Pigeon doesn't have trait data will be dropped

#using the Migration column to subset the avonet data to include only residents (1)
#and partial migrants (2), Migration = 3 are long-distance migrants
avonet <- subset(avonet, Migration != 3)

#cleaning column names
avonet_col_names <- colnames(avonet) %>% 
  str_replace("Species2", "scientific_name") %>%
  str_replace("Hand-Wing.Index", "HWI") %>% 
  str_replace("Trophic.Niche", "diet") %>%
  str_to_lower()
setnames(avonet, avonet_col_names)

#creating a list of resident species
resident_species <- unique(avonet$scientific_name)

#filter ebird data to include only resident species
dat <- dat %>% 
  filter(scientific_name %in% resident_species)

#write a csv file with final data frame, filtered avonet data and unique locations 
write.csv(dat, file = "data/ebird_elev_residents_WG.csv", row.names = FALSE)
write.csv(locations, file = "data/unique_chklst_locations_WG.csv", row.names = FALSE)
write.csv(avonet, file = "data/avonet_WG.csv", row.names = FALSE)
