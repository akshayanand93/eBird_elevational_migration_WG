#clear environment and perform garbage collection
rm(list = ls())
gc()

#load required_libraries
library(data.table)
library(raster)
library(ggplot2)
library(tidyverse)
library(stringr)
library(ggpubr)
library(boot)
library(cowplot)

theme_set(theme_pubclean())
mycolors <- c(summer = "#C5283D", winter = "#1C3144", monsoon = "#228B22")

#Process trait and eBrid data

#Process functional trait data
#import trait data
residentlist_ft <- read.csv("Taiwan_Breeding_Bird_Trait.csv")

#temperature range
residentlist_ft$TempRange_Global <- residentlist_ft$TempMax_Global - residentlist_ft$TempMin_Global

#dietary diversity
residentlist_ft$DietDiv <- apply(dplyr::select(residentlist_ft, starts_with("Diet_")), 1, 
                                 function(x) (-1) * sum((x/100) * log(x/100), na.rm = T))

#fruits and nectar in diet
residentlist_ft$Diet_fruit_nect <- residentlist_ft$Diet_fruit + residentlist_ft$Diet_nect

#log body mass
residentlist_ft$BodyMassLog <- log10(residentlist_ft$BodyMassMean)

#nest structure (closed nest as the reference group)
residentlist_ft$NestStr <- 1
residentlist_ft$NestStr[apply(dplyr::select(residentlist_ft, starts_with("NestStr"))[, 4:7], 1, sum) > 0] <- 0
residentlist_ft$NestStr <- factor(residentlist_ft$NestStr)

#standardize trait values
residentlist_ft.scaled <- cbind(dplyr::select(residentlist_ft, starts_with("Scientific"), NestStr), 
                                scale(dplyr::select(residentlist_ft, Diet_inv, TempMin_Global, TempRange_Global:BodyMassLog)))

#save processed trait data
saveRDS(residentlist_ft.scaled, file = "output/residentlist_ft_scaled.rds")

#Process eBird data
#import ebird data
ebird <- fread("ebd_TW_relFeb-2019_ebirdp.txt")

#select complete checklists for resident birds in breeding and non-breeding seasons after 2000
ebirdym <- subset(ebird, year(ebird$'OBSERVATION DATE') >= 2000 &
                    month(ebird$'OBSERVATION DATE') %in% c(1, 4:6, 11:12) &
                    ebird$"ALL SPECIES REPORTED" == 1 &
                    ebird$"SCIENTIFIC NAME" %in% residentlist_ft.scaled$ScientificName_eBird)

#remove duplicate records
dataselect <- ebirdym[!duplicated(ebirdym[, c("SCIENTIFIC NAME", "OBSERVATION DATE", "LATITUDE", "LONGITUDE", "GROUP IDENTIFIER")]), 
                      c("SCIENTIFIC NAME", "OBSERVATION DATE", "OBSERVER ID", "LATITUDE", "LONGITUDE", "SAMPLING EVENT IDENTIFIER", "PROTOCOL TYPE", 
                        "PROTOCOL CODE", "DURATION MINUTES", "EFFORT DISTANCE KM", "EFFORT AREA HA", "NUMBER OBSERVERS", "TIME OBSERVATIONS STARTED")]

#remove records with extreme sampling effort
dataselect.final <- dataselect[!which((dataselect$`EFFORT DISTANCE KM` > 1.00) | (dataselect$`EFFORT AREA HA` > 1.00) | 
                                        (dataselect$`DURATION MINUTES` > 300) | (dataselect$`NUMBER OBSERVERS` > 10)), ]

#obtain the elevation information of observation records 
#import DEM data
DEM <- raster("twdtm_asterV2_30m/tif file/twdtm_asterV2_30m.tif")

#import Taiwan main island boundary
TW <- shapefile("Taiwan_diss_without_island/Taiwan_diss.shp")

#mask out values outside Taiwan
DEM.mask <- mask(DEM, TW)

#calculate the mean elevation within 100 m of each sampling event
sites <- unique(dataselect.final[, c("LONGITUDE","LATITUDE")])
sites$elevation <- raster::extract(DEM.mask, sites, buffer = 100, fun = mean)
dataselect.elev <- merge(dataselect.final, sites, by = c("LONGITUDE", "LATITUDE"), all.x = T)

#remove the records outside Taiwan main island
dataselect.elev.final <- dataselect.elev[complete.cases(dataselect.elev$elevation), ]

#estimate unit-area sampling effort for each elevation band
DEM.values <- getValues(DEM.mask)
DEM.values <- DEM.values[complete.cases(DEM.values)]
Area <- hist(DEM.values, breaks = c(-Inf, 500, 1000, 1500, 2000, 2500, Inf), plot = F)$counts

#save processed data
save(dataselect.elev.final, DEM.values, Area, file = "output/dataselect.elev.final_twmar.RData")


#Resample sampling events under different sampling efforts  

#import data
load(file = "output/dataselect.elev.final_twmar.RData")

#obtain sample event IDs in each season and in each elevation band 
samEvents <- dataselect.elev.final[!duplicated(dataselect.elev.final$"SAMPLING EVENT IDENTIFIER"), ]
samEvents$elev_level <- cut(samEvents$elevation, breaks = c(-Inf, 500, 1000, 1500, 2000, 2500, Inf), labels = 1:6)
samEvents.B <- subset(samEvents, lubridate::month(samEvents$"OBSERVATION DATE") %in% 4:6)
samEvents.W <- subset(samEvents, lubridate::month(samEvents$"OBSERVATION DATE") %in% c(1, 11, 12))

ID.B <- lapply(1:6, function(x) samEvents.B$"SAMPLING EVENT IDENTIFIER"[samEvents.B$elev_level == x])
ID.W <- lapply(1:6, function(x) samEvents.W$"SAMPLING EVENT IDENTIFIER"[samEvents.W$elev_level == x])

#estimate sampling efforts
efforts.B <- summary(samEvents.B$elev_level) / (Area * 0.03 * 0.03)
efforts.W <- summary(samEvents.W$elev_level)  / (Area * 0.03 * 0.03)

uniSpe <- as.data.frame(unique(dataselect.elev.final[, "SCIENTIFIC NAME"]))

#obtain 3 levels of sampling effort (q1, q2, q3 of the number of sampling events for each elevation band)
qEffort <- quantile(c(efforts.B, efforts.W), c(0.25, 0.50, 0.75))

#resample for elevations of lower boundary (5th percentile), center (median), and upper boundary (95th percentile) for each species
system.time(for (qt in c(0.05, 0.50, 0.95)) {
  # resample using q1, q2 and q3 efforts
  for (i in 1:3) {
    samplingN <- round((Area*0.03*0.03) * qEffort[i]) 
    
    # resample sampling events for each elevation band 
    set.seed(56789)
    sampleID.B <- lapply(1:1000, function(y) {unlist(lapply(1:6, function(x) sample(ID.B[[x]], samplingN[x], replace=T)))})
    set.seed(56789)
    sampleID.W <- lapply(1:1000, function(y) {unlist(lapply(1:6, function(x) sample(ID.W[[x]], samplingN[x], replace=T)))})
    
    #calculate the elevation percentiles (5, 50, 95) for each bird species in the two seasons (using parallel processing)
    ptm <- proc.time()
    rs <- mclapply(1:1000, function(y) { 
      m <- t(sapply(uniSpe[, 1], function(x){
        occs.B <- subset(dataselect.elev.final, `SCIENTIFIC NAME`==x & `SAMPLING EVENT IDENTIFIER` %in% sampleID.B[[y]], select="elevation")
        occs.W <- subset(dataselect.elev.final, `SCIENTIFIC NAME`==x & `SAMPLING EVENT IDENTIFIER` %in% sampleID.W[[y]], select="elevation")
        B.n <- nrow(occs.B)
        W.n <- nrow(occs.W)
        B.elev <- if (B.n > 0) quantile(occs.B$elevation, qt) else NA
        W.elev <- if (W.n > 0) quantile(occs.W$elevation, qt) else NA
        return(c(B.elev = B.elev, B.n = B.n, W.elev = W.elev, W.n = W.n))
      }))
    }, mc.cores = 1)
    proc.time() - ptm
    save(rs, file = paste0("output/taiwan_resampled_", sprintf("%02d", qt*100), ".q", i, ".RData"))
  }
}
)

#calculate the difference in the elevation between the two seasons (non-breeding - breeding)
#import processed ebird data
load("output/dataselect.elev.final_twmar.RData")
ebirddata <- dataselect.elev.final
ebirddata$obdate <- as.Date(ebirddata$'OBSERVATION DATE')
ebirdset <- ebirddata[, c("LONGITUDE","LATITUDE","SCIENTIFIC NAME","obdate","elevation")]

#compute median breeding elevation
ebirdset_m <- subset(ebirdset , month(ebirdset$obdate) %in% c(4:6))
ebirdset_breedele <- aggregate(ebirdset_m[,5], list(ebirdset_m$`SCIENTIFIC NAME`), median)
ebirdset_breedele$spcies <- ebirdset_breedele$Group.1
saveRDS(ebirdset_breedele, file = "output/ebirdset_breedele.rds")

#calculate and save migration data
system.time(for (minSample in c(30, 100)) {
  for (ds in c("05.q1", "05.q2", "05.q3", "95.q1", "95.q2", "95.q3", "50.q1", "50.q2", "50.q3")) {
    load(paste("output/taiwan_resampled_", ds, ".RData", sep = ""))
    
    B.elev <- sapply(1:1000, function(x) rs[[x]][, 1])
    B.n <-    sapply(1:1000, function(x) rs[[x]][, 2])
    W.elev <- sapply(1:1000, function(x) rs[[x]][, 3])
    W.n <-    sapply(1:1000, function(x) rs[[x]][, 4])
    
    diff <- W.elev - B.elev
    dimnames(diff) <- list(uniSpe[,1], 1:1000)
    if (minSample == 30) diff[W.n < 30 | B.n < 30] <- NA
    if (minSample == 100) diff[W.n < 100 | B.n < 100] <- NA
    
    diff.sel <- diff[apply(diff, 1, FUN = function(x) sum(is.na(x))) < 1000, ]
    
    med.CI <- apply(diff.sel, 1, FUN = function(x) quantile((x), c(0.025, .5, .975), na.rm = TRUE))
    med.CI <- t(med.CI)
    med.CI  <- as.data.frame(med.CI)
    med.CI$Group <- with(med.CI, ifelse(`2.5%` > 0 & `97.5%` > 0, 2, ifelse(`2.5%` < 0 & `97.5%` < 0, 1, 0)))
    med.CI$LCI_diff <- as.numeric(med.CI$`2.5%`)
    med.CI$Median_diff <- as.numeric(med.CI$`50%`)
    med.CI$UCI_diff <- as.numeric(med.CI$`97.5%`)
    med.CI$Species <- rownames(med.CI)
    med.CI <- med.CI[, -c(1,2,3)]
    
    #adding median breeding elevation
    dimnames(B.elev) <- list(uniSpe[,1], 1:1000)
    if (minSample == 30) B.elev[W.n < 30 | B.n < 30] <- NA
    if (minSample == 100) B.elev[W.n < 100 | B.n < 100] <- NA
    
    B.sel <- B.elev[apply(B.elev, 1, FUN = function(x) sum(is.na(x))) < 1000, ]
    
    med.S.CI <- apply(B.sel, 1, FUN = function(x) quantile((x), c(0.025, .5, .975), na.rm = TRUE))
    med.S.CI <- t(med.S.CI)
    med.S.CI  <- as.data.frame(med.S.CI)
    med.S.CI$LCI_S <- as.numeric(med.S.CI$`2.5%`)
    med.S.CI$Median_S <- as.numeric(med.S.CI$`50%`)
    med.S.CI$UCI_S <- as.numeric(med.S.CI$`97.5%`)
    med.S.CI$Species <- rownames(med.S.CI)
    med.S.CI <- med.S.CI[, -c(1,2,3)]
    
    #adding median breeding elevation
    dimnames(W.elev) <- list(uniSpe[,1], 1:1000)
    if (minSample == 30) W.elev[W.n < 30 | B.n < 30] <- NA
    if (minSample == 100) W.elev[W.n < 100 | B.n < 100] <- NA
    
    W.sel <- W.elev[apply(W.elev, 1, FUN = function(x) sum(is.na(x))) < 1000, ]
    
    med.W.CI <- apply(W.sel, 1, FUN = function(x) quantile((x), c(0.025, .5, .975), na.rm = TRUE))
    med.W.CI <- t(med.W.CI)
    med.W.CI  <- as.data.frame(med.W.CI)
    med.W.CI$LCI_W <- as.numeric(med.W.CI$`2.5%`)
    med.W.CI$Median_W <- as.numeric(med.W.CI$`50%`)
    med.W.CI$UCI_W <- as.numeric(med.W.CI$`97.5%`)
    med.W.CI$Species <- rownames(med.W.CI)
    med.W.CI <- med.W.CI[, -c(1,2,3)]
    
    ts <- left_join(med.CI, med.S.CI, by = "Species") %>% 
      left_join(., med.W.CI, by = "Species", )
    ts <- ts[c("Species", "Group", "Median_S", "LCI_S", "UCI_S", 
               "Median_W", "LCI_W", "UCI_W", "Median_diff", "LCI_diff", "UCI_diff")]
    write.csv(ts, file = paste("output/elev_migration_taiwan_", ds, ".sam", minSample, ".csv", sep = ""), row.names = F)
  }
}
)

#------------------------------------------------------------------------------#

##using weighted bootstrapped quantile estimates to calculate median elevation quantiles
#load in ebird data
dat <- dataselect.elev.final

#tidy up column names
new_col_names <- colnames(dat) %>% 
  str_replace("elevation", "elev") %>%
  str_replace(" ", "_") %>% 
  str_to_lower()
setnames(dat, new_col_names)

#additional filtering
dat_filt <- dat %>% 
  #adding year month day column
  mutate(date = as_date(observation_date),
         year = year(observation_date),
         month = month(observation_date),
         day = day(observation_date),
         DoY = yday(observation_date))

#removing duplicates
dat_filt <- dat_filt[!duplicated(dat_filt$`sampling_event identifier`)]

#adding a season column
seasons = function(x){
  if(x %in% 4:6) return("summer")
  if(x %in% c(1, 11, 12)) return("winter")
  
}
dat_filt$season <- sapply(dat_filt$month, seasons)

#rounding elevation values
dat_filt$elev <- round(dat_filt$elev, digits = 2)

#setting limits
minElev          <- min(dat_filt$elev)
maxElev          <- max(dat_filt$elev)
seasonList       <- unique(dat_filt$season)
speciesList      <- unique(dat_filt$scientific_name)

#calculate cumulative minutes sampled for each season/elevation combination
samplingEffort <- dat_filt %>%
  dplyr::group_by(season, elev) %>%
  dplyr::summarise(cum_min = sum(duration_minutes, na.rm = F))

#plot
samplingEffort %>%
  ggplot() +
  geom_col(aes(x=elev, y = cum_min), color = "black") +
  facet_wrap(~season)

#GAM of sampling effort.
samplingEffort_gam <- lapply(seasonList, function(s) {
  df <- dplyr::filter(samplingEffort, season == s)
  m <- mgcv::gam(cum_min~s(elev), data=df, family = "poisson")
  pred_df <- data.frame(elev = seq(from = minElev, to = maxElev, by = 0.01))
  pred <- predict(m, pred_df)
  out <- data.frame(pred_df, effort = pred, season = s)
  return(out)
}) %>%
  bind_rows()

#plot
samplingEffort_gam %>%
  ggplot() +
  geom_col(aes(x=elev, y = effort), color = "black") +
  facet_wrap(~season)

#round elev column in GAM df
samplingEffort_gam$elev <- round(samplingEffort_gam$elev, 2)

#create data frame
dat_effort <- dat_filt %>%
  left_join(samplingEffort_gam, by = c("elev", "season"))%>%
  dplyr::mutate(inverse_cum_min = 1/effort)

#calculate median elevation with CI
boot_fun_50 <- function(df) {
  quantile(sample(df$elev, size = round(nrow(df) * .8), prob = df$inverse_cum_min, replace = T), 0.5)
}


set.seed(42)

system.time(out <- lapply(speciesList, function(species) {
  out <- lapply(seasonList, function(s) {
    
    writeLines(paste(species, s, sep = " - "))
    df <- dplyr::filter(dat_effort, scientific_name == species, season == s)
    if(nrow(df) < 10) {
      writeLines("skipping");
      return(NULL)
    } else {
      
      nreps <- 20000
      out_rep <- replicate(nreps, boot_fun_50(df))
      df95 <- quantile(out_rep, c(.025, .5, .975))
      
      data.table(season = s, scientific_name = species, t(data.frame(df95)))
    }
    
    
  }) %>%
    data.table::rbindlist()
}) %>%
  data.table::rbindlist()
)

#sample size
sample_size_spp <- dat_filt %>% 
  group_by(scientific_name, season) %>% 
  summarise(sample = n())

#plotting the two methods
#subset and clean weighted bootstrap estimate (wbe) data
wbe_summer <- subset(out, season == "summer")
wbe_summer <- wbe_summer[, -1]
colnames(wbe_summer) <- c("Species", "LCI_S_ECQ", "Median_S_ECQ", "UCI_S_ECQ")

wbe_winter <- subset(out, season == "winter")
wbe_winter <- wbe_winter[, -1]
colnames(wbe_winter) <- c("Species", "LCI_W_ECQ", "Median_W_ECQ", "UCI_W_ECQ")

#read in taiwan data
median_taiwan_q1_100 <- read.csv("output/elev_migration_taiwan_50.q1.sam100.csv")
median_taiwan_q1_30 <- read.csv("output/elev_migration_taiwan_50.q1.sam30.csv")
median_taiwan_q2_100 <- read.csv("output/elev_migration_taiwan_50.q2.sam100.csv")
median_taiwan_q2_30 <- read.csv("output/elev_migration_taiwan_50.q2.sam30.csv")
median_taiwan_q3_100 <- read.csv("output/elev_migration_taiwan_50.q3.sam100.csv")
median_taiwan_q3_30 <- read.csv("output/elev_migration_taiwan_50.q3.sam30.csv")

#select useful columns
median_taiwan_q1_100 <- median_taiwan_q1_100[, c(1, 3, 6)]
colnames(median_taiwan_q1_100) <- c("Species", "Median_S_q1_100", "Median_W_q1_100")
median_taiwan_q1_30 <- median_taiwan_q1_30[, c(1, 3, 6)]
colnames(median_taiwan_q1_30) <- c("Species", "Median_S_q1_30", "Median_W_q1_30")
median_taiwan_q2_100 <- median_taiwan_q2_100[, c(1, 3, 6)]
colnames(median_taiwan_q2_100) <- c("Species", "Median_S_q2_100", "Median_W_q2_100")
median_taiwan_q2_30 <- median_taiwan_q2_30[, c(1, 3, 6)]
colnames(median_taiwan_q2_30) <- c("Species", "Median_S_q2_30", "Median_W_q2_30")
median_taiwan_q3_100 <- median_taiwan_q3_100[, c(1, 3, 6)]
colnames(median_taiwan_q3_100) <- c("Species", "Median_S_q3_100", "Median_W_q3_100")
median_taiwan_q3_30 <- median_taiwan_q3_30[, c(1, 3, 6)]
colnames(median_taiwan_q3_30) <- c("Species", "Median_S_q3_30", "Median_W_q3_30")
median_taiwan_q1_100_summer <- median_taiwan_q1_100[c(1, 2)]
median_taiwan_q1_100_winter <- median_taiwan_q1_100[c(1, 3)]
median_taiwan_q1_30_summer <- median_taiwan_q1_30[c(1, 2)]
median_taiwan_q1_30_winter <- median_taiwan_q1_30[c(1, 3)]
median_taiwan_q2_100_summer <- median_taiwan_q2_100[c(1, 2)]
median_taiwan_q2_100_winter <- median_taiwan_q2_100[c(1, 3)]
median_taiwan_q2_30_summer <- median_taiwan_q2_30[c(1, 2)]
median_taiwan_q2_30_winter <- median_taiwan_q2_30[c(1, 3)]
median_taiwan_q3_100_summer <- median_taiwan_q3_100[c(1, 2)]
median_taiwan_q3_100_winter <- median_taiwan_q3_100[c(1, 3)]
median_taiwan_q3_30_summer <- median_taiwan_q3_30[c(1, 2)]
median_taiwan_q3_30_winter <- median_taiwan_q3_30[c(1, 3)]

#join
df_plot_s <- wbe_summer %>% 
  left_join(median_taiwan_q1_100_summer, by = "Species") %>% 
  left_join(., median_taiwan_q1_30_summer, by = "Species") %>% 
  left_join(., median_taiwan_q2_100_summer, by = "Species") %>% 
  left_join(., median_taiwan_q2_30_summer, by = "Species") %>% 
  left_join(., median_taiwan_q3_100_summer, by = "Species") %>% 
  left_join(., median_taiwan_q3_30_summer, by = "Species")

df_plot_w <- wbe_winter %>% 
  left_join(median_taiwan_q1_100_winter, by = "Species") %>% 
  left_join(., median_taiwan_q1_30_winter, by = "Species") %>% 
  left_join(., median_taiwan_q2_100_winter, by = "Species") %>% 
  left_join(., median_taiwan_q2_30_winter, by = "Species") %>% 
  left_join(., median_taiwan_q3_100_winter, by = "Species") %>% 
  left_join(., median_taiwan_q3_30_winter, by = "Species")

#log transform for comparison
df_plot_s <- log(df_plot_s[, c(2:10)])
df_plot_w <- log(df_plot_w[, c(2:10)])

#scatter plot of Tsai et al. (2021) vs weighted bootstrap estimates
#summer
#first quantile sampling effort with 100 occurrence records threshold vs our method
scatter_s_q1_100 <- ggplot(df_plot_s) +
  aes(x = Median_S_q1_100, y = Median_S_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(title = "Summer", x = "Tsai et al Median Quantile (Q1, 100)") +
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#first quantile sampling effort with 30 occurrence records threshold vs our method
scatter_s_q1_30 <- ggplot(df_plot_s) +
  aes(x = Median_S_q1_30, y = Median_S_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q1, 30)")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#second quantile sampling effort with 100 occurrence records threshold vs our method
scatter_s_q2_100 <- ggplot(df_plot_s) +
  aes(x = Median_S_q2_100, y = Median_S_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q2, 100)", 
       y = "Weighted Bootstrap Median Quantile")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8))

#second quantile sampling effort with 30 occurrence records threshold vs our method
scatter_s_q2_30 <- ggplot(df_plot_s) +
  aes(x = Median_S_q2_30, y = Median_S_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q2, 30)")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#third quantile sampling effort with 100 occurrence records threshold vs our method
scatter_s_q3_100 <- ggplot(df_plot_s) +
  aes(x = Median_S_q3_100, y = Median_S_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q3, 100)") +
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#third quantile sampling effort with 30 occurrence records threshold vs our method
scatter_s_q3_30 <- ggplot(df_plot_s) +
  aes(x = Median_S_q3_30, y = Median_S_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q3, 30)")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())


#winter
#first quantile sampling effort with 100 occurrence records threshold vs our method
scatter_w_q1_100 <- ggplot(df_plot_w) +
  aes(x = Median_W_q1_100, y = Median_W_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(title = "Winter", x = "Tsai et al Median Quantile (Q1, 100)") +
  theme_classic()+
  theme(plot.title = element_text(size = 8, face = "bold"), 
        axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#first quantile sampling effort with 30 occurrence records threshold vs our method
scatter_w_q1_30 <- ggplot(df_plot_w) +
  aes(x = Median_W_q1_30, y = Median_W_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q1, 30)")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#second quantile sampling effort with 100 occurrence records threshold vs our method
scatter_w_q2_100 <- ggplot(df_plot_w) +
  aes(x = Median_W_q2_100, y = Median_W_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q2, 100)", 
       y = "Weighted Bootstrap Median Quantile")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8))

#second quantile sampling effort with 30 occurrence records threshold vs our method
scatter_w_q2_30 <- ggplot(df_plot_w) +
  aes(x = Median_W_q2_30, y = Median_W_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q2, 30)")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#third quantile sampling effort with 100 occurrence records threshold vs our method
scatter_w_q3_100 <- ggplot(df_plot_w) +
  aes(x = Median_W_q3_100, y = Median_W_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q3, 100)") +
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#third quantile sampling effort with 30 occurrence records threshold vs our method
scatter_w_q3_30 <- ggplot(df_plot_w) +
  aes(x = Median_W_q3_30, y = Median_W_ECQ) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  scale_y_continuous(breaks = seq(0, 8, by = 1), limits = c(2, 8)) +
  labs(x = "Tsai et al Median Quantile (Q3, 30)")+
  theme_classic()+
  theme(axis.title = element_text(size = 8), 
        axis.text = element_text(size = 8), 
        axis.title.y = element_blank())

#create panel
winter_panel <- plot_grid(scatter_w_q1_100, scatter_w_q1_30, scatter_w_q2_100, 
                          scatter_w_q2_30, scatter_w_q3_100, scatter_w_q3_30, 
                          ncol = 2, nrow = 3, align = "vh")

summer_panel <- plot_grid(scatter_s_q1_100, scatter_s_q1_30, scatter_s_q2_100, 
                          scatter_s_q2_30, scatter_s_q3_100, scatter_s_q3_30, 
                          ncol = 2, nrow = 3, align = "vh")

#save plots
ggsave("figs/method_compare_s_S1.2.tif", summer_panel, device = "tif", height = 230, width = 168, units = "mm", dpi = 600)
ggsave("figs/method_compare_w_S1.3.tif", winter_panel, device = "tif", height = 230, width = 168, units = "mm", dpi = 600)
