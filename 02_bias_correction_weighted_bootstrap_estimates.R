library(data.table)
library(tidyverse)
library(ggpubr)
library(mgcv)
library(boot)

theme_set(theme_classic())
mycolors <- c(summer = "#C7E02F", winter = "#471164FF", monsoon = "#24868EFF")

dat <- read_csv("data/ebird_elev_residents_WG.csv")
min_elev <- min(dat$elev)
max_elev <- max(dat$elev)
season_list <- unique(dat$season)
species_list <- unique(dat$common_name)


#calculate cumulative minutes sampled for each season/elevation combination
sampling_effort <- dat %>%
  dplyr::group_by(season, elev) %>%
  dplyr::summarise(cum_min = sum(duration_minutes, na.rm = F))

#plot.
sampling_effort_raw <- sampling_effort %>%
  ggplot() +
  geom_col(aes(x=elev, y = cum_min), color = "black") +
  facet_wrap(~season) +
  labs(x = "Elevation (m)", y = "Effort (Minutes)", title = "Raw Data")+
  theme(axis.title = element_text(size = 18), title = element_text(size = 18),
        strip.text = (element_text(size = 18)), 
        axis.text = element_text(size = 12))


#fit GAM for each season and predict sampling effort based on elevation
sampling_effort_gam <- lapply(season_list, function(s) {
  #filter data for the current season
  df <- dplyr::filter(sampling_effort, season == s)
  #fit GAM of cumulative minutes sampled ~ elevation
  mod <- mgcv::gam(cum_min~s(elev), data=df, family = "poisson")
  #create data frame with a sequence of elevation values from min to max for prediction
  pred_df <- data.frame(elev = seq(from = min_elev, to = max_elev, by = 0.01))
  #predict sampling effort per elevation using the GAM
  pred <- predict(mod, pred_df)
  #output data frame
  out <- data.frame(pred_df, effort = pred, season = s)
  return(out)
}) %>%
  #bind rows of output to get all seasons in one data frame
  bind_rows()

#plot
sampling_effort_gam <- samplingEffort_gam %>%
  ggplot() +
  geom_col(aes(x=elev, y = effort), color = "black") +
  facet_wrap(~season) +
  labs(x = "Elevation (m)", y = "Effort (Minutes)", title = "GAM")+
  theme(axis.title = element_text(size = 18), title = element_text(size = 18),
        strip.text = (element_text(size = 18)), 
        axis.text = element_text(size = 12))

#round elevation values in GAM data frame
sampling_effort_gam$elev <- round(sampling_effort_gam$elev, 2)

#unite data with filtered ebird data and calculate the inverse of cumulative minutes sampled
dat_effort <- dat %>%
  left_join(sampling_effort_gam, by = c("elev", "season")) %>%
  dplyr::mutate(
    inverse_cum_min = 1/effort)

#calculate sample size per season
sample_size_season <- dat %>% 
  group_by(scientific_name, season) %>% 
  summarise(n = n()) %>% 
  ungroup()
sample_size_season <- sample_size_season %>% 
  pivot_wider(names_from = season, values_from = n)

max(sample_size_season$monsoon)#11726
max(sample_size_season$summer)#11771
max(sample_size_season$winter)#22882

#our data has a max sample size in the winter of 22882, to standardize sampling 
#for each species we will resample each species/season 20000 times

#subsample test species (Indian Blackbird) known to be a elevational migrant (anecdotal)
inbl_raw <- dat_effort %>%
  dplyr::filter(common_name == "Indian Blackbird")

#resample weighted by sampling effort
inbl_boot <- inbl_raw %>%
  group_by(season) %>%
  dplyr::sample_n(size = 20000, replace = T, weight = inverse_cum_min)

#extract median elevation
#raw
inbl_median_raw <- inbl_raw %>% 
  group_by(common_name, season) %>% 
  summarise(med = median(elev))
#resampled
inbl_median_boot <- inbl_boot %>% 
  group_by(common_name, season) %>% 
  summarise(med = median(elev))

#plot - test species raw observation data.
inbl_raw_plot <- inbl_raw %>%
  ggplot() +
  geom_histogram(aes(elev, color = season, fill = season), alpha = 0.7) +
  facet_wrap(~season)+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  labs(title = "Indian Blackbird Raw Data", x = "Elevation (m)", y = "Count")+
  theme_pubclean(base_size = 16)+
  theme(legend.position = "none")

#plot - test species bootstrapped data data.
inbl_boot_plot <- inbl_boot %>%
  ggplot() +
  geom_histogram(aes(elev, color = season, fill = season), alpha = 0.7) +
  facet_wrap(~season)+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors)+
  labs(title = "Indian Blackbird Corrected", x = "Elevation (m)", y = "Count")+
  theme_pubclean(base_size = 16)+
  theme(legend.position = "none")
#method seems to be working well by smoothing out sampling for the test species

#calculate weighted bootstrap median elevation estimates per species/season combination
weighted_boot_fun_50 <- function(df) {
  quantile(sample(df$elev, size = round(nrow(df) * .8), prob = df$inverse_cum_min, replace = T), 0.5)
}


set.seed(42)

out_50 <- lapply(species_list, function(species) {
  out_50 <- lapply(season_list, function(s) {
    
    writeLines(paste(species, s, sep = " - "))
    df <- dplyr::filter(dat_effort, common_name == species, season == s)
    if(nrow(df) < 5) {
      writeLines("skipping");
      return(NULL)
    } else {
      
      nreps <- 20000
      out_rep <- replicate(nreps, weighted_boot_fun_50(df))
      df95 <- quantile(out_rep, c(.025, .5, .975))
      
      data.table(season = s, common_name = species, t(data.frame(df95)))
    }
    
    
  }) %>%
    data.table::rbindlist()
}) %>%
  data.table::rbindlist()

#save data
if(!dir.exists("output")) dir.create("output")
fwrite(out_50, file = file.path("output", "median_quantile_estimates.csv"), row.names = F)

#plot
med_30 <- out_50 %>%
  dplyr::arrange(desc(`50%`)) %>%
  dplyr::filter( common_name %in% c("Indian Blackbird", 
                                    "Nilgiri Pipit", 
                                    "Cinereous Tit", 
                                    "Chestnut-headed Bee-eater", 
                                    species_list[1:6])) %>%
  ggplot() +
  aes(y = reorder(common_name, `50%`), yend = common_name, color = season) +
  geom_segment(aes(x = `2.5%`, xend = `97.5%`), linewidth = 1, alpha = 0.7) +
  geom_point(aes(x = `50%`), size = 3, alpha = 0.7) +
  scale_color_manual(values = mycolors) +
  labs(x = "Elevation (m)", title = "Quantifying Elevational Migration") +
  theme_pubclean(base_size = 16) +
  theme(axis.title.y = element_blank(), 
        axis.text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position =  "right")
med_30