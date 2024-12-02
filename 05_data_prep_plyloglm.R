#clear environment and perform garbage collection
rm(list = ls())
gc()

#load required libraries
library(tidyverse)
library(data.table)
library(phytools)
library(ape)
library(rtrees)
library(car)
library(reshape2)
library(ggpubr)
library(cowplot)
library(phylolm)

#read in data
out_50 <- read_csv("output/median_quantile_estimates.csv")
dat <- read_csv("data/ebird_elev_residents_WG.csv")

#creating seasonal (s, m, w) lower CI (lci), median (med), and upper CI (uci) columns
mig_df <- out_50 %>% 
  mutate(s_lci = ifelse(season == "summer" & `2.5%` > 0, `2.5%`, NA), 
         s_med = ifelse(season == "summer" & `50%` > 0, `50%`, NA), 
         s_uci = ifelse(season == "summer" & `97.5%` > 0, `97.5%`, NA),
         m_lci = ifelse(season == "monsoon" & `2.5%` > 0, `2.5%`, NA), 
         m_med = ifelse(season == "monsoon" & `50%` > 0, `50%`, NA), 
         m_uci = ifelse(season == "monsoon" & `97.5%` > 0, `97.5%`, NA), 
         w_lci = ifelse(season == "winter" & `2.5%` > 0, `2.5%`, NA), 
         w_med = ifelse(season == "winter" & `50%` > 0, `50%`, NA), 
         w_uci = ifelse(season == "winter" & `97.5%` > 0, `97.5%`, NA))

#removing NA's and summarizing df
mig_df <- mig_df %>%
  group_by(common_name) %>%
  summarise_all(~ na.omit(.)[1]) %>% 
  dplyr::select(common_name, s_lci, s_med, s_uci, m_lci, m_med, m_uci, w_lci, 
         w_med, w_uci)

#calculating migration values and classifying migrants if shifts are significant in one or more seasons
mig_df <- mig_df %>% 
  mutate(s_m_mig = m_med - s_med, #monsoon distance 
         m_w_mig = w_med - m_med, #winter distance
         w_s_mig = s_med - w_med, #summer distance
         migrant = ifelse((s_uci < m_lci | m_uci < s_lci) | 
                            (m_uci < w_lci | w_uci < m_lci) | 
                            (w_uci < s_lci | s_uci < w_lci), 
                          "migrant", "non-migrant"))

#adding scientific names
taxonomy <- read_csv("data/ebird_taxonomy_v2022.csv")
resident_species <- unique(dat$scientific_name)
#filter
taxonomy <- taxonomy %>% 
  dplyr::filter(SCI_NAME %in% resident_species) %>% 
  dplyr::select(PRIMARY_COM_NAME, SCI_NAME)
#cleaning column names
new_col_names <- colnames(taxonomy) %>%
  str_replace("SCI_NAME", "scientific_name") %>% 
  str_replace("PRIMARY_COM_NAME", "common_name") %>% 
  str_to_lower()
setnames(taxonomy, new_col_names)
#join
mig_df <- left_join(taxonomy, mig_df, by = "common_name")

#save data
if(!dir.exists("output")) dir.create("output")
fwrite(mig_df, file = file.path("output", "elev_migration_values.csv"), row.names = F)

#adding trait data
avonet <- read_csv("data/avonet_WG.csv")

#join
mig_trait <- left_join(mig_df, avonet, by = "scientific_name")

#select useful columns
mig_trait <- mig_trait %>% 
  dplyr::select(-"migration")

#reorder columns
col_order <- c("common_name", "scientific_name", "migrant", "s_lci", "s_med", "s_uci", 
               "m_lci", "m_med", "m_uci", "w_lci", "w_med", "w_uci", 
               "s_m_mig", "m_w_mig", "w_s_mig", "hwi", "mass", "diet")
mig_trait <- mig_trait[, col_order]

#add "_" to scientific names
mig_trait$scientific_name <- str_replace(mig_trait$scientific_name, " ", "_")

#save data
fwrite(mig_trait, file = file.path("output", "elev_mig_trait.csv"), row.names = F)

#adding environmental tolerances
climatic_vars_wg <- read.csv("output/eBird_climatic_vars_wg.csv")

#join with data
elev_mig_wg <- left_join(mig_trait, climatic_vars_wg, by = "scientific_name")

#check NA's
which(is.na(elev_mig_wg), arr.ind = TRUE)#no NA's

#create column to denote significant shifts
elev_mig_wg <- elev_mig_wg %>% 
  mutate(s_mig = ifelse(w_uci < s_lci | s_uci < w_lci, w_s_mig, 0), #summer
         m_mig = ifelse(s_uci < m_lci | m_uci < s_lci, s_m_mig, 0), #monsoon
         w_mig = ifelse(m_uci < w_lci | w_uci < m_lci, m_w_mig, 0)) #winter

#classify seasonal migrants
elev_mig_wg <- elev_mig_wg %>% 
  mutate(s_migrant = ifelse(s_mig == 0, 0, 1), #summer
         m_migrant = ifelse(m_mig == 0, 0, 1), #monsoon
         w_migrant = ifelse(w_mig == 0, 0, 1)) #winter

#calculating absolute migration distance
elev_mig_wg <- elev_mig_wg %>% 
  mutate(m_dist_abs = abs(s_m_mig), #monsoon
         w_dist_abs = abs(m_w_mig), #winter
         s_dist_abs = abs(w_s_mig)) #summer

#save data
fwrite(elev_mig_wg, file.path("output", "elev_mig_climate_wg.csv"), row.names = FALSE)

#sub-setting for model
elev_mig_mod_wg <- elev_mig_wg %>% 
  dplyr::select(common_name, scientific_name, migrant, s_m_mig, m_w_mig, 
                w_s_mig, hwi, mass, diet, precip_breadth, temp_breadth, wind_breadth, 
                precip_mean, temp_mean, wind_mean, precip_med, temp_med, wind_med, 
                s_dist_abs, m_dist_abs, w_dist_abs, s_mig, m_mig, w_mig, s_migrant, 
                m_migrant, w_migrant)


#scaling climatic and trait variables
elev_mig_mod_wg[c(7, 8, 10:18)] <- lapply(elev_mig_mod_wg[c(7, 8, 10:18)], 
                                       function(x) c(scale(x)))

#converting diet to factor and setting reference level
elev_mig_mod_wg$diet <- as.factor(elev_mig_mod_wg$diet)
elev_mig_mod_wg$diet <- factor(relevel(elev_mig_mod_wg$diet, ref = "Omnivore"))

#histograms of predictors
hist(elev_mig_wg$hwi, breaks = 50)
hist(elev_mig_wg$mass, breaks = 50)
hist(elev_mig_wg$precip_mean, breaks = 50)
hist(elev_mig_wg$precip_med, breaks = 50)
hist(elev_mig_wg$precip_breadth, breaks = 50)
hist(elev_mig_wg$temp_mean, breaks = 50)
hist(elev_mig_wg$temp_med, breaks = 50)
hist(elev_mig_wg$temp_breadth, breaks = 50)
hist(elev_mig_wg$wind_mean, breaks = 50)
hist(elev_mig_wg$wind_med, breaks = 50)
hist(elev_mig_wg$wind_breadth, breaks = 50)

#------------------------------------------------------------------------------#

#testing within season variation in direction of movement
#filter significant migrants per season
#monsoon
m_slope_df <- elev_mig_wg %>% 
  filter(s_uci < m_lci | m_uci < s_lci) %>% 
  dplyr::select(scientific_name, s_m_mig)

#winter
w_slope_df <- elev_mig_wg %>% 
  filter(m_uci < w_lci | w_uci < m_lci) %>% 
  dplyr::select(scientific_name, m_w_mig)

#summer
s_slope_df <- elev_mig_wg %>% 
  filter(w_uci < s_lci | s_uci < w_lci) %>% 
  dplyr::select(scientific_name, w_s_mig)

#categorize slope
m_slope_df$slope <- factor(ifelse(m_slope_df$s_m_mig > 0, "upslope", "downslope"))
w_slope_df$slope <- factor(ifelse(w_slope_df$m_w_mig > 0, "upslope", "downslope"))
s_slope_df$slope <- factor(ifelse(s_slope_df$w_s_mig > 0, "upslope", "downslope"))

#running a chi-squared test
#monsoon
m_slope_table <- table(m_slope_df$slope)
#test
m_chi_test <- chisq.test(m_slope_table)#highly significant
#check assumptions
m_chi_test$expected#all good

#winter
w_slope_table <- table(w_slope_df$slope)
#test
w_chi_test <- chisq.test(w_slope_table)#non significant
#check assumptions
w_chi_test$expected#all good

#summer
s_slope_table <- table(s_slope_df$slope)
#test
s_chi_test <- chisq.test(s_slope_table)#highly significant
#check assumptions
s_chi_test$expected#all good


#------------------------------------------------------------------------------#

##binomial model to test if traits explain migration status (migrant vs non-migrant)
#getting phylo tree
resident_species <- unique(elev_mig_mod_wg$scientific_name)

#create data frame to extract tree
resident_list <- sp_list_df(sp_list = resident_species, taxon = "bird")
#The following genus are not in birdtree classification database: 
#Yungipicus, Rubigula, Montecincla, Pterorhinus, Sholicola

#five genus missing from bird tree, write file to manually input genus and family
fwrite(resident_list, file.path("data", "resident_list.csv"), row.names = FALSE)
#family of the following species was changed - Yungipicus_nanus -> Picidae, Rubigula_gularis -> Pycnonotidae,
#Montecincla_cachinnans - Montecincla_fairbanki - Montecincla_meridionalis - Pterorhinus_delesserti -> Timaliidae
#Sholicola_major - Sholicola_albiventris -> Muscicapidae

#read in complete species list
resident_list <- read_csv("data/resident_list.csv")

#get tree
resident_tree <- get_tree(sp_list = resident_list, taxon = "bird", scenario = "at_basal_node")

#consensus tree
set.seed(42)
cons_resident_tree <- ls.consensus(resident_tree)

#save the consensus tree
write.tree(cons_resident_tree, "output/cons_resident_tree_wg.tree")

#convert response (migrant) variable to binary
elev_mig_mod_wg$migrant <- ifelse(elev_mig_mod_wg$migrant == "non-migrant", 0, 1)

#creating the model data frame
elev_mig_binom_df_wg <- as.data.frame(elev_mig_mod_wg[, c(2, 3, 7:18, 25:27)])
rownames(elev_mig_binom_df_wg) <- elev_mig_binom_df_wg$scientific_name
elev_mig_binom_df_wg <- elev_mig_binom_df_wg[, c(2, 15, 16, 17, 3:14)]
elev_mig_binom_df_wg$migrant <- factor(elev_mig_binom_df_wg$migrant, levels = c(0, 1))
elev_mig_binom_df_wg$s_migrant <- factor(elev_mig_binom_df_wg$s_migrant)
elev_mig_binom_df_wg$m_migrant <- factor(elev_mig_binom_df_wg$m_migrant)
elev_mig_binom_df_wg$w_migrant <- factor(elev_mig_binom_df_wg$w_migrant)

#running seasonal binomial models to understand seasonal migration status vs traits
#summer
#fitting models with mean and median to assess overall fit
s_binom_mod_mean <- phyloglm(s_migrant ~ hwi + mass + diet + temp_breadth + precip_breadth +
                              wind_breadth + temp_mean + precip_mean + wind_mean, 
                            data = elev_mig_binom_df_wg, phy = cons_resident_tree, 
                            method = "logistic_MPLE")
summary(s_binom_mod_mean) 

s_binom_mod_med <- phyloglm(s_migrant ~ hwi + mass + diet + temp_breadth + precip_breadth +
                              wind_breadth + temp_med + precip_med + wind_med, 
                            data = elev_mig_binom_df_wg, phy = cons_resident_tree, 
                            method = "logistic_MPLE")

summary(s_binom_mod_med) 
#the model with the mean environmental tolerances performs marginally better, however
#we will move forward with the median envt. tolerances as that gives us a clearer estimate 
#of typical environmental values a species experiences

#model diagnostics
vif(glm(s_migrant ~ hwi + mass + diet + temp_breadth + precip_breadth +
          wind_breadth + temp_med + precip_med + wind_med, 
        data = elev_mig_binom_df_wg, family = binomial(link = "logit")))#nothing above 5

#backward model selection
s_binom_mod <- phyloglmstep(s_migrant ~ hwi + mass + diet + temp_breadth + precip_breadth +
                              wind_breadth + temp_med + precip_med + wind_med, 
                            data = elev_mig_binom_df_wg, phy = cons_resident_tree, 
                            method = "logistic_MPLE", direction = "backward")
summary(s_binom_mod)#temp_breadth & precip_med +ve significant, precip_breadth -ve near-significant

#save results
summary_s_binom_mod <- summary(s_binom_mod)
s_binom_mod_table <- as.data.frame(cbind(summary_s_binom_mod$coefficients, 
                                          confint(s_binom_mod)))
s_binom_mod_table$Predictor <- row.names(s_binom_mod_table)
row.names(s_binom_mod_table) <- NULL
s_binom_mod_table <- s_binom_mod_table[, c(7, 1, 5, 6, 2, 3, 4)]
colnames(s_binom_mod_table) <- c("Predictor", "Estimate", "LCI", "UCI","SE", 
                                  "Z value", "P value")

#changing predictor names for easier interpretation
#create function
modify_predictor_names <- function(name) {
  name <- gsub("precip", "Precipitation", name)
  name <- gsub("temp", "Temperature", name)
  name <- gsub("wind", "Wind Speed", name)
  name <- gsub("med", "Median", name)
  name <- gsub("mean", "Mean", name)  
  name <- gsub("breadth", "Tolerance Range", name)  
  name <- gsub("hwi", "HWI", name, ignore.case = TRUE)  
  name <- gsub("mass", "Body Mass", name, ignore.case = TRUE)  
  name <- gsub("diet", "", name, ignore.case = TRUE)
  name <- gsub("_", " ", name)
  name <- gsub("\\b([a-z])", "\\U\\1", name, perl = TRUE)
  
  return(name)
  
}

#apply function
s_binom_mod_table$Predictor <- sapply(s_binom_mod_table$Predictor, modify_predictor_names)

#monsoon
m_binom_mod <- phyloglmstep(m_migrant ~ hwi + mass + diet + temp_breadth + precip_breadth +
                              wind_breadth + temp_med + precip_med + wind_med, 
                            data = elev_mig_binom_df_wg, phy = cons_resident_tree, 
                            method = "logistic_MPLE", direction = "backward")
summary(m_binom_mod)#frugivore & precip_breadth -ve significant, temp_breadth +ve significant

#save results
summary_m_binom_mod <- summary(m_binom_mod)
m_binom_mod_table <- as.data.frame(cbind(summary_m_binom_mod$coefficients, 
                                         confint(m_binom_mod)))
m_binom_mod_table$Predictor <- row.names(m_binom_mod_table)
row.names(m_binom_mod_table) <- NULL
m_binom_mod_table <- m_binom_mod_table[, c(7, 1, 5, 6, 2, 3, 4)]
colnames(m_binom_mod_table) <- c("Predictor", "Estimate", "LCI", "UCI","SE", 
                                 "Z value", "P value")

#changing predictor names for easier interpretation
m_binom_mod_table$Predictor <- sapply(m_binom_mod_table$Predictor, modify_predictor_names)

#winter
w_binom_mod <- phyloglmstep(w_migrant ~ hwi + mass + diet + temp_breadth + precip_breadth +
                              wind_breadth + temp_med + precip_med + wind_med, 
                            data = elev_mig_binom_df_wg, phy = cons_resident_tree, 
                            method = "logistic_MPLE", direction = "backward")
summary(w_binom_mod)#temp_breadth & precip_med +ve significant

#save results
summary_w_binom_mod <- summary(w_binom_mod)
w_binom_mod_table <- as.data.frame(cbind(summary_w_binom_mod$coefficients, 
                                         confint(w_binom_mod)))
w_binom_mod_table$Predictor <- row.names(w_binom_mod_table)
row.names(w_binom_mod_table) <- NULL
w_binom_mod_table <- w_binom_mod_table[, c(7, 1, 5, 6, 2, 3, 4)]
colnames(w_binom_mod_table) <- c("Predictor", "Estimate", "LCI", "UCI","SE", 
                                 "Z value", "P value")

#changing predictor names for easier interpretation
w_binom_mod_table$Predictor <- sapply(w_binom_mod_table$Predictor, modify_predictor_names)

#save all seasons in a data frame
binom_mod_results <- rbind(s_binom_mod_table, m_binom_mod_table, w_binom_mod_table)

#round numeric values to 3 decimal places
binom_mod_results[, 2:7] <- lapply(binom_mod_results[, 2:7], function(x) {
  if (is.numeric(x)) {
    round(x, 3)
  } else {
    x
  }
})

#write binomial model results
write.csv(binom_mod_results, "output/table_2_binom_mod_results.csv", row.names = FALSE)

#plotting results
#summer
#separate predictor
s_binom_mod_table$Predictor <- c("(Intercept)","Temperature\nTolerance\nRange", 
                                 "Precipitation\nTolerance\nRange", 
                                 "Temperature\nMedian", 
                                 "Precipitation\nMedian")

#plot
s_binom_mod_plot <- s_binom_mod_table %>% 
  ggplot(aes(y = reorder(Predictor, -Estimate), x = Estimate)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, height = 0.4), linewidth = 0.3) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey20") +
  geom_text(aes(x = UCI, label = ifelse(`P value` < 0.05, "*", "")), 
            vjust = 0.5, hjust = -0.5, color = "black") +
  theme_pubclean(base_size = 9) +
  theme(text = element_text(family = "sans"), 
        axis.title = element_blank(),
        plot.background = element_rect(color = "grey20"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

s_binom_mod_plot

#monsoon
#separate predictor names
m_binom_mod_table$Predictor <- c(
  "(Intercept)",
  "HWI",
  "Aquatic\nPredator",
  "Frugivore",
  "Granivore",
  "Invertivore",
  "Nectarivore",
  "Vertivore",
  "Temperature\nTolerance\nRange",
  "Precipitation\nTolerance\nRange",
  "Wind Speed\nTolerance\nRange",
  "Precipitation\nMedian",
  "Wind Speed\nMedian"
)

#plot
m_binom_mod_plot <- m_binom_mod_table %>% 
  ggplot(aes(y = reorder(Predictor, -Estimate), x = Estimate)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, height = 0.4), linewidth = 0.3) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey20") +
  geom_text(aes(x = UCI, label = ifelse(`P value` < 0.05, "*", "")), 
            vjust = 0.5, hjust = -0.5, color = "black") +
  scale_y_discrete(limits = c(
    "Temperature\nTolerance\nRange",
    "Precipitation\nMedian",
    "Wind Speed\nTolerance\nRange",
    "HWI",
    "Invertivore",
    "Granivore",
    "Nectarivore",
    "Frugivore",
    "Aquatic\nPredator",
    "Vertivore",
    "Precipitation\nTolerance\nRange",
    "Wind Speed\nMedian", 
    "(Intercept)"
  )) +
  theme_pubclean(base_size = 9) +
  theme(text = element_text(family = "sans"), 
        axis.title = element_blank(),
        plot.background = element_rect(color = "grey20"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

m_binom_mod_plot

#winter
#separate predictor
w_binom_mod_table$Predictor <- c("(Intercept)", "Temperature\nTolerance\nRange", 
                                 "Precipitation\nMedian")

#plot
w_binom_mod_plot <- w_binom_mod_table %>% 
  ggplot(aes(y = reorder(Predictor, -Estimate), x = Estimate)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, height = 0.4), linewidth = 0.3) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey20") +
  geom_text(aes(x = UCI, label = ifelse(`P value` < 0.05, "*", "")), 
            vjust = 0.5, hjust = -0.5, color = "black") +
  theme_pubclean(base_size = 9) +
  theme(text = element_text(family = "sans"), 
        axis.title = element_blank(),
        plot.background = element_rect(color = "grey20"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

w_binom_mod_plot

wg_binom_mod_plot <- plot_grid(
  s_binom_mod_plot, NULL, m_binom_mod_plot, NULL, w_binom_mod_plot, 
  labels = c("Summer", "", "Monsoon", "", "Winter"), 
  label_fontface = "plain",
  label_fontfamily = "sans",
  label_size = 8,
  nrow = 5, 
  align = "v",
  rel_heights = c(1.05, 0.025, 2.3, 0.025, 0.85), 
  hjust = -0.5, 
  vjust = 1.5)

wg_binom_mod_plot <- ggdraw() +
  draw_plot(wg_binom_mod_plot, 0, 0, 1, 1) + 
  draw_label(
    "Model Estimates ± CI", 
    x = 0.5, 
    y = 0.02, 
    vjust = 0.5, 
    hjust = 0.5, 
    size = 8, 
    fontfamily = "sans")

wg_binom_mod_plot

#save plot
if(!dir.exists("figs")) dir.create("figs")
ggsave("figs/binom_model_estimates_wg.tif", 
       wg_binom_mod_plot, 
       device = "tif", 
       height = 200, 
       width = 110, 
       units = "mm", 
       dpi = 800)

#save workspace image
save.image("05_data_prep_phyloglm.RData")