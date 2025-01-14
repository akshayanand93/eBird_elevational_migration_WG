#clear environment and perform garbage collection
rm(list = ls())
gc()

#load required libraries
library(tidyverse)
library(cowplot)
library(ggridges)
library(ggpubr)
library(stringr)
library(grid)

#read in data
dat <- read_csv("data/ebird_elev_residents_WG.csv")
mig_df <- read_csv("output/elev_migration_values.csv")
elev_mig_wg <- read_csv("output/elev_mig_climate_wg.csv")

##plotting elevational migration for figure 2
#create column to denote significant shifts
elev_mig_wg <- elev_mig_wg %>% 
  mutate(m_mig = ifelse(s_uci < m_lci | m_uci < s_lci, s_m_mig, 0), 
         w_mig = ifelse(m_uci < w_lci | w_uci < m_lci, m_w_mig, 0), 
         s_mig = ifelse(w_uci < s_lci | s_uci < w_lci, w_s_mig, 0))

#select useful columns
dist_df <- elev_mig_wg %>% 
  dplyr::select(common_name, scientific_name, s_mig, m_mig, w_mig)

#filter for significant migrants
#summer
s_dist <- dist_df %>% 
  filter(s_mig != 0) %>% 
  mutate(slope = ifelse(s_mig > 0, "upslope", "downslope"))
s_dist$scientific_name <- str_replace(s_dist$scientific_name, "_", " ")

#monsoon
m_dist <- dist_df %>% 
  filter(m_mig != 0) %>% 
  mutate(slope = ifelse(m_mig > 0, "upslope", "downslope"))
m_dist$scientific_name <- str_replace(m_dist$scientific_name, "_", " ")

#winter
w_dist <- dist_df %>% 
  filter(w_mig != 0) %>% 
  mutate(slope = ifelse(w_mig > 0, "upslope", "downslope"))
w_dist$scientific_name <- str_replace(w_dist$scientific_name, "_", " ")

#plotting for each season
#summer
s_dist_plot <- ggplot(s_dist) +
  aes(x = reorder(scientific_name, -s_mig), xend = scientific_name) +
  geom_point(aes(y = s_mig), shape = ifelse(s_dist$s_mig > 0, 24, 25), 
             color = "black", fill = "black") +
  geom_segment(aes(y = 0, yend = s_mig), color = "black") +
  scale_y_continuous(breaks = seq(-100, 600, 100)) +
  theme_pubclean(base_size = 8) +
  labs(y = "Distance (m)") +
  theme(text = element_text(family = "sans"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, face = "italic", hjust = 1),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

#monsoon
m_dist_plot <- ggplot(m_dist) +
  aes(x = reorder(scientific_name, -m_mig), xend = scientific_name) +
  geom_point(aes(y = m_mig), shape = ifelse(m_dist$m_mig > 0, 24, 25), 
             color = "black", fill = "black") +
  geom_segment(aes(y = 0, yend = m_mig), color = "black") +
  scale_y_continuous(breaks = seq(-400, 400, 100)) +
  theme_pubclean(base_size = 8) +
  labs(y = "Distance (m)") +
  theme(text = element_text(family = "sans"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, face = "italic", hjust = 1), 
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

m_dist_plot


#winter
w_dist_plot <- ggplot(w_dist) +
  aes(x = reorder(scientific_name, -w_mig), xend = scientific_name) +
  geom_point(aes(y = w_mig), shape = ifelse(w_dist$w_mig > 0, 24, 25), 
             color = "black", fill = "black") +
  geom_segment(aes(y = 0, yend = w_mig), color = "black") +
  scale_y_continuous(breaks = seq(-500, 400, 100)) +
  theme_pubclean(base_size = 8) +
  labs(y = "Distance (m)") +
  theme(text = element_text(family = "sans"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, face = "italic", hjust = 1), 
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank())

w_dist_plot

#create the panel
dist_plot <- plot_grid(s_dist_plot, m_dist_plot, w_dist_plot, 
                       labels = c("Summer", "Monsoon", "Winter"), 
                       label_size = 8,
                       label_fontface = "plain",
                       label_fontfamily = "sans",
                       nrow = 3,
                       align = "vh", 
                       rel_heights = c(1.2, 0.9, 0.9),
                       hjust = -0.1, 
                       vjust = c(1.3, -0.1, -0.1))

dist_plot

ggsave("figs/figure_1_dist.png", dist_plot, height = 230, width = 168, dpi = 800, 
       device = "png", units = "mm")

##Supplementary figs
#cleaning season col in dat
dat$season <- str_to_sentence(dat$season)
dat$season <- factor(dat$season, levels = c("Summer", "Monsoon", "Winter"))

#plotting elevational distribution of occurence records
occ_elev_distr <- ggplot(dat, aes(x = elev, y = season)) +
  geom_density_ridges(scale = 1, quantile_lines = TRUE, quantiles = 0.50, alpha = 0.7) +
  labs(x = "Elevation (m)", y = "Density of Occurence Records") +
  scale_y_discrete(limits = c("Winter", "Monsoon", "Summer")) +
  theme_pubclean(base_size = 8) +
  theme(legend.position = "none")

#plotting number of checklists per season
#summarize data
chklist_season <- dat %>% 
  group_by(season) %>% 
  summarise(n_chklst = n_distinct(sampling_id))

#plot
chklst_season_plot <- chklist_season %>% 
  ggplot(aes(x = reorder(season, -n_chklst), y = n_chklst)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  labs(x = "Season",
       y = "No. of Checklists") +
  theme_pubclean(base_size = 8) +
  theme(legend.position = "none")

#create the panel
sam_effort_plot <- plot_grid(occ_elev_distr, chklst_season_plot, 
                             labels = c("(a)", "(b)"), 
                             label_fontface = "plain",
                             label_size = 8,
                             ncol = 2, align = "hv")
#save figs
ggsave("figs/sam-effort-S1.1-a-b.tif", sam_effort_plot, device = "tif", 
       height = 84, width = 168, units = "mm", dpi = 800)

#plotting execution time for tsai vs our method
#create data frame to store execution time values
exec_time <- data.frame("method" = c("Tsai et al.", "Weighted Bootstrap"), 
                        "ex_time" = c(30546.28, 1124.62))

#convert to minutes
exec_time$ex_time <- exec_time$ex_time/60

#plot
exec_time_plot <- ggplot(exec_time) +
  geom_bar(aes(x = method, y = ex_time), stat = "identity", alpha = 0.7) +
  labs(x = "Method", y = "Execution Time (Minutes)") +
  theme_pubclean(base_size = 8)

#save
ggsave("figs/method-exec-time-S1.4.tif", exec_time_plot, device = "tif", 
       height = 79, width = 79, units = "mm", dpi = 800)


#creating figure S1.5
#create a column to denote significant shifts
elev_mig_wg <- elev_mig_wg %>% 
  mutate(m_sig = ifelse(m_mig == 0, "non_sig", "sig"), 
         w_sig = ifelse(w_mig == 0, "non_sig", "sig"), 
         s_sig = ifelse(s_mig == 0, "non_sig", "sig"))

#monsoon plot
m_dir_plot <- ggplot(elev_mig_wg) +
  aes(x = reorder(scientific_name, -s_med), xend = scientific_name, color = m_sig, fill = m_sig) +
  geom_point(aes(y = m_med), shape = 
               ifelse(elev_mig_wg$s_m_mig < 0, 25, 
                      ifelse(elev_mig_wg$s_m_mig > 0, 24, 1)), size = 0.5) +
  geom_linerange(aes(ymin = m_med, ymax = s_med), linewidth = 0.2) +
  scale_color_manual(values = c("non_sig" = "black", "sig" = "#E69F00")) +
  scale_fill_manual(values = c("non_sig" = "black", "sig" = "#E69F00")) +
  theme_pubclean(base_size = 8) +
  labs(y = "Elevation (m)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

m_dir_plot

#winter plot
w_dir_plot <- ggplot(elev_mig_wg) +
  aes(x = reorder(scientific_name, -s_med), xend = scientific_name, color = w_sig, fill = w_sig) +
  geom_point(aes(y = w_med), shape = 
               ifelse(elev_mig_wg$m_w_mig < 0, 25, 
                      ifelse(elev_mig_wg$m_w_mig > 0, 24, 1)), size = 0.5) +
  geom_linerange(aes(ymin = w_med, ymax = m_med), linewidth = 0.2) +
  scale_color_manual(values = c("non_sig" = "black", "sig" = "#E69F00")) +
  scale_fill_manual(values = c("non_sig" = "black", "sig" = "#E69F00")) +
  theme_pubclean(base_size = 8) +
  labs(y = "Elevation (m)", x = "Species") +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

w_dir_plot

#summer plot
s_dir_plot <- ggplot(elev_mig_wg) +
  aes(x = reorder(scientific_name, -s_med), xend = scientific_name, color = s_sig, fill = s_sig) +
  geom_point(aes(y = s_med), shape = 
               ifelse(elev_mig_wg$w_s_mig < 0, 25, 
                      ifelse(elev_mig_wg$w_s_mig > 0, 24, 1)), size = 0.5) +
  geom_linerange(aes(ymin = s_med, ymax = w_med), linewidth = 0.2) +
  scale_color_manual(values = c("non_sig" = "black", "sig" = "#E69F00")) +
  scale_fill_manual(values = c("non_sig" = "black", "sig" = "#E69F00")) +
  theme_pubclean(base_size = 8) +
  labs(y = "Elevation (m)") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

s_dir_plot

#arranging the panels for figure
fig_s1_5_plot <- plot_grid(s_dir_plot, m_dir_plot, w_dir_plot, 
                           labels = c("Summer", "Monsoon", "Winter"), 
                           label_size = 8, 
                           nrow = 3, 
                           align = "vh",
                           vjust = 1.5,
                           hjust = -0.3,
                           rel_heights = c(1, 1, 1))

fig_s1_5_plot

ggsave("figs/figure_s1.5.tif", fig_s1_5_plot, height = 200, width = 168, dpi = 800, 
       device = "tif", units = "mm")

#create figure s1.6
#summarize data
diet_mass <- elev_mig_wg %>% 
  group_by(diet) %>% 
  summarise(mean_mass = mean(mass), 
            n = n(), 
            sd = sd(mass), 
            se = sd/sqrt(n))

#plot
diet_mass_plot <- ggplot(diet_mass) +
  geom_bar(aes(x = reorder(diet, -mean_mass), y = mean_mass), 
           stat = "identity", alpha = 0.7) +
  geom_errorbar(aes(x = reorder(diet, -mean_mass), ymin = mean_mass - se, 
                    ymax = mean_mass + se), width = 0.1) +
  scale_y_continuous(breaks = seq(0, 300, 50)) +
  theme_pubclean(base_size = 8) +
  labs(y = "Mean Body Mass (g)") +
  theme(axis.title.x = element_blank())

diet_mass_plot

ggsave("figs/figure_s1.6_diet_mass.tif", diet_mass_plot, device = "tif", 
       dpi = 800, width = 168, height = 84, units = "mm")

#creating table S1.1
table_s1.1 <- mig_df[, c(2, 1, 3:14)]

table_s1.1 <- table_s1.1 %>% 
  mutate(m_sig = ifelse(s_uci < m_lci | m_uci < s_lci, 1, 0), 
         w_sig = ifelse(m_uci < w_lci | w_uci < m_lci, 1, 0), 
         s_sig = ifelse(w_uci < s_lci | s_uci < w_lci, 1, 0))

table_s1.1[, 3:14] <- lapply(table_s1.1[, 3:14], function(x) {
  if (is.numeric(x)) {
    round(x, 3)
  } else {
    x
  }
})

write.csv(table_s1.1, "output/table_s1.1_migration_values.csv", row.names = FALSE)

#summarizing endemics in the dataset
endemics <- read.csv("data/India-Checklist_v8_2.csv") %>% 
  filter(Endemic..Western.Ghats. == "X") %>% 
  pull(Scientific.Name)

endemics <- str_replace(endemics, " ", "_")

elev_mig_wg_endemics <- elev_mig_wg %>% 
  filter(scientific_name %in% endemics) %>% 
  filter(m_mig != 0 | s_mig != 0 | w_mig != 0)

#save workspace image
save.image("07_figs.RData")
