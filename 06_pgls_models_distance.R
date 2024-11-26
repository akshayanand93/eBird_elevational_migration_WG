library(tidyverse)
library(MuMIn)
library(caper)
library(nlme)

#read in data
elev_mig_mod_wg <- read.csv("output/elev_mig_climate_wg.csv")
cons_resident_tree <- read.tree("output/cons_resident_tree_wg.tree")

#sub-setting for model
elev_mig_mod_wg <- elev_mig_mod_wg %>% 
  dplyr::select(common_name, scientific_name, migrant, s_m_mig, m_w_mig, 
                w_s_mig, hwi, mass, diet, precip_breadth, temp_breadth, wind_breadth, 
                precip_mean, temp_mean, wind_mean, precip_med, temp_med, wind_med, 
                s_dist_abs, m_dist_abs, w_dist_abs, s_mig, m_mig, w_mig, s_migrant, 
                m_migrant, w_migrant)


#scaling climatic and trait variables
elev_mig_mod_wg[c(7, 8, 10:18)] <- lapply(elev_mig_mod_wg[c(7, 8, 10:18)], 
                                          function(x) c(scale(x)))


#converting diet to factor and setting reference level
elev_mig_mod_wg$diet <- factor(relevel(elev_mig_mod_wg$diet, ref = "Invertivore"))

#running seasonal PGLS models of migration distance ~ traits
#convert to data frame
elev_mig_wg_pgls <- as.data.frame(elev_mig_mod_wg[, -1])

#build comparative data with tree and predictors
elev_mig_wg_pgls <- comparative.data(phy = cons_resident_tree, 
                                     data = elev_mig_wg_pgls, 
                                     names.col = "scientific_name")


#monsoonal shifts
#global model with mean envt tolerances
m_mod_wg_mean <- pgls(m_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                        wind_breadth + precip_mean + temp_mean + wind_mean, 
                      data = elev_mig_wg_pgls, lambda = "ML")#precip mean +ve significant
summary(m_mod_wg_mean)
plot(m_mod_wg_mean)

#global model with median envt tolerances
m_mod_wg_med <- pgls(m_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                       wind_breadth + precip_med + temp_med + wind_med, 
                     data = elev_mig_wg_pgls, lambda = "ML")
summary(m_mod_wg_med)
plot(m_mod)

#computing AIC scores
AIC(m_mod_wg_mean)#1871.603
AIC(m_mod_wg_med)#1862.392

#the model using median Western Ghats envt tolerances fits the data better based on residual
#std.error, adjusted R-squared and AIC. Moving forward to model selection with this model

#model selection
options(na.action = "na.fail")
m_mod_wg_drg <- dredge(m_mod_wg_med)
options(na.action = "na.omit")

#examining sum of weights to assess importance of predictor variables
m_mod_wg_sw <- sw(m_mod_wg_drg)

#model averaging for models with delta AICc < 2
m_mod_wg_avg <- model.avg(m_mod_wg_drg, subset = delta <= 2, fit = TRUE)
summary(m_mod_wg_avg)#11 models, precip_med & wind_med +ve significant


#model diagnostics
#examining VIF, vif() does not take caper::pgls, fitting a nlme::gls
vif(nlme::gls(m_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                wind_breadth + precip_med + temp_med + wind_med, data = elev_mig_mod_wg))
#all GVIFs are lower than 5, no strong signal for multicollinearity

m_predicted <- predict(m_mod_wg_avg)

m_resid <- elev_mig_mod$m_dist_abs - m_predicted

#qqplots
qqnorm(m_resid)
qqline(m_resid)

#residuals vs fitted
plot(m_predicted, m_resid)
abline(h = 0, col = "red")

#best fit model
m_mod_wg_best <- get.models(m_mod_wg_drg, subset = 1)[[1]]
summary(m_mod_wg_best)

#winter shifts
#global model with mean envt tolerances
w_mod_wg_mean <- pgls(w_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                        wind_breadth + precip_mean + temp_mean + wind_mean, data = elev_mig_wg_pgls, 
                      lambda = "ML")#mass +ve significant

summary(w_mod_wg_mean)
plot(w_mod_wg_mean)

#global model with median envt tolerances
w_mod_wg_med <- pgls(w_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                       wind_breadth + precip_med + temp_med + wind_med, data = elev_mig_wg_pgls, 
                     lambda = "ML")#mass & precip_med +ve significant

summary(w_mod_wg_med)
plot(w_mod_wg_med)

#computing AIC scores
AIC(w_mod_wg_mean)#1896.313
AIC(w_mod_wg_med)#1893.027

#model with median envt tolerances fits better based on adjusted R-squared,
#residual std.error and AIC. Moving forward with this model

#model selection
options(na.action = "na.fail")
w_mod_wg_drg <- dredge(w_mod_wg_med)
options(na.action = "na.omit")


#examining sum of weights to assess importance of predictor variables
w_mod_wg_sw <- sw(w_mod_wg_drg)


#model averaging for models with delta AICc < 2
w_mod_wg_avg <- model.avg(w_mod_wg_drg, subset = delta <= 2, fit = TRUE)
summary(w_mod_wg_avg)#9 models, precip_med & wind_med +ve significant


#model diagnostics
#examining vif
vif(nlme::gls(w_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                wind_breadth + precip_med + temp_med + wind_med, data = elev_mig_mod_wg))#nothing above 5

w_predicted <- predict(w_mod_wg_avg)

w_resid <- elev_mig_mod$w_dist_abs - w_predicted

#qqplots
qqnorm(w_resid)
qqline(w_resid)

#residuals vs fitted
plot(w_predicted, w_resid)
abline(h = 0, col = "red")

#best fit model
w_mod_wg_best <- get.models(w_mod_wg_drg, subset = 1)[[1]]
summary(w_mod_wg_best)


#summer shifts
#global model with mean envt tolerances
s_mod_wg_mean <- pgls(s_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                        wind_breadth + precip_mean + temp_mean + wind_mean, 
                      data = elev_mig_wg_pgls, lambda = "ML")#no significant predictors

summary(s_mod_wg_mean)
plot(s_mod_wg_mean)

#global model with median envt tolerances
s_mod_wg_med <- pgls(s_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                       wind_breadth + precip_med + temp_med + wind_med, 
                     data = elev_mig_wg_pgls, lambda = "ML")

summary(s_mod_wg_med)
plot(s_mod_wg_med)

#computing AIC
AIC(s_mod_wg_mean)#1954.261
AIC(s_mod_wg_med)#1955.869

#the model with mean envt tolerances is performing slightly better, but using the 
#median model for consistency

#model selection
options(na.action = "na.fail")
s_mod_wg_drg <- dredge(s_mod_wg_med)
options(na.action = "na.omit")


#examining sum of weights to assess importance of predictor variables
s_mod_wg_sw <- sw(s_mod_wg_drg)

#model averaging for models with delta AICc < 2
s_mod_wg_avg <- model.avg(s_mod_wg_drg, subset = delta <= 2, fit = TRUE)
summary(s_mod_wg_avg)#8 models, precip_med & temp_breadth +ve significant


#model diagnostics
#examining vif
vif(nlme::gls(s_dist_abs ~ hwi + mass + diet + precip_breadth + temp_breadth + 
                wind_breadth + precip_med + temp_med + wind_med, data = elev_mig_mod_wg))#nothing above 5

s_predicted <- predict(s_mod_wg_avg)

s_resid <- elev_mig_mod$s_dist_abs - s_predicted

#qqplots
qqnorm(s_resid)
qqline(s_resid)

#residuals vs fitted
plot(s_predicted, s_resid)
abline(h = 0, col = "red")

#best fit model
s_mod_wg_best <- get.models(s_mod_wg_drg, subset = 1)[[1]]
summary(s_mod_wg_best)#precip_med & temp_breadth +ve significant


#saving results and plotting
summary_m_mod_wg <- summary(m_mod_wg_avg)
m_mod_wg_table <- as.data.frame(cbind(summary_m_mod_wg$coefmat.full, 
                                      confint(m_mod_wg_avg)))
m_mod_wg_table$Predictor <- row.names(m_mod_wg_table)
row.names(m_mod_wg_table) <- NULL
m_mod_wg_table <- m_mod_wg_table[, c(7, 1, 5, 6, 2, 3, 4)]
colnames(m_mod_wg_table) <- c("Predictor", "Estimate", "LCI", "UCI","SE", 
                              "Z.value", "P.value")

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
m_mod_wg_table$Predictor <- sapply(m_mod_wg_table$Predictor, modify_predictor_names)

#winter
summary_w_mod_wg <- summary(w_mod_wg_avg)
w_mod_wg_table <- as.data.frame(cbind(summary_w_mod_wg$coefmat.full, 
                                      confint(w_mod_wg_avg)))
w_mod_wg_table$Predictor <- row.names(w_mod_wg_table)
row.names(w_mod_wg_table) <- NULL
w_mod_wg_table <- w_mod_wg_table[, c(7, 1, 5, 6, 2, 3, 4)]
colnames(w_mod_wg_table) <- c("Predictor", "Estimate", "LCI", "UCI","SE", 
                              "Z.value", "P.value")

#changing predictor names
w_mod_wg_table$Predictor <- sapply(w_mod_wg_table$Predictor, modify_predictor_names)

#summer
summary_s_mod_wg <- summary(s_mod_wg_avg)
s_mod_wg_table <- as.data.frame(cbind(summary_s_mod_wg$coefmat.full, 
                                      confint(s_mod_wg_avg)))
s_mod_wg_table$Predictor <- row.names(s_mod_wg_table)
row.names(s_mod_wg_table) <- NULL
s_mod_wg_table <- s_mod_wg_table[, c(7, 1, 5, 6, 2, 3, 4)]
colnames(s_mod_wg_table) <- c("Predictor", "Estimate", "LCI", "UCI","SE", 
                              "Z.value", "P.value")


#changing predictor names
s_mod_wg_table$Predictor <- sapply(s_mod_wg_table$Predictor, modify_predictor_names)

#save
dist_mod_results_wg <- rbind(s_mod_wg_table, m_mod_wg_table, w_mod_wg_table)
write.csv(dist_mod_results, "output/distance_model_results.csv", row.names = FALSE)

#plotting model estimates
#monsoon
m_mod_wg_plot <- m_mod_wg_table %>% 
  ggplot() +
  aes(y = Predictor, yend = Predictor) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, height = 0.5), linewidth = 0.3) +
  geom_point(aes(x = Estimate)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey20") +
  geom_text(aes(x = UCI, label = ifelse(P.value < 0.05, "*", "")), 
            vjust = 0.5, hjust = -0.5) +
  theme_pubclean(base_size = 8) +
  theme(axis.title = element_blank(), 
        plot.background = element_rect(color = "grey70"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
m_mod_wg_plot

#winter
w_mod_wg_plot <- w_mod_wg_table %>% 
  ggplot() +
  aes(y = Predictor, yend = Predictor) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, height = 0.5), linewidth = 0.3) +
  geom_point(aes(x = Estimate)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey20") +
  geom_text(aes(x = UCI, label = ifelse(P.value < 0.05, "*", "")), 
            vjust = 0.5, hjust = -0.5) +
  theme_pubclean(base_size = 8) +
  theme(axis.title = element_blank(), 
        plot.background = element_rect(color = "grey70"), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
w_mod_wg_plot

#summer
s_mod_wg_plot <- s_mod_wg_table %>% 
  ggplot() +
  aes(y = Predictor, yend = Predictor) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, height = 0.5), linewidth = 0.3) +
  geom_point(aes(x = Estimate)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey20") +
  geom_text(aes(x = UCI, label = ifelse(P.value < 0.05, "*", "")), 
            vjust = 0.5, hjust = -0.5) +
  theme_pubclean(base_size = 8) +
  theme(axis.title = element_blank(), 
        plot.background = element_rect(color = "grey70"), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
s_mod_wg_plot

wg_dist_mod_plot <- plot_grid(
  s_mod_wg_plot, NULL, m_mod_wg_plot, NULL, w_mod_wg_plot, 
  labels = c("(a)", "", "(b)", "", "(c)"), 
  label_fontface = "plain",
  label_size = 8,
  nrow = 5, 
  align = "v",
  rel_heights = c(1, 0.025, 1, 0.025, 1), 
  hjust = -0.5, 
  vjust = 1.5 
)

wg_dist_mod_plot <- ggdraw() +
  draw_plot(wg_mod_plot, 0, 0.025, 1, 0.975) + 
  draw_label(
    "Model Estimate (m) ± CI", 
    x = 0.5, 
    y = 0.02, 
    vjust = 1, 
    hjust = 0.5, 
    size = 8 
  ) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

#creating model selection table
s_mod_sel <- s_mod_wg_drg %>% 
  filter(delta <= 2)
#get r-squared
s_r_sq <- sapply(get.models(s_mod_wg_drg, subset = delta <= 2), 
                 function(m) summary(m)$r.squared)
s_mod_sel$r_sq <- s_r_sq
s_mod_sel$season <- "summer"
s_mod_sel <- as.data.frame(s_mod_sel)

m_mod_sel <- m_mod_wg_drg %>% 
  filter(delta <= 2)
m_r_sq <- sapply(get.models(m_mod_wg_drg, subset = delta <= 2), 
                 function(m) summary(m)$r.squared)
m_mod_sel$r_sq <- m_r_sq
m_mod_sel$season <- "monsoon"
m_mod_sel <- as.data.frame(m_mod_sel)

w_mod_sel <- w_mod_wg_drg %>% 
  filter(delta <= 2)
w_r_sq <- sapply(get.models(w_mod_wg_drg, subset = delta <= 2), 
                 function(m) summary(m)$r.squared)
w_mod_sel$r_sq <- w_r_sq
w_mod_sel$season <- "winter"
w_mod_sel <- as.data.frame(w_mod_sel)

table_mod_sel <- rbind(s_mod_sel, m_mod_sel, w_mod_sel)

#round all numeric columns to 3 decimal places
table_mod_sel[] <- lapply(table_mod_sel, function(x) {
  if (is.numeric(x)) {
    x <- round(x, 3)
  }
  return(x)
})

#replace NAs with "X"
table_mod_sel[] <- lapply(table_mod_sel, function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  x[is.na(x)] <- "X"
  return(x)
})

#write file
write.csv(table_mod_sel, "output/table_mod_sel.csv", row.names = FALSE)