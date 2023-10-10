####WoolMixMod Project: Biometric Modeling####
#This script describes the biometric analyses in Price and Wolfhagen (submitted). "Wool They, Won't They: Untangling the Zooarchaeological Evidence for Intensive Wool Production in Northern Mesopotamia c. 4500-1500 cal. BC
#Standard faunal measurements (following von den Driesch, 1976) were collected from 2357 sheep bones from 32 assemblages in northern Mesopotamia
#These bones were analyzed using a multilevel mixture model of LSI values that viewed assemblages as a mixture of immature, adult female, and adult male specimens (following Wolfhagen 2023)
#The multilevel structure is thus: element portions within assemblages, assemblages within periods
#See GitHub page (http://www.github.com/wolfhagenj/woolmixmod) for more details about the project, references, scripts, and datasets

####0. Importing Libraries####
library("data.table")
library("cmdstanr")
library("ggplot2")
library("ggpubr")
library("rnaturalearth")
library("rnaturalearthdata")
library("rstan")
library("sf")
library("zoolog")
library("ggdist")
library("Cairo")

####1. Importing and Wrangling Data####

#Calculate LSI (log_e) values using the Clutton-Brock, et al. (1990) animal
sheep_standard_animal <- data.table(zoolog::referencesDatabase$`Ovis aries`$Clutton)[EL %in% c("Humerus", "Radius", "Metacarpus", "Metatarsus", "Tibia", "Astragalus", "Calcaneus", "P1 ant", "P1 post", "P2 ant", "P2 post")]
sheep_standard_animal[EL %in% "Humerus", Element := "humerus"]
sheep_standard_animal[EL %in% "Radius", Element := "radius bone"]
sheep_standard_animal[EL %in% "Metacarpus", Element := "fused metacarpal bones 3 and 4"]
sheep_standard_animal[EL %in% "Tibia", Element := "tibia"]
sheep_standard_animal[EL %in% "Astragalus", Element := "talus"]
sheep_standard_animal[EL %in% "Calcaneus", Element := "calcaneus"]
sheep_standard_animal[EL %in% "Metatarsus", Element := "fused metatarsal bones 3 and 4"]
sheep_standard_animal[EL %in% c("P1 ant", "P1 post"), Element := "phalanx 1"]
sheep_standard_animal[EL %in% c("P2 ant", "P2 post"), Element := "phalanx 2"]

#organizing measurement names to include element portion
sheep_standard_animal[Element %in% "humerus" & Measure %in% "Bd", Measurement := "Hum_Bd"]
sheep_standard_animal[Element %in% "humerus" & Measure %in% "BT", Measurement := "Hum_BT"]
sheep_standard_animal[Element %in% "radius bone" & Measure %in% "Bp", Measurement := "Rad_Bp"]
sheep_standard_animal[Element %in% "radius bone" & Measure %in% "Bd", Measurement := "Rad_Bd"]
sheep_standard_animal[Element %in% "fused metacarpal bones 3 and 4" & Measure %in% "Bd", Measurement := "Mtc_Bd"]
sheep_standard_animal[Element %in% "tibia" & Measure %in% "Bd", Measurement := "Tib_Bd"]
sheep_standard_animal[Element %in% "talus" & Measure %in% "Bd", Measurement := "Ast_Bd"]
sheep_standard_animal[Element %in% "calcaneus" & Measure %in% "GB", Measurement := "Cal_GB"]
sheep_standard_animal[Element %in% "fused metatarsal bones 3 and 4" & Measure %in% "Bd", Measurement := "Mtt_Bd"]
sheep_standard_animal[Element %in% "phalanx 1" & Measure %in% "GL", Measurement := "PH1_Glpe"]
sheep_standard_animal[Element %in% "phalanx 1" & Measure %in% "Bp", Measurement := "PH1_Bp"]
sheep_standard_animal[Element %in% "phalanx 2" & Measure %in% "GL", Measurement := "PH2_GL"]
sheep_standard_animal[Element %in% "phalanx 2" & Measure %in% "Bp", Measurement := "PH2_Bp"]

#average the measurements of the anterior and posterior phalanges (PH1, PH2)
sheep_standard_animal <- rbind(sheep_standard_animal[grepl("P2", EL) == F & grepl("P1", EL) == F],
                               sheep_standard_animal[grepl("P1", EL), .(EL = "P1 Average", Standard = mean(Standard)), .(TAX, Measure, Element, Measurement)][, .(TAX, EL, Measure, Standard, Element, Measurement)],
                               sheep_standard_animal[grepl("P2", EL), .(EL = "P2 Average", Standard = mean(Standard)), .(TAX, Measure, Element, Measurement)][, .(TAX, EL, Measure, Standard, Element, Measurement)])

#Archaeological data (measurement data and demographic observations for sites)
#NOTE: three assemblages that were in the original paper are missing from this dataset
#Thus there are gaps in the "Specimen" and "Site_No" data values
wool_mixmod_data <- fread("./Data/woolmixmod_metric_data.csv")

#match the standard measurement (from reference animal) to the archaeological measurements
wool_mixmod_data <- wool_mixmod_data[sheep_standard_animal[!is.na(Measurement), .(Measurement, Reference_value = Standard)], on = c("Measurement")][!is.na(Specimen)]
#calculate the LSI_e values (for visualization, raw measurements are used for modeling)
wool_mixmod_data[, "LSI_ln" := log(Measurement_value / Reference_value)]

#Demographic data (observations of unfused phalanges and sexed pelves)
#NOTE: three assemblages that were in the original paper are missing from this dataset
#Thus there are gaps in the "Specimen" and "Site_No" data values
wool_demographic_observations <- fread("./Data/woolmixmod_demographic_data.csv")

#Sidebar: visualizing the sites#
world <- ne_countries(scale = "medium", returnclass = "sf") #for visualization/examining spatial trends
#
ggplot(data = world) + geom_sf(color = "black", fill = "grey80") +
  coord_sf(xlim = c(25.5, 46), ylim = c(32, 45), expand = FALSE) +
  geom_rect(xmin = 36, xmax = 45.25, ymin = 34, ymax = 38.5, fill = NA, linetype = "dashed", col = "black") +
  geom_point(data = wool_demographic_observations, aes(x = Longitude, y = Latitude, col = wool_mixmod_data[, .(LSI = mean(LSI_ln)), Site_No][order(Site_No), LSI])) + scale_color_continuous(name = "LSI")
ggplot(data = world) + geom_sf(color = "black", fill = "grey80") +
  coord_sf(xlim = c(36, 45.25), ylim = c(34, 38.5), expand = FALSE) +
  geom_point(data = wool_demographic_observations, aes(x = Longitude, y = Latitude, col = wool_mixmod_data[, .(LSI = mean(LSI_ln)), Site_No][order(Site_No), LSI])) + scale_color_continuous(name = "LSI")

#Assigning sites to analytical periods (metric data and demographic data)
wool_mixmod_data[Period %in% c("PPNB", "PN"), Period_No := 1]
wool_mixmod_data[Period %in% "Halaf", Period_No := 2]
wool_mixmod_data[Period %in% "LC 1-3", Period_No := 3]
wool_mixmod_data[Period %in% "LC 4-5", Period_No := 4]
wool_mixmod_data[Period %in% "EBA", Period_No := 5]
wool_mixmod_data[Period %in% "MBA", Period_No := 6]
#
wool_demographic_observations[Period %in% c("PPNB", "PN"), Period_No := 1]
wool_demographic_observations[Period %in% "Halaf", Period_No := 2]
wool_demographic_observations[Period %in% "LC 1-3", Period_No := 3]
wool_demographic_observations[Period %in% "LC 4-5", Period_No := 4]
wool_demographic_observations[Period %in% "EBA", Period_No := 5]
wool_demographic_observations[Period %in% "MBA", Period_No := 6]

#Troubleshooting
wool_mixmod_data[Specimen %in% "Animal.7752", Measurement_value := NA] #had been 29.1 originally, removing because it is implausibly large for an Astragalus Bd
#Removing the specimen
wool_mixmod_data <- wool_mixmod_data[!is.na(Period_No) & !is.na(Measurement_value)]

#Ensuring that numerical categories do not have any missing numbers (i.e., go from 1-max)
wool_mixmod_data[, Site_No := as.numeric(as.factor(Site_No))]
wool_mixmod_data[, Specimen_No := as.numeric(as.factor(Specimen_No))]
wool_mixmod_data[, Element_Portion := as.numeric(as.factor(Element_Portion))]
wool_mixmod_data[, Measurement_Set := as.numeric(as.factor(Measurement_Set))]
#
wool_demographic_observations[, Site_No := as.numeric(as.factor(Site_No))]

#Sidebar: Evaluating the dataset--how many assemblages are small (N_Specimens < 20), medium (20 <= N_Specimens < 50), or large (N_Specimens >= 50)?#
dcast(wool_mixmod_data[, .N, .(Specimen_No, Site_No, Site, Period, Period_No)][, .N, .(Site_No, Site, Period, Period_No)][, .(Size = ifelse(N < 20, "Small", ifelse(N < 50, "Medium", "Large"))), .(Site_No, Site, Period, Period_No)][, .N, .(Size = factor(Size, levels = c("Small", "Medium", "Large")), Period, Period_No)], Period_No + Period ~ Size, value.var = "N", fill = 0)

####2. Fitting the Bayesian model####
#Setting up the data for the Stan model
#Model is adapted from Wolfhagen (2023) to group sites by periods
#See associated Stan model (LSI_mixture_model_multiperiod.stan) to see details of the model

wool_multi_period_standata <- list(
  #Sample sizes
  N_Periods = wool_mixmod_data[, .N, Period_No][, .N], #added for this model
  N_Sites = wool_mixmod_data[, .N, Site_No][, .N],
  N_Specimens = wool_mixmod_data[, .N, Specimen_No][, .N],
  N_Measurements = wool_mixmod_data[, .N],
  N_Element_Portions = wool_mixmod_data[, .N, Element_Portion][, .N],
  N_Dimensions = wool_mixmod_data[, .N, Measurement_Set][, .N],
  #Site observations
  Period = wool_mixmod_data[, .N, .(Site_No, Period_No)][order(Site_No), Period_No], #added for this model
  #Specimen observations
  Site = wool_mixmod_data[, .N, .(Specimen_No, Element_Portion, Site_No, Immature)][order(Specimen_No), Site_No],
  Element_Portion = wool_mixmod_data[, .N, .(Specimen_No, Element_Portion, Site_No, Immature)][order(Specimen_No), Element_Portion],
  Immature = wool_mixmod_data[, .N, .(Specimen_No, Element_Portion, Site_No, Immature)][order(Specimen_No), Immature],
  Immature_Proportion = as.matrix(dcast(wool_mixmod_data[, .(Immature_Proportion = mean(Immature)), .(Site_No, Element_Portion)], Site_No ~ Element_Portion, value.var = "Immature_Proportion", fill = 0))[, 2:(wool_mixmod_data[, .N, Element_Portion][, .N] + 1)],
  #Measurement observations
  Measurement_obs = wool_mixmod_data[, Measurement_value],
  Measurement_sd = wool_mixmod_data[, Measurement_value * 0.01], #Calculate measurement error for observed measurements and reference data (1% based on data from Breslawksi and Byers 2015)
  Reference_obs = wool_mixmod_data[, .N, .(Measurement_Set, Reference_value)][order(Measurement_Set), Reference_value],
  Reference_sd = wool_mixmod_data[, .N, .(Measurement_Set, Reference_value)][order(Measurement_Set), Reference_value * 0.01],
  Dimension = wool_mixmod_data[, Measurement_Set],
  Specimen = wool_mixmod_data[, Specimen_No],
  #Demographic observations
  N_Immature_obs = wool_demographic_observations[, .N],
  Immature_obs_site = wool_demographic_observations[, Site_No],
  Immature_obs = wool_demographic_observations[, N_Unfused],
  Immature_obs_n = wool_demographic_observations[, N_Ageable],
  N_Female_obs = wool_demographic_observations[, .N],
  Female_obs_site = wool_demographic_observations[, Site_No],
  Female_obs = wool_demographic_observations[, N_Female],
  Female_obs_n = wool_demographic_observations[, N_Sexable],
  #Prior distributions for hyper-parameters
  prior_theta_raw_1 = c(-0.5, 1.5),
  prior_theta_raw_2 = c(0, 1.5),
  prior_mu_female = c(0, 0.2),
  prior_logdelta_immature = c(-3.5, 0.25),
  prior_logdelta_male = c(-2.7, 0.2),
  prior_logsigma_immature = c(-3.05, 0.1),
  prior_logsigma_female = c(-3.1, 0.1),
  prior_logsigma_male = c(-3.1, 0.1)
)
#make sure the model is compiled
LSI_multi_period_model <- cmdstan_model("./Scripts/LSI_mixture_model_multiperiod.stan")
#fit the data, producing posterior samples of the parameters
wool_multi_period_samples <- LSI_multi_period_model$sample(
  data = wool_multi_period_standata,
  chains = 4,
  parallel_chains = 4,
  refresh = 250,
  adapt_delta = 0.95,
  seed = 668160008, #this seed ensures strict reproducibility, remove to get conceptual reproducibility/apply to new datasets
  max_treedepth = 15
)

#evaluate model fit
wool_multi_period_samples$summary(c("grand_theta_raw",
                                    "grand_mu_female", "grand_logdelta_immature", "grand_logdelta_male",
                                    "grand_logsigma_immature", "grand_logsigma_female", "grand_logsigma_male"))
#save the posterior distributions
wool_multi_period_stanfit <- rstan::read_stan_csv(wool_multi_period_samples$output_files())
wool_multi_period_post <- extract(wool_multi_period_stanfit)

####3. Posterior Analysis####
#collect the posterior data at both the site-level and period-level
wool_metric_period_posterior_data <- data.table(Iteration = rep(1:4000, 6), Period_No = rep(1:6, each = 4000),
                                         Period = rep(factor(c("Neolithic", "Halaf", "LC 1-3", "LC 4-5", "EBA", "MBA"), levels = c("Neolithic", "Halaf", "LC 1-3", "LC 4-5", "EBA", "MBA")), each = 4000),
                                         p_immature = c(wool_multi_period_post$period_p_immature),
                                         sex_ratio = c(wool_multi_period_post$period_theta_female),
                                         p_male = c(wool_multi_period_post$period_theta[, 3, ]),
                                         mu_female = c(wool_multi_period_post$period_mu_female),
                                         mu_male = c(wool_multi_period_post$period_mu_male),
                                         delta_male = c(wool_multi_period_post$period_mu_male) - c(wool_multi_period_post$period_mu_female))
#
wool_metric_site_posterior_data <- data.table(Iteration = rep(1:4000, wool_demographic_observations[, .N]), Site_No = rep(wool_demographic_observations[order(Site_No), Site_No], each = 4000),
                                       Site = rep(wool_demographic_observations[order(Site_No), Site], each = 4000),
                                       Site_No = rep(wool_demographic_observations[order(Site_No), Site_No], each = 4000),
                                       Period = rep(wool_demographic_observations[order(Site_No), factor(ifelse(Period %in% c("PPNB", "PN"), "Neolithic", Period), levels = c("Neolithic", "Halaf", "LC 1-3", "LC 4-5", "EBA", "MBA"))], each = 4000),
                                       Period_No = rep(wool_demographic_observations[order(Site_No), as.numeric(factor(ifelse(Period %in% c("PPNB", "PN"), "Neolithic", Period), levels = c("Neolithic", "Halaf", "LC 1-3", "LC 4-5", "EBA", "MBA")))], each = 4000),
                                       p_immature = c(wool_multi_period_post$site_p_immature),
                                       sex_ratio = c(wool_multi_period_post$site_theta_female),
                                       p_male = c(wool_multi_period_post$site_theta[, 3, ]),
                                       mu_female = c(wool_multi_period_post$site_mu_female),
                                       mu_male = c(wool_multi_period_post$site_mu_male),
                                       delta_male = c(wool_multi_period_post$site_mu_male) - c(wool_multi_period_post$site_mu_female))
#Specimen-level results (membership probabilities for each specimen)
wool_metric_specimen_posterior_data <- wool_mixmod_data[, .N, .(Specimen_No, Specimen, Site_No, Site, Period_No, Period, Element_Portion, Element)][order(Specimen_No), .(Specimen_No, Specimen, Site_No, Site, Period_No, Period, Element_Portion, Element, LSI = round(colMeans(wool_multi_period_post$LSI), 3), p_immature = round(colMeans(wool_multi_period_post$specimen_prob[, , 1]), 2), p_female = round(colMeans(wool_multi_period_post$specimen_prob[, , 2]), 2), p_male = round(colMeans(wool_multi_period_post$specimen_prob[, , 3]), 2))]

#Figure 1: LSI values for sheep measurements in the biometric dataset (in LSIe scale), grouped by analytical period
site_lsiplot <-
  ggplot(wool_mixmod_data[, .(Site, Site_No, Period = factor(ifelse(Period %in% c("PPNB", "PN"), "Neolithic", Period), levels = c("Neolithic", "Halaf", "LC 1-3", "LC 4-5", "EBA", "MBA")), Period_No, LSI_ln, Plot_Order = paste0(Period_No, ".", Site_No))]) + aes(y = LSI_ln, x = Plot_Order) +
  stat_halfeye(adjust=1.25, aes(fill = Period), orientation = "vertical", justification = .11) +
  stat_dots(side="left", aes(fill = Period), binwidth=0.002, justification = 1.11, col = "black") +
  geom_boxplot(width = .2, aes(fill = Period), alpha = 0.3, orientation = "vertical", outlier.shape = NA)  +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_y_continuous(name = "LSI Value") +
  scale_x_discrete(name = "", breaks = wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio, p_male, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))][, .N, .(Plot_Order, Site)][order(Plot_Order), Plot_Order], labels = wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio, p_male, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))][, .N, .(Plot_Order, Site)][order(Plot_Order), Site]) +
  #  coord_cartesian(ylim = c(-0.1, 0.3)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
Cairo(width = 8, height = 5, units = "in", file = "./Outputs/Figure 1 - LSI plot.png", type = "png", dpi = 600)
site_lsiplot
dev.off()

#Figure 2: Histogram of LSIe values for sheep measurements from the Late Chalcolithic 1-3 site Tell Surezha with mixture modeling results.
result_histogram_function <- function(mixmod_data, posterior, site_no) {
  site_name <- mixmod_data[Site_No %in% site_no, .N, Site][, Site]
  
  #get the mixture model parameters (posterior means)
  immature_parameters <- c(mean(posterior$site_theta[, 1, site_no]), mean(posterior$site_mu_immature[, site_no]), mean(posterior$site_sigma_immature[, site_no]))
  female_parameters <- c(mean(posterior$site_theta[, 2, site_no]), mean(posterior$site_mu_female[, site_no]), mean(posterior$site_sigma_female[, site_no]))
  male_parameters <- c(mean(posterior$site_theta[, 3, site_no]), mean(posterior$site_mu_male[, site_no]), mean(posterior$site_sigma_male[, site_no]))
  fitted_normal_mix_values <- data.table(
    Group = rep(factor(c("Immature", "Female", "Male", "Total"), levels = c("Immature", "Female", "Male", "Total")), each = 1000),
    X_Value = rep(seq(-0.40, 0.40, length.out = 1000), 4),
    Fit_Value = c(immature_parameters[1] * dnorm(seq(-0.40, 0.40, length.out = 1000), immature_parameters[2], immature_parameters[3]),
                  female_parameters[1] * dnorm(seq(-0.40, 0.40, length.out = 1000), female_parameters[2], female_parameters[3]),
                  male_parameters[1] * dnorm(seq(-0.40, 0.40, length.out = 1000), male_parameters[2], male_parameters[3]),
                  immature_parameters[1] * dnorm(seq(-0.40, 0.40, length.out = 1000), immature_parameters[2], immature_parameters[3]) + female_parameters[1] * dnorm(seq(-0.40, 0.40, length.out = 1000), female_parameters[2], female_parameters[3]) + male_parameters[1] * dnorm(seq(-0.1, 0.15, length.out = 1000), male_parameters[2], male_parameters[3]))
  )
  
  #calculate the most likely identity for the specimens (based on posterior mean probabilities)
  mixmod_data[, Group := factor(ifelse(p_immature >= p_female & p_immature >= p_male, "Immature", ifelse(p_female >= p_male, "Female", "Male")), levels = c("Immature", "Female", "Male"))]
  
  #update the fitted values to match the scale of the histogram
  basic_histogram <- ggplot() +
    geom_histogram(data = mixmod_data[, .(Group, LSI)], binwidth = 0.03, aes(x = LSI), position = "stack", alpha = 0.8) +
    theme_classic()
  fitted_normal_mix_values[, Scaled_Fit := .(Fit_Value * (max(layer_data(basic_histogram)$count) / max(Fit_Value)))]
  
  #draw up the histogram
  fitted_histogram <- ggplot() +
    geom_histogram(data = mixmod_data, binwidth = 0.03, aes(x = LSI), position = "stack", alpha = 0.3) +
    geom_line(data = fitted_normal_mix_values[Group %in% c("Immature", "Female", "Male")], aes(x = X_Value, y = Scaled_Fit, color = Group), lwd = 1, linetype = "dashed") +
    scale_color_manual(name = "Group", values = c("black", "blue", "red"), labels = c("Immature", "Female", "Male")) +
    labs(x = "LSI Value", y = "Frequency", title = site_name) +
    coord_cartesian(xlim = c(-0.3, 0.3)) +
    theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
  
  fitted_histogram
}
site_plot_no <- 25
surezha_histogram <- result_histogram_function(wool_metric_specimen_posterior_data[Site_No %in% site_plot_no], posterior = wool_multi_period_post, site_no = site_plot_no)
#
Cairo(width = 8, height = 5, units = "in", file = "./Outputs/Figure 2 - example result histogram (Tell Surezha).png", type = "png", dpi = 600)
surezha_histogram
dev.off()

#Figure 3: Posterior distributions of modeled average female body size (in LSIe scale) for each period
period_sizeplot <- ggplot(wool_metric_period_posterior_data[, .(Iteration, Period, Period_No, p_immature, sex_ratio, mu_female, mu_male, delta_male)]) + aes(y = mu_female, x = Period) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_y_continuous(name = "Female Average LSI Value") +
  scale_x_discrete(name = "", breaks = c("Neolithic", "Halaf", "LC 1-3", "LC 4-5", "EBA", "MBA")) +
  coord_cartesian(ylim = c(-0.1, 0.3)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
Cairo(width = 8, height = 5, units = "in", file = "./Outputs/Figure 3 - female size (period).png", type = "png", dpi = 600)
period_sizeplot
dev.off()

#Figure 4: Posterior distributions of modeled average female body size (in LSIe scale) for each site, grouped by analytical period
site_sizeplot <- ggplot(wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))]) + aes(y = mu_female, x = Plot_Order) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_y_continuous(name = "Female Average LSI Value") +
  scale_x_discrete(name = "", breaks = wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio, p_male, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))][, .N, .(Plot_Order, Site)][order(Plot_Order), Plot_Order], labels = wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio, p_male, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))][, .N, .(Plot_Order, Site)][order(Plot_Order), Site]) +
  coord_cartesian(ylim = c(-0.1, 0.3)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
Cairo(width = 8, height = 5, units = "in", file = "./Outputs/Figure 4 - female size (site).png", type = "png", dpi = 600)
site_sizeplot
dev.off()

#Figure 5: Posterior distributions of modeled average adult sex ratios (% males out of adult females and adult males) for each site, grouped by analytical period
site_sexratioplot <- ggplot(wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio = 1 - sex_ratio, p_male, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))]) + aes(y = sex_ratio, x = Plot_Order) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_y_continuous(name = "Sex Ratio (% Male)") +
  scale_x_discrete(name = "", breaks = wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio, p_male, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))][, .N, .(Plot_Order, Site)][order(Plot_Order), Plot_Order], labels = wool_metric_site_posterior_data[, .(Iteration, Site, Site_No, Period, Period_No, p_immature, sex_ratio, p_male, mu_female, mu_male, delta_male, Plot_Order = paste0(Period_No, ".", Site_No))][, .N, .(Plot_Order, Site)][order(Plot_Order), Site]) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
Cairo(width = 8, height = 5, units = "in", file = "./Outputs/Figure 5 - sex ratio (site).png", type = "png", dpi = 600)
site_sexratioplot
dev.off()

#Figure 6: Posterior distributions of modeled average adult sex ratios (% males out of adult females and adult males) for each period
period_sexratioplot <- ggplot(wool_metric_period_posterior_data[, .(Iteration, Period, Period_No, p_immature, sex_ratio = 1 - sex_ratio, p_male, mu_female, mu_male, delta_male)]) + aes(y = sex_ratio, x = Period) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  scale_y_continuous(name = "Sex Ratio (% Male)") +
  scale_x_discrete(name = "", breaks = c("Neolithic", "Halaf", "LC 1-3", "LC 4-5", "EBA", "MBA")) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
Cairo(width = 8, height = 5, units = "in", file = "./Outputs/Figure 6 - sex ratio (period).png", type = "png", dpi = 600)
period_sexratioplot
dev.off()

####4. Supplemental Table Files####

#Supplemental Table 2: metric data
write.csv(wool_mixmod_data, file = "./Outputs/Supplemental Table 2 - Wool Biometric Model Metric Data.csv", row.names = F)

#Supplemental Table 3: assemblage info
write.csv(wool_demographic_observations, file = "./Outputs/Supplemental Table 3 - Wool Biometric Model Demographic Data.csv", row.names = F)

####5. Additional Result Tables (online-only)####

#Result Table 1: metric results (site-level)
write.csv(wool_metric_site_posterior_data, file = "./Outputs/Result Table 1 - Posterior Metric Data (site).csv", row.names = F)

#Result Table 2: metric results (period-level)
write.csv(wool_metric_period_posterior_data, file = "./Outputs/Result Table 2 - Posterior Metric Data (period).csv", row.names = F)

#Result Table 3: metric results (average membership probabilities for each specimen)
write.csv(wool_metric_specimen_posterior_data, file = "./Outputs/Result Table 3 - Posterior Metric Data (specimen).csv", row.names = F)
