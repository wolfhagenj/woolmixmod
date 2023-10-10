####WoolMixMod Project: Age Modeling####
#This script describes the age modeling analyses in Price and Wolfhagen (submitted). "Wool They, Won't They: Untangling the Zooarchaeological Evidence for Intensive Wool Production in Northern Mesopotamia c. 4500-1500 cal. BC
#Mandible ages from 77 assemblages were collected with the following data:
#The number of caprine (sheep and goat) mandibles aged in Payne Stages G, H, or I (or Habermehl's M3++ or M3+++): N_GHI
#The total number of caprine (sheep and goat) mandibles with age data: N_Total

#These data are then modeled using a Binomial model that allows sites to have individual "theta_site" values
#N_GHI ~ Binomial(N_total, theta_site)

#Site-specific values vary around period-specific means ("period_mu") and standard deviations ("period_sigma")
#logit(theta_site) ~ Normal(period_mu, period_sigma)

#The multilevel structure is thus: assemblages within periods
#See GitHub page (http://www.github.com/wolfhagenj/woolmixmod) for more details about the project, references, scripts, and datasets

####0. Import Libraries####
library("data.table")
library("cmdstanr")
library("rstan")
library("ggdist")
library("ggplot2")
library("ggpubr")
library("Cairo")

####1. Importing and Wrangling Data####
wool_mandible_data <- fread("./Data/wool_mandible_data.csv")

#Put in the period codes
wool_mandible_data[, Period := factor(Period, levels = c("Neolithic", "Halaf-Ubaid", "LC 1-3", "LC 4-5", "EBA", "MBA"))]
wool_mandible_data[, Period_No := as.numeric(Period)]

####2. Fitting the Bayesian model####
#Setting up the data for the Stan model
#See associated Stan model (GHI_multiperiod_normal_binomial.stan) to see details of the model
wool_mandible_age_standata <- list(
  N_Observations = wool_mandible_data[, .N],
  N_Sites = wool_mandible_data[, .N, Site_No][, .N],
  N_Periods = wool_mandible_data[, .N, Period_No][, .N],
  Site = wool_mandible_data[order(Site_No), Site_No],
  Period = wool_mandible_data[, .N, .(Site_No, Period_No)][order(Site_No), Period_No],
  N_GHI = wool_mandible_data[order(Site_No), N_GHI],
  N_Total_Mandibles = wool_mandible_data[order(Site_No), N_Total]
)

#
mandible_binomial_model <- cmdstan_model("./Scripts/GHI_multiperiod_normal_binomial.stan")
wool_mandible_age_samples <- mandible_binomial_model$sample(
  data = wool_mandible_age_standata,
  chains = 4,
  parallel_chains = 4,
  refresh = 250,
  adapt_delta = 0.99,
  seed = 495625296,
  max_treedepth = 15
)
wool_mandible_age_samples$summary(c("theta_average", "theta_period"))
wool_age_stanfit <- rstan::read_stan_csv(wool_mandible_age_samples$output_files())
wool_age_post <- extract(wool_age_stanfit)

####3. Posterior Analysis####
wool_age_site_posterior_data <- data.table(Iteration = rep(1:4000, wool_mandible_data[, .N, Site_No][, .N]),
                                      Period = rep(wool_mandible_data[, .N, .(Site_No, Period, Period_No)][order(Site_No), Period], each = 4000),
                                      Period_No = rep(wool_mandible_data[, .N, .(Site_No, Period, Period_No)][order(Site_No), Period_No], each = 4000),
                                      Site = rep(wool_mandible_data[, .N, .(Site_No, Site_Name = paste(Period, Site, sep = "."))][order(Site_No), Site_Name], each = 4000),
                                      Site_No = rep(wool_mandible_data[, .N, .(Site_No, Site_Name = paste(Period, Site, sep = "."))][order(Site_No), Site_No], each = 4000),
                                      theta = c(wool_age_post$theta_site))
wool_age_period_posterior_data <- data.table(Iteration = rep(1:4000, wool_mandible_data[, .N, Period][, .N]),
                                        Period = rep(wool_mandible_data[, .N, .(Period, Period_No)][order(Period_No), Period], each = 4000),
                                        Period_No = rep(wool_mandible_data[, .N, .(Period, Period_No)][order(Period_No), Period_No], each = 4000),
                                        theta_mean = c(wool_age_post$theta_period),
                                        theta_random = c(wool_age_post$theta_period_random))

#Table 2: Empirical and modeled frequencies of assemblages with at least 40% old (Payne GHI) caprine mandibles
empirical_frequencies <- wool_mandible_data[, .(`N Sites` = .N, `N Old` = sum(N_GHI / N_Total >= 0.4), `Empirical Percentage` = paste0(100 * round(mean(N_GHI / N_Total >= 0.4), 2), "%")), .(Period)][order(Period), .(Period, `Empirical Percentage` = paste0(`Empirical Percentage`, " (", `N Old`, " / ", `N Sites`, ")"))]
modeled_frequencies <- wool_age_site_posterior_data[, .(.N, theta = mean(theta), pct_over40 = mean(theta > 0.4)), .(Iteration, Period, Period_No)][, .(`Modeled Percentage` = paste0(100 * round(mean(pct_over40), 2), "%"), `Confidence Interval (Modeled Percentage)` = paste0(100 * round(quantile(pct_over40, 0.025), 2), "-", 100 * round(quantile(pct_over40, 0.975), 2), "%")), .(Period, N)][order(Period), .(Period, `Modeled Percentage` = paste0(`Modeled Percentage`, " (", `Confidence Interval (Modeled Percentage)`, ")"), `Expected Percentage (Model)` = paste0(100 * round(colMeans(wool_age_post$theta_period_random > 0.4), 2), "%"))]
table_2 <- empirical_frequencies[modeled_frequencies, on = "Period"]
write.csv(table_2, file = "./Outputs/Table 2 - Proportions of Old Sites.csv", row.names = F)

#names to produce for site-level plots
ageplot_names <- data.table(
  plot_order = wool_age_site_posterior_data[, .(Iteration, Period, Period_No, Site, Site_No, theta, plot_order = paste0(Period_No, Site_No))][, .N, .(plot_order, Site)][order(plot_order), plot_order],
  Display_Name = wool_mandible_data[, .N, .(Site_No, Site_Name = paste(Period, Site, sep = "."), Site)][order(Site_No), Site]
)

#Figure 7: Posterior distributions of the average proportion of old (% Payne GHI) mandibles per site, grouped by analytical period
site_ageplot_neochalco <- ggplot(wool_age_site_posterior_data[Period_No %in% 1:4][, .(Iteration, Period, Period_No, Site, Site_No, theta, plot_order = paste0(Period_No, Site_No))]) + aes(y = theta, x = plot_order) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  geom_hline(yintercept = 0.4, linetype = "dashed") +
  scale_fill_manual(name = "Period", drop = F, values = c("Neolithic" = "#F8766D", "Halaf-Ubaid" = "#B79F00", "LC 1-3" = "#00BA38", "LC 4-5" = "#00BFC4", "EBA" = "#619CFF", "MBA" = "#F564E3")) +
  scale_color_manual(name = "Period", drop = F, values = c("Neolithic" = "#F8766D", "Halaf-Ubaid" = "#B79F00", "LC 1-3" = "#00BA38", "LC 4-5" = "#00BFC4", "EBA" = "#619CFF", "MBA" = "#F564E3")) +
  scale_x_discrete(name = "", breaks = ageplot_names[order(plot_order), plot_order], labels = ageplot_names[order(plot_order), Display_Name]) +
  scale_y_continuous(name = "% Old Adult (>= GHI)", labels = scales::percent, limits = c(0, 1.02)) + coord_cartesian(expand = FALSE) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
site_ageplot_ba <- ggplot(wool_age_site_posterior_data[Period_No %in% 5:6][, .(Iteration, Period, Period_No, Site, Site_No, theta, plot_order = paste0(Period_No, Site_No))]) + aes(y = theta, x = plot_order) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  geom_hline(yintercept = 0.4, linetype = "dashed") +
  scale_fill_manual(name = "Period", drop = F, values = c("Neolithic" = "#F8766D", "Halaf-Ubaid" = "#B79F00", "LC 1-3" = "#00BA38", "LC 4-5" = "#00BFC4", "EBA" = "#619CFF", "MBA" = "#F564E3")) +
  scale_color_manual(name = "Period", drop = F, values = c("Neolithic" = "#F8766D", "Halaf-Ubaid" = "#B79F00", "LC 1-3" = "#00BA38", "LC 4-5" = "#00BFC4", "EBA" = "#619CFF", "MBA" = "#F564E3")) +
  scale_x_discrete(name = "", breaks = ageplot_names[order(plot_order), plot_order], labels = ageplot_names[order(plot_order), Display_Name]) +
  scale_y_continuous(name = "% Old Adult (>= GHI)", labels = scales::percent, limits = c(0, 1.02)) + coord_cartesian(expand = FALSE) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
site_ageplot_sites <- ggarrange(site_ageplot_neochalco, site_ageplot_ba, ncol = 1, nrow = 2, common.legend = T, legend = "bottom")
Cairo(width = 8, height = 10, units = "in", file = "./Outputs/Figure 7 - age plot (site).png", type = "png", dpi = 600)
site_ageplot_sites
dev.off()

#Figure 8: Posterior distributions of the average proportion (top) and expected proportion (bottom) of old (% Payne GHI) mandibles by analytical period
period_ageplot <- ggplot(wool_age_period_posterior_data[, .(theta = theta_mean), .(Iteration, Period, Period_No, plot_order = paste0(Period_No, Period))]) + aes(y = theta, x = Period) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  geom_hline(yintercept = 0.4, linetype = "dashed") +
  labs(title = "Average %GHI") +
  scale_y_continuous(name = "% Old Adult (>= GHI)", labels = scales::percent, limits = c(0, 1.02)) + coord_cartesian(expand = FALSE) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
period_random_ageplot <- ggplot(wool_age_period_posterior_data[, .(theta = theta_random), .(Iteration, Period, Period_No)]) + aes(y = theta, x = Period) +
  stat_slab(normalize = "groups", aes(fill = Period, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(1, 0.95, 0.80), labels = scales::percent_format())))) +
  stat_pointinterval(.width = c(0.8, 0.95), aes(color = Period), position = position_dodge(width = 0.2, preserve = "single")) +
  scale_fill_ramp_discrete(name = "Interval", range = c(0.8, 0.3), na.translate = F) +
  geom_hline(yintercept = 0.4, linetype = "dashed") +
  labs(title = "Expected %GHI Distribution") +
  scale_y_continuous(name = "% Old Adult (>= GHI)", labels = scales::percent, limits = c(0, 1.02)) + coord_cartesian(expand = FALSE) +
  theme_classic() + theme(legend.text = element_text(size = 10), legend.position = "bottom", legend.justification = c("center"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6), axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 0.95), axis.title = element_text(size = 12), axis.text.y = element_text(angle = 0, size = 10, hjust = 0.5))
#
period_theta_plot <- ggarrange(period_ageplot, period_random_ageplot, ncol = 1, nrow = 2, common.legend = T, legend = "bottom")
Cairo(width = 8, height = 10, units = "in", file = "./Outputs/Figure 8 - age plot (period - average and expected).png", type = "png", dpi = 600)
period_theta_plot
dev.off()

####4. Supplemental Table File####

#Supplemental Table 1: mandible data
write.csv(wool_mandible_data, file = "./Outputs/Supplemental Table 1 - Wool Aging Model Mandible Data.csv", row.names = F)

####5. Additional Result Table Files####

#Result Table 4: aging results (site-level)
write.csv(wool_age_site_posterior_data, file = "./Outputs/Result Table 4 - Posterior Aging Data (site).csv", row.names = F)

#Result Table 5: aging results (period-level)
write.csv(wool_age_period_posterior_data, file = "./Outputs/Result Table 5 - Posterior Aging Data (period).csv", row.names = F)
