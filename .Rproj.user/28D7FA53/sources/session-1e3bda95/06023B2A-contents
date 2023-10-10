#Prior predictive modeling: age model (Binomial-Normal)
#The goal with this prior predictive check is to evaluate the consistency
#of the priors for the model to explore the parameter space effectively

####0. Importing Libraries####
library("boot")
library("data.table")
library("doParallel")
library("ggplot2")
library("ggdist")
library("parallel")

#setting resources to do simulations with multicore processing (on Windows)
num_cores <- detectCores(logical = T)
cl <- makeCluster(num_cores - 4) #creating a virtual cluster
registerDoParallel(cl) #registering the cluster

####1. Functions####

#samples from a normal distribution, but only accepts positive values
rposnorm <- function(n, mu, sigma) {
  #samples from a defined normal distribution until it finds a positive value
  accept = FALSE
  while(!accept) {
    value <- rnorm(n, mean = mu, sd = sigma)
    if(sum(value > 0) == length(value)) accept <- TRUE
  }
  value
}

#creates a correlation matrix using the sampled marginal correlation values
rho_matrix_maker <- function(lkjcorr, K) {
  Rho_matrix <- matrix(NA, ncol = K, nrow = K)
  Rho_matrix[upper.tri(Rho_matrix)] <- lkjcorr
  diag(Rho_matrix) <- 1
  for(i in 1:K) {
    for(j in 1:K) {
      if(is.na(Rho_matrix[i,j])) {
        Rho_matrix[i,j] <- Rho_matrix[j,i]
      }
    }
  }
  Rho_matrix
}

#
agemodel_prior_predictive_check <- function(grand_prior = agemodel_grandprior_list, observation_data = agemodel_observation_data) {
  #Observation data#
  N_Observations <- observation_data$N_Observations
  N_Sites <- observation_data$N_Sites
  N_Periods <- observation_data$N_Periods
  Period <- observation_data$Period
  Site <- observation_data$Site
  N_Total_Mandibles <- observation_data$N_Total_Mandibles
  
  #Hyper-parameters#
  mu_alpha <- rnorm(1, grand_prior$prior_mu_alpha[1], grand_prior$prior_mu_alpha[2])
  logmu_sigma <- rnorm(1, grand_prior$prior_logmu_sigma[1], grand_prior$prior_logmu_sigma[2])
#  sigma_period <- rposnorm(2, 0, 1)
  sigma_period <- c(rposnorm(1, 0, 0.5), rposnorm(1, 0, 0.1))
  #number of random correlations: ((Number of parameters)^2 - (Number of parameters)) / 2
  Rho_period <- ggdist::rlkjcorr_marginal((length(sigma_period)^2 - length(sigma_period)) / 2, K = 2, eta = 2)
  z_period <- matrix(rnorm(2 * N_Periods, 0, 1), ncol = N_Periods, nrow = 2)
  z_site <- rnorm(N_Sites, 0, 1)
  
  #calculating output variables (untransformed)
  v_period <- t(diag(sigma_period) %*% rho_matrix_maker(Rho_period, 2) %*% diag(sigma_period) %*% z_period)
  #
  period_mu <- mu_alpha + v_period[, 1]
  period_sigma <- exp(logmu_sigma + v_period[, 2])
  
  #calculating site-level thetas
  theta_site <- rep(NA, N_Sites)
  for(i in 1:N_Sites) {
    theta_site[i] <- boot::inv.logit(period_mu[Period[i]] + period_sigma[Period[i]] * z_site[i])
  }
  
  #the actual observation
  N_GHI <- rep(NA, N_Observations)
  for(i in 1:N_Observations) {
    N_GHI[i] <- rbinom(n = 1, size = N_Total_Mandibles[i], prob = theta_site[Site[i]])
  }
  
  #output
  output <- list(Overall_Parameters = data.table(mu_alpha, mu_sigma = exp(logmu_sigma)),
                 Period_Parameters = data.table(period_mu, period_sigma, theta_period = boot::inv.logit(period_mu)),
                 Site_Parameters = data.table(theta_site, N_GHI, N_Total_Mandibles))
  output
}

####3. Running the prior ####

#parameters about the run (N_Observations, N_Sites, N_Regions, Region, Site, N_Total_Mandibles)
agemodel_grandprior_list <- list(
  prior_mu_alpha = c(0, 1),
  prior_logmu_sigma = c(-1, 0.5)
)
agemodel_observation_data <- list(
  N_Observations = 50,
  N_Sites = 50,
  N_Periods = 5,
  Period = rep(c(1:5), each = 10),
  Site = 1:50,
  N_Total_Mandibles = rpois(50, 30)
)

agemodel_prior_predictive_check()

clusterExport(cl, list('data.table', 'rbindlist', 'rposnorm', 'rho_matrix_maker', 'rlkjcorr_marginal', 'inv.logit', 'agemodel_prior_predictive_check', 'agemodel_grandprior_list', 'agemodel_observation_data'))

agemodel_predicted_assemblages <- parLapply(cl = cl, 1:1000, fun = function(x) agemodel_prior_predictive_check())
beepr::beep()

prior_predictive_perioddata <- rbindlist(lapply(1:length(agemodel_predicted_assemblages), function(x) data.table(Iteration = x, agemodel_predicted_assemblages[[x]]$Period_Parameters)))
prior_predictive_sitedata <- rbindlist(lapply(1:length(agemodel_predicted_assemblages), function(x) data.table(Iteration = x, Period = agemodel_observation_data$Period, agemodel_predicted_assemblages[[x]]$Site_Parameters)))

####4. Visualizations####
#
hist(prior_predictive_perioddata[, theta_period], main = "Histogram of Simulated Theta Values (Period-level)", xlab = "Theta (Period)")
#
hist(prior_predictive_sitedata[, theta_site], main = "Histogram of Simulated Theta Values (Site-level)", xlab = "Theta (Site)")

#
hist(prior_predictive_perioddata[, max(theta_period) - min(theta_period), Iteration][, V1], main = "Histogram of Ranges of Period-Level Theta Values", xlab = "Range of Theta Values (Period-Level)")

#
hist(rbindlist(lapply(agemodel_predicted_assemblages, function(x) x$Overall_Parameters))[, mu_sigma], main = "Histogram of Simulated Average Variability of Regional Theta Values", xlab = "mu_sigma (Variability of Regional Theta Values)")

#Examinations of the spread of theta values within a region#
prior_predictive_sitedata[, .(Site_Theta_sd = sd(theta_site), Site_Theta_range = max(theta_site) - min(theta_site)), .(Iteration, Period)][, quantile(Site_Theta_range, c(0.025, 0.50, 0.975))]
prior_predictive_sitedata[, .(Site_Theta_sd = sd(theta_site), Site_Theta_range = max(theta_site) - min(theta_site)), .(Iteration, Period)][, mean(Site_Theta_range)]
hist(prior_predictive_sitedata[, .(Site_Theta_sd = sd(theta_site), Site_Theta_range = max(theta_site) - min(theta_site)), .(Iteration, Period)][, Site_Theta_range], main = "Range of Simulated Site Theta Values per Region", xlab = "Range of Theta Values")
hist(prior_predictive_sitedata[, .(Site_Theta_sd = sd(theta_site), Site_Theta_range = max(theta_site) - min(theta_site)), .(Iteration, Period)][, Site_Theta_sd], main = "Standard Deviation of Simulated Site Theta Values per Region", xlab = "Standard Deviation of Theta Values")
