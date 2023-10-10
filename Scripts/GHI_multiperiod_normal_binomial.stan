//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower = 1> N_Observations;
  int<lower = 1> N_Sites;
  int<lower = 1> N_Periods;
  int<lower = 1, upper = N_Sites> Site[N_Observations];
  int<lower = 1, upper = N_Periods> Period[N_Sites];
  int<lower = 0> N_GHI[N_Observations];
  int<lower = 0> N_Total_Mandibles[N_Observations];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // parametric set up is thus:
  // logit(theta[Site]) ~ Normal(mu[Period], sigma[Period])
  // N_GHI[Observation] ~ Binomial(N_Total_Mandibles[Observation], theta[Site])
  // parameters that describe the overall 
  real mu_alpha; // overall mean (mu_mu)
  real logmu_sigma; // overall sigma (logmu_sigma)
  // variation across periods
  vector<lower = 0>[2] sigma_period; // variation in mu_alpha and logsigma_alpha
  matrix[2, N_Periods] z_period;
  cholesky_factor_corr[2] L_Rho_period;
  // site-specific thetas
  vector[N_Sites] z_site;
}

transformed parameters {
  // period-level parameters
  matrix[N_Periods, 2] v_period;
  matrix[2, 2] Rho_period;
  vector[N_Periods] period_mu;
  vector[N_Periods] period_sigma;
  // site-specific thetas
  vector[N_Sites] theta_site;
  //calculating the offsets
  v_period = (diag_pre_multiply(sigma_period, L_Rho_period) * z_period)';
  Rho_period = L_Rho_period * L_Rho_period';
  period_mu = mu_alpha + col(v_period, 1);
  period_sigma = exp(logmu_sigma + col(v_period, 2));
  for(i in 1:N_Sites) {
    theta_site[i] = inv_logit(period_mu[Period[i]] + period_sigma[Period[i]] * z_site[i]);
  }
}

model {
  //prior distributions
  mu_alpha ~ normal(0, 1); //reflects prior knowledge
  logmu_sigma ~ normal(-1, 0.5); // trying to approximate a half-normal(0, 1)
  // parameters for inter-period variation
  sigma_period[1] ~ normal(0, 0.5);
  sigma_period[2] ~ normal(0, 0.1);
  L_Rho_period ~ lkj_corr_cholesky(2);
  to_vector(z_period) ~ normal(0, 1);
  // site-level variation
  z_site ~ normal(0, 1);
  for(i in 1:N_Observations) {
    N_GHI[i] ~ binomial(N_Total_Mandibles[i], theta_site[Site[i]]);
  }
}

generated quantities {
  real theta_average;
  vector[N_Periods] theta_period;
  vector[N_Periods] theta_period_random;
  theta_average = inv_logit(mu_alpha);
  theta_period = inv_logit(period_mu);
  for(i in 1:N_Periods) {
    theta_period_random[i] = inv_logit(normal_rng(period_mu[i], period_sigma[i]));
  }
}
