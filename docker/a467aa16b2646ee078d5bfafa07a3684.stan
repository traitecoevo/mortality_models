
      data {
        int<lower=1> n_obs;
        int<lower=1> n_spp;
        int<lower=1> spp[n_obs];
        int<lower=0, upper=1> y[n_obs];
        vector[n_obs] census_length;
        vector[n_obs] growth_dt;
        vector[n_spp] rho_c;
        
        // Held out data
        int<lower=1> n_obs_heldout;
        int<lower=1> n_spp_heldout;
        int<lower=1> spp_heldout[n_obs_heldout];
        int<lower=0, upper=1> y_heldout[n_obs_heldout];
        vector[n_obs_heldout] census_length_heldout;
        vector[n_obs_heldout] growth_dt_heldout;
        vector[n_spp_heldout] rho_c_heldout;
      }
      
      parameters { 
        
      // Mortality model parameters
      real raw_log_a0[n_spp];
      real mu_log_a0;
      real<lower=0> sigma_log_a0;

      real raw_log_b0[n_spp];
      real mu_log_b0;
      real<lower=0> sigma_log_b0;

      real raw_log_c0[n_spp];
      real mu_log_c0;
      real<lower=0> sigma_log_c0;

      
      
      
      }
  
      model {
        
    // Declaring mortality parameters
    real alpha;
    real a0[n_spp];

    real beta;
    real b0[n_spp];

    real gamma;
    real c0[n_spp];
  
    real cumulative_hazard;

    // Calculating species random effects
    for (s in 1:n_spp) {
      a0[s] <- exp(raw_log_a0[s] * sigma_log_a0 + mu_log_a0); // e.g. implies lognormal(mu_log_a0, sigma_log_a0)
      b0[s] <- exp(raw_log_b0[s] * sigma_log_b0 + mu_log_b0);
      c0[s] <- exp(raw_log_c0[s] * sigma_log_c0 + mu_log_c0);
    }

    for (i in 1:n_obs) {
      // Calculating mortality parameters
      alpha <- a0[spp[i]]; 
      beta <- b0[spp[i]]; 
      gamma <- c0[spp[i]];

    // Likelihood for hazard model
    cumulative_hazard <- -census_length[i] * (alpha * exp(-beta * growth_dt[i]) + gamma);

      if (y[i] == 0) {
        increment_log_prob(cumulative_hazard);
      } else {
        increment_log_prob(log1m_exp(cumulative_hazard));
      }
    }

    // Priors

    //Mortality model priors
    raw_log_a0 ~ normal(0,1);
    mu_log_a0 ~ normal(-0.61, 0.57);
    sigma_log_a0 ~ cauchy(0, 2.5);

    raw_log_b0 ~ normal(0, 1);
    mu_log_b0 ~ normal(0, 5);
    sigma_log_b0 ~ cauchy(0, 2.5);

    raw_log_c0 ~ normal(0, 1);
    mu_log_c0 ~ normal(-3.86, 0.2);
    sigma_log_c0 ~ cauchy(0, 2.5);

    
    
    
      }
      
      generated quantities {
        
    // Declaring fitted parameters
    real a0[n_spp];
    real b0[n_spp];
    real c0[n_spp];

    real alpha_fit;
    real beta_fit;
    real gamma_fit;

    real cumulative_hazard_fit;
    real log_lik_fit;
    real sum_log_lik_fit;

    // Declaring heldout parameters
    real alpha_heldout;
    real beta_heldout;
    real gamma_heldout;

    real cumulative_hazard_heldout;
    real log_lik_heldout[n_obs_heldout];
    real sum_log_lik_heldout;

    // Initialization of summed log likelihoods
    sum_log_lik_fit <- 0;
    sum_log_lik_heldout <- 0;

    // recalulate species random effects
    for (s in 1:n_spp) {
      a0[s] <- exp(raw_log_a0[s] * sigma_log_a0 + mu_log_a0);
      b0[s] <- exp(raw_log_b0[s] * sigma_log_b0 + mu_log_b0);
      c0[s] <- exp(raw_log_c0[s] * sigma_log_c0 + mu_log_c0);
    }


    // log likelihood for fitted model
    for (i in 1:n_obs) {
      alpha_fit <- a0[spp[i]]; 
      beta_fit <- b0[spp[i]]; 
      gamma_fit <- c0[spp[i]];

      cumulative_hazard_fit <- -census_length[i] * (alpha_fit * exp(-beta_fit * growth_dt[i]) + gamma_fit);

      if (y[i] == 0) {
        log_lik_fit <- cumulative_hazard_fit;
      }
      else {
        log_lik_fit <- log1m_exp(cumulative_hazard_fit);
      }
      sum_log_lik_fit <- sum_log_lik_fit + log_lik_fit;
    }

    // log likelihood for held out data
    for (j in 1:n_obs_heldout) {
      alpha_heldout <- a0[spp_heldout[j]]; 
      beta_heldout <- b0[spp_heldout[j]]; 
      gamma_heldout <- c0[spp_heldout[j]];

      cumulative_hazard_heldout <- -census_length_heldout[j] * (alpha_heldout * exp(-beta_heldout * growth_dt_heldout[j]) + gamma_heldout);

      if (y_heldout[j] == 0) {
        log_lik_heldout[j] <- cumulative_hazard_heldout;
      }
      else {
        log_lik_heldout[j] <- log1m_exp(cumulative_hazard_heldout);
      }
      sum_log_lik_heldout <- sum_log_lik_heldout + log_lik_heldout[j];
    }
      }
