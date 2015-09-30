
      data {
        
        int<lower=1> n_obs;
        int<lower=1> n_spp;
        int<lower=1> spp[n_obs];
        int<lower=0, upper=1> y[n_obs];
        vector[n_obs] census_length;
        vector[n_obs] growth_dt;
        
        // Held out data
        int<lower=1> n_obs_heldout;
        int<lower=1> spp_heldout[n_obs_heldout];
        int<lower=0, upper=1> y_heldout[n_obs_heldout];
        vector[n_obs_heldout] census_length_heldout;
        vector[n_obs_heldout] growth_dt_heldout;
      }

      parameters {
        
      // Mortality model parameters
      real raw_log_alpha[n_spp];
      real mu_log_alpha;
      real<lower=0> sigma_log_alpha;
      
      real raw_log_beta[n_spp];
      real mu_log_beta;
      real<lower=0> sigma_log_beta;
      
      real raw_log_gamma[n_spp];
      real mu_log_gamma;
      real<lower=0> sigma_log_gamma;
      }

      model {
        
      // Declaring mortality parameters
      real alpha[n_spp];
      real beta[n_spp];
      real gamma[n_spp];
      real cumulative_hazard;
      
      // Calculating species random effects
      for (s in 1:n_spp) {
        alpha[s] <- exp(raw_log_alpha[s] * sigma_log_alpha + mu_log_alpha); // e.g. implies lognormal(mu_log_alpha, sigma_log_alpha)
        beta[s] <- exp(raw_log_beta[s] * sigma_log_beta + mu_log_beta);
        gamma[s] <- exp(raw_log_gamma[s] * sigma_log_gamma + mu_log_gamma);
      }
      
      for (i in 1:n_obs) {
        // Likelihood for hazard model
        cumulative_hazard <- -census_length[i] * (alpha[spp[i]] * exp(-beta[spp[i]] * growth_dt[i]) + gamma[spp[i]]);
        
        if (y[i] == 0) {
          increment_log_prob(cumulative_hazard);
        } else {
          increment_log_prob(log1m_exp(cumulative_hazard));
        }
      }
      
      // Priors
      
      //Mortality model priors
      raw_log_alpha ~ normal(0,1);
      mu_log_alpha ~ normal(-0.61, 0.57);
      sigma_log_alpha ~ cauchy(0, 2.5);
      
      raw_log_beta ~ normal(0, 1);
      mu_log_beta ~ normal(0, 5);
      sigma_log_beta ~ cauchy(0, 2.5);
      
      raw_log_gamma ~ normal(0, 1);
      mu_log_gamma ~ normal(-3.86, 0.2);
      sigma_log_gamma ~ cauchy(0, 2.5);
      }

      generated quantities {
        
      // Declaring fitted parameters
      real alpha[n_spp];
      real beta[n_spp];
      real gamma[n_spp];
      
      real cumulative_hazard_fit;
      real log_lik_fit;
      real sum_log_lik_fit;
      
      real cumulative_hazard_heldout;
      real log_lik_heldout[n_obs_heldout];
      real sum_log_lik_heldout;
      
      // Initialization of summed log likelihoods
      sum_log_lik_fit <- 0;
      sum_log_lik_heldout <- 0;
      
      // recalulate species random effects
      for (s in 1:n_spp) {
        alpha[s] <- exp(raw_log_alpha[s] * sigma_log_alpha + mu_log_alpha);
        beta[s] <- exp(raw_log_beta[s] * sigma_log_beta + mu_log_beta);
        gamma[s] <- exp(raw_log_gamma[s] * sigma_log_gamma + mu_log_gamma);
      }
      
      
      // log likelihood for fitted model
      for (i in 1:n_obs) {
        
        cumulative_hazard_fit <- -census_length[i] * (alpha[spp[i]] * exp(-beta[spp[i]] * growth_dt[i]) + gamma[spp[i]]);
        
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
        
        cumulative_hazard_heldout <- -census_length_heldout[j] * (alpha[spp_heldout[j]] * exp(-beta[spp_heldout[j]] * growth_dt_heldout[j]) + gamma[spp_heldout[j]]);
        
        if (y_heldout[j] == 0) {
          log_lik_heldout[j] <- cumulative_hazard_heldout;
        }
        else {
          log_lik_heldout[j] <- log1m_exp(cumulative_hazard_heldout);
        }
        sum_log_lik_heldout <- sum_log_lik_heldout + log_lik_heldout[j];
    }
      }
