
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
      real mu_log_a0;
      real mu_log_b0;
      real mu_log_c0;
      
      
      real b1;
      real c1;
      }

      model {
        
      // Declaring mortality parameters
      real alpha;
      real beta;
      real gamma;
      real cumulative_hazard;
      
      // Calculating species random effects
      
      for (i in 1:n_obs) {
        // Calculating mortality parameters
        alpha <- exp(mu_log_a0); 
        beta <- exp(mu_log_b0) * pow(rho_c[spp[i]], b1); 
        gamma <- exp(mu_log_c0) * pow(rho_c[spp[i]], c1);
        
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
      mu_log_a0 ~ normal(-0.61, 0.57);
      mu_log_b0 ~ normal(0, 5);
      mu_log_c0 ~ normal(-3.86, 0.2);
      
      
      b1 ~ normal(0,5);
      c1 ~ normal(0,5);
      }

      generated quantities {
        
      // Declaring fitted parameters
      real alpha[n_spp];
      real beta[n_spp];
      real gamma[n_spp];
      
      real cumulative_hazard_fit;
      real log_lik_fit;
      real sum_log_lik_fit;
      
      // Declaring heldout parameters
      
      real cumulative_hazard_heldout;
      real log_lik_heldout[n_obs_heldout];
      real sum_log_lik_heldout;
      
      // Initialization of summed log likelihoods
      sum_log_lik_fit <- 0;
      sum_log_lik_heldout <- 0;
      
      // log likelihood for fitted model
      for (s in 1:n_spp) {
        alpha[s] <- exp(mu_log_a0); 
        beta[s] <- exp(mu_log_b0) * pow(rho_c[s], b1); 
        gamma[s] <- exp(mu_log_c0) * pow(rho_c[s], c1);
      }
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
      
      // log likelihood for heldout data
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
