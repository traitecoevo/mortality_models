
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
      real raw_log_hazard[n_spp];
      real mu_log_hazard;
      real<lower=0> sigma_log_hazard;
      }

      model {
        
      // Declaring mortality parameters
      real hazard[n_spp];
      real cumulative_hazard;
      
      // Calculating species random effects
      for (s in 1:n_spp) {
        hazard[s] <- exp(raw_log_hazard[s] * sigma_log_hazard + mu_log_hazard); // e.g. implies lognormal(mu_log_hazard, sigma_log_hazard)
      }
      
      for (i in 1:n_obs) {
        // Likelihood for hazard model
        cumulative_hazard <- -census_length[i] * hazard[spp[i]];
        
        if (y[i] == 0) {
          increment_log_prob(cumulative_hazard);
        } else {
          increment_log_prob(log1m_exp(cumulative_hazard));
        }
      }
      
      // Priors
      
      //Mortality model priors
      raw_log_hazard ~ normal(0, 1);
      mu_log_hazard ~ normal(-3.86, 0.2);
      sigma_log_hazard ~ cauchy(0, 2.5);
      }

      generated quantities {
        
      // Declaring fitted parameters
      real hazard[n_spp];
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
        hazard[s] <- exp(raw_log_hazard[s] * sigma_log_hazard + mu_log_hazard);
      }
      
      
      // log likelihood for fitted model
      for (i in 1:n_obs) {
        
        cumulative_hazard_fit <- -census_length[i] * hazard[spp[i]];
        
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
        
        cumulative_hazard_heldout <- -census_length_heldout[j] * hazard[spp_heldout[j]];
        
        if (y_heldout[j] == 0) {
          log_lik_heldout[j] <- cumulative_hazard_heldout;
        }
        else {
          log_lik_heldout[j] <- log1m_exp(cumulative_hazard_heldout);
        }
        sum_log_lik_heldout <- sum_log_lik_heldout + log_lik_heldout[j];
    }
      }
