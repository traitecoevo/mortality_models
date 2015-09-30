
      data {
        
        int<lower=1> n_obs;
        int<lower=0, upper=1> y[n_obs];
        vector[n_obs] census_length;
        
        // Held out data
        int<lower=1> n_obs_heldout;
        int<lower=0, upper=1> y_heldout[n_obs_heldout];
        vector[n_obs_heldout] census_length_heldout;
      
      }

      parameters {
        
      // Mortality model parameters
      real log_hazard;
      }

      model {
        
      // Declaring mortality parameters
      real hazard;
      real cumulative_hazard;
      
      // Calculate fixed effects on normal scale
      hazard <- exp(log_hazard);
      
      for (i in 1:n_obs) {
        // Likelihood for hazard model
        cumulative_hazard <- -census_length[i] * hazard;
        
        if (y[i] == 0) {
          increment_log_prob(cumulative_hazard);
        } else {
          increment_log_prob(log1m_exp(cumulative_hazard));
        }
      }
      
      // Priors
      
      //Mortality model priors
      log_hazard ~ normal(-3.86, 0.2);
      }

      generated quantities {
        
      // Declaring fitted parameters
      real hazard;
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
      
      // recalulate fixed effects
      hazard <- exp(log_hazard);
      
      // log likelihood for fitted model
      for (i in 1:n_obs) {
        
        cumulative_hazard_fit <- -census_length[i] * hazard;
        
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
        
        cumulative_hazard_heldout <- -census_length_heldout[j] * hazard;
        
        if (y_heldout[j] == 0) {
          log_lik_heldout[j] <- cumulative_hazard_heldout;
        }
        else {
          log_lik_heldout[j] <- log1m_exp(cumulative_hazard_heldout);
        }
        sum_log_lik_heldout <- sum_log_lik_heldout + log_lik_heldout[j];
      }
      }
