
      data {
        
        int<lower=1> n_obs;
        int<lower=0, upper=1> y[n_obs];
        vector[n_obs] census_length;
        vector[n_obs] growth_dt;
        
        // Held out data
        int<lower=1> n_obs_heldout;
        int<lower=0, upper=1> y_heldout[n_obs_heldout];
        vector[n_obs_heldout] census_length_heldout;
        vector[n_obs_heldout] growth_dt_heldout;
      }

      parameters {
        
      // Mortality model parameters
      real log_alpha;
      real log_beta;
      }

      model {
        
      // Declaring mortality parameters
      real alpha;
      real beta;
      real cumulative_hazard;
      
      alpha <- exp(log_alpha);
      beta <- exp(log_beta);
      
      for (i in 1:n_obs) {
        // Likelihood for hazard model
        cumulative_hazard <- -census_length[i] * (alpha * exp(-beta * growth_dt[i]));
        
        if (y[i] == 0) {
          increment_log_prob(cumulative_hazard);
        } else {
          increment_log_prob(log1m_exp(cumulative_hazard));
        }
      }
      
      // Priors
      
      //Mortality model priors
      log_alpha ~ normal(-0.61, 0.57);
      log_beta ~ normal(0, 5);
      }

      generated quantities {
        
      // Declaring fitted parameters
      real alpha;
      real beta;
      
      real cumulative_hazard_fit;
      real log_lik_fit;
      real sum_log_lik_fit;
      
      real cumulative_hazard_heldout;
      real log_lik_heldout[n_obs_heldout];
      real sum_log_lik_heldout;
      
      // Initialization of summed log likelihoods
      sum_log_lik_fit <- 0;
      sum_log_lik_heldout <- 0;
      
      alpha <- exp(log_alpha);
      beta <- exp(log_beta);
      
      
      // log likelihood for fitted model
      for (i in 1:n_obs) {
        
        cumulative_hazard_fit <- -census_length[i] * (alpha * exp(-beta * growth_dt[i]));
        
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
        
        cumulative_hazard_heldout <- -census_length_heldout[j] * (alpha * exp(-beta * growth_dt_heldout[j]));
        
        if (y_heldout[j] == 0) {
          log_lik_heldout[j] <- cumulative_hazard_heldout;
        }
        else {
          log_lik_heldout[j] <- log1m_exp(cumulative_hazard_heldout);
        }
        sum_log_lik_heldout <- sum_log_lik_heldout + log_lik_heldout[j];
    }
      }
