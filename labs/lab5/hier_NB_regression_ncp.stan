functions {
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    if (gamma_rate >= exp(20.79))
      return -9;
      
    return poisson_rng(gamma_rate);
  }
}

data {
  int<lower=1> N;                     
  int<lower=0> complaints[N];              
  vector<lower=0>[N] traps;                
  
  // 'exposure'
  vector[N] log_sq_foot;  
  
  // building-level data
  int<lower=1> K; # numero di predittori
  int<lower=1> J; # numero di palazzi
  
  // Questa è in sostanza una lista di chiavi che puntano al 
  // building di riferimento.
  int<lower=1, upper=J> building_idx[N];
  matrix[J,K] building_data;
}
parameters {
  real<lower=0> inv_phi;   
  real beta;               
  
  vector[J] mu_raw;
  real<lower=0> sigma_mu; # varianza delle medie
  real alpha;             
  vector[K] zeta;    
}
transformed parameters {
  real phi = inv(inv_phi);
  
  // non centered parametrization
  vector[J] mu = alpha + building_data * zeta 
                  + sigma_mu * mu_raw;
}
model {
  // prior per i parametri e likelihood
  
  
  alpha ~ normal(log(4),1);
  zeta ~ normal(0,1);
  beta ~ normal(-0.25,1);
  inv_phi ~ normal(0,1);
  
  // new parameter for mu exploration
  mu_raw ~ normal(0,1);
  
  // Anche se non può essere negativa, essendoci un valore 
  // minimo se mettiamo la normale la taglia. Altrimenti va 
  // bene una chi^2...
  sigma_mu ~ normal(0,1);
  complaints ~ neg_binomial_2_log(mu[building_idx] + 
               beta * traps + log_sq_foot, phi);
} 

generated quantities {
  int y_rep[N];
  for (n in 1:N) {
    real eta_n = mu[building_idx[n]] + beta * traps[n] + log_sq_foot[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
  }
}
