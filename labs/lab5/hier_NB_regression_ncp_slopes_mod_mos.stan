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
  int<lower=1> K; // numero di predittori
  int<lower=1> J; // numero di palazzi
  
  // Questa è in sostanza una lista di chiavi che puntano al 
  // building di riferimento.
  int<lower=1, upper=J> building_idx[N];
  matrix[J,K] building_data;
  
  // temporal information
  int<lower=1> M; // number of time intervals (12 if month)
  int<lower=1, upper=M> mo_idx[M]; // matching between data and months
  
}

parameters {
  real<lower=0> inv_phi;   
  real beta;               
  
  vector[J] mu_raw;
  vector[J] kappa_raw; // anche qui parametrizzazione non centrata
  real<lower=0> sigma_mu; // varianza delle medie
  real<lower=0> sigma_kappa;
  real alpha;             
  vector[K] zeta;    
  vector[K] gamma; 
  
  real<lower=0, upper=1> rho_raw; // parametrizzazione non centrata
  vector[M] mo_raw;
  real<lower=0> sigma_mo;
}

transformed parameters {
  real phi = inv(inv_phi);
  vector[J] mu; 
  vector[J] kappa;
  
  // autoregressive process
  real rho = 2.0*rho_raw - 1.0;
  vector[M] mo = sigma_mo * mo_raw;
  mo[1] /= sqrt(1-rho^2);
  for(m in 2:M){
    mo[m] += rho * mo[m-1];
  }
    
  // non centered parametrization
  mu = alpha + building_data * zeta + sigma_mu * mu_raw;
  kappa = beta + building_data * gamma + sigma_kappa * kappa_raw;
}

model {
  // prior per i parametri e likelihood
  
  alpha ~ normal(log(4),1);
  zeta ~ normal(0,1);
  beta ~ normal(-0.25,1);
  gamma ~ normal(0,1);
  inv_phi ~ normal(0,1);
  
  // new parameter for mu exploration
  mu_raw ~ normal(0,1);
  kappa_raw ~ normal(0,1);
  
  // Anche se non può essere negativa, essendoci un valore 
  // minimo se mettiamo la normale la taglia. Altrimenti va 
  // bene una chi^2...
  sigma_mu ~ normal(0,1);
  complaints ~ neg_binomial_2_log(mu[building_idx] + 
               beta * traps + mo[mo_idx] + log_sq_foot, phi);
               
  // new priors
  rho_raw ~ beta(10,5);
  mo_raw ~ normal(0,1);
  sigma_mo ~ normal(0,1);
} 

generated quantities {
  int y_rep[N];
  for (n in 1:N) {
    real eta_n = mu[building_idx[n]] + beta * traps[n] + mo[mo_idx[n]] + log_sq_foot[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
  }
}
