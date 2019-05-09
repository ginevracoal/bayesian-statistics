// Poisson model for our data

functions {
  /*
  * Alternative to poisson_log_rng() that 
  * avoids potential numerical problems during warmup
  */
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi/exp(eta));
    if (gamma_rate >= exp(20.79))
      return -9;
    // perché tende ad una poisson
    return poisson_rng(gamma_rate);
  }
}

data {
  //input data
  
  // lower è il valore minimo della variabile
  int<lower=1>N;
  int<lower=0> complaints[N];
  vector<lower=0>[N] traps;
  vector<lower=0>[N] live_in_super;
  vector<lower=0>[N] log_sq_foot;
}

parameters {
  //parameters of the model
  real alpha;
  real beta;
  real beta_super;
  real inv_phi;
}

transformed parameters{
  real phi = inv(inv_phi);
}

model {
  //model
  // se non scriviamo niente di default prende una uniforme
  beta ~ normal(-0.25,1);
  alpha ~ normal(log(4),1);
  beta_super ~ normal(-0.5,1);
  inv_phi ~ normal(0,1);
  complaints ~ neg_binomial_2(alpha + beta * traps + beta_super * live_in_super + log_sq_foot, phi);
} 

generated quantities {
  // sample predicted values from the model for posterior 
  // predictive checks
  int y_rep[N];

  for (n in 1:N) {
    real eta_n = alpha + beta * traps[n] + beta_super * live_in_super[n] + log_sq_foot[n];
    // si chiama safe perché se il valore è troppo 
    // grande lo riportiamo ad un valore sicuro
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
  }
}
