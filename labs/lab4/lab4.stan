// write a Poisson model for our data

functions {
  /*
  * Alternative to poisson_log_rng() that 
  * avoids potential numerical problems during warmup
  */
  int poisson_log_safe_rng(real eta) {
    real pois_rate = exp(eta);
    if (pois_rate >= exp(20.79))
      return -9;
    return poisson_rng(pois_rate);
  }
}
data {
  //input data
  
  // cosa significa lower?
  int<lower=1>N;
  int<lower=0> complaints[N];
  vector<lower=0>[N] traps;
}
parameters {
  //parameters of the model
  real alpha;
  real beta;
}
model {
  //model
  // se non scriviamo niente di default prende una uniforme
  
} 

generated quantities {
  // sample predicted values from the model for posterior predictive checks
  int y_rep[N];

  for (n in 1:N) {
    real eta_n = alpha + beta * traps[n];
    y_rep[n] = poisson_log_safe_rng(eta_n);
  }
}
