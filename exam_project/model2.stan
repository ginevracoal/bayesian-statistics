data {
  int<lower=1> N; // number of days
  real<lower=0> avg_mort[N]; // daily deaths
  // predictors
  vector[N] SO2; 
  vector[N] TSP;
  vector[N] mean_temp;
}

parameters {
  real alpha; 
  real beta; 
  real gamma;
  real delta;
  real<lower = 0> sigma; // sd error
}

transformed parameters{}

model {
  // Priors
  alpha ~ normal(0,30);
  beta ~ normal(0,30);
  gamma ~ normal(0,30);
  delta ~ normal(0,30);
  sigma ~ cauchy(0,30);
  
  // Likelihood
  for (n in 1:N){
    avg_mort[n] ~ normal(alpha + beta * SO2[n] +
    gamma * TSP[n] + delta * mean_temp[n], sigma );
  }
  //for(n in 1:N){
  //  target +=  normal_lpdf(y[n] | mu, sigma);
} 

generated quantities {
  // sample predicted values from the model for posterior 
  // predictive checks
  real y_rep[N];
  real mu;
  vector[N] log_lik;
  
  for (n in 1:N) {
    mu = alpha + beta * SO2[n] + gamma * TSP[n] + delta * mean_temp[n];
    y_rep[n] = normal_rng(mu, sigma);
    log_lik[n] = normal_lpdf(avg_mort[n]| mu, sigma);
  }
}
