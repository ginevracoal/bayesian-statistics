data {
  int<lower=1> N; // number of days
  real<lower=0> avg_mort[N]; // avg deaths
  // predictors
  vector[N] SO2; 
  vector[N] TSP;
}

parameters {
  real alpha; 
  real beta;
  real gamma;
  real<lower = 0> sigma; // sd error
}

transformed parameters{}

model {
  // Priors
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  gamma ~ normal(0,10);
  sigma ~ cauchy(0,30);
  
  // Likelihood
  for (n in 1:N){ 
    avg_mort[n] ~ normal(alpha + beta * SO2[n] + gamma * TSP[n],
    sigma);
  }
  //for(n in 1:N){
  //  target +=  normal_lpdf(tot_mort[n] | mu, sigma);
} 

generated quantities {
  // sample predicted values from the model for posterior 
  // predictive checks
  real y_rep[N];
  real mu;
  
  for (n in 1:N) {
    mu = alpha + beta * SO2[n] + gamma * TSP[n];
    y_rep[n] = normal_rng(mu, sigma);
  }
}