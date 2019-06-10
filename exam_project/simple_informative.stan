data {
  int<lower=1> N; // number of days
  int<lower=0> tot_mort[N]; // daily deaths
  vector[N] SO2; // predictor
}

parameters {
  real alpha; // intercept
  real beta; // slope
  real<lower = 0> sigma; // sd error
}

transformed parameters{}

model {
  // Priors
  
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  sigma ~ cauchy(0,30);
  
  // Likelihood
  for (n in 1:N){ 
    tot_mort[n] ~ normal(alpha + beta * SO2[n], sigma);
  }
  //for(n in 1:N){
  //  target +=  normal_lpdf(y[n] | mu, sigma);
} 

generated quantities {
  // sample predicted values from the model for posterior 
  // predictive checks
  real y_rep[N];
  real mu;
  
  for (n in 1:N) {
    mu = alpha + beta * SO2[n];
    y_rep[n] = normal_rng(mu, sigma);
    //y_rep[n] =  mu + sigma * inv_Phi(SO[n]);  
  }
}
