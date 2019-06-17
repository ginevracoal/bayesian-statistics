data {
  int<lower=1> N; // number of observations
  int<lower=1> D; // number of weekdays
  //int<lower=1> K; // number of predictors
  real<lower=0> tot_mort[N]; // daily deaths
  // predictors
  vector[N] SO2; 
  vector[N] TSP;
  vector[N] mean_temp;
  int<lower=1, upper=D> day_of_week[N];
  // weekday data with D variables
  //int<lower=1, upper=D> weekday_idx[N];
  //matrix[D,K] weekday_data;
}
parameters {
  vector[D] alpha; 
  real beta;
  real gamma;
  real delta;
  vector<lower=0>[D] sigma; // sd error
}
transformed parameters{}
model {
  // Priors
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  gamma ~ normal(0,10);
  delta ~ normal(0,10);
  sigma ~ cauchy(0,10);
  
  // Likelihood
  tot_mort ~ normal(alpha[day_of_week] + 
                    beta * SO2 + gamma * TSP +
                    delta * mean_temp, 
                    sigma[day_of_week]);
  //for(n in 1:N){
  //  target +=  normal_lpdf(tot_mort[n] | mu, sigma);
} 
generated quantities {
  // sample predicted values from the model for posterior 
  // predictive checks
  real y_rep[N];
  vector[N] log_lik;
  
  for (n in 1:N) {
    real mu = alpha[day_of_week[n]] + beta * SO2[n] + gamma * TSP[n] + delta * mean_temp[n];
    y_rep[n] = normal_rng(mu, sigma[day_of_week[n]]);
    log_lik[n] = normal_lpdf(tot_mort[n]| mu, sigma[day_of_week[n]]);
  }
}
