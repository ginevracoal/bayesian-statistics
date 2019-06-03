data{
  int N;          // number of voters
  int vote[N];    // vote: 0 (Clinton), 1 (Bush)
  int income[N];  // 1-5 income scale
}
parameters{
  real alpha;    // intercept
  real beta;     // income coefficient
}
model{
  for (n in 1:N){
    target += bernoulli_lpmf(vote[n] | Phi(alpha+income[n]*beta));
                // likelihood          
  }
  
  target += normal_lpdf(alpha | 0, 100); // intercept prior
  target += normal_lpdf(beta | 0, 100); // income prior
}

