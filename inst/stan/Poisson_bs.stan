data {
  int<lower=0> J;  // number of gene
  int<lower=0> S;  // number of samples
  int<lower=0> N;  // J * S
  int<lower=0> jj[N]; // index of j
  int<lower=0> ss[N]; // index of samples
  int<lower=0> y[N];
  real offset[J];
}
parameters {
  real<lower=0, upper= 100000> bs[S];
}


model {
  for (n in 1:N) {
    target +=   poisson_lpmf(y[n] | bs[ss[n]] * offset[jj[n]]);
  }
}
