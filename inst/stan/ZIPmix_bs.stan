data {
  int<lower=0> J;  // number of gene
  int<lower=0> S;  // number of samples
  int<lower=0> N;  // J * S
  int<lower=0> jj[N]; // index of j
  int<lower=0> ss[N]; // index of samples
  int<lower=0> y[N];
  real offset[J];
  real<lower=0, upper=1> zero_p[J];
}
parameters {
  real<lower=0, upper= 100000> bs[S];
}


model {
  for (n in 1:N) {
    if (y[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | zero_p[jj[n]]),
                            bernoulli_lpmf(0 | zero_p[jj[n]])
                              + poisson_lpmf(y[n] | bs[ss[n]] * offset[jj[n]]));
      else
        target += bernoulli_lpmf(0 | zero_p[jj[n]])
        + poisson_lpmf(y[n] |  bs[ss[n]]* offset[jj[n]]);
  }
}

