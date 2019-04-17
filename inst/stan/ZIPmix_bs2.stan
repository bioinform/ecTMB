data {
  int<lower=0> J;  // number of gene
  int<lower=0> S;  // number of samples
  int<lower=0> N;  // J * S
  int<lower=0> jj[N]; // index of j
  int<lower=0> ss[N]; // index of samples
  int<lower=0> y[N];
  real offset[J];
  real<lower=0, upper=1> zero_p[J];
  real<lower=0> mu1;
  real<lower=0> mu2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0, upper=1> lambda1;
  real<lower=0, upper=1> lambda2;
}
parameters {
  real<lower=0, upper= 100000> bs[S];
}


model {
  for (n in 1:N) {
    if (y[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | zero_p[jj[n]]),
                            bernoulli_lpmf(0 | zero_p[jj[n]])
                              + poisson_lpmf(y[n] | bs[ss[n]] * offset[jj[n]]))
      + log_sum_exp(log(lambda1) + normal_lpdf(log(bs[ss[n]]) | mu1, sigma1),
                     log(lambda2) + normal_lpdf(log(bs[ss[n]]) | mu2, sigma2));
      else
        target += bernoulli_lpmf(0 | zero_p[jj[n]])
        + poisson_lpmf(y[n] |  bs[ss[n]]* offset[jj[n]])
        + log_sum_exp(log(lambda1) + normal_lpdf(log(bs[ss[n]]) | mu1, sigma1),
                      log(lambda2) + normal_lpdf(log(bs[ss[n]]) | mu2, sigma2));
  }
}

