data {
  int<lower=0> I;  // number of patients
  int<lower=0> J;  // number of gene
  int N;
  int<lower=1,upper=I> ii[N]; // patient for observation n
  int<lower=1,upper=J> jj[N]; // gene for observation n
  int<lower=0> y[N];
  real offset[N];
}
parameters {
  real<lower=0, upper=1> zero_p[J];
  real<lower=0> lambda[J];
}

// transformed parameters {
//   real mu = lambda[jj[N]] * offset;
// }
model {
  for (n in 1:N) {
    if (y[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | zero_p[jj[n]]),
                            bernoulli_lpmf(0 | zero_p[jj[n]])
                              + poisson_lpmf(y[n] | lambda[jj[n]] * offset[n]));
    else
      target += bernoulli_lpmf(0 | zero_p[jj[n]])
      + poisson_lpmf(y[n] | lambda[jj[n]] * offset[n]);
  }
}
