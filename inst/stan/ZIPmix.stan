data {
  int<lower=0> I;  // number of patients
  int<lower=0> J;  // number of gene
  int N;
  int<lower=1,upper=I> ii[N]; // patient for observation n
  int<lower=1,upper=J> jj[N]; // gene for observation n
  int<lower=0> y[N];
  real<lower=0> a_zip;
  real offset[N];
  real<lower=0> sigma;
}
parameters {
  real<lower=0, upper=1> zero_p;
  real<lower=0> a_opt;
}


model {
  for (n in 1:N) {
    if (y[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | zero_p),
                            bernoulli_lpmf(0 | zero_p)
                              + poisson_lpmf(y[n] | a_opt * offset[n]))
                              + normal_lpdf(log(a_opt) | log(a_zip), sigma);
    else
      target += bernoulli_lpmf(0 | zero_p)
      + poisson_lpmf(y[n] |  a_opt* offset[n])
      + normal_lpdf(log(a_opt) | log(a_zip), sigma);
  }
}
