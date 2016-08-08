# stan models for Indonesian fatality model
# expecting datalist consiting of Ndat(number), x (mmi), and y(fatality ratio)

# bernoullie-log-normal truncated [0, 1]
# note that lognormal in stan takes meanlog and sdlog
model_bernlnorm = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real c; // parameter of linear predictor
  real d; // c+dx
  real<lower=0> s; // sigma lognormal constant cov
}
model {
  a ~ cauchy(0, 2.5);
  b ~ cauchy(0, 2.5);
  c ~ cauchy(0, 2.5);
  d ~ cauchy(0, 2.5);
  s ~ cauchy(0, 2.5); // half cauchy

  for (n in 1:Ndata) {
  (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
  if (y[n] > 0)
     y[n] ~ lognormal(a + b*x[n], s) T[,1.0];
  }
}
"

# same as bernoullie-log-normal but added generated quantities for prediction
model_bernlnorm_add = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate

  int<lower=0> Nnew ;
  vector[Nnew] xnew ; // xnew
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real c; // parameter of linear predictor
  real d; // c+dx
  real<lower=0> s; // sigma lognormal constant cov
}
model {
  a ~ cauchy(0, 2.5);
  b ~ cauchy(0, 2.5);
  c ~ cauchy(0, 2.5);
  d ~ cauchy(0, 2.5);
  s ~ cauchy(0, 2.5); // half cauchy

  for (n in 1:Ndata) {
  (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
  if (y[n] > 0)
     y[n] ~ lognormal(a + b*x[n], s);
  }
}
generated quantities{
  vector[Nnew] ynew ;
  for (n in 1:Nnew) {
    ynew[n] <- bernoulli_rng(1.0 - inv_logit(c + d*xnew[n]));
    if (ynew[n] > 0)
      ynew[n] <- lognormal_rng(a + b*xnew[n], s);
   }
}
"

# bernoullie-gamma truncated [0, 1]
# took so long
# note that gamma in stan takes shape(alpha) and rate(beta)
# mean = shape/rate, variance = shape/(rate**2)
# mean = link(a+b*x) (either exp or inv_logit)
model_berngamma_logit = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real c; // parameter of linear predictor
  real d; // c+dx
  real<lower=0> s; // sigma lognormal constant cov
}
model {
  a ~ cauchy(0, 2.5);
  b ~ cauchy(0, 2.5);
  c ~ cauchy(0, 2.5);
  d ~ cauchy(0, 2.5);
  s ~ cauchy(0, 2.5); // shape, half cauchy

  for (n in 1:Ndata) {
  (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
  if (y[n] > 0)
    y[n] ~ gamma(s, s / inv_logit(a + b*x[n])); // shape, rate
    /*
    y[n] ~ gamma(s, s / inv_logit(a + b*x[n]))  T[,1.0]; // shape, rate
    */
  }
}
"

# bernoullie-gamma truncated [0, 1]
# note that gamma in stan takes shape(alpha) and rate(beta)
# mean = shape/rate, variance = shape/(rate**2)
# mean = link(a+b*x) (either exp or inv_logit)
model_berngamma_log = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real c; // parameter of linear predictor
  real d; // c+dx
  real<lower=0> s; // shape, half cauchy
}
model {
  a ~ cauchy(0, 2.5);
  b ~ cauchy(0, 2.5);
  c ~ cauchy(0, 2.5);
  d ~ cauchy(0, 2.5);
  s ~ cauchy(0, 2.5); // half cauchy

  for (n in 1:Ndata) {
  (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
  if (y[n] > 0)
    y[n] ~ gamma(s, s / exp(a + b*x[n])); // shape, rate
    /*
    y[n] ~ gamma(s, s / exp(a + b*x[n])) T[,1.0]; // shape, rate
    */
  }
}
"

# same as b-gamma but added generated quantities for prediction
model_berngamma_log_add = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate

  int<lower=0> Nnew ;
  vector[Nnew] xnew ; // xnew
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real c; // parameter of linear predictor
  real d; // c+dx
  real<lower=0> s; // shape, half cauchy
}
model {
  a ~ cauchy(0, 2.5);
  b ~ cauchy(0, 2.5);
  c ~ cauchy(0, 2.5);
  d ~ cauchy(0, 2.5);
  s ~ cauchy(0, 2.5); // half cauchy

  for (n in 1:Ndata) {
  (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
  if (y[n] > 0)
    y[n] ~ gamma(s, s / exp(a + b*x[n])); // shape, rate
    /*
    y[n] ~ gamma(s, s / exp(a + b*x[n])) T[,1.0]; // shape, rate
    */
  }
}
generated quantities{
  vector[Nnew] ynew ;
  for (n in 1:Nnew) {
    ynew[n] <- bernoulli_rng(1.0 - inv_logit(c + d*xnew[n]));
    if (ynew[n] > 0)
      ynew[n] <- gamma_rng(s, s / exp(a + b*xnew[n]));
   }
}
"

# apply binonmial-gamma
model_binomialgamma_log_add = "
data {
  int<lower=0> Ndata ;
  vector<lower=0>[Ndata] x ; // MMI
  vector<lower=0, upper=1>[Ndata] y ; // fatality rate
  vector<lower=0>[Ndata] pop ; // MMI
  vector<lower=0>[Ndata] fat ; // MMI

  int<lower=0> Nnew ;
  vector[Nnew] xnew ; // xnew
}
parameters {
  real a; // parameter of linear predictor
  real b; // a+bx
  real c; // parameter of linear predictor
  real d; // c+dx
  real<lower=0> s; // shape, half cauchy
}
model {
  a ~ cauchy(0, 2.5);
  b ~ cauchy(0, 2.5);
  c ~ cauchy(0, 2.5);
  d ~ cauchy(0, 2.5);
  s ~ cauchy(0, 2.5); // half cauchy

  for (n in 1:Ndata) {
  (fat[n] == 0) ~ binomial_logit(pop[n], c + d * x[n]);
  if (fat[n] > 0)
    y[n] ~ gamma(s, s / exp(a + b*x[n])); // shape, rate
  }
}
generated quantities{
  vector[Nnew] ynew ;
  for (n in 1:Nnew) {
    ynew[n] <- binomial_rng(1, 1.0 - inv_logit(c + d*xnew[n]));
    if (ynew[n] > 0)
      ynew[n] <- gamma_rng(s, s / exp(a + b*xnew[n]));
   }
}
"
