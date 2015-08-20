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

# bernoullie-gamma truncated [0, 1]
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
    s ~ cauchy(0, 2.5); // half cauchy

    for (n in 1:Ndata) {
    (y[n] == 0) ~ bernoulli_logit(c + d * x[n]);
    if (y[n] > 0)
       y[n] ~ gamma(s, s / inv_logit(a + b*x[n]))  T[,1.0]; // shape, rate
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
       y[n] ~ gamma(s, s / exp(a + b*x[n]))  T[,1.0]; // shape, rate
    }    
}
"
