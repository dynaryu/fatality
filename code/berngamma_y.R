# regression bern-gamma for vector
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
require('qmap')

set.seed(123)

n <- 100
prob <- 0.8 # non-zero probability
shape <- 0.6
rate <- 2.0
scale <- 1/rate

y <- rberngamma(n, prob, scale=scale, shape=shape)

dataList = list(
    y = y,
    N = length(y)
)

model_string = "
data {
  int<lower=0> N;
  vector<lower=0>[N] y; // 
}
parameters {
  real<lower=0> shape; // shape
  real<lower=0> rate; // rate
  real<lower=0, upper=1> theta; // non-zero probability
}
model {
    shape ~ cauchy(0, 2.5); # half-cauchy
    rate ~ cauchy(0, 2.5); # half-cauchy
    theta ~ cauchy(0, 2.5); # half-cauchy

    for (n in 1:N) {
    (y[n] == 0) ~ bernoulli(1-theta);
    if (y[n] > 0)
       y[n] ~ gamma(shape, rate); //
    }    
}
"

stanDso <- stan_model( model_code = model_string)

nChains = 4
iterSteps = 5000
burnInSteps = 500

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso ,
                         data = dataList ,
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps ,
                         warmup = burnInSteps ,
                         #init = initsList , # optional
                         thin = 1 )

sampled <- extract(stanFit) # list
paramNames <- names(sampled)
paramNames <- paramNames[-length(paramNames)]

print(stanFit, digits_summary=3, pars=paramNames, probs=c(0.025, 0.5, 0.975))

traceplot(stanFit, pars=paramNames)

# try with truncated gamma
y1 <- y[y<1.0]

dataList = list(
    y = y1,
    N = length(y1)
)

model_string_tr = "
data {
  int<lower=0> N;
  vector<lower=0, upper=1>[N] y; // 
}
parameters {
  real<lower=0> shape; // shape
  real<lower=0> rate; // rate
  real<lower=0, upper=1> theta; // non-zero probability
}
model {
    shape ~ cauchy(0, 2.5); # half-cauchy
    rate ~ cauchy(0, 2.5); # half-cauchy
    theta ~ cauchy(0, 2.5); # half-cauchy

    for (n in 1:N) {
    (y[n] == 0) ~ bernoulli(1-theta);
    if (y[n] > 0)
       y[n] ~ gamma(shape, rate) T[,1]; //
    }    
}
"
stanDso <- stan_model( model_code = model_string_tr)

nChains = 4
iterSteps = 5000
burnInSteps = 500

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso ,
                         data = dataList ,
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps ,
                         warmup = burnInSteps ,
                         #init = initsList , # optional
                         thin = 1 )

print(stanFit, digits_summary=3, pars=paramNames, probs=c(0.025, 0.5, 0.975))

traceplot(stanFit, pars=paramNames)
