# regression bern-log-normal
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
require(rstan)
require('qmap')

set.seed(123)

n <- 100
prob <- 0.8 # non-zero probability
meanlog <- 0.0
sdlog <- 1.0

y <- rbernlnorm(n, prob, meanlog, sdlog)

dataList = list(
    y = y,
    N = length(y)
)

model_string = "
data {
  int<lower=0> N;
  vector<lower=0>[N] y; // fatality rate
}
parameters {
  real mu; // meanlog
  real<lower=0> sigma; // sdlog
  real<lower=0, upper=1> theta; // non-zero probability
}
model {
    mu ~ cauchy(0, 2.5);
    sigma ~ cauchy(0, 2.5); # half-cauchy
    theta ~ cauchy(0, 2.5); # half-cauchy

    for (n in 1:N) {
    (y[n] == 0) ~ bernoulli(1-theta);
    if (y[n] > 0)
       y[n] ~ lognormal(mu, sigma);
    }    
}
"

stanDso <- stan_model( model_code = model_string)

nChains = 4
iterSteps = 10000
burnInSteps =1000

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso ,
                         data = dataList ,
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps ,
                         warmup = burnInSteps ,
                         #init = initsList , # optional
                         thin = 1 )

print(stanFit, digits_summary=3, pars=c('mu','sigma','theta'), probs=c(0.025, 0.5, 0.975))

traceplot(stanFit, pars=c('mu','sigma','theta'))

# increase npts to 500
y <- rbernlnorm(500, prob, meanlog, sdlog)

dataList = list(
    y = y,
    N = length(y)
)

# Get MC sample of posterior:
stanFit <- sampling( object=stanDso ,
                         data = dataList ,
                         #pars = parameters , # optional
                         chains = nChains ,
                         iter = iterSteps ,
                         warmup = burnInSteps ,
                         #init = initsList , # optional
                         thin = 1 )

print(stanFit, digits_summary=3, pars=c('mu','sigma','theta'), probs=c(0.025, 0.5, 0.975))

traceplot(stanFit, pars=c('mu','sigma','theta'))
