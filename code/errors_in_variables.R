# copied from http://mbjoseph.github.io/blog/2013/11/27/measure/

library(rstan)
library(modeest)

n.reps <- 1
n.repeated <- 10
n <- 50

# true covariate values
x <- runif(n, -3, 3)
y <- x + rnorm(n)  # alpha=0, beta=1, sdy=1

# random subset to perform repeat covariate measurements
which.repeated <- sample(n, n.repeated)
xsd <- 1  # measurement error
xerr <- rnorm(n + (n.repeated * (n.reps - 1)), 0, xsd)

# indx assigns measurements to sample units
indx <- c(1:n, rep(which.repeated, each = n.reps - 1))
indx <- sort(indx)
nobs <- length(indx)
xobs <- x[indx] + xerr
plot(x[indx], xobs,
    xlab = "True covariate value",
    ylab = "Observed covariate value")
abline(0, 1, lty = 2)
segments(x0 = x[indx], x1 = x[indx],
    y0 = x[indx], y1 = xobs, col = "red")
abline(v = x[which.repeated], col = "green", lty = 3)

# write the .stan file
model_string = "
data{
  int n;
  int nobs;
  real xobs[nobs];
  real y[n];
  int indx[nobs];
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigmay;
  real<lower=0> sigmax;
  real x[n];
}

model {
  // priors
  alpha ~ normal(0,100);
  beta ~ normal(0,100);
  sigmay ~ uniform(0,1000);
  sigmax ~ uniform(0,1000);
  
  // model structure  
  for (i in 1:nobs){
    xobs[i] ~ normal(x[indx[i]], sigmax);
  }
  for (i in 1:n){
    y[i] ~ normal(alpha + beta*x[i], sigmay);
  }
}
"

stan_d <- c("y", "xobs", "nobs", "n", "indx")
chains <- 3
iter <- 1000
thin <- 1

mod1 <- stan(model_code = model_string, data = stan_d,
    chains = chains, iter = iter,
    thin = thin)

posteriors <- extract(mod1)

# highest density interval helper function (thanks to Joe Mihaljevic)
HDI <- function(values, percent = 0.95) {
    sorted <- sort(values)
    index <- floor(percent * length(sorted))
    nCI <- length(sorted) - index
    width <- rep(0, nCI)
    for (i in 1:nCI) {
        width[i] <- sorted[i + index] - sorted[i]
    }
    HDImin <- sorted[which.min(width)]
    HDImax <- sorted[which.min(width) + index]
    HDIlim <- c(HDImin, HDImax)
    return(HDIlim)
}

# comparing estimated true x values to actual x values
Xd <- array(dim = c(n, 3))
for (i in 1:n) {
    Xd[i, 1:2] <- HDI(posteriors$x[, i])
    Xd[i, 3] <- mlv(posteriors$x[, i], method = "shorth")$M
}

lims <- c(min(Xd), max(Xd))
plot(x, Xd[, 3], xlab = "True covariate value",
    ylab = "Estimated covariate value",
    col = "purple", pch = 19, ylim = lims)
abline(0, 1, lty = 2)
segments(x0 = x, x1 = x, y0 = Xd[, 1], y1 = Xd[, 2], col = "purple")