library(rstan)
library(parallel)

args <- commandArgs(TRUE)
if (length(args)==0){
  s <- 1
} else {
  s <- as.numeric(args[1])
}
print(paste("Core number",s, "started at", date()))


rstan_schools <- function(y, sigma){
  schools_code <- '
    data {
    int<lower=0> J; // number of schools 
    real y[J]; // estimated treatment effects
    real<lower=0> sigma[J]; // s.e. of effect estimates 
  }
  parameters {
    real mu; 
    real<lower=0> tau;
    real eta[J];
  }
  transformed parameters {
    real theta[J];
    for (j in 1:J)
      theta[j] <- mu + tau * eta[j];
  }
  model {
    eta ~ normal(0, 1);
    y ~ normal(theta, sigma);
  }'
  
  
  schools_dat <- list(J = length(y), 
                      y = y,
                      sigma = sigma)

  fit <- stan(model_code = schools_code, data = schools_dat, iter = 1000, chains = 4)
  return(fit)
}

schools.sim <- function(i,mu0,A,k,sigma){
  out = NULL
  try({
  ##sim true theta
  theta <- rnorm(k, mu0, sqrt(A))
  y <- rnorm(k, theta, sigma)

  ##fit model
  fit <- rstan_schools(y,sigma)
  fit.means <- apply(fit, c(3), mean)
  fit.sds <- apply(fit, c(3), sd)
  fit.25s <- apply(fit, c(3), function(x){quantile(x, 0.025)})
  fit.975s <- apply(fit, c(3), function(x){quantile(x, 0.975)})
  mu.mean <- fit.means[grep("mu", names(fit.means))]
  theta.mean <- fit.means[grep("theta", names(fit.means))]
  theta.25 <- fit.25s[grep("theta", names(fit.means))]
  theta.975 <- fit.975s[grep("theta", names(fit.means))]
  tau.mean <- fit.means[grep("tau", names(fit.means))]

  postmean <- mu0 * (sigma^2 / (sigma^2 + A)) + y * (A / (sigma^2 + A))
  postsd <- sqrt(sigma^2 * (A / (sigma^2 + A)))

  ##calc coverage
  coverage.ind <- as.numeric((theta.25 < theta) & (theta.975 > theta))
  coverage.RB <- pnorm(theta.975, postmean, postsd) - pnorm(theta.25, postmean, postsd)
out <- list(coverage.ind = coverage.ind, coverage.RB = coverage.RB)
})
  return(out)
}

##true values
mu0 <- 8.168
A <- 117.71
k <- 8
sigma <- c(9, 10, 10, 11, 11, 15, 16, 18)
##nits <- 1000

## out <- mclapply(1:nits, schools.sim, mu0 = mu0, A = A, k = k, sigma = sigma, mc.cores = 3)
out <- schools.sim(s, mu0 = mu0, A = A, k = k, sigma = sigma)
save(out, file=paste("output/mcmcout",s,".RData", sep = ""))
