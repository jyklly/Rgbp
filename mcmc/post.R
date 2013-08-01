library(Rgbp)

nits <- 1000
results <- lapply(1:nits, function(s){load(file=paste("output/mcmcout",s,".RData", sep = ""));return(out)})
ind <- as.numeric(apply(do.call("rbind",lapply(results,function(x){x$coverage.ind})),2,mean))
rb <- apply(do.call("rbind",lapply(results,function(x){x$coverage.RB})),2,mean)

y <- c(-1, 8, 18, 7, 1, 28, -3, 12)
sigma <- c(9, 10, 10, 11, 11, 15, 16, 18)
k <- length(y)
m1 <- gbp(y,sigma)
pdf("mcmccompare.pdf")
coverage(m1, nsim = 1000)
lines(rb, lty = 1, col = "green", lwd = 2)
lines(ind, lty = 2, col = "darkgreen", lwd = 2)
legend("bottomright", c("MCMC: RB", "MCMC: Simple"), lty = c(1,2), col = c("green", "darkgreen"))
dev.off()

post.approx.rgbp <- as.matrix(cbind(y,sigma,gbp(y, sigma, Alpha=0.5)$post.mean, gbp(y, sigma, Alpha=0.95)$post.sd, gbp(y, sigma, Alpha=0.95)$post.intv.low, gbp(y, sigma, Alpha=0.5)$post.intv.low, gbp(y, sigma, Alpha=0)$post.intv.upp, gbp(y, sigma, Alpha=0.5)$post.intv.upp, gbp(y, sigma, Alpha=0.95)$post.intv.upp))

post.approx.rgbp <- as.data.frame(apply(post.approx.rgbp, 2, as.numeric))
names(post.approx.rgbp) <- c("y", "sigma", "post.mean", "post.sd", "post.2.5%", "post.25%","post.50%", "post.75%", "post.97.5%")
post.approx.rgbp

##mcmc

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

library(rstan)
m2 <- rstan_schools(y, sigma)

