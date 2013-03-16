\name{Rgbp}
\alias{Rgbp-package}
\alias{Rgbp}
\docType{package}
\title{Bayesian Hierarchical Modeling using Stein's Harmonic Prior}
\description{
Bayesian Hierarchical modeling for Gaussian (GRIMM), Binomial (BRIMM) and Poisson (PRIMM) data \cr using generalized Stein's harmonic prior for good frequentist repeated sampling property.
}
\details{
  \tabular{ll}{
    Package: \tab Rgbp\cr
    Type: \tab Package\cr
    Version: \tab 1.0.0\cr
    Date: \tab 2013-03-16\cr
    License: \tab GPL-2\cr
  }
}
\author{
  Joseph Kelly, Carl Morris, and Hyungsuk Tak

  Maintainer: Joseph Kelly <kelly2@fas.harvard.edu>
}
\references{
  Morris, C. and Tang, R. (2011). Estimating Random Effects via Adjustment for Density Maximization. \emph{Statistical Science}. \bold{26}. 271-287.

  Morris, C. and Lysy, M. (2012). Shrinkage Estimation in Multilevel Normal Models. \emph{Statistical Science}. \bold{27}. 115-134.
}

\examples{

  # Loading datasets
  data(schools)

  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)
  x1 <- rep(c(-1, 0, 1), 6)

  y <- schools$y
  se <- schools$se
  x2 <- rep(c(-1, 0, 1, 2), 2)
  
  ####################################################################################
  # If we do not have any covariate and do not know a mean of the prior distribution #
  ####################################################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se)
  g
  summary(g)
  plot(g)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(g, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of A and of a regression coefficient
  ### (intercept), not using estimated values as true ones
  coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, model = "br")
  b
  summary(b)
  plot(b)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  ### (intercept), not using estimated values as true ones
  coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, model = "pr")
  p
  summary(p)
  plot(p)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  ### (intercept), not using estimated values as true ones
  coverage(p, A.or.r = 60, reg.coef = -5, nsim = 10)  

  ##################################################################################
  # If we have one covariate and do not know a mean of the prior distribution yet, #
  ##################################################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se, x2, model = "gr")
  g
  summary(g)
  plot(g)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(g, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of A, of regression coefficients, and
  ### of covariate, not using estimated values as true ones
  coverage(g, A.or.r = 9, reg.coef = c(10, 1), covariates = x2, nsim = 10)  
  ### two values of reg.coef are for beta0 and beta1

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, x1, model = "br")
  b
  summary(b)
  plot(b)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r, of regression coefficients,
  # and of covariate, not using estimated values as true ones
  coverage(b, A.or.r = 60, reg.coef = c(-1, 0), covariates = x1, nsim = 10)  
  ### two values of reg.coef are for beta0 and beta1

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, x1, model = "pr")
  p
  summary(p)
  plot(p)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r, of regression coefficients,
  ### and of covariate, not using estimated values as true ones
  coverage(p, A.or.r = 60, reg.coef = c(-2, 0), covariates = x1, nsim = 10)  
  ### two values of reg.coef are for beta0 and beta1

  ################################################
  # If we know a mean of the prior distribution, #
  ################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se, mean.PriorDist = 8)
  g
  summary(g)
  plot(g)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(g, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of A and of 2nd level mean 
  ### as true ones, not using estimated values as true ones
  coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  #########
  # BRIMM #
  #########
 
  b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")
  b
  summary(b)
  plot(b)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and of 2nd level mean 
  ### as true ones, not using estimated values as true ones
  coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  #########
  # PRIMM #
  #########

  p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")
  p
  summary(p)
  plot(p)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and of 2nd level mean 
  ### as true ones, not using estimated values as true ones
  coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

}

\keyword{package}