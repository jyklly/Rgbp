\name{Rgbp}
\alias{Rgbp-package}
\alias{Rgbp}
\docType{package}
\title{
Bayesian Hierarchical Modeling using Stein's Harmonic Prior
}
\description{
Bayesian Hierarchical modeling for Gaussian (GRIMM), Binomial (BRIMM) and Poisson (PRIMM) data assuming generalized Stein's harmonic priors.
}
\details{
  \tabular{ll}{
    Package: \tab Rgbp\cr
    Type: \tab Package\cr
    Version: \tab 1.0.0\cr
    Date: \tab 2012-10-31\cr
    License: \tab GPL-2\cr
  }
}
\author{
  Joseph Kelly, Carl Morris, and Hyungsuk Tak

  Maintainer: Joseph Kelly <kelly2@fas.harvard.edu>
}
\references{
  Morris, C. and Tang, R. (2011). \emph{Estimating Random Effects via Adjustment for Density Maximization}. Statistical Science. 26. 271-287.
  Morris, C. and Lysy, M. (2012). \emph{Shrinkage Estimation in Multilevel Normal Models}. Statistical Science. 27. 115-134.
}

\examples{
 
  # Loading datasets
  data(baseball)
  data(schools)

  z <- baseball$Hits
  n <- baseball$At.Bats
  x1 <- rep(c(-1, 0, 1), 6)

  y <- schools$y
  se <- schools$se
  x2 <- rep(c(-1, 0, 1, 2), 2)

  # If we know a prior mean, the second-level mean parameter,
  # GRIMM
  g <- gbp(y, se, mean.PriorDist = 8)
  g
  summary(g)
  plot(g)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(g, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of A and of 2nd level mean 
  # as true ones, not using estimated values as true ones
  coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  # BRIMM
  b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")
  b
  summary(b)
  plot(b)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of 2nd level mean 
  # as true ones, not using estimated values as true ones
  coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  # PRIMM
  p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")
  p
  summary(p)
  plot(p)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of 2nd level mean 
  # as true ones, not using estimated values as true ones
  coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  # If we do not know a prior mean (2nd-level mean) and do not have any covariate,
  # GRIMM
  g2 <- gbp(y, se)
  g2
  summary(g2)
  plot(g2)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(g2, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of A and of a regression coefficient
  # (intercept), not using estimated values as true ones
  coverage(g2, A.or.r = 9, reg.coef = 10, nsim = 10)  

  # BRIMM
  b2 <- gbp(z, n, model = "br")
  b2
  summary(b2)
  plot(b2)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b2, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  # (intercept), not using estimated values as true ones
  coverage(b2, A.or.r = 60, reg.coef = -1, nsim = 10)  

  # PRIMM
  p2 <- gbp(z, n, model = "pr")
  p2
  summary(p2)
  plot(p2)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p2, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  # (intercept), not using estimated values as true ones
  coverage(p2, A.or.r = 60, reg.coef = -5, nsim = 10)  

  # If we have one covariate,
  # GRIMM
  g3 <- gbp(y, se, x2, model = "gr")
  g3
  summary(g3)
  plot(g3)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(g3, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of A, of regression coefficients, and
  # of covariate, not using estimated values as true ones
  coverage(g3, A.or.r = 9, reg.coef = c(10, 1), covariates = x2, nsim = 10)  

  # BRIMM
  b3 <- gbp(z, n, x1, model = "br")
  b3
  summary(b3)
  plot(b3)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b2, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r, of regression coefficients,
  # and of covariate, not using estimated values as true ones
  coverage(b2, A.or.r = 60, reg.coef = c(-1, 0), covariates = x1, nsim = 10)  

  # PRIMM
  p3 <- gbp(z, n, x1, model = "pr")
  p3
  summary(p3)
  plot(p3)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p3, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r, of regression coefficients,
  # and of covariate, not using estimated values as true ones
  coverage(p3, A.or.r = 60, reg.coef = c(-2, 0), covariates = x1, nsim = 10)  
}

\keyword{package}