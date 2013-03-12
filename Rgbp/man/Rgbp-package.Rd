\name{Rgbp}
\alias{Rgbp-package}
\alias{Rgbp}
\docType{package}
\title{
Bayesian Hierarchical Modeling with Steins Harmonic Prior
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
  data(schools)
  y <- schools$y
  se <- schools$se
  g <- gbp(y, se, mu0 = 8)
  print(g)
  summary(g)
  plot(g)
  coverage(g, nsim = 100)  
}
\keyword{ package }
