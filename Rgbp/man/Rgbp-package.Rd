\name{Rgbp}
\alias{Rgbp-package}
\alias{Rgbp}
\docType{package}
\title{
Bayesian Hierarchical Modeling with Stein's and Generalized Stein's Harmonic Priors
}
\description{
Bayesian Hierarchical modeling for Gaussian (Stein's), Binomial (Generalized Stein's), and Poisson (Generalized Stein's)Data
}
\details{
  \tabular{ll}{
    Package: \tab Rgbp\cr
    Type: \tab Package\cr
    Version: \tab 0.1-0\cr
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
library(Rgbp)
# If the prior mean is known in advance (the mean of the second level Hierarchy, for example, 0.8),
g <- gbp(y,se,mu=0.8)
b <- gbp(z,n,mu=0.8,model="br")
p <- gbp(z,n,mu=0.8,model="pr")

# If the prior mean is not known in advance (the mean of the second level Hierarchy, for example, 0.8), the regression will be an intercept model.
g<-gbp(y,se)
b<-gbp(z,n,model="br")
p<-gbp(z,n,model="pr")

# If you have any covariate(s),
g<-gbp(y,se,x)
b<-gbp(z,n,x,model="br")
p<-gbp(z,n,x,model="pr")
}
\keyword{ package }
