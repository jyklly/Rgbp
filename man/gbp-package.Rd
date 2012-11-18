\name{gbp-package}
\alias{gbp-package}
\alias{gbp}
\docType{package}
\title{
Empirical Bayesian Hierarchical modeling
}
\description{
Empirical Bayesian Hierarchical modeling for Gaussian, Binomial, and Poisson Data
}
\details{
\tabular{ll}{
Package: \tab gbp\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2012-10-31\cr
License: \tab GPL-2\cr
}
 An overview of how to use the package, including the most important functions
}
\author{
Joseph Kelly, Hyung Suk Tak, and Carl Morris
Maintainer: Joseph Kelly <kelly2@fas.harvard.edu>
}
\references{
Morris, C. and Tang, R. (2011), "Estimating Random Effects via Adjustment for Density Maximization", Statistical Science 26 271-287
Morris, C. and Lysy, M. (2012), "Shrinkage Estimation in Multilevel Normal Models", Statistical Science 27 115-134
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(sn)
source("../R/gr.R")
source("../R/gbp.R")
source("../R/bp.R")
load("../data/baseball.RData")
load("../data/schools.RData")

z<-baseball$Hits
n<-baseball$At.Bats

b <- gbp(z, n, model="br")
p <- gbp(z, n, model="pr")
g <- gbp(schools$y,schools$se)
b
summary(b)
plot(b)
p
summary(p)
plot(p)
g
summary(g)
plot(g)

# if we have any covariates
p <- 2
X <- matrix(rnorm(length(z)*p),nrow=length(z),ncol=p)
X2 <- matrix(rnorm(length(schools$y)*p),nrow=length(schools$y),ncol=p)
b2 <- gbp(z, n, X, model="br")
p2 <- gbp(z, n, X, model="pr")
g2 <- gbp(schools$y, schools$se, X2, model="gr")
b2
summary(b2)
plot(b2)
p2
summary(p2)
plot(p2)
g2
summary(g2)
plot(g2)
}
