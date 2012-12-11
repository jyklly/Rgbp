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
If you know the prior mean in advance (the mean of the second level Hierarchy, for example, 0.8),
g<-gbp(y,se,mu=0.8)
b<-gbp(z,n,mu=0.8,model="br")
p<-gbp(z,n,mu=0.8,model="pr")

If you do not know the prior mean in advance (the mean of the second level Hierarchy, for example, 0.8), the regression will be an intercept model.
g<-gbp(y,se)
b<-gbp(z,n,model="br")
p<-gbp(z,n,model="pr")

If you have any covariate(s),
g<-gbp(y,se,x)
b<-gbp(z,n,x,model="br")
p<-gbp(z,n,x,model="pr")
}
\author{
Joseph Kelly, Hyung Suk Tak, and Carl Morris

Maintainer: Joseph Kelly <kelly2@fas.harvard.edu>
}
\references{
Morris, C. and Tang, R. (2011). \emph{Estimating Random Effects via Adjustment for Density Maximization}. Statistical Science. 26. 271-287.
Morris, C. and Lysy, M. (2012). \emph{Shrinkage Estimation in Multilevel Normal Models}. Statistical Science. 27. 115-134.
}
\keyword{ package }
\examples{
data(baseball)
data(schools)

z<-baseball$Hits
n<-baseball$At.Bats
y<-schools$y
se<-schools$se

# If we know a prior mean, the second-level mean parameter,
# GRIMM when prior mean is known
g <- gbp(y,se,mu=8.75)
g
summary(g)
plot(g)

# BRIMM when prior mean is known
b <- gbp(z, n, mu=0.265, model="br")
b
summary(b)
plot(b)

# PRIMM when prior mean is known
p <- gbp(z, n, mu=0.265, model="pr")
p
summary(p)
plot(p)

# If we do not know a prior mean and do not have any covariate,
# GRIMM without covariate when prior mean is unknown
g <- gbp(y,se)
g
summary(g)
plot(g)

# BRIMM without covariate when prior mean is unknown
b <- gbp(z, n, model="br")
b
summary(b)
plot(b)

# PRIMM without covariate when prior mean is unknown
p <- gbp(z, n, model="pr")
p
summary(p)
plot(p)

# If we have two covariates
m <- 2
X <- matrix(rnorm(length(z)*m),nrow=length(z),ncol=m)
X2 <- matrix(rnorm(length(y)*m),nrow=length(y),ncol=m)

# GRIMM with two covariates
g2 <- gbp(schools$y, schools$se, X2, model="gr")
g2
summary(g2)
plot(g2)

# BRIMM with two covariates
b2 <- gbp(z, n, X, model="br")
b2
summary(b2)
plot(b2)

# BRIMM with two covariates
p2 <- gbp(z, n, X, model="pr")
p2
summary(p2)
plot(p2)
}
