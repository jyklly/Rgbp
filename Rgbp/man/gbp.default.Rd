\name{gbp}
\alias{gbp}
\alias{gbp.default}
\title{Fitting Bayesian Hierarchical Models}
\description{
  \code{gbp} is used to fit Bayesian hierarchical models for Gaussian (GRIMM), Binomial (BRIMM), 
  and Poisson (PRIMM) data using generalized Stein's harmonic prior for good frequntist repeated sampling
  perperty.
}
\usage{\method{gbp}{default}(x, y, covariates, mean.PriorDist, model = "gr", intercept = TRUE, Alpha = 0.95)}
\arguments{
  \item{x}{
    a (\emph{k} by 1) vector of the sample mean for GRIMM and of the number of successful 
    trials for BRIMM and PRIMM,\cr where \emph{k} is the number of groups (or units) in a dataset.
  }
  \item{y}{
    a (\emph{k} by 1) vector composed of the standard errors of all groups for GRIMM and 
    of the total numbers of trials for BRIMM and PRIMM.
  }
  \item{covariates}{
    (optional) a (\emph{k} by \emph{t}) matrix of covariate values, where \emph{t} is the number of 
    covariates (\emph{t} \eqn{\ge} 1). To fit an intercept, do not include a vector of ones but rather use 
    the intercept argument.
  }
  \item{mean.PriorDist}{
    (optional) a numeric value for the second-level mean parameter, \emph{i.e.} the mean of prior
    distribution, if you know this value a priori.
  }
  \item{model}{
    a character string indicating which hierarchical model to fit. "gr" for Gaussian, "br" for Binomial, 
    and "pr" for Poisson. Default is "gr"
  }
  \item{intercept}{
    \code{TRUE} or \code{FALSE} flag indicating whether an intercept should be
    included in the regressuin. Default is \code{TRUE}.
  }
  \item{Alpha}{
    a \code{float} between 0 and 1 to estimate 100*Alpha\% intervals. Default is 0.95.
  }
}

\details{
  \code{gbp} fits Bayesian hierarchical models using generalized Stein's harmonic prior for good frequentist
  repeated sampling property.

  For Normal data, \code{gbp} constructs a two-level Normal-Normal multilevel model 
  (\eqn{\sigma^{2}}{\sigma^2} is assumed to be known):
  \deqn{(y[j]|\theta[j]) ~ N(\theta[j], \sigma^{2})}{(y|\theta) ~ N(\theta, \sigma^2)}
  \deqn{(\theta[j]|A, \mu[0j]) ~ N(\mu[0j], A)}{(\theta|A, \mu0) ~ N(\mu0, A)}
  \deqn{\mu[0j] = x[j]'\beta}{\mu0 = x'\beta}

  For Poisson data, \code{gbp} builds a two-level Poisson-Gamma multilevel model:
  \deqn{(z[j]|\lambda[j]) ~ Pois(n[j]\lambda[j])}{(z|\lambda) ~ Pois(n\lambda)}
  \deqn{(\lambda[j]|r, \mu[0j]) ~ Gamma(r\mu[0j], r) ~ Gamma[\mu[0j], \mu[0j] / r)}
  {(\lambda|r, \mu0) ~ Gamma(r\mu0, scale = 1 / r) ~ Gamma[\mu0, \mu0 / r]}
  \deqn{log(\mu[0j]) = x[j]'\beta}{log(\mu0) = x'\beta}

  For Binomial data, \code{gbp} sets a two-level Binomial-Beta multilevel model:
  \deqn{(z[j]|p[j]) ~ Bin(n[j], p[j])}{(z|p) ~ Bin(n, p)}
  \deqn{(p[j]|r, \mu[0j]) ~ Beta(r\mu[0j], r\mu[0j])}{(p|r, \mu0) ~ Beta(r\mu0, r(1 - \mu0))}
  \deqn{logit(\mu[0j]) = x[j]'\beta}{logit(\mu0) = x'\beta}


  the distribution of hierarchy is Normal with mean \eqn{\theta[j]}{\theta} and variance
  \eqn{sigma^[2]}{sigma^2} assumed to be known. Asfor the 2nd hierarchy, a conjugate prior distribution on the parameter 
  of  1st level is 
  \eqn{p(\theta[j] | \mu[0j], A (or r)) = p[\mu[0j], A]}{p(\theta | \mu0, A (or r)) = p[\mu0, A]} for GRIMM,
  \eqn{p[\mu[0j], \mu[0j] / r]}{p[\mu0, \mu0 / r]} for PRIMM, 
  and \eqn{p[\mu[0j], \mu[0j] * (1 - \mu[0j]) / (r + 1)]}{p[\mu0, \mu0 * (1 - \mu0) / (r + 1)]} for BRIMM,
  where \eqn{g(\mu[0j]) = x[j]'\beta}{g(\mu0) = x'\beta}.
by maximizing the adjusted
  posterior using the quasi-Newton method "BFGS" via \code{optim}. By
  adjusting the posterior greater coverage is achieved for simiarly
  sized intervals (Morris, Tang, 2011). The function returns an object
  of class \code{gbp} which provides the functions \code{plot} and \code{summary}.
}

\value{
  An object of class \code{gbp} comprises of:
  \item{sample.mean}{sample mean of each group}
  \item{se}{if GRIMM, standard error of each group}
  \item{n}{if BRIMM and PRIMM, total number of trials of each group}
  \item{prior.mean}{numeric if entered, NA if not entered}
  \item{prior.mean.hat}{estimate of prior mean if prior mean is not assigned a priori}
  \item{shrinkage}{shrinkage estimate of each group}
  \item{sd.shrinkage}{standard error of shrinkage estimate}
  \item{post.mean}{posterior mean of each group}
  \item{post.sd}{posterior standard deviation of each group}
  \item{post.intv.low}{lower bound of 100*Alpha\% posterior interval}
  \item{post.intv.upp}{upper bound of 100*Alpha\% posterior interval}
  \item{model}{"gr" for GRIMM, "br" for BRIMM, and "pr" for PRIMM}
  \item{X}{covariate vector or matrix if any. NA if not}
  \item{beta.new}{regression coefficient estimates}
  \item{beta.var}{estimated variance matrix of regression coefficient}
  \item{intercept}{whether T or F}
  \item{a.new}{alpha estimate}
  \item{a.var}{variance of alpha estimate}
}

\references{
  Morris, C. and Tang, R. (2011). Estimating Random Effects via
  Adjustment for Density Maximization. \emph{Statistical Science}. \bold{26}. 271-287.

  Morris, C. and Lysy, M. (2012). Shrinkage Estimation in Multilevel Normal Models. 
  \emph{Statistical Science}. \bold{27}. 115-134.
}

\author{
  Joseph Kelly, Carl Morris, and Hyungsuk Tak
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


\keyword{methods}
