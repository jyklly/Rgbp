\name{gbp.default}
\alias{gbp}
\alias{gbp.default}
\title{
  Bayesian Hierarchical Modeling with (Generalized) Stein's Harmonic Priors
}
\description{
  Bayesian Hierarchical modeling for Gaussian (Stein's), Binomial (Generalized Stein's), 
  and Poisson (Generalized Stein's)Data
}
\usage{
  \method{gbp}{default}(x, y, covariates, mu0, model = "gr", intercept = T, Alpha = 0.95, ...)
}
\arguments{
  \item{x}{
	The first argument of gbp is a vector of the sample mean for GRIMM and of the number of successful 
    trials for BRIMM and PRIMM.
  }
  \item{y}{
	The second argument of gbp is a vector composed of the standard errors of all groups for GRIMM and 
    of the total numbers of trials for BRIMM and PRIMM.
  }
  \item{covariates}{
	Any covariate, if any. The intercept part will be automatically generated. The right format for 
    one covariate is a vector and for the two covariate is k by 2 matrix.
  }
  \item{mu0}{
	A prior mean, the second-level mean parameter, if you know this value a priori.
  }
  \item{model}{
	"gr" for GRIMM, "br" for BRIMM, and "pr" for PRIMM. Default is "gr"
  }
  \item{intercept}{
	User can determine whether the regression part includes the intercept term or not. Default is True.
  }
  \item{Alpha}{
	This function will estimate 100*Alpha\% intervals. Default is 95\%
  }
  \item{...}{
	Any argument for future improvement
  }
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
\author{
  Joseph Kelly and Hyungsuk Tak
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
  g <- gbp(y, se, mu0 = 8)
  g
  summary(g)
  plot(g)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(g, nsim = 100)  

  # when we want to simulate psuedo datasets based on different values of A and of 2nd level mean 
  # as true ones, not using estimated values as true ones
  coverage(g, y = 9, mu0 = 5, nsim = 10)  

  # BRIMM
  b <- gbp(z, n, mu0 = 0.265, model = "br")
  b
  summary(b)
  plot(b)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of 2nd level mean 
  # as true ones, not using estimated values as true ones
  coverage(b, y = 60, mu0 = 0.3, model = "br", nsim = 10)  

  # PRIMM
  p <- gbp(z, n, mu0 = 0.265, model = "pr")
  p
  summary(p)
  plot(p)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of 2nd level mean 
  # as true ones, not using estimated values as true ones
  coverage(p, y = 60, mu0 = 0.3, model = "br", nsim = 10)  

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
  coverage(g2, y = 9, beta = 10, nsim = 10)  

  # BRIMM
  b2 <- gbp(z, n, model = "br")
  b2
  summary(b2)
  plot(b2)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b2, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  # (intercept), not using estimated values as true ones
  coverage(b2, y = 60, beta = -1, nsim = 10)  

  # PRIMM
  p2 <- gbp(z, n, model = "pr")
  p2
  summary(p2)
  plot(p2)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p2, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  # (intercept), not using estimated values as true ones
  coverage(p2, y = 60, beta = -5, nsim = 10)  

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
  coverage(g3, y = 9, beta = c(10, 1), X = x2, nsim = 10)  

  # BRIMM
  b3 <- gbp(z, n, x1, model = "br")
  b3
  summary(b3)
  plot(b3)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b2, nsim = 10)  

  # when we want to simulate psuedo datasets based on different values of r, of regression coefficients,
  # and of covariate, not using estimated values as true ones
  coverage(b2, y = 60, beta = c(-1, 0), X = x1, nsim = 10)  

  # PRIMM
  p3 <- gbp(z, n, x1, model = "pr")
  p3
  summary(p3)
  plot(p3)

  # when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p3, nsim = 100)  

  # when we want to simulate psuedo datasets based on different values of r, of regression coefficients,
  # and of covariate, not using estimated values as true ones
  coverage(p3, y = 60, beta = c(-2, 0), X = x1, nsim = 10)  
}

\keyword{methods}
