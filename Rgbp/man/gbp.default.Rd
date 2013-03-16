\name{gbp}			
\alias{gbp}
\alias{gbp.default}
\title{Fitting Bayesian Hierarchical Models}
\description{
  \code{gbp} is used to fit Bayesian hierarchical models for Gaussian (GRIMM), Binomial (BRIMM), 
  and Poisson (PRIMM) data using generalized Stein's harmonic prior for
  good frequentist repeated sampling property.
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
    included in the regression. Default is \code{TRUE}.
  }
  \item{Alpha}{
    a \code{float} between 0 and 1 to estimate 100*Alpha\% intervals. Default is 0.95.
  }
}

\details{
  \code{gbp} fits a Bayesian hierarchical model using generalized
  Stein's harmonic prior which enables good frequentist repeated sampling properties. The first-level is a distribution of observed data (likelihood) and the second-level is the prior distribution on the first-level parameter. 

  To be specific, for Normal data, \code{gbp} constructs a two-level Normal-Normal multilevel model 
  (\eqn{\sigma[j]^{2}}{\sigma_j^2} is assumed to be known and subscript \emph{j} indicates \emph{j}-th group 
  (or unit) in a dataset):
  \deqn{(y[j] | \theta[j]) ~ N(\theta[j], \sigma[j]^{2})}{(y_j | \theta_j) ~ N(\theta_j, \sigma_j^2)}
  \deqn{(\theta[j] | A, \mu[0j]) ~ N(\mu[0j], A)}{(\theta_j | A, \mu0_j) ~ N(\mu0_j, A)}
  \deqn{\mu[0j] = x[j]'\beta}{\mu0_j = x_j'\beta}
  for \eqn{j = 1, \ldots, k}.

  For Poisson data, \code{gbp} builds a two-level Poisson-Gamma multilevel model 
  (a square bracket below indicates [mean, variance] of distribution):
  \deqn{(z[j] | \theta[j]) ~ Pois(n[j]\theta[j])}{(z_j | \theta_j) ~ Pois(n_j\theta_j)}
  \deqn{(\theta[j] | r, \mu[0j]) ~ Gamma(r\mu[0j], scale = 1 / r) ~ Gamma[\mu[0j], \mu[0j] / r]}{(\theta_j | r, \mu0_j) ~ Gamma(r\mu0_j, scale = 1 / r) ~ Gamma[\mu0_j, \mu0_j / r]}
  \deqn{log(\mu[0j]) = x[j]'\beta}{log(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}.

  For Binomial data, \code{gbp} sets a two-level Binomial-Beta multilevel model:
  \deqn{(z[j] | \theta[j]) ~ Bin(n[j], \theta[j])}{(z_j | \theta_j) ~ Bin(n_j, \theta_j)}
  \deqn{(\theta[j] | r, \mu[0j]) ~ Beta(r\mu[0j], r\mu[0j]) ~ Beta[\mu[0j], mu[0j](1 - mu[0j]) / (r + 1)]}{(\theta_j | r, \mu0_j) ~ Beta(r\mu0_j, r(1 - \mu0_j)) ~ Beta[\mu0_j, \mu0_j(1 - \mu0_j) / (r + 1)]}
  \deqn{logit(\mu[0j]) = x[j]'\beta}{logit(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}.

  Theoretically, generalized Stein's harmonic prior is Uniform on the second level variance component
  (variance of the prior distribution), \emph{i.e.}, dA for GRIMM and d(1 / r) 
  (=\eqn{\frac{1}{r^{2}}}{dr / r^2}) for BRIMM and PRIMM, leading to proper posterior distributions. 

  Under this setting, the argument \code{x} in \code{gbp} is a (\emph{k} by 1) vector of the sample mean
  (\emph{y}) for GRIMM and of the number of successful trials (\emph{z}) for BRIMM and PRIMM, where \emph{k} 
  is the number of groups (or units) in a dataset.

  The argument \code{y} in \code{gbp} is a (\emph{k} by 1) vector composed of the standard errors 
  (\eqn{\sigma}) of all groups for GRIMM and the total numbers of trials (\emph{n}) for BRIMM and PRIMM.

  As for two optional arguments, \code{covariates} and \code{mean.PriorDist}, there are three feasible 
  combinations of them to run \code{gbp}. The first situation is when we do not have any covariate and do not 
  know mean of the prior distribution (\eqn{\mu[0]}{\mu0}) a priori. In this case, assigning none of two
  optional arguments, such as "\code{gbp(z, n, model = "br")}", will lead to a correct model. \code{gbp} 
  will automatically fit a regression with only an intercept term to estimate common mean of the prior
  distribution (exchangeability).

  The second situation is when we have some covariates (a \emph{k} by \emph{t} matrix, where \emph{t} 
  \eqn{\ge} 1) and do not know mean of the prior distribution (\eqn{\mu[0]}{\mu0}) a priori. In this case, 
  assigning a \emph{k} by \emph{t} matrix (each column corresponds to one covariate), \emph{X}, such as 
  "\code{gbp(z, n, X, model = "pr")}", will lead to a right model. Default of \code{gbp} is to fit a
  regression with an intercept term to estimate mean of the prior distribution. Double exchangeability will
  hold in this case.

  The last case is when we know mean of the prior distribution (\eqn{\mu[0]}{\mu0}) a priori. Now, we do
  not need to estimate regression coefficients at all because we know a true value of \eqn{\mu[0]}{\mu0}.
  Designating this value into the argument of \code{gbp} like 
  "\code{gbp(y, se, mean.PriorDist = 3, model = "gr")}" is enough to account for it. For reference, 
  \code{mean.PriorDist} has the stronger priority than \code{covariates}, which means that when both
  arguments are designated, \code{gbp} will fit a hierarchical model with known mean of prior distribution.

  When it comes to estimating hyper-parameters, (\emph{A} or \emph{r} and \eqn{\beta}, \code{gbp} uses 
  a mixture model with \eqn{\theta} integrated out, searching for a value that maximizes marginal posterior
  distribution of \eqn{\alpha} (= log(\emph{A}) for GRIMM and -log(\emph{r}) for BRIMM and PRIMM). 
  And then, \code{gbp} looks for estimated \eqn{\beta} that maximizes its likelihood given previously 
  estimated \eqn{\alpha}. \code{optim} is used for maximizing functions. 

  \code{gbp} returns an object of class \code{gbp} which provides the functions \code{plot}, \code{print},
  and \code{summary}.
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
  \item{X}{a covariate vector or matrix if any. NA if not}
  \item{beta.new}{regression coefficient estimates}
  \item{beta.var}{estimated variance matrix of regression coefficient}
  \item{intercept}{whether TRUE or FALSE}
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
  data(schools)

  # baseball data where z is Hits and n is at bats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)
  x1 <- rep(c(-1, 0, 1), 6)

  y <- schools$y
  se <- schools$se
  x2 <- rep(c(-1, 0, 1, 2), 2)
  
  ####################################################################################
  # If we do not have any covariates and do not know a mean of the prior distribution #
  ####################################################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se)
  g
  summary(g)
  plot(g)

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(g, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of A and of a regression coefficient
  ### (intercept), not using estimated values as true ones
  coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, model = "br")
  b
  summary(b)
  plot(b)

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of r and of a regression coefficient
  ### (intercept), not using estimated values as true ones
  coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, model = "pr")
  p
  summary(p)
  plot(p)

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of r and of a regression coefficient
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

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of r, of regression coefficients,
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

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of r, of regression coefficients,
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

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(g, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of A and of 2nd level mean 
  ### as true ones, not using estimated values as true ones
  coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  #########
  # BRIMM #
  #########
 
  b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")
  b
  summary(b)
  plot(b)

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of r and of 2nd level mean 
  ### as true ones, not using estimated values as true ones
  coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  #########
  # PRIMM #
  #########

  p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")
  p
  summary(p)
  plot(p)

  ### when we want to simulated pseudo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  ### when we want to simulate pseudo datasets based on different values of r and of 2nd level mean 
  ### as true ones, not using estimated values as true ones
  coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  
}


\keyword{methods}

%%  LocalWords:  gbp BRIMM PRIMM
