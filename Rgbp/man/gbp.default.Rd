\name{gbp}			
\alias{gbp}
\alias{gbp.default}
\title{Fitting Bayesian Hierarchical Models}
\description{
  \code{gbp} fits a Bayesian hierarchical model using the Uniform prior distribution on the second level variance component (variance of the prior distribution), which enables good frequentist repeated sampling properties. 
}
\usage{\method{gbp}{default}(x, w, covariates, mean.PriorDist, model = "gr", intercept = TRUE, Alpha = 0.95)}
\arguments{
  \item{x}{
    a (\emph{k} by 1) vector of \emph{k} groups' sample means for Gaussian or of each group's number of successful trials for Binomial and Poisson data, where \emph{k} is the number of groups (or units) in a dataset.
  }
  \item{w}{
    a (\emph{k} by 1) vector composed of the standard errors of all groups for Gaussian or 
    of each group's total number of trials for Binomial and Poisson data.
  }
  \item{covariates}{
    (optional) a matrix of covariate(s) each column of which corresponds to one covariate.
  }
  \item{mean.PriorDist}{
    (optional) a numeric value for the second-level mean parameter, \emph{i.e.} the mean of prior
    distribution, if you know this value a priori.
  }
  \item{model}{
    a character string indicating which hierarchical model to fit. "gr" for Gaussian data, "br" for Binomial, 
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
  \code{gbp} fits a hierarchical model whose first-level hierarchy is a distribution of observed data and second-level is a conjugate prior distribution on the first-level parameter. To be specific, for Normal data, \code{gbp} constructs a two-level Normal-Normal multilevel model
  (\eqn{\sigma^{2}_{j}}{\sigma_j^2} is assumed to be known or to be accurately estimated, and subscript \emph{j} indicates \emph{j}-th group 
  (or unit) in a dataset):
  \deqn{(y_{j}~ |~ \theta_{j})~ \sim ~ indep~ N(\theta_{j}, \sigma^{2}_{j})}{(y_j | \theta_j) ~ indep N(\theta_j, \sigma_j^2)}
  \deqn{(\theta_{j}~ |~ A,~ \mu_{0j})~ \sim ~  indep~ N(\mu_{0j},~ A)}{(\theta_j | A, \mu0_j) ~ indep N(\mu0_j, A)}
  \deqn{\mu_{0j}~ =~ x_{j}'\beta}{\mu0_j = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For Poisson data, \code{gbp} builds a two-level Poisson-Gamma multilevel model (for reference, a square bracket below indicates [mean, variance] of distribution and a constant multiplied to the notation representing Gamma distribution (Gam) is a scale. Also, for consistent notation, \eqn{y_{j}=\frac{z_{j}}{n_{j}}}{y_j = z_j / n_j}):
  \deqn{(z_{j}~ |~ \theta_{j})~ \sim ~ indep~ Pois(n_{j}\theta_{j})}{(z_j | \theta_j) ~ indep Pois(n_j\theta_j)}
  \deqn{(\theta_{j}~ |~ r,~ \mu_{0j})~ \sim ~ indep~ \frac{1}{r}Gam(r\mu_{0j})~ \sim ~ indep~ Gam[\mu_{0j}, \mu_{0j} / r]}{(\theta_j | r, \mu0_j) ~ indep Gam(r\mu0_j) / r ~ indep Gam[\mu0_j, \mu0_j / r]}
  \deqn{log(\mu_{0j})~ =~ x_{j}'\beta}{log(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For Binomial data, \code{gbp} sets a two-level Binomial-Beta multilevel model (for reference, a square bracket below indicates [mean, variance] of distribution and \eqn{y_{j} = \frac{z_{j}}{n_{j}}}{y_j = z_j / n_j}):
  \deqn{(z_{j}~ |~ \theta_{j})~ \sim ~ indep~ Bin(n_{j}, \theta_{j})}{(z_j | \theta_j) ~ indep Bin(n_j, \theta_j)}
  \deqn{(\theta_{j}~ |~ r, \mu_{0j})~ \sim ~ indep~ Beta(r\mu_{0j},~ r(1-\mu_{0j}))~ \sim ~ indep~ Beta[\mu_{0j},~ \mu_{0j}(1 - \mu_{0j})~ /~ (r + 1)]}{(\theta_j | r, \mu0_j) ~ indep Beta(r\mu0_j, r(1 - \mu0_j)) ~ indep Beta[\mu0_j, \mu0_j(1 - \mu0_j) / (r + 1)]}
  \deqn{logit(\mu_{0j})~ =~ x_{j}'\beta}{logit(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For reference, based on the above notations, the Uniform prior distribution on the second level variance component (variance of the prior distribution) is \emph{dA} for Gaussian and \eqn{d(\frac{1}{r})}{d(1 / r)} 
  (= \eqn{\frac{dr}{r^{2}}}{dr / r^2}) for Binomial and Poisson data.

  Under this setting, the argument \code{x} in \code{gbp} is a (\emph{k} by 1) vector of \emph{k} groups' sample means (\eqn{y_{j}'s}{y_j's} in the description of Normal-Normal model above) for Gaussian or of each group's number of successful trials (\eqn{z_{j}'s}{z_j's}) for Binomial and Poisson data, where \emph{k} is the number of groups (or units) in a dataset.

  The argument \code{w} is a (\emph{k} by 1) vector composed of the standard errors (\eqn{\sigma_{j}'s}{\sigma_j's} or accurately estimated \eqn{s_{j}'s}{s_j's}) of all groups for Gaussian or of each group's total number of trials (\eqn{n_{j}'s}{n_j's}) for Binomial and Poisson data.

  As for two optional arguments, \code{covariates} and \code{mean.PriorDist}, there are three feasible 
  combinations of them to run \code{gbp}. The first situation is when we do not have any covariate and do not 
  know a mean of the prior distribution (\eqn{\mu_{0}}{\mu0}) a priori. In this case, assigning none of two
  optional arguments, such as "\code{gbp(z, n, model = "br")}", will lead to a correct model. \code{gbp} 
  will automatically fit a regression with only an intercept term to estimate a common mean of the prior
  distribution (exchangeability).

  The second situation is when we have covariate(s) and do not know a mean of the prior distribution (\eqn{\mu_{0}}{\mu0}) a priori. In this case, assigning a matrix, \emph{X}, each column of which corresponds to one covariate, such as "\code{gbp(z, n, X, model = "pr")}", will lead to a correct model. Default of \code{gbp} is to fit a regression including an intercept term to estimate a mean of the prior distribution. Double exchangeability will hold in this case.

  The last case is when we know a mean of the prior distribution (\eqn{\mu_{0}}{\mu0}) a priori. Now, we do
  not need to estimate regression coefficients at all because we know a true value of \eqn{\mu_{0}}{\mu0} (strong assumption).
  Designating this value into the argument of \code{gbp} like 
  "\code{gbp(y, se, mean.PriorDist = 3)}" is enough to account for it. For reference, 
  \code{mean.PriorDist} has a stronger priority than \code{covariates}, which means that when both
  arguments are designated, \code{gbp} will fit a hierarchical model using the known mean of prior distribution, \code{mean.PriorDist}.

  \code{gbp} returns an object of class "\code{gbp}" which provides three relevant functions \code{plot}, \code{print}, and \code{summary}.
}

\value{
  An object of class \code{gbp} comprises of:
  \item{sample.mean}{sample mean of each group (or unit)}
  \item{se}{if Gaussian data, standard error of sample mean in each group (or unit)}
  \item{n}{if Binomial and Poisson data, total number of trials of each group (or unit)}
  \item{prior.mean}{numeric if entered, NA if not entered}
  \item{prior.mean.hat}{estimate of prior mean by a regression if prior mean is not assigned a priori}
  \item{shrinkage}{shrinkage estimate of each group}
  \item{sd.shrinkage}{standard deviation of shrinkage estimate}
  \item{post.mean}{posterior mean of each group}
  \item{post.sd}{posterior standard deviation of each group}
  \item{post.intv.low}{lower bound of 100*Alpha\% posterior interval}
  \item{post.intv.upp}{upper bound of 100*Alpha\% posterior interval}
  \item{model}{"gr" for Gaussian, "br" for Binomial, and "pr" for Poisson data}
  \item{X}{a covariate vector or matrix if designated. NA if not}
  \item{beta.new}{regression coefficient estimates}
  \item{beta.var}{estimated variance matrix of regression coefficient}
  \item{intercept}{whether TRUE or FALSE}
  \item{a.new}{a posterior mode of \eqn{\alpha} defined as log(\emph{A}) for Gaussian or -log(\emph{r}) for Binomial and Poisson data. Practical meaning of estimating \eqn{\alpha} will appear in \code{summary(gbp.object)}.}
  \item{a.var}{variance of alpha estimate}
}

\references{
  Morris, C. and Lysy, M. (2012). Shrinkage Estimation in Multilevel Normal Models. 
  \emph{Statistical Science}. \bold{27}. 115-134.
}

\author{
  Joseph Kelly, Carl Morris, and Hyungsuk Tak
}
\examples{

  # Loading datasets
  data(schools)
  y <- schools$y
  se <- schools$se

  # Arbitrary covariate for schools data
  x2 <- rep(c(-1, 0, 1, 2), 2)
  
  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)

  # Arbitrary covariate for baseball data
  x1 <- rep(c(-1, 0, 1), 6)

  ########################################################
  # Gaussioan Regression Interactive Multilevel Modeling #
  ########################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    g <- gbp(y, se)
    g
    summary(g)
    plot(g)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of A and 
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    gcv <- coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gr")
    g
    summary(g)
    plot(g)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  
 
    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of A,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    gcv <- coverage(g, A.or.r = 9, reg.coef = c(10, 1), covariates = x2, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8)
    g
    summary(g)
    plot(g)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of A and
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  #######################################################
  # Binomial Regression Interactive Multilevel Modeling #
  #######################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "br")
    b
    summary(b)
    plot(b)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of r and 
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "br")
    b
    summary(b)
    plot(b)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    bcv <- coverage(b, A.or.r = 60, reg.coef = c(-1, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")
    b
    summary(b)
    plot(b)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of r and
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  ######################################################
  # Poisson Regression Interactive Multilevel Modeling #
  ######################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    p <- gbp(z, n, model = "pr")
    p
    summary(p)
    plot(p)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of r and 
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, reg.coef = -5, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    p <- gbp(z, n, x1, model = "pr")
    p
    summary(p)
    plot(p)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    pcv <- coverage(p, A.or.r = 60, reg.coef = c(-2, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")
    p
    summary(p)
    plot(p)

    ### when we want to simulate psuedo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate psuedo datasets based on different values of r and
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

}


\keyword{methods}
