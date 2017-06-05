\name{gbp}			
\alias{gbp}
\title{Fitting Gaussian, Poisson, and Binomial Hierarchical Models}
\description{
  \code{gbp} fits Bayesian hierarchical models using the Uniform distribution on the second level variance component (variance of the prior distribution), which enables good frequentist repeated sampling properties. 
}

\usage{gbp(y, se.or.n, covariates, mean.PriorDist, model, intercept, 
           confidence.lvl, n.AR, n.AR.factor, trial.scale, save.result,
           normal.CI, t, u)}

\arguments{
  \item{y}{
    a (\emph{k} by 1) vector of \emph{k} groups' sample means for Gaussian or of each group's number of successful trials for Binomial and Poisson data, where \emph{k} is the number of groups (or units) in a dataset.
  }
  \item{se.or.n}{
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
    a character string indicating which hierarchical model to fit. "gaussian" for Gaussian data, 
    "poisson" for Poisson, and "binomial" for Binomial. Default is "gaussian"
  }
  \item{intercept}{
    \code{TRUE} or \code{FALSE} flag indicating whether an intercept should be
    included in the regression. Default is \code{TRUE}.
  }
  \item{confidence.lvl}{
    a \code{float} between 0 and 1 to estimate 100*confidence.lvl\% intervals. Default is 0.95.
  }

  \item{n.AR}{
    Only for Binomial model. If \code{n.AR = 1000}, all the results will be based on 1000 posterior samples using the acceptance-rejection sampling. Default is 0.
  }
  \item{n.AR.factor}{
    Only for Binomial model. If \code{n.AR = 1000} and \code{n.AR.factor = 4}, then the acceptance-rejection method will sample 4 * 1000 trial samples, and accept or reject them until the method gets n.AR posterior samples. Default is 4.
  }

  \item{trial.scale}{
    A scale used in the trial distribution of \eqn{\alpha}. If resultant weight has too mamy 0's, scale should be smaller than before. If resultant weight has too few 0's, scale should be larger than before. If there are relatively huge weights, scale should be larger than before. Default is 1.3.
  }

  \item{save.result}{
    Only for Binomial model with the acceptance-rejection sampling. If \code{save.result} is TRUE, all the results of weights and posterior samples will be saved in the gbp object. Default is TRUE.
  }

  \item{normal.CI}{
    Only applicable for Gaussian data. If \code{TRUE} then a Normal
    approximation is used to construct intervals for the level 1 group
    means. If \code{FALSE} (default) then a Skew-Normal distribution is
    used. Setting the value to \code{TRUE} may result in speed
    improvements but may lead to intervals that under cover. 
  }

  \item{t}{
    Non-negative constant to determine the hyper-prior distribution of r for the Binomial model with the acceptance-rejection method. If t is positive, then the hyper-prior distribution of r is proper, otherwise improper. \eqn{dr/(t + r)^{u+1}}.
  }

  \item{u}{
    A positive constant to determine the hyper-prior distribution of r for the Binomial model with the acceptance-rejection method. \eqn{dr/(t + r)^{u+1}}.
  }

}

\details{
  \code{gbp} fits a hierarchical model whose first-level hierarchy is a distribution of observed data and second-level is a conjugate prior distribution on the first-level parameter. To be specific, for Normal data, \code{gbp} constructs a two-level Normal-Normal multilevel model. \eqn{V_{j}}{V_j} (=\eqn{\sigma^{2}/n_{j}}{\sigma^2/n_j}) is assumed to be known or to be accurately estimated, and subscript \emph{j} indicates \emph{j}-th group 
  (or unit) in a dataset.
  \deqn{(y_{j} | \theta_{j}) \stackrel{ind}{\sim} N(\theta_{j}, \sigma^{2}_{j})}{(y_j | \theta_j) ~ indep N(\theta_j, \sigma_j^2)}
  \deqn{(\theta_{j} |\mu_{0j} , A) \stackrel{ind}{\sim} N(\mu_{0j}, A)}{(\theta_j | \mu0_j, A) ~ indep N(\mu0_j, A)}
  \deqn{\mu_{0j} = x_{j}'\beta}{\mu0_j = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For Poisson data, \code{gbp} builds a two-level Poisson-Gamma multilevel model. A square bracket below indicates [mean, variance] of distribution, a constant multiplied to the notation representing Gamma distribution (Gam) is a scale, and \eqn{y_{j}=\frac{z_{j}}{n_{j}}}{y_j = z_j / n_j}.
  \deqn{(z_{j} | \theta_{j}) \stackrel{ind}{\sim} Pois(n_{j}\theta_{j})}{(z_j | \theta_j) ~ indep Pois(n_j\theta_j)}
  \deqn{(\theta_{j} | r, \mu_{0j}) \stackrel{ind}{\sim} \frac{1}{r}Gam(r\mu_{0j})\stackrel{ind}{\sim}Gam[\mu_{0j}, \mu_{0j} / r] }{(\theta_j | r, \mu0_j) ~ indep Gam(r\mu0_j) / r ~ indep Gam[\mu0_j, \mu0_j / r]}
  \deqn{log(\mu_{0j}) = x_{j}'\beta}{log(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For Binomial data, \code{gbp} sets a two-level Binomial-Beta multilevel model. A square bracket below indicates [mean, variance] of distribution and \eqn{y_{j} = \frac{z_{j}}{n_{j}}}{y_j = z_j / n_j}.
  \deqn{(z_{j} | \theta_{j}) \stackrel{ind}{\sim} Bin(n_{j}, \theta_{j})}{(z_j | \theta_j) ~ indep Bin(n_j, \theta_j)}
  \deqn{(\theta_{j} | r, \mu_{0j}) \stackrel{ind}{\sim} Beta(r\mu_{0j}, r(1-\mu_{0j})) \stackrel{ind}{\sim} Beta[\mu_{0j}, \mu_{0j}(1 - \mu_{0j}) / (r + 1)]}{(\theta_j | r, \mu0_j) ~ indep Beta(r\mu0_j, r(1 - \mu0_j)) ~ indep Beta[\mu0_j, \mu0_j(1 - \mu0_j) / (r + 1)]}
  \deqn{logit(\mu_{0j}) = x_{j}'\beta}{logit(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For reference, based on the above notations, the Uniform prior distribution on the second level variance component (variance of the prior distribution) is \emph{dA} for Gaussian and \eqn{d(\frac{1}{r})}{d(1 / r)} 
  (= \eqn{\frac{dr}{r^{2}}}{dr / r^2}) for Binomial and Poisson data. The second level variance component can be interpreted as variation among the first-level parameters (\eqn{\theta_{j}}{\theta_j}) or variance of ensemble information.

  Under this setting, the argument \code{y} in \code{gbp} is a (\emph{k} by 1) vector of \emph{k} groups' sample means (\eqn{y_{j}'s}{y_j's} in the description of Normal-Normal model above) for Gaussian or of each group's number of successful trials (\eqn{z_{j}'s}{z_j's}) for Binomial and Poisson data, where \emph{k} is the number of groups (or units) in a dataset.

  The argument \code{se.or.n} is a (\emph{k} by 1) vector composed of the standard errors (\eqn{V_{j}'s}{V_j's}) of all groups for Gaussian or of each group's total number of trials (\eqn{n_{j}'s}{n_j's}) for Binomial and Poisson data.

  As for two optional arguments, \code{covariates} and \code{mean.PriorDist}, there are three feasible 
  combinations of them to run \code{gbp}. The first situation is when we do not have any covariate and do not 
  know a mean of the prior distribution (\eqn{\mu_{0}}{\mu0}) a priori. In this case, assigning none of two
  optional arguments, such as "\code{gbp(z, n, model = "binomial")}", will lead to a correct model. \code{gbp} 
  will automatically fit a regression with only an intercept term to estimate a common mean of the prior
  distribution (exchangeability).

  The second situation is when we have covariate(s) and do not know a mean of the prior distribution (\eqn{\mu_{0}}{\mu0}) a priori. In this case, assigning a matrix, \emph{X}, each column of which corresponds to one covariate, such as "\code{gbp(z, n, X, model = "poisson")}", will lead to a correct model. Default of \code{gbp} is to fit a regression including an intercept term to estimate a mean of the prior distribution. Double exchangeability will hold in this case.

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
  \item{shrinkage}{shrinkage estimate of each group (adjusted posterior mean)}
  \item{sd.shrinkage}{posterior standard deviation of shrinkage}
  \item{post.mean}{posterior mean of each group}
  \item{post.sd}{posterior standard deviation of each group}
  \item{post.intv.low}{lower bound of 100*confidence.lvl\% posterior interval (quantile of posterior distribution)}
  \item{post.intv.upp}{upper bound of 100*confidence.lvl\% posterior interval (quantile of posterior distribution)}
  \item{model}{"gaussian" for Gaussian, "poisson" for Poisson, and "binomial" for Binomial data}
  \item{X}{a covariate vector or matrix if designated. NA if not}
  \item{beta.new}{regression coefficient estimates}
  \item{beta.var}{estimated variance matrix of regression coefficient}
  \item{intercept}{whether TRUE or FALSE}
  \item{a.new}{a posterior mode of \eqn{\alpha} defined as log(\emph{A}) for Gaussian or log(\eqn{\frac{1}{r}}{1 / r}) for Binomial and Poisson data. Practical meaning (variation of ensemble information) of estimating \eqn{\alpha} will appear in \code{summary(gbp.object)}.}
  \item{a.var}{posterior variance of \eqn{\alpha}}
  \item{confidence.lvl}{confidence level based on which confidence interval is constructed }
  \item{weight}{(only for Binomial model) weights for acceptance-rejection method}
  \item{p}{(only for Binomial model) posterior samples of p based on the acceptance-rejection method, if this method was used. This is a (k by nsim) matrix. Each row contains posteiror samples of each random effect.}
  \item{alpha}{(only for Binomial model) posterior samples of alpha based on the acceptance-rejection method, if this method was used}
  \item{beta}{(only for Binomial model) posterior samples of beta based on the acceptance-rejection method, if this method was used}
  \item{accept.rate}{(only for Binomial model) the acceptance rate of the acceptance-rejection method, if this method was used}
  \item{n.AR}{(Only for Binomial model) If \code{n.AR = 1000}, all the results will be based on 1000 posterior samples using the acceptance-rejection sampling. Default is 0.
  }
  \item{n.AR.factor}{(only for Binomial model) If \code{n.AR = 1000} and \code{n.AR.factor = 4}, then the acceptance-rejection method will sample 4 * 1000 trial posteior samples, and accept or reject them. Default is 4.}
}

\references{
  Tak, H., Kelly, J., and Morris, C. (2017) Rgbp: An R Package for Gaussian, Poisson, and Binomial Random Effects Models with Frequency Coverage Evaluations. \emph{Journal of Statistical Software}. \bold{78}, 5, 1--33.

  Morris, C. and Lysy, M. (2012). Shrinkage Estimation in Multilevel Normal Models. 
  \emph{Statistical Science}. \bold{27}, 1, 115--134.
}

\author{
  Hyungsuk Tak, Joseph Kelly, and Carl Morris
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

  # One covariate: 1 if a player is an outfielder and 0 otherwise
  x1 <- c(1,  1,  1,  1,  1,  0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0)

  ################################################################
  # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
  ################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    g <- gbp(y, se, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A
    ### and a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    gcv <- coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  
 
    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A
    ### and regression coefficients, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    gcv <- coverage(g, A.or.r = 9, reg.coef = c(10, 1), nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A and
    ### 2nd level mean as true ones, not using estimated values as true ones.
    coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  ###############################################################
  # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
  ###############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r
    ### and a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r
    ### and regression coefficients, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    bcv <- coverage(b, A.or.r = 60, reg.coef = c(-1, 0), nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and
    ### 2nd level mean as true ones, not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  ##############################################################
  # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
  ##############################################################

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "poisson")
    p
    print(p, sort = FALSE)
    summary(p)
    plot(p)
    plot(p, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and
    ### 2nd level mean as true ones, not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

}


\keyword{methods}
