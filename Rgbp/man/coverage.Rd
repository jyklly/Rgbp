\name{coverage}
\alias{coverage}
\title{Estimating Coverage Probability}
\description{\code{coverage} estimates Rao-Blackwellized and simple unbiased coverage probabilities.}
\usage{coverage(gbp.object, A.or.r, reg.coef, mean.PriorDist, nsim = 10)}

\arguments{
  \item{gbp.object}{
    a resultant object of \code{gbp} function. 
  }
  \item{A.or.r}{
    (optional) a given true numeric value of \emph{A} for Gaussian data or of \emph{r} for Binomial and Poisson data. If not designated, the estimated value in the \code{gbp.object} object will be considered as a true value.
  }
  \item{reg.coef}{
    (optional) a given true (\emph{m} by 1) vector for regression coefficients, \eqn{\beta}, where \emph{m} is the number of regression coefficients including an intercept. If not designated, the estimated value in the \code{gbp.object} object will be considered as a true value.
  }
  \item{mean.PriorDist}{
    (optional) a given true numeric value for the mean of (second-level) prior distribution. If not designated, the previously known value in the \code{gbp.object} object will be considered as a known prior mean.
  }
  \item{nsim}{
    number of datasets to be generated. Default is 10.
  }
}

\details{
  As for the argument \code{gbp.object}, if the result of \code{gbp} is designated to 
  \code{b}, for example \cr "\code{b <- gbp(z, n, model = "binomial")}", the argument \code{gbp.object} indicates this \code{b}.

  Data generating process is based on a second-level hierarchical model. The first-level hierarchy is 
  a distribution of observed data and the second-level is a conjugate prior distribution 
  on the first-level parameter.
 
  To be specific, for Normal data, \code{gbp} constructs a two-level Normal-Normal 2-level model. \eqn{\sigma^{2}_{j}}{\sigma_j^2} below is assumed to be known or to be accurately estimated (\eqn{s^{2}_{j}}{s^2}) and subscript \emph{j} indicates \emph{j}-th group in a dataset.
  \deqn{(y_{j} | \theta_{j}) \stackrel{ind}{\sim} N(\theta_{j}, \sigma^{2}_{j})}{(y_j | \theta_j) ~ indep N(\theta_j, \sigma_j^2)}
  \deqn{(\theta_{j} |\mu_{0j} , A) \stackrel{ind}{\sim} N(\mu_{0j}, A)}{(\theta_j | \mu0_j, A) ~ indep N(\mu0_j, A)}
  \deqn{\mu_{0j} = x_{j}'\beta}{\mu0_j = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For Poisson data, \code{gbp} builds a two-level Poisson-Gamma multi-level model. A square bracket below indicates [mean, variance] of distribution and a constant multiplied to the notation representing Gamma distribution (Gam) is a scale. Also, for consistent notation, \eqn{y_{j}=\frac{z_{j}}{n_{j}}}{y_j = z_j / n_j} and \eqn{n_{j}}{n_j} can be interpreted as \emph{j}-th group's exposure only in this Poisson-Gamma hierarchical model.
  \deqn{(z_{j} | \theta_{j}) \stackrel{ind}{\sim} Pois(n_{j}\theta_{j})}{(z_j | \theta_j) ~ indep Pois(n_j\theta_j)}
  \deqn{(\theta_{j} | r, \mu_{0j}) \stackrel{ind}{\sim} \frac{1}{r}Gam(r\mu_{0j})\stackrel{ind}{\sim}\_Gam[\mu_{0j}, \mu_{0j} / r] }{(\theta_j | r, \mu0_j) ~ indep Gam(r\mu0_j) / r ~ indep Gam[\mu0_j, \mu0_j / r]}
  \deqn{log(\mu_{0j}) = x_{j}'\beta}{log(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  For Binomial data, \code{gbp} sets a two-level Binomial-Beta multi-level model. For reference, a square bracket below indicates [mean, variance] of distribution and \eqn{y_{j} = \frac{z_{j}}{n_{j}}}{y_j = z_j / n_j}.
  \deqn{(z_{j} | \theta_{j}) \stackrel{ind}{\sim} \_Bin(n_{j}, \theta_{j})}{(z_j | \theta_j) ~ indep Bin(n_j, \theta_j)}
  \deqn{(\theta_{j} | r, \mu_{0j}) \stackrel{ind}{\sim} Beta(r\mu_{0j}, r(1-\mu_{0j})) \stackrel{ind}{\sim} \_Beta[\mu_{0j}, \mu_{0j}(1 - \mu_{0j}) / (r + 1)]}{(\theta_j | r, \mu0_j) ~ indep Beta(r\mu0_j, r(1 - \mu0_j)) ~ indep Beta[\mu0_j, \mu0_j(1 - \mu0_j) / (r + 1)]}
  \deqn{logit(\mu_{0j}) = x_{j}'\beta}{logit(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}, where \emph{k} is the number of groups (units) in a dataset.

  From now on, the subscript \emph{(i)} means \emph{i}-th simulation and the subscript \emph{j} indicates \emph{j}-th group. So, notations with a subscript \emph{(i)} are (\emph{k} by 1) vectors, for example \eqn{\theta_{(i)}' = (\theta_{(i)1}, \theta_{(i)2}, \ldots, \theta_{(i)k})}{\theta_i' = (\thate_(i)1, \theta_(i)2, ..., \theta_(i)k)}.

  Pseudo-data generating process starts from the second-level hierarchy to the first-level. \code{coverage} first generates true parameters (\eqn{\theta_{(i)}}{\theta_(i)}) for \emph{k} groups at the second-level and then moves onto the first-level to simulate pseudo-data sets, \eqn{y_{(i)}}{y_(i)} for Gaussian or \eqn{z_{(i)}}{z_(i)} for Binomial and Poisson data, given previously generated true parameters (\eqn{\theta_{(i)}}{\theta_(i)}). 

  So, in order to generate pseudo-datasets, \code{coverage} needs parameters of prior distribution,  
  (\emph{A} (or \emph{r}) and \eqn{\beta} (\code{reg.coef})) 
  or (\emph{A} (or \emph{r}) and \eqn{\mu_{0}}{\mu0}). From here, we have four options to run \code{coverage}.

  First, if any values related to the prior distribution are not designated like 
  \code{coverage(b, nsim = 10)}, then \code{coverage} will regard estimated values (or known prior mean, \eqn{\mu_{0}}{\mu0}) in \code{b} (\code{gbp.object}) as given true values when it generates lots of pseudo-datasets. After sampling \eqn{\theta_{(i)}}{\theta_(i)} from the prior distribution determined by these estimated values (or known prior mean) in \code{b} (\code{gbp.object}), \code{coverage} creates an \emph{i}-th pseudo-dataset based on \eqn{\theta_{(i)}}{\theta_(i)} just sampled.

  Second, \code{coverage} allows us to try different true values in generating datasets. Suppose \code{gbp.object} is based on the model with a known prior mean, \eqn{\mu_{0}}{\mu0}. Then, we can try either different \code{A.or.r} or \code{mean.PriorDist}. For example, \code{coverage(b, A.or.r = 20, nsim = 10)}, \cr
\code{coverage(b, mean.PriorDist = 0.5, nsim = 10)}, or \cr
\code{coverage(b, A.or.r = 20, mean.PriorDist = 0.5, nsim = 10)}. Note that we cannot set \code{reg.coef} because the second-level mean (prior mean) is known in \code{gbp.object} to begin with.

  Suppose \code{gbp.object} is based on the model with an unknown prior mean. In this case, \code{gbp.object} has the estimation result of regression model, linear regression for Normal-Normal, log-linear regression for Poisson-Gamma, or logistic regression for Binomial-Beta, (only intercept term if there is no covariate) to estimate the unknown prior mean. Then, we can try some options: one or two of (\code{A.or.r}, \code{mean.PriorDist}, \code{reg.coef}). For example, \code{coverage(b, A.or.r = 20, nsim = 10)}, \code{coverage(b, mean.PriorDist = 0.5, nsim = 10)}, or \code{coverage(b, reg.coef = 0.1, nsim = 10)} with no covariate where \code{reg.coef} is a designated intercept term. Estimates in \code{gbp.object} will be used for undesignated values. Also, we can try appropriate combinations of two arguments. For example, \code{coverage(b, A.or.r = 20, mean.PriorDist = 0.5, nsim = 10)} and \cr
\code{coverage(b, A.or.r = 20, reg.coef = 0.1, nsim = 10)}. If we have one covariate, a 2 by 1 vector should be designated for \code{reg.coef}, one for an intercept term and the other for a regression coefficient of the covariate. Note that the two arguments, \code{mean.PriorDist} and \code{reg.coef}, cannot be assigned together because we do not need \code{reg.coef} given \code{mean.PriorDist}. 

  The simple unbiased estimator of coverage probability in \emph{j}-th group is a sample mean of indicators over all simulated datasets. The \emph{j}-th indicator in \emph{i}-th simulation is 1 if the estimated interval of the \emph{j}-th group on \emph{i}-th simulated dataset contains a true parameter 
  \eqn{\theta_{(i)j}}{\theta_(i)j} that generated the observed value of the \emph{j}-th group in the 
  \emph{i}-th dataset.

  Rao-Blackwellized unbiased estimator for group \emph{j} is a conditional expectation of the simple unbiased estimator given a sufficient statistic, \eqn{y_{j}}{y_j} for Gaussian or \eqn{z_{j}}{z_j} for Binomial and Poisson data.
}

\value{
  \item{coverageRB}{
    Rao-Blackwellized unbiased coverage estimate for each group averaged over all simulations.
  }
  \item{coverageS}{
    Simple unbiased coverage estimate for each group averaged over all simulations.
  }
  \item{average.coverageRB}{
    Overall Rao-Blackwellized unbiased coverage estimate across all the groups and simulations.
  }
  \item{average.coverageS}{
    Overall simple unbiased coverage estimate across all the groups and simulations.
  }
  \item{se.coverageRB}{
    Standard error of Rao-Blackwellized unbiased coverage estimate for each group.
  }
  \item{se.coverageS}{
    Standard error of simple unbiased coverage estimate for each group.
  }
  \item{raw.resultRB}{
    all the Rao-Blackwellized unbiased coverage estimates for every group and for every simulation.
  }
  \item{raw.resultS}{
    all the simple unbiased coverage estimates for every group and for every simulation.
  }
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
  
  #################################################################
  # Gaussian Regression Interactive Multi-level Modeling (GRIMM) #
  #################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    g <- gbp(y, se, model = "gaussian")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverageS, gcv$average.coverageRB, gcv$average.coverageS,
    ### gcv$minimum.coverageRB, gcv$raw.resultRB, gcv$raw.resultS

    ### gcv <- coverage(g, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, nsim = 100)
    ### gcv <- coverage(g, reg.coef = 10, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, reg.coef = 10, nsim = 100)

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gaussian")
 
    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  
 
    ### gcv$coverageRB, gcv$coverageS, gcv$average.coverageRB, gcv$average.coverageS,
    ### gcv$minimum.coverageRB, gcv$raw.resultRB, gcv$raw.resultS

    ### gcv <- coverage(g, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 200, nsim = 100)
    ### gcv <- coverage(g, reg.coef = c(10, 2), nsim = 100)
    ### gcv <- coverage(g, A.or.r = 200, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 200, reg.coef = c(10, 2), nsim = 100)

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8, model = "gaussian")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverageS, gcv$average.coverageRB, gcv$average.coverageS,
    ### gcv$minimum.coverageRB, gcv$raw.resultRB, gcv$raw.resultS

    ### gcv <- coverage(g, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, mean.PriorDist = 3, nsim = 100)

  ################################################################
  # Binomial Regression Interactive Multi-level Modeling (BRIMM) #
  ################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "binomial")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverageS, bcv$average.coverageRB, bcv$average.coverageS,
    ### bcv$minimum.coverageRB, bcv$raw.resultRB, bcv$raw.resultS

    ### bcv <- coverage(b, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, nsim = 100)
    ### bcv <- coverage(b, reg.coef = -1.5, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, reg.coef = -1.5, nsim = 100)

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "binomial")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverageS, bcv$average.coverageRB, bcv$average.coverageS,
    ### bcv$minimum.coverageRB, bcv$raw.resultRB, bcv$raw.resultS

    ### bcv <- coverage(b, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, nsim = 100)
    ### bcv <- coverage(b, reg.coef = c(-1.5, 0), nsim = 100)
    ### bcv <- coverage(b, A.or.r = 40, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 40, reg.coef = c(-1.5, 0), nsim = 100)

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "binomial")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverageS, bcv$average.coverageRB, bcv$average.coverageS,
    ### bcv$minimum.coverageRB, bcv$raw.resultRB, bcv$raw.resultS

    ### bcv <- coverage(b, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 40, mean.PriorDist = 0.2, nsim = 100)

  ###############################################################
  # Poisson Regression Interactive Multi-level Modeling (PRIMM) #
  ###############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    p <- gbp(z, n, model = "poisson")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverageS, pcv$average.coverageRB, pcv$average.coverageS,
    ### pcv$minimum.coverageRB, pcv$raw.resultRB, pcv$raw.resultS

    ### pcv <- coverage(p, mean.PriorDist = 0.265, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 120, nsim = 100)
    ### pcv <- coverage(p, reg.coef = -2, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 120, mean.PriorDist = 0.265, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 120, reg.coef = -2, nsim = 100)

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    p <- gbp(z, n, x1, model = "poisson")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverageS, pcv$average.coverageRB, pcv$average.coverageS,
    ### pcv$minimum.coverageRB, pcv$raw.resultRB, pcv$raw.resultS

    ### pcv <- coverage(p, mean.PriorDist = 0.265, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 170, nsim = 100)
    ### pcv <- coverage(p, reg.coef = c(-2, 1), nsim = 100)
    ### pcv <- coverage(p, A.or.r = 170, mean.PriorDist = 0.265, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 170, reg.coef = c(-2, 1), nsim = 100)

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "poisson")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverageS, pcv$average.coverageRB, pcv$average.coverageS,
    ### pcv$minimum.coverageRB, pcv$raw.resultRB, pcv$raw.resultS

    ### pcv <- coverage(p, mean.PriorDist = 0.265, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 150, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 150, mean.PriorDist = 0.265, nsim = 100)

}

\references{
Christiansen, C. and Morris, C. (1997). Hierarchical Poisson Regression Modeling. \emph{Journal of the American Statistical Association}. \bold{92}. 618-632.
}

\author{Joseph Kelly, Hyungsuk Tak, and Carl Morris}

\keyword{methods}