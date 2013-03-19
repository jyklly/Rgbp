\name{coverage}
\alias{coverage}
\title{Estimating Coverage Probability}
\description{\code{coverage} estimates Rao-Blackwellized and unbiased coverage probabilities.}
\usage{coverage(gbp.object, A.or.r, reg.coef, covariates, mean.PriorDist, nsim = 10)}

\arguments{
  \item{gbp.object}{
    a resultant object of \code{gbp} function. 
  }
  \item{A.or.r}{
    (optional) a numeric value of \emph{A} for GRIMM or of \emph{r} for BRIMM (and PRIMM).\cr
    Designating this argument should come with other arguments, for example, \cr
    (\code{A.or.r}, \code{reg.coef}, \code{covariates} (if any)) or (\code{A.or.r}, \code{mean.PriorDist}).
  }
  \item{reg.coef}{
    (optional) a (\emph{m} by 1) vector for regression coefficients, \eqn{\beta}, where\cr
    \emph{m} is the number of regression coefficients including an intercept.
  }
  \item{covariates}{
    (optional) a (\emph{k} by \emph{t}) matrix of covariates without a column of ones for an intercept, where \cr
    \emph{k} is the number of groups (or units) in a dataset and \cr
    \emph{t} is the number of covariates (t \eqn{\ge} 1).\cr
    If \code{gbp} fits an intercept in the regression, \emph{t = (m - 1)}.
  }
  \item{mean.PriorDist}{
    (optional) a numeric value for the mean of (second-level) prior distribution.
  }
  \item{nsim}{
    number of simulations (datasets to be generated). Default is 10.
  }
}

\details{
  As for the argument \code{gbp.object}, if the result of \code{gbp} is designated to 
  \code{b}, for example \cr "\code{b <- gbp(z, n, model = "br")}", the argument \code{gbp.object} indicates this \code{b}.

  Data generating process is based on a second-level hierarchical model. The first-level is 
  a distribution of observed data (Likelihood) and the second-level is a conjugate prior distribution 
  on the first-level parameter. Covariates appear in the second-level because covariates are obtainable
  before we observe data.
 
  To be specific, for Normal data, \code{gbp} constructs a two-level Normal-Normal multilevel model 
  (\eqn{\sigma^{2}_{j}}{\sigma_j^2} is assumed to be known and subscript \emph{j} indicates \emph{j}-th group 
  (or unit) in a dataset):
  \deqn{(y_{j}~ |~ \theta_{j})~ \sim ~ indep~ N(\theta_{j}, \sigma^{2}_{j})}{(y_j | \theta_j) ~ indep N(\theta_j, \sigma_j^2)}
  \deqn{(\theta_{j}~ |~ A,~ \mu_{0j})~ \sim ~  indep~ N(\mu_{0j},~ A)}{(\theta_j | A, \mu0_j) ~ indep N(\mu0_j, A)}
  \deqn{\mu_{0j}~ =~ x_{j}'\beta}{\mu0_j = x_j'\beta}
  for \eqn{j = 1, \ldots, k}.

  For Poisson data, \code{gbp} builds a two-level Poisson-Gamma multilevel model 
  (a square bracket below indicates [mean, variance] of distribution):
  \deqn{(z_{j}~ |~ \theta_{j})~ \sim ~ indep~ Pois(n_{j}\theta_{j})}{(z_j | \theta_j) ~ indep Pois(n_j\theta_j)}
  \deqn{(\theta_{j}~ |~ r,~ \mu_{0j})~ \sim ~ indep~ Gam(r\mu_{0j}, scale = 1 / r)~ \sim ~ indep~ Gam[\mu_{0j}, \mu_{0j} / r]}{(\theta_j | r, \mu0_j) ~ indep Gam(r\mu0_j, scale = 1 / r) ~ indep Gam[\mu0_j, \mu0_j / r]}
  \deqn{log(\mu_{0j})~ =~ x_{j}'\beta}{log(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}.

  For Binomial data, \code{gbp} sets a two-level Binomial-Beta multilevel model:
  \deqn{(z_{j}~ |~ \theta_{j})~ \sim ~ indep~ Bin(n_{j}, \theta_{j})}{(z_j | \theta_j) ~ indep Bin(n_j, \theta_j)}
  \deqn{(\theta_{j}~ |~ r, \mu_{0j})~ \sim ~ indep~ Beta(r\mu_{0j},~ r(1-\mu_{0j}))~ \sim ~ indep~ Beta[\mu_{0j},~ \mu_{0j}(1 - \mu_{0j})~ /~ (r + 1)]}{(\theta_j | r, \mu0_j) ~ indep Beta(r\mu0_j, r(1 - \mu0_j)) ~ indep Beta[\mu0_j, \mu0_j(1 - \mu0_j) / (r + 1)]}
  \deqn{logit(\mu_{0j})~ =~ x_{j}'\beta}{logit(\mu0_j) = x_j'\beta}
  for \eqn{j = 1, \ldots, k}.

  From now on, the subscript \emph{i} means \emph{i}-th simulation and \emph{j} indicates \emph{j}-th group 
  (or unit).

  Pseudo-data generating process starts from the second-level to the first-level hiearchy. \code{coverage} first generates true parameters (\eqn{\theta_{i}}{\theta_i}) for \emph{k} groups (units) at the second-level and then moves onto the first-level to simulate pseudo-data sets, \eqn{y_{i}}{y_i} for GRIMM or \eqn{z_{i}}{z_i} for BRIMM and PRIMM, given previously generated true parametes (\eqn{\theta_{i}}{\theta_i}). 

  So, in order to generate pseudo-datasets, \code{coverage} needs parameters of prior distribution,  
  \emph{A} (or \emph{r}), \eqn{\beta} (\code{reg.coef}), and \emph{X} (\code{covariates}) (if any), 
  or \emph{A} (or \emph{r}) and \eqn{\mu_{0}}{\mu0} (\code{mean.PriorDist}) to start the pseudo-data generating process described above. From here, we have four options to run \code{gbp}.

  First, if any values related to the prior distribution are not designated like 
  \code{coverage(b, nsim = 10)}, then \code{coverage} will regard estimated values in \code{b} 
  (\code{gbp.object}) as given true values when it generates bunch of pseudo datasets. After sampling a (\emph{k} by 1) vector 
  \eqn{\theta_{i}}{\theta_i} from the prior
  distribution determined by those estimated values in \code{b} (\code{gbp.object}), \code{coverage}
  creates an \emph{i}-th pseudo-dataset based on \eqn{\theta_{i}}{\theta_i} just sampled.

  Second, \code{coverage} allows us to designate different true values in generating datasets, for example\cr
  \code{coverage(b, A.or.r = 15, reg.coef = 3, nsim = 100)} assuming we do not have any covariates and 
  do not know a mean of the prior distribution. One value designated in \code{reg.coef} will be used to calculate the mean of second-level distribution by \eqn{g(\mu_{0}) = \beta_{0} = 3}{g(\mu0) = \beta0 = 3}. Then, \code{coverage} samples a (\emph{k} by 1) vector
  \eqn{\theta_{i}}{\theta_i} from the prior distribution determined by designated values, \code{A.or.r} 
  and \code{reg.coef} (only intercept term). Sampling \emph{i}-th pseudo-data is based on 
  \eqn{\theta_{i}}{\theta_i} just sampled.

  Third, \code{coverage} enables us to designate different true values in generating datasets like\cr
  \code{coverage(b, A.or.r = 15, reg.coef = c(3, -1), covariates = X, nsim = 100)} when we have one covariate
  (can be more than one but \code{reg.coef} should reflect on the number of regression coefficients including
  an intercept term) without any knowledge about the mean of prior distribution a priori. For reference, a covariate matrix, \code{X}, should not include a column of ones for an intercept term in the regression and the mean of prior distribution will be set as 
  \eqn{g(\mu_{0j}) = x_{j}'\beta}{g(\mu0_j) = x_j'\beta}, where \eqn{x_{j}}{x_j} is (1, \emph{j}-th row of \code{X})'. Then, 
  \code{coverage} samples a (\emph{k} by 1) vector \eqn{\theta_{i}}{\theta_i} from the prior distribution 
  determined by designated values, \code{A.or.r}, \code{reg.coef}, and \code{covariates}. 
  Sampling \emph{i}-th pseudo-data is based on \eqn{\theta_{i}}{\theta_i} just sampled.

  Lastly, \code{coverage} provides us a way to designate different true values in generating datasets like
  \code{coverage(b, A.or.r = 15, mean.PriorDist = 0.45, nsim = 100)} when we know the mean of prior 
  distribution a priori. Then, \code{coverage} samples a (\emph{k} by 1) vector \eqn{\theta_{i}}{\theta_i} 
  from the prior distribution determined by designated values, \code{A.or.r} and \code{mean.PriorDist}.
  The \emph{i}-th Pseudo-datasets are generated based on \eqn{\theta_{i}}{\theta_i} just sampled.

  The unbiased estimator of coverage probability in \emph{j}-th group (or unit) is a sample mean of indicators
  over all simulated datasets. The \emph{j}-th indicator in \emph{i}-th simulation is 1 if the estimated 
  interval of the \emph{j}-th group on \emph{i}-th simulated dataset contains a true parameter 
  \eqn{\theta_{ij}}{\theta_ij} that generated the observed value of the \emph{j}-th group in the 
  \emph{i}-th dataset.

  Rao-Blackwellized estimator is an expectation of the unbiased estimator described above given a sufficient
  statistic, \emph{y}.
}

\value{
  \item{coverageRB}{
    estimated Rao-Blackwellized coverage probability for each group (or unit) averaged over all simulations.
  }
  \item{coverage10}{
    estimated unbiased coverage probability for each group (or unit) averaged over all simulations.
  }
  \item{average.coverageRB}{
    average value over \code{coverageRB}.
  }
  \item{average.coverage10}{
    average value over \code{coverage10}.
  }
  \item{minimum.coverageRB}{
    minimum value among \code{coverageRB}.
  }
  \item{minimum.coverage10}{
    minimum value among \code{coverage10}.
  }
  \item{raw.resultRB}{
    all the Rao-Blackwellized coverage probabilities for every group and for every simulation.
  }
  \item{raw.result10}{
    all the unbiased coverage probabilities for every group and for every simulation.
  }
}

\examples{

  # Loading datasets
  data(schools)

  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)
  x1 <- rep(c(-1, 0, 1), 6)

  y <- schools$y
  se <- schools$se
  x2 <- rep(c(-1, 0, 1, 2), 2)
  
  ####################################################################################
  # If we do not have any covariate and do not know a mean of the prior distribution #
  ####################################################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se)

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(g, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of A and 
  ### of a regression coefficient (intercept), 
  ### not using estimated values as true ones.
  coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, model = "br")

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(b, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and 
  ### of a regression coefficient (intercept), 
  ### not using estimated values as true ones.
  coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, model = "pr")

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(p, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and 
  ### of a regression coefficient (intercept), 
  ### not using estimated values as true ones.
  coverage(p, A.or.r = 60, reg.coef = -5, nsim = 10)  

  ##################################################################################
  # If we have one covariate and do not know a mean of the prior distribution yet, #
  ##################################################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se, x2, model = "gr")

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(g, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of A,
  ### of regression coefficients, and of covariate, not using estimated values 
  ### as true ones.
  coverage(g, A.or.r = 9, reg.coef = c(10, 1), covariates = x2, nsim = 10)  
  ### two values of reg.coef are for beta0 and beta1

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, x1, model = "br")

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(b, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r,
  ### of regression coefficients, and of covariate, not using estimated values 
  ### as true ones.
  coverage(b, A.or.r = 60, reg.coef = c(-1, 0), covariates = x1, nsim = 10)  
  ### two values of reg.coef are for beta0 and beta1

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, x1, model = "pr")

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(p, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r,
  ### of regression coefficients, and of covariate, not using estimated values 
  ### as true ones.
  coverage(p, A.or.r = 60, reg.coef = c(-2, 0), covariates = x1, nsim = 10)  
  ### two values of reg.coef are for beta0 and beta1

  ################################################
  # If we know a mean of the prior distribution, #
  ################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se, mean.PriorDist = 8)

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(g, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of A and
  ### of 2nd level mean as true ones, not using estimated values as true ones.
  coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  #########
  # BRIMM #
  #########
 
  b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(b, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and
  ### of 2nd level mean as true ones, not using estimated values as true ones.
  coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  #########
  # PRIMM #
  #########

  p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")

  ### when we want to simulated psuedo datasets considering the estimated values 
  ### as true ones.
  coverage(p, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and
  ### of 2nd level mean as true ones, not using estimated values as true ones.
  coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

}

\references{
Christiansen, C. and Morris, C. (1997). Hierarchical Poisson Regression Modeling. \emph{Journal of the American Statistical Association}. \bold{92}. 618-632.
}

\author{Joseph Kelly, Carl Morris, and Hyungsuk Tak}

\keyword{methods}