\name{Rgbp}
\alias{Rgbp-package}
\alias{Rgbp}
\docType{package}
\title{Hierarchical Modeling and Frequency Method Checking on Overdispersed Gaussian, Poisson, and Binomial Data}
\description{
Bayesian-frequentist reconciliation via approximate Bayesian hierarchical modeling and frequency method checking for over-dispersed Gaussian, Binomial, and Poisson data.
}
\details{
  \tabular{ll}{
    Package: \tab Rgbp\cr
    Type: \tab Package\cr
    Version: \tab 1.0.0\cr
    Date: \tab 2013-03-16\cr
    License: \tab GPL-2\cr
    Main functions: \tab \code{\link{gbp}}, \code{\link{coverage}}\cr
  }

\code{Rgbp} is an R package that utilizes approximate Bayesian machinery to provide a method of estimating two-level hierarchical models for Gaussian, Poisson, and Binomial data in a fast and computationally efficient manner. The main products of this package are point and interval estimates for the true parameters, whose good frequency properties can be validated via its repeated sampling procedure called frequency method checking.  It is found that such Bayesian-frequentist reconciliation allows \code{Rgbp} to have attributes desirable from both perspectives, working well in small samples and yielding good coverage probabilities for its interval estimates.
}
\author{
  Hyungsuk Tak, Joseph Kelly, and Carl Morris

  Maintainer: Joseph Kelly <josephkelly@google.com>
}
\references{
  Morris, C. and Lysy, M. (2012). Shrinkage Estimation in Multilevel Normal Models. \emph{Statistical Science}. \bold{27}. 115-134.
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
    ### gcv <- coverage(g, nsim = 10)  
    ### more details in ?coverage

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
    ### gcv <- coverage(g, nsim = 10)  
    ### more details in ?coverage 

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
    ### gcv <- coverage(g, nsim = 10)  
    ### more details in ?coverage 

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
    ### bcv <- coverage(b, nsim = 10)  
    ### more details in ?coverage 

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
    ### bcv <- coverage(b, nsim = 10)  
    ### more details in ?coverage 

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
    ### bcv <- coverage(b, nsim = 10)  
    ### more details in ?coverage 

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
    ### pcv <- coverage(p, nsim = 10)  
    ### more details in ?coverage 

}

\keyword{package}