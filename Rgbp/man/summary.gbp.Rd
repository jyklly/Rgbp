\name{summary.gbp}
\alias{summary.gbp}
\title{Summarizing Estimation Result}
\description{\code{summary.gbp} prepares a summary of estimation result saved in the object defined as "gbp" class creating "summary.gbp" class}
\usage{\method{summary}{gbp}(object, ...)}
\arguments{
  \item{object}{
    a resultant object of \code{gbp} function. 
  }
  \item{...}{
    further arguments passed to other methods.
  }
}
\value{
  \code{summary.gbp} prepares below contents:
  \item{main}{
    a table to be displayed by \code{summary(gbp.object)}. \code{\link{print.summary.gbp}}.
  }
  \item{sec.var}{
    a vector containing an estimation result of the second-level variance component. \code{\link{print.summary.gbp}}.
  }
  \item{reg}{
    a vector composed of a summary of regression fit (if fitted). \cr \code{\link{print.summary.gbp}}.
  }
}

\author{Hyungsuk Tak, Joseph Kelly, and Carl Morris}

\examples{

  data(hospital)

  z <- hospital$d
  n <- hospital$n
  y <- hospital$y
  se <- hospital$se
  
  ###################################################################################
  # We do not have any covariates and do not know a mean of the prior distribution. #
  ###################################################################################

    ###############################################################
    # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
    ###############################################################

    g <- gbp(y, se, model = "gaussian")
    summary(g)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "binomial")
    summary(b)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, mean.PriorDist = 0.03, model = "poisson")
    summary(p)

}
\keyword{method}