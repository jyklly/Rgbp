\name{print.gbp}
\alias{print.gbp}
\title{Displaying 'gbp' Class}
\description{
\code{print.gbp} enables users to see a compact group-level (unit-level) estimation result of \code{gbp} function.
}
\usage{\method{print}{gbp}(x, sort = TRUE, ...)}
\arguments{
  \item{x}{
    a resultant object of \code{gbp} function. 
  }
  \item{sort}{
    \code{TRUE} or \code{FALSE} flag. If \code{TRUE}, the result will appear by the order of \code{se} for Gaussian, or of \code{n} for Binomial and Poisson data. If \code{FALSE}, it will do by the order of data input. Default is \code{TRUE}.
  }
  \item{...}{
    further arguments passed to other methods.
  }
}
\details{
  As for the argument \code{x}, if the result of \code{gbp} is designated to 
  \code{b} like "\code{b <- gbp(z, n, model = "binomial")}", the argument \code{x} is supposed to be \code{b}.

  We do not need to type "\code{print(b, sort = TRUE)}" but "\code{b}" itself is enough to call \code{print(b, sort = TRUE)}. But if we want to see the result NOT sorted by the order of \code{se} for Gaussian, or of \code{n} for Binomial and Poisson data, \code{print(b, sort = FALSE)} will show the result by the order of data input.
}
\value{
  \code{print(gbp.object)} will display:
  \item{obs.mean}{sample mean of each group}
  \item{se}{if Gaussian data, standard error of each group}
  \item{n}{if Binomial or Poisson data, total number of trials of each group}
  \item{X}{a covariate vector or matrix if designated. NA if not}
  \item{prior.mean}{numeric if entered, NA if not entered}
  \item{prior.mean.hat}{estimate of prior mean by a regression if prior mean is not assigned a priori. The variable name on the display will be "prior.mean"}
  \item{prior.mean.AR}{the posterior mean(s) of the expected random effects, if the acceptance-rejection method is used for the binomial model. The variable name on the display will be "prior.mean".}
  \item{shrinkage}{shrinkage estimate of each group (adjusted posterior mean)}
  \item{shrinkage.AR}{the posterior mean of the shrinkage factor, if the acceptance-rejection method is used for the binomial model. The variable name on the display will be "shrinkage".}
  \item{low.intv}{lower bound of 100*Alpha\% posterior interval}
  \item{post.mean}{posterior mean of each group}
  \item{upp.intv}{upper bound of 100*Alpha\% posterior interval}
  \item{post.sd}{posterior standard deviation of each group}
}
\author{
  Hyungsuk Tak, Joseph Kelly, and Carl Morris
}

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
    g
    print(g, sort = FALSE)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "binomial")
    b
    print(b, sort = FALSE)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, mean.PriorDist = 0.03, model = "poisson")
    p
    print(p, sort = FALSE)

}
\keyword{methods}