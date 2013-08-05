\name{plot.gbp}
\alias{plot.gbp}
\title{Drawing Shrinkage and Posterior Interval Plots}
\description{
\code{plot(gbp.object)} draws shrinkage and posterior interval plots
}
\usage{\method{plot}{gbp}(x, sort = TRUE, ...)}
\arguments{
  \item{x}{
    a resultant object of \code{gbp} function. 
  }
  \item{sort}{
    \code{TRUE} or \code{FALSE} flag. If \code{TRUE}, the interval plot (second plot) will be drawn by the order of \code{se} for Gaussian, or of \code{n} for Binomial and Poisson data. If \code{FALSE}, it will be by the order of data input. Default is \code{TRUE}.
  }
  \item{...}{
    further arguments passed to other methods.
  }
}

\details{
  As for the argument \code{x}, if the result of \code{gbp} is designated to 
  \code{b} like "\code{b <- gbp(z, n, model = "binomial")}", the argument \code{x} is supposed to be \code{b}.

  This function produces two plots containing information about the prior, sample, and posterior means.

  The first plot is a shrinkage plot representing sample means (black circle) on the
  upper horizontal line and prior (blue line) and posterior means (red
  circle) on the lower horizontal line. The aim of this plot is to get a sense of the
  magnitude of the shrinkage and to observe if any change in ordering of the
  groups has occurred. Crossovers (changes of order) are noted by a
  black square as indicated in the legend. If the points plotted have
  the same value then a sunflower plot is produced where each petal
  (line protruding from the point) represent the count of points with
  that value. The plot also aims to incorporate the uncertainty and the lengths of the violet and green lines are
  proportional to the standard error and the posterior standard
  deviation respectively.

  The final plot shows interval estimates of all the groups (units) in a dataset. Two short horizontal ticks at both ends of each black vertical line indicate 97.5\% and 2.5\% quantiles of a posterior distribution for each group (Normal for Gaussian, Beta for Binomial, and Gamma for Poisson). Red dots (posterior mean) are between black circles (sample mean) and blue line(s) (prior mean) as a result of shrinkage (regression toward the mean).

  If we want to see the interval plot (the second plot) NOT sorted by the order of \code{se} for Gaussian, or of \code{n} for Binomial and Poisson data, \code{plot(b, sort = FALSE)} will show this plot by the order of data input.
}

\value{Two plots described in \emph{details} will be displayed.}

\author{Joseph Kelly, Hyungsuk Tak, and Carl Morris}

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
    plot(g)
    plot(g, sort = FALSE)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "binomial")
    plot(b)
    plot(b, sort = FALSE)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, model = "poisson")
    plot(p)
    plot(p, sort = FALSE)

}

\keyword{methods}
