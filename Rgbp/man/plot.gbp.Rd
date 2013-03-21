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
a}

\value{Two plots described in \emph{details} will be displayed.}

\author{Joseph Kelly, Carl Morris, and Hyungsuk Tak}

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

    g <- gbp(y, se)
    plot(g)
    plot(g, sort = FALSE)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "br")
    plot(b)
    plot(b, sort = FALSE)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, model = "pr")
    plot(p)
    plot(p, sort = FALSE)

}

\keyword{methods}
