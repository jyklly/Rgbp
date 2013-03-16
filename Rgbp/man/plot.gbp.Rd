\name{plot.gbp}
\alias{plot.gbp}
\title{
Shrinkage and Posterior Interval Plots
}
\description{
This function plots two graphs, Shrinkage and Posterior Interval Plots, based on the object defined as \code{gbp} class.
}
\usage{
\method{plot}{gbp}(object, ...)
}
\arguments{
  \item{object}{
    The result of gbp function, which has \code{gbp} class
  }
  \item{...}{
    further arguments passed to other methods.
  }
}

\value{
  Two plots containing shrinkage and posterior interval plots will be displayed.
}

\author{
Hyungsuk Tak
}

\examples{
  data(schools)
  y <- schools$y
  se <- schools$se
  g <- gbp(y, se, mu0 = 8)
  plot(g)
}
\keyword{methods}
