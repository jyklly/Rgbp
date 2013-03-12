\name{plot.gbp}
\alias{plot.gbp}
\title{
Shrinkage and Posterior Interval Plots
}
\description{
This function plots two graphs, Shrinkage and Posterior Interval Plots, based on the object defined as 'gbp' class.
}
\usage{
\method{plot}{gbp}(x, ...)
}
\arguments{
  \item{x}{
    The result of gbp function, which has 'gbp' class
  }
  \item{...}{
    Optional arguments that may be included in future versions.
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
