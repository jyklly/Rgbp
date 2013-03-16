\name{summary.gbp}
\alias{summary.gbp}

\title{
  Summary of Estimation
}
\description{
  It retuns the summary of the object defined as "gbp" class
}
\usage{
  \method{summary}{gbp}(object, ...)
}
\arguments{
  \item{object}{
    The resultant object of "gbp" function 
  }
  \item{...}{
    further arguments passed to other methods.
  }
}
\value{
  \item{main.summary}{summary of result of the print.gbp function}
  \item{second.level.variance.summary}{summary of the second level variance component}
  \item{regression.summary}{summary of the regression, if any}
}
\author{
  Hyungsuk Tak
}
\examples{
  data(schools)
  y <- schools$y
  se <- schools$se
  g <- gbp(y, se, mean.PriorDist = 8)
  summary(g)
}
\keyword{method}