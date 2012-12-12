\name{coverage.gbp}
\alias{coverage.gbp}
\title{
Estimating coverage probability from the object defined as "gbp" class
}
\description{
This function estimates Rao-Blackwellized and unbiased coverage probability based on Gaussian Quadrature
}
\usage{
\method{coverage}{gbp}(x, nsim, detail=F, ...)
}
\arguments{
  \item{x}{
The resultant object of "gbp" function
}
  \item{nsim}{
The number of simulations we want
}
  \item{detail}{
If TRUE, this function will give us detailed results
}
  \item{...}{
Anything else
}
}
\value{
  \item{coverage}{Estimated Rao-Blackwellized coverage probability for each group}
  \item{coverage.1or0}{Estimated Unbiased coverage probability, proportion of intervals including true values, for each group}
  \item{average.coverage}{Mean value of coverage}
  \item{average.coverage.1or0}{Mean value of coverage.1or0}
  \item{minimum.coverage}{Minimum value of coverage}
  \item{minimum.coverage.1or0}{Mean value of coverage.1or0}
  \item{raw.result}{(if detail=T) All the coverage probability of each group for every simulation}
  \item{raw.result.1or0}{(if detail=T) All the coverage probability of each group for every simulation}
}
\author{
Hyung Suk Tak
}

\keyword{methods}