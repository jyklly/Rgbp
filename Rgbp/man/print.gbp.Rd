\name{print.gbp}
\alias{print.gbp}
\title{
Print object defined as "gbp" class
}
\description{
This function returns the resulting object of "gbp" function
}
\usage{
  \method{print}{gbp}(x, ...)
}
\arguments{
  \item{x}{
    The resultant object of "gbp" function
  }
  \item{...}{
    further arguments passed to other methods.
  }
}
\value{
  \item{sample.mean}{sample mean}
  \item{se}{standard error if model="gr"}
  \item{n}{number of total trials if model="br" or "pr"}
  \item{X}{if we put covariates}
  \item{prior.mean}{prior mean}
  \item{shrinkage}{shrinkage estimates}
  \item{se.shrinkage}{standard error of shrinkage estimate}
  \item{post.intv.low}{lower bound for 95\% posterior intervals}
  \item{post.mean}{posterior mean}
  \item{post.intv.upp}{upper bound for 95\% posterior intervals}
  \item{post.se}{square root of posterior variance}
}
\author{
  Hyungsuk Tak
}

\examples{
  data(schools)
  y <- schools$y
  se <- schools$se
  g <- gbp(y, se, mean.PriorDist = 8)
  print(g)
}
\keyword{methods}