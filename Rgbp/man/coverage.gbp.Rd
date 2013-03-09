\name{coverage}
\alias{coverage}
\title{
Estimating Coverage Probability
}
\description{
This function estimates Rao-Blackwellized and unbiased coverage probabilities. Unbiased coverage probability is calculated based on counting how many times the interval covers the true parameter value among all simulations. Rao-Blackellized coverage probability is expected value of unbiased coverage probability given sufficient statistic over all simulations.
}
\usage{
coverage(x, y, beta, X, mu0, nsim = 10, ...)
}
\arguments{
  \item{x}{
The resultant object of "gbp" function.  For example, if you ran the gbp function and saved the result in "b" like "b <- gbp(n, z, model = "br")", x here means b. If you did NOT assign y (r for BR and PR, A for GR), beta (regression coefficients), X (covariate), or mu0 (2nd level prior mean) in the argument of the coverage function, this function will consider the estimated values in the object x as true parameter values to generate datasets.
}
  \item{y}{
This argument is useful and should be assigned when you want to try a different true value of r for BR (and PR) or A for GR. For reference, if you want to try a different true value instead of treating its estimated value as a true value, then arguments, (y, beta, X(if any)) or (y, mu0), should be filled in with object x. If you assign this value, the coverage function will consider the assigned value as a true value to generate datasets, ignoring estimated r or A value in the object x. 
}
  \item{beta}{
abc
}
  \item{X}{
abc
}
  \item{mu0}{
abc
}
  \item{nsim}{
abc default is 10
}
  \item{...}{
Any argument for future improvement
}
}
\value{
  \item{coverageRB}{Estimated Rao-Blackwellized coverage probability for each group}
  \item{coverage10}{Estimated Unbiased coverage probability, proportion of intervals including true values, for each group}
  \item{average.coverageRB}{Mean value of coverage}
  \item{average.coverage10}{Mean value of coverage.1or0}
  \item{minimum.coverageRB}{Minimum value of coverage}
  \item{minimum.coverage10}{Mean value of coverage.1or0}
  \item{raw.resultRB}{All the coverage probability of each group for every simulation}
  \item{raw.result10}{All the coverage probability of each group for every simulation}
}
\author{
Hyungsuk Tak
}

\keyword{methods}