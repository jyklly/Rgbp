\name{coverage.plot}
\alias{coverage.plot}
\title{Drawing the coverage plot}
\description{In a case where users closed the default coverage plot that the \code{coverage} function generated, the function \code{coverage.plot} redraws the coverage plot using the coverage object.}
\usage{coverage.plot(cov)}

\arguments{
  \item{cov}{
    Any saved result of the \code{coverage} function. 
  }
}

\details{
  It is possible that a user want to redraw the coverage plot for any reasons. If the result of the \code{coverage} function was saved into a variable, this \code{coverage.plot} redraw the coverage plot using the saved result.
}

\value{The coverage plot will be displayed again.}

\examples{

  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)

  b <- gbp(z, n, model = "binomial")
  cov <- coverage(b, nsim = 10)  
  coverage.plot(cov)

}

\author{Hyungsuk Tak, Joseph Kelly, and Carl Morris}

\keyword{methods}