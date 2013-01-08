\name{schools}
\alias{schools}
\docType{data}

\title{
Eight Schools Data
}

\description{
The batting averages of 18 major league players through their first 45 official at bats of the 1970 season. Batting averages are published weekly in the New York Times, and by April 26, 1970 Clemente had batted 45 times.
}

\usage{data(schools)}

\format{
  A data set of 18 players with 10 covariates.
  \describe{
    \item{\code{y}}{Each player's first name}
    \item{\code{se}}{Each player's last name}
  }
}

\source{
Efron, B. and Morris, C. (1975). \emph{Data Analysis Using Stein's Estimator and its Generalizations}. JASA. 70. 311-319.
}
\references{
Efron, B. and Morris, C. (1975). \emph{Data Analysis Using Stein's Estimator and its Generalizations}. JASA. 70. 311-319.
}
\examples{
data(schools)
}

\keyword{datasets}