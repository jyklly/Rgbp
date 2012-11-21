\name{baseball}
\alias{baseball}
\docType{data}

\title{
Baseball Data
}

\description{
The batting averages of 18 major league players through their first 45 official at bats of the 1970 season. Batting averages are published weekly in the New York Times, and by April 26, 1970 Clemente had batted 45 times.
}

\usage{data(baseball)}

\format{
  A data set of 18 players with 10 covariates.
  \describe{
    \item{\code{FirstName}}{Each player's first name}
    \item{\code{LastName}}{Each player's last name}
    \item{\code{At.Bats}}{The number of times (45) each batter (player) faces a pitcher}
    \item{\code{Hits}}{Each player's number of hits among 45 trials}
    \item{\code{BattingAverage}}{Batting averages among 45 at bats}
    \item{\code{RemainingAt.Bats}}{The number of times each batter (player) faces a pitcher after 45 AB until the end of season}
    \item{\code{RemainingAverage}}{Batting averages after 45 AB until the end of season}
    \item{\code{SeasonAt.Bats}}{The number of times each batter (player) faces a pitcher during the whole season}
    \item{\code{SeasonHits}}{Each player's number of hits during the whole season}
    \item{\code{SeasonAverage}}{Batting averages during the whole season}
  }
}

\source{
Efron, B. and Morris, C. (1975). \emph{Data Analysis Using Stein's Estimator and its Generalizations}. JASA. 70. 311-319.
}
\references{
Efron, B. and Morris, C. (1975). \emph{Data Analysis Using Stein's Estimator and its Generalizations}. JASA. 70. 311-319.
}

\examples{
data(baseball)
}

\keyword{datasets}