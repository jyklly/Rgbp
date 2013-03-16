\name{baseball}

\docType{data}
\alias{baseball}

\title{
  Baseball Data
}

\description{
  Batting averages of 18 major league players through their first 45 official at bats of the 1970 season.\cr
  These batting averages were published weekly in the New York Times, and by April 26, 1970.
}

\usage{data(baseball)}

\format{
  A data set of 18 players with 10 covariates:
  \describe{
    \item{\code{FirstName}}{each player's first name}
    \item{\code{LastName}}{each player's last name}
    \item{\code{At.Bats}}{number of times batted}
    \item{\code{Hits}}{each player's number of hits among 45 at bats}
    \item{\code{BattingAverage}}{batting averages among 45 at bats}
    \item{\code{RemainingAt.Bats}}{number of times batted after 45 at bats until the end of season}
    \item{\code{RemainingAverage}}{batting averages after 45 at bats until the end of season}
    \item{\code{SeasonAt.Bats}}{number of times batted over the whole season}
    \item{\code{SeasonHits}}{each player's number of hits over the whole season}
    \item{\code{SeasonAverage}}{batting averages over the whole season}
  }
}

\source{
Efron, B. and Morris, C. (1975). Data Analysis Using Stein's Estimator and its Generalizations. \emph{Journal of the American Statistical Association}. \bold{70}. 311-319.
}

\examples{

  data(baseball)
  z <- baseball$Hits
  n <- baseball$At.Bats

  ########################################################################################
  # Since we do not have any covariate and do not know a mean of the prior distribution, #
  ########################################################################################

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, model = "br")
  b
  summary(b)
  plot(b)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(b, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  ### (intercept), not using estimated values as true ones
  coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, model = "pr")
  p
  summary(p)
  plot(p)

  ### when we want to simulated psuedo datasets considering the estimated values as true ones
  coverage(p, nsim = 10)  

  ### when we want to simulate psuedo datasets based on different values of r and of a regression coefficient
  ### (intercept), not using estimated values as true ones
  coverage(p, A.or.r = 60, reg.coef = -5, nsim = 10)  

}

\keyword{datasets}
