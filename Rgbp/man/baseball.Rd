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
  A data set of 18 players with 10 covariates.
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

  # number of successful trials 
  z <- baseball$Hits  

  # total number of trials
  n <- baseball$At.Bats  
  
  # running BRIMM without any covariate and saving resultant object to 'b'
  b <- gbp(z, n, model = "br")  
  
  # printing the result
  b  

  # printing summary of the result
  summary(b)

  # visualizing shrinkage and 95\% interval plots
  plot(b)

  # estimating coverage probability letting 'coverage' regard the estimates 
  # of 'r' (proportional to inverse of the variance of prior distribution) and 'beta' (regression coefficients) 
  # in 'b' as true values when 'coverage' generates simulated datasets.
  coverage(b, nsim = 10)

  # estimating coverage probability letting 'coverage' regard the assigned values for
  # 'y' (= 'r') and 'beta' in the argument of 'coverage' as true values 
  # when 'coverage' generates simulated datasets.
  coverage(b, y = 60, beta = 0, nsim = 10)   
}

\keyword{datasets}