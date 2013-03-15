\name{hospital}

\docType{data}
\alias{hospital}

\title{
31 Hospital Data
}

\description{
1992 medical profiling evaluation of 31 New York hospitals. We are to consider these as Normally-distributed indices of successful outcome rates for patients at these 31 hospitals following coronary artery bypass graft (CABG) surgeries. The indices are centered (y) so that the New York statewide average outcome over all hospitals lies near 0. Larger estimates (y) indicate hospitals that performed better for these surgeries.
}

\usage{data(hospital)}

\format{
  A data set of 31 hospitals
  \describe{
    \item{\code{y}}{Values obtained by a variance stabilizing transformation of the unbiased success rate estimates, d/n, assuming Binomial data. Details in the reference}
    \item{\code{sd}}{Approximated standard error of y, sqrt(var(y))}
    \item{\code{d}}{The number of deaths}
    \item{\code{n}}{The number of patients}
  }
}

\source{
  Morris, C. and Lysy, M. (2012). Shrinkage Estimation in Multilevel Normal Models. 
  \emph{Statistical Science}. \bold{27}. 115-134.
}
\examples{
data(hospital)
}

\keyword{datasets}