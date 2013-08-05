\name{hospital}

\docType{data}
\alias{hospital}

\title{
Thirty-one Hospital Data
}

\description{
Medical profiling evaluation of 31 New York hospitals in 1992. We are to consider these as Normally-distributed indices of successful outcome rates for patients at these 31 hospitals following Coronary Artery Bypass Graft (CABG) surgeries. The indices are centered so that the New York statewide average outcome over all hospitals lies near 0. Larger estimates of \code{y} indicate hospitals that performed better for these surgeries.
}

\usage{data(hospital)}

\format{
  A dataset of 31 hospitals comprises of:
  \describe{
    \item{\code{y}}{values obtained through a variance stabilizing transformation of the unbiased death rate estimates, \code{d} / \code{n}, assuming Binomial data. Details in the reference.}
    \item{\code{se}}{approximated standard error of y.}
    \item{\code{d}}{the number of deaths within a month of CABG surgeries in each hospital}
    \item{\code{n}}{total number of patients receiving CABG surgeries (case load) in each hospital}
  }
}

\source{
Morris, C. and Lysy, M. (2012). Shrinkage Estimation in Multilevel Normal Models. \emph{Statistical Science}. \bold{27}. 115-134.
}

\examples{
  data(hospital)
}

\keyword{datasets}