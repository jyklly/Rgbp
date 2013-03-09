\name{print.summary.gbp}
\alias{print.summary.gbp}
\title{
Print object defined as "summary.gbp" class
}
\description{
This function returns the resulting object of "gbp.summary" function}
\usage{
\method{print}{summary.gbp}(x,...)
}
\arguments{
  \item{x}{
	The result of gbp function, which has \code{gbp.summary} class
}
  \item{...}{
	Any argument for future purpose
}
}
\value{
   \item{Main summary}{The estimation results for units with minimum, median, and maximum values of se or n (if all se's or n's are the same, then units with minimum, median, and maximum values of sample mean) and overall mean of the estimation result
}
  \item{Second-level Variance Component Estimation Summary}{we estimate alpha and transform it into r:
    \describe{
    \item{r.hat}{(or A.hat if Grimm) second-level variance component of BRIMM or PRIMM. The amount of prior information}
    \item{se.r}{(or se.A if Grimm) standard error of r.hat based on Delta method (Normal assumed)}
	\item{low.r.hat}{(or low.A.hat if Grimm) transformed lower bound of 95.835\% posterior interval of r}
	\item{upp.r.hat}{(or upp.A.hat if Grimm) transformed upper bound of 95.835\% posterior interval of r}
    }
}
  \item{Regression Summary}{if we did not assign a known prior mean (2nd level mean):
    \describe{
    \item{estimate}{coefficient estimates}
    \item{se}{standard error of estimates}
	\item{z.val}{estimate/se}
	\item{p.val}{two-sided p-values}
    }
}
}  
\author{
Hyungsuk Tak
}

\keyword{methods}
