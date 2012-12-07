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
   \item{Main summary}{minimum, 1st quantile, median, mean, 3rd quantile, maximum values of sample mean, se (or n), prior mean, shrinkage, standard error of shrinkage, posterior estimates(lower and upper bounds, mean and standard deviation)
}
  \item{Second-level Variance Component Estimation Summary}{we estimate alpha and transform it into r:
    \describe{
    \item{r.hat}{(or A.hat if Grimm) second-level variance component of BRIMM or PRIMM. The amount of prior information}
    \item{se.r}{(or se.A if Grimm) standard error of r.hat based on Delta method (Normal assumed)}
	\item{low.r.hat}{(or low.A.hat if Grimm) transformed lower bound of 95.835\% posterior interval of r}
	\item{upp.r.hat}{(or upp.A.hat if Grimm) transformed upper bound of 95.835\% posterior interval of r}
    }
}
  \item{Regression Summary}{if there is any covariate:
    \describe{
    \item{estimate}{coefficient estimates}
    \item{se}{standard error of estimates}
	\item{z.val}{estimate/se}
	\item{p.val}{two-sided p-values}
    }
}
}  
\author{
Hyung Suk Tak
}

\keyword{methods}
