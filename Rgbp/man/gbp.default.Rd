\name{gbp.default}
\alias{gbp}
\alias{gbp.default}
\title{
Empirical Bayesian Hierarchical modeling
}
\description{
Empirical Bayesian Hierarchical modeling for Gaussian, Binomial, and Poisson Data
}
\usage{
\method{gbp}{default}(x, y, X, mu, model="gr", Alpha=0.95, intercept=T, ...)
}
\arguments{
  \item{x}{
	The first argument of gbp is a vector of the sample mean for GRIMM and of the number of successful trials for BRIMM and PRIMM.
}
  \item{y}{
	The second argument of gbp is a vector of the standard error for GRIMM and of the total number of trials for BRIMM and PRIMM.
}
  \item{X}{
	Any covariate, if any. The intercept part will be automatically generated. The right format for one covariate is a vector and for the two covariate is k by 2 matrix.
}
  \item{mu}{
	A prior mean, the second-level mean parameter, if known.
}
  \item{model}{
	"gr" for GRIMM, "br" for BRIMM, and "pr" for PRIMM. Default is "gr"
}
  \item{Alpha}{
	This function will estimate the lower and upper bounds of 100*Alpha\% posterior interval. Default is 0.95
}
  \item{intercept}{
	User can determine whether the regression part includes the intercept term or not. Default is T.
}
  \item{...}{
	Any argument for future purpose
}
}
\value{
	An object of class \code{gbp} comprises of:
	\item{sample.mean}{sample mean of each group}
	\item{se}{if GRIMM, standard error of each group}
	\item{n}{if BRIMM and PRIMM, total number of trials of each group}
	\item{prior.mean}{numeric if entered, NA if not entered}
	\item{shrinkage}{shrinkage estimate of each group}
	\item{se.shrinkage}{standard error of shrinkage estimate}
	\item{post.intv.low}{lower bound of 100*Alpha\% posterior interval}
	\item{post.mean}{posterior mean of each group}
	\item{post.intv.upp}{upper bound of 100*Alpha\% posterior interval}
	\item{post.sd}{posterior standard deviation of each group}
	\item{model}{"gr" for GRIMM, "br" for BRIMM, and "pr" for PRIMM}
	\item{x}{covariate vector or matrix if any. NA if not}
	\item{beta.new}{regression coefficient estimates}
	\item{beta.hess}{hessian matrix of regression coefficient estimates}
	\item{intercept}{whether T or F}
	\item{a.new}{alpha estimate}
	\item{a.var}{variance of alpha estimate}
}
\author{
Joseph Kelly and Hyung Suk Tak
}
\keyword{methods}
