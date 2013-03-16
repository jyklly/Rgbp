\name{print.gbp}
\alias{print.gbp}
\title{Displying "gbp" Class}
\description{
\code{print.gbp} enables users to see a compact group-level (unit-level) estimation result of \code{gbp} function.
}
\usage{\method{print}{gbp}(x, ...)}
\arguments{
  \item{x}{
    a resultant object of \code{gbp} function. 
  }
  \item{...}{
    further arguments passed to other methods.
  }
}
\details{
  As for the argument \code{x}, if the result of \code{gbp} is designated to 
  \code{b} like "\code{b <- gbp(z, n, model = "br")}", the argument \code{x} is supposed to be \code{b}.

  Users do not need to type \code{print(b)} but \code{b} itself is enough to call \code{print.gbp}.
}
\value{
  \code{print(gbp.object)} will display
  \item{sample.mean}{sample mean of each group}
  \item{se}{if GRIMM, standard error of each group}
  \item{n}{if BRIMM and PRIMM, total number of trials of each group}
  \item{X}{a covariate vector or matrix if designated. NA if not}
  \item{prior.mean}{numeric if entered, NA if not entered}
  \item{prior.mean.hat}{estimate of prior mean by a regression if prior mean is not assigned a priori}
  \item{shrinkage}{shrinkage estimate of each group}
  \item{sd.shrinkage}{standard deviation of shrinkage estimate}
  \item{post.intv.low}{lower bound of 100*Alpha\% posterior interval}
  \item{post.mean}{posterior mean of each group}
  \item{post.intv.upp}{upper bound of 100*Alpha\% posterior interval}
  \item{post.sd}{posterior standard deviation of each group}
}
\author{
  Joseph Kelly, Carl Morris, and Hyungsuk Tak
}

\examples{

  data(hospital)

  z <- hospital$d
  n <- hospital$n
  y <- hospital$y
  se <- hospital$se
  
  ###################################################################################
  # We do not have any covariates and do not know a mean of the prior distribution. #
  ###################################################################################

  #########
  # GRIMM #
  ######### 

  g <- gbp(y, se)
  g

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, model = "br")
  b

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, model = "pr")
  p
}
\keyword{methods}