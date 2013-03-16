\name{plot.gbp}
\alias{plot.gbp}
\title{Shrinkage and Posterior Interval Plots}
\description{
\code{plot(gbp.object)} draws shrinkage and posterior interval plots
}
\usage{\method{plot}{gbp}(x, ...)}
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

  The overall window popping up as a result of \code{plot(b)} has three parts. The first part (column) of this window is for a legend describing symbols used in the plots. In the legend, a black circle represents sample mean, a red dot does posterior mean, a blue line does prior mean, a violet line (additional explanation needed).

  The second column is about the shrinkage plot and  it has two horizontal lines; the observed sample means are on the upper line and the posterior means are on the lower line. (additional explanation is expected of Joey)

  The final plot shows interval estimates of all the groups (units) in a dataset. Two short horizontal ticks at both ends of each black vertical line indicate 97.5\% and 2.5\% quantiles of a posterior distribution for each group (Normal for GRIMM, Gamma for PRIMM, and Beta for BRIMM). Red dots (posterior mean) are between black circles (sample mean) and blue line(s) (prior mean) as a result of shrinkage (regression toward the mean).


}
\value{Two plots described in \emph{details} will be displayed.}

\author{Joseph Kelly, Carl Morris, and Hyungsuk Tak}

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
  plot(g)

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, model = "br")
  plot(b)

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, model = "pr")
  plot(p)

}

\keyword{methods}
