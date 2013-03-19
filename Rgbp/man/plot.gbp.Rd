\name{plot.gbp}
\alias{plot.gbp}
\title{Drawing Shrinkage and Posterior Interval Plots}
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

  This function produces two plots containing information about the
  prior, sample and posterior means.

  The first plot is a shrinkage plot representing sample means on the
  upper horizontal line and prior and posterior means on the lower
  horizontal line. The aim of this plot is to get a sense of the
  magnitude of the shrinkage and to observe if any change in ordering of the
  groups has occurred. Cross-overs (changes of order) are noted by a red
  hollow circle as indicated in the legend. 


  In the legend, a black circle
  represents sample mean, a red represents posterior mean and a blue
  line the prior  mean. The lengths of the violet and green are
  proportional to the standard error and the posterior standard
  deviation respectively.

  The second column is the shrinkage plot represented by two horizontal lines; the upper line being the observed sample means and the lower line the posterior means. In this plot, the length of violet line is relative length to max(se) (i.e. se of each group / max(se) for GRIMM or n of each group / max(n) for BRIMM) and the length of darkgreen line is relative length to max(post.sd) (i.e. posterior sd of each group / max(posterior sd)). 

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
