\name{summary.gbp}
\alias{summary.gbp}
\title{Summarizing Estimation Result from "gbp" Object}
\description{\code{summary.gbp} prepares the summary of the object defined as "gbp" class creating "summary.gbp" class}
\usage{\method{summary}{gbp}(object, ...)}
\arguments{
  \item{object}{
    a resultant object of \code{gbp} function. 
  }
  \item{...}{
    further arguments passed to other methods.
  }
}
\value{
  \code{summary.gbp} prepares below contents:
  \item{main}{
    a table to be displayed by \code{print(gbp.object)}. \code{\link{print.summary.gbp}}.
  }
  \item{sec.var}{
    second-level variance component estimation summary. \code{\link{print.summary.gbp}}.
  }
  \item{reg}{
    regression summary (if fitted). \code{\link{print.summary.gbp}}.
  }
}

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
  summary(g)

  #########
  # BRIMM #
  ######### 

  b <- gbp(z, n, model = "br")
  summary(b)

  #########
  # PRIMM #
  ######### 

  p <- gbp(z, n, model = "pr")
  summary(p)
}
\keyword{method}