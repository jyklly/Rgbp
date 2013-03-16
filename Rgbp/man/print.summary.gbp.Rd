\name{print.summary.gbp}
\alias{print.summary.gbp}
\title{The Way to Disply "summary.gbp" Class}
\description{
\code{summary(gbp.object)} enables users to see a compact summary of estimation result.
}
\usage{\method{print}{summary.gbp}(x, ...)}
\arguments{
  \item{x}{
    a resultant object of \code{gbp} function. 
  }
  \item{...}{
    further arguments passed to other methods.
  }
}

\details{
  The summary has three parts depending on the model fitted by \code{gbp}; \code{Main Summary}, \code{Second-level Variance Component Estimation Summary}, and \code{Regression Summary} (if fitted).

  A display of \code{Main Summary} changes depending first on whether all the groups (units) has the same standard error for GRIMM (the same total number of trials for BRIMM and PRIMM). If they are not the same, 
\code{Main Summary} lists groups (units) with minimun, median, and maximum values of the standard error for GRIMM (of the total number of trials for BRIMM and PRIMM). And the last line is about the overall average for all the groups (unites) within each column. Note it is not the average over groups (units) listed above.

  If groups (units) have the same standard error for GRIMM (the same total number of trials for BRIMM and PRIMM), \code{Main Summary} lists groups (units) with minimun, median, and maximum values of the sample mean. And the last row shows the overall average for all the groups (units) within each column.

  For reference, if there are several units with median values, they will show up with numbering.

  The second part is about the \code{Second-level Variance Component Estimation Summary}. To be specific, it shows estimate of \eqn{\alpha} defined as log(\emph{A}) for GRIMM and -log(\emph{r}) for BRIMM and PRIMM and its standard deviation. It is actually a posterior mode.

  The last part depends on whether \code{gbp} fitted a regression or not. For reference, \code{gbp} fits a regression if the second-level mean (\code{mean.PriorDist}) was not designated. In case, \code{gbp.object} includes a regression result, \code{summary(gbp.object)} will display the result of resgression fit.

}
\value{
  \code{summary(gbp.object)} shows compact summary of estimation result such as:
  \item{Main summary}{
    \describe{
      \item{Group w/ min(se or n)}{an estimation result of a group (unit) with the minimum standard error for GRIMM or the minimum total number of trials for BRIMM and PRIMM.}
      \item{Group w/ min(sample.mean)}{appears instead of \code{Group w/ min(se or n)} when all the groups (units) have the same standard error for GRIMM or the same total number of trials for BRIMM and PRIMM.}
      \item{Group w/ median(se or n)}{an estimation result of group(s) (unit(s)) with the median standard error for GRIMM or the median total number of trials for BRIMM and PRIMM.}
      \item{Group w/ median(sample.mean)}{appears instead of \code{Group w/ median(se or n)} when all the groups (units) have the same standard error for GRIMM or the same total number of trials for BRIMM and PRIMM.}
      \item{Group w/ max(se or n)}{an estimation result of a group (unit) with the meximum standard error for GRIMM or the maximum total number of trials for BRIMM and PRIMM.}
      \item{Group w/ max(sample.mean)}{appears instead of \code{Group w/ max(se or n)} when all the groups (units) have the same standard error for GRIMM or the same total number of trials for BRIMM and PRIMM.}
      \item{Mean over all groups}{the overall average for all the groups (units) within each column.}
    }
  }
  \item{Second-level Variance Component Estimation Summary}{
    \describe{
      \item{alpha.hat}{a posterior mode of \eqn{\alpha} defined as log(\emph{A}) for GRIMM and -log(\emph{r}) for BRIMM and PRIMM.}
      \item{alpha.hat.sd}{standard deviation of alpha.hat.}
    }
  }
  \item{Regression Summary (if fitted)}{
    \describe{
      \item{estimate}{regression coefficient estimates.}
      \item{se}{estimated standard error of regression coefficients.}
      \item{z.val}{estimate / se.}
      \item{p.val}{two-sided p-values.}
    }
  }
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

\keyword{methods}
