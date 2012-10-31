\name{pb}
\alias{pb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
main function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pb(z, n, x = NA, prior.mean = NA, model = "br", intercept = T, result.order = F, graph = T, graph.order = F, graph.legend.position = 2, grp.name = NA, eps = 0.0001, r.alpha = 0.0833)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
%%     ~~Describe \code{z} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{prior.mean}{
%%     ~~Describe \code{prior.mean} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{intercept}{
%%     ~~Describe \code{intercept} here~~
}
  \item{result.order}{
%%     ~~Describe \code{result.order} here~~
}
  \item{graph}{
%%     ~~Describe \code{graph} here~~
}
  \item{graph.order}{
%%     ~~Describe \code{graph.order} here~~
}
  \item{graph.legend.position}{
%%     ~~Describe \code{graph.legend.position} here~~
}
  \item{grp.name}{
%%     ~~Describe \code{grp.name} here~~
}
  \item{eps}{
%%     ~~Describe \code{eps} here~~
}
  \item{r.alpha}{
%%     ~~Describe \code{r.alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z, n, x = NA, prior.mean = NA, model = "br", intercept = T, 
    result.order = F, graph = T, graph.order = F, graph.legend.position = 2, 
    grp.name = NA, eps = 0.0001, r.alpha = 0.0833) 
{
    y <- sample.mean <- as.vector(z/n)
    dif <- 1
    n.iter <- 0
    a.ini <- ini.r.a(model, prior.mean, y)$a.ini
    a.res <- alpha.est.prior.known(dif, eps, a.ini, model, prior.mean, 
        n.iter)
    B.res <- shrink.est(a.res, model, grp.name, sample.mean, 
        result.order)
    r.res <- r.est(a.res, r.alpha)
    post.res <- switch(model, br = br.post.est(B.res, prior.mean, 
        y, grp.name, result.order), pr1 = pr1.post.est(B.res, 
        prior.mean, y, grp.name, result.order))
    if (graph == T) {
        post.graph(post.res, graph.legend.position, graph.order, 
            prior.mean)
    }
    options(scipen = 3)
    output <- list(post.est = post.res$post.table, shrinkage = B.res$B.table, 
        r.hat = r.res$r.table, n.iterations = a.res$n.iter)
    return(output)
    options(digits = 7)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
