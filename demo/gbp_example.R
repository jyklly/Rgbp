library(sn)
source("../R/gr.R")
source("../R/gbp.R")
source("../R/bp.R")
source("../R/print.gbp.R")
source("../R/plot.gbp.R")
source("../R/summary.gbp.R")
load("../data/baseball.RData")
load("../data/schools.RData")

z<-baseball$Hits
n<-baseball$At.Bats

b <- gbp(z, n, model="br")
p <- gbp(z, n, model="pr")
g <- gbp(schools$y, schools$se)

b
summary(b)
plot(b)

p
summary(p)
plot(p)

g
summary(g)
plot(g)

#covariates
p <- 2
X <- matrix(rnorm(length(z)*p),nrow=length(z),ncol=p)
X2 <- matrix(rnorm(length(schools$y)*p),nrow=length(schools$y),ncol=p)

b2 <- gbp(z, n, X, model="br")
p2 <- gbp(z, n, X, model="pr")
g2 <- gbp(schools$y, schools$se, X2, model="gr")

b2
summary(b2)
plot(b2)

p2
summary(p2)
plot(p2)

g2
summary(g2)
plot(g2)