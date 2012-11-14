library(sn)
source("../R/gr.R")
source("../R/gbp.R")
source("../R/bp.R")
load("../data/baseball.RData")
load("../data/schools.RData")


br <- gbp(baseball$Hits, baseball$At.Bats, model="br")
pr1 <- gbp(baseball$Hits, baseball$At.Bats, model="pr1")
m1 <- gbp(schools$y,schools$se)

#covariates
p <- 1
X <- matrix(rnorm(length(baseball$hits*p)),nrow=length(baseball$hits),ncol=p)
br <- gbp(baseball$Hits, baseball$At.Bats,X=X, model="br")
pr1 <- gbp(baseball$Hits, baseball$At.Bats,X=X, model="pr1")
X <- matrix(rnorm(length(schools$y*p)),nrow=length(schools$y),ncol=p)
m1 <- gbp(schools$y,schools$se,X=X)
