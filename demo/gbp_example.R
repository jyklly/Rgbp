library(sn)
source("../R/gr.R")
source("../R/gbp.R")
source("../R/bp.R")
load("../data/baseball.RData")
load("../data/schools.RData")


br <- gbp(baseball$Hits, baseball$At.Bats, model="br")
pr1 <- gbp(baseball$Hits, baseball$At.Bats, model="pr")
m1 <- gbp(schools$y,schools$se)

#covariates
p <- 2
X <- matrix(rnorm(length(baseball$Hits)*p),nrow=length(baseball$Hits),ncol=p)
X
br <- gbp(baseball$Hits, baseball$At.Bats,X=X, model="br")
pr1 <- gbp(baseball$Hits, baseball$At.Bats,X=X, model="pr")
X <- matrix(rnorm(length(schools$y)*p),nrow=length(schools$y),ncol=p)
X
m1 <- gbp(schools$y,schools$se,X=X,model="gr")
