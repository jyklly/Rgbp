library(sn)
source("../R/gr.R")
source("../R/gbp.R")
source("../R/bp.R")
load("../data/baseball.RData")
load("../data/schools.RData")


br <- bp(baseball$Hits, baseball$At.Bats, model="br")
pr1 <- bp(baseball$Hits, baseball$At.Bats, model="pr1")

source("../R/gbp_shrinkage_graph.R")
m1 <- gr(schools$y,schools$se,mu=8)
