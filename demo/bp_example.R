source("../R/bp.R")
load("../data/baseball.RData")
br <- bp(baseball$Hits, baseball$At.Bats, model="br")
pr1 <- bp(baseball$Hits, baseball$At.Bats, model="pr1")
