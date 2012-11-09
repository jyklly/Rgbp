#library(gbp)
library(sn)
source("../R/gr.R")
source("../R/gbp_shrinkage_graph.R")
load("../data/schools.RData")
m1 <- gr(schools$y,schools$se,mu=8)
m1
m1 <- gr(schools$y,schools$se)
m1
gbp.shrinkage.graph.R(m1$y,m1$se,8,m1$theta,m1$shat)
