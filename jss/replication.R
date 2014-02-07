################################################################################
###  Load Package
################################################################################
library("Rgbp")


################################################################################
###  5.1. 31 Hospitals: Known Second-level mean
################################################################################
z <- c( 3,   2,   5,  11,   9,  12,  12,   4,  10,  13,  14,   7,  12,
       11,  13,  22,  15,  11,  14,  11,  16,  14,   9,  15,  13,  35,
       26,  25,  20,  35,  27)
n <- c(67,  68, 210, 256, 269, 274, 278, 295, 347, 349, 358, 396, 431,
      441, 477, 484, 494, 501, 505, 540, 563, 593, 602, 629, 636, 729,
      849, 914, 940, 1193, 1340)
## or 
data("hospital")
z <- hospital$d # variable name for the number of deaths is d in dataset
n <- hospital$n

p <- gbp(z, n, mean.PriorDist = 0.03, model = "poisson")
p
summary(p)
plot(p)
dev.off()
pcv <- coverage(p, nsim = 1000)
dev.off()
pcv$coverageRB
pcv$coverageS
pcv$se.coverageRB
pcv$se.coverageS


################################################################################
###  5.2. 8 Schools: Unknown Second-level Mean with No Covariates
################################################################################
y <- c(12, -3, 28, 7, 1, 8, 18, -1)
se <- c(18, 16, 15, 11, 11, 10, 10, 9)
## or
data("schools")
y <- schools$y
se <- schools$se

g <- gbp(y, se, model = "gaussian")
g
summary(g)
plot(g)
dev.off()
gcv <- coverage(g, nsim = 1000)
dev.off()
gcv$coverageRB
gcv$se.coverageRB


################################################################################
###  5.3. 18 Baseball Players: Unknown Second-level Mean and One Covariate
################################################################################
z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10, 9, 8, 7)
n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)
x <- c( 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0)
##or
data("baseball")
z <- baseball$Hits
n <- baseball$At.Bats
x <- ifelse(baseball$Position == "fielder", 1, 0)

b <- gbp(z, n, x, model = "binomial")
b
summary(b)
plot(b)
dev.off()
bcv <- coverage(b, nsim = 1000)
dev.off()
bcv$coverageRB
bcv$se.coverageRB
