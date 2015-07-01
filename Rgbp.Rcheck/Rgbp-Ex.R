pkgname <- "Rgbp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "Rgbp-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('Rgbp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Rgbp-package")
### * Rgbp-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Rgbp
### Title: Hierarchical Modeling and Frequency Method Checking on
###   Overdispersed Gaussian, Poisson, and Binomial Data
### Aliases: Rgbp-package Rgbp
### Keywords: package

### ** Examples


  # Loading datasets
  data(schools)
  y <- schools$y
  se <- schools$se

  # Arbitrary covariate for schools data
  x2 <- rep(c(-1, 0, 1, 2), 2)
  
  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)

  # One covariate: 1 if a player is an outfielder and 0 otherwise
  x1 <- c(1,  1,  1,  1,  1,  0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0)

  ################################################################
  # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
  ################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    g <- gbp(y, se, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    ### gcv <- coverage(g, nsim = 10)  
    ### more details in ?coverage

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    ### gcv <- coverage(g, nsim = 10)  
    ### more details in ?coverage 

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    ### gcv <- coverage(g, nsim = 10)  
    ### more details in ?coverage 

  ###############################################################
  # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
  ###############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    ### bcv <- coverage(b, nsim = 10)  
    ### more details in ?coverage 

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    ### bcv <- coverage(b, nsim = 10)  
    ### more details in ?coverage 

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    ### bcv <- coverage(b, nsim = 10)  
    ### more details in ?coverage 

  ##############################################################
  # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
  ##############################################################

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "poisson")
    p
    print(p, sort = FALSE)
    summary(p)
    plot(p)
    plot(p, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    ### pcv <- coverage(p, nsim = 10)  
    ### more details in ?coverage 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Rgbp-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("baseball")
### * baseball

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: baseball
### Title: Baseball Data
### Aliases: baseball
### Keywords: datasets

### ** Examples

  data(baseball)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("baseball", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coverage")
### * coverage

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coverage
### Title: Estimating Coverage Probability
### Aliases: coverage
### Keywords: methods

### ** Examples


  # Loading datasets
  data(schools)
  y <- schools$y
  se <- schools$se

  # Arbitrary covariate for schools data
  x2 <- rep(c(-1, 0, 1, 2), 2)

  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)

  # One covariate: 1 if a player is an outfielder and 0 otherwise
  x1 <- c(1,  1,  1,  1,  1,  0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0)
  
  #################################################################
  # Gaussian Regression Interactive Multi-level Modeling (GRIMM) #
  #################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    g <- gbp(y, se, model = "gaussian")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverageS, gcv$average.coverageRB, gcv$average.coverageS,
    ### gcv$minimum.coverageRB, gcv$raw.resultRB, gcv$raw.resultS

    ### gcv <- coverage(g, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, nsim = 100)
    ### gcv <- coverage(g, reg.coef = 10, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, reg.coef = 10, nsim = 100)

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gaussian")
 
    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  
 
    ### gcv$coverageRB, gcv$coverageS, gcv$average.coverageRB, gcv$average.coverageS,
    ### gcv$minimum.coverageRB, gcv$raw.resultRB, gcv$raw.resultS

    ### gcv <- coverage(g, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 200, nsim = 100)
    ### gcv <- coverage(g, reg.coef = c(10, 2), nsim = 100)
    ### gcv <- coverage(g, A.or.r = 200, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 200, reg.coef = c(10, 2), nsim = 100)

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8, model = "gaussian")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverageS, gcv$average.coverageRB, gcv$average.coverageS,
    ### gcv$minimum.coverageRB, gcv$raw.resultRB, gcv$raw.resultS

    ### gcv <- coverage(g, mean.PriorDist = 3, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, nsim = 100)
    ### gcv <- coverage(g, A.or.r = 150, mean.PriorDist = 3, nsim = 100)

  ################################################################
  # Binomial Regression Interactive Multi-level Modeling (BRIMM) #
  ################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "binomial")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverageS, bcv$average.coverageRB, bcv$average.coverageS,
    ### bcv$minimum.coverageRB, bcv$raw.resultRB, bcv$raw.resultS

    ### bcv <- coverage(b, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, nsim = 100)
    ### bcv <- coverage(b, reg.coef = -1.5, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, reg.coef = -1.5, nsim = 100)

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "binomial")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverageS, bcv$average.coverageRB, bcv$average.coverageS,
    ### bcv$minimum.coverageRB, bcv$raw.resultRB, bcv$raw.resultS

    ### bcv <- coverage(b, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, nsim = 100)
    ### bcv <- coverage(b, reg.coef = c(-1.5, 0), nsim = 100)
    ### bcv <- coverage(b, A.or.r = 40, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 40, reg.coef = c(-1.5, 0), nsim = 100)

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "binomial")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverageS, bcv$average.coverageRB, bcv$average.coverageS,
    ### bcv$minimum.coverageRB, bcv$raw.resultRB, bcv$raw.resultS

    ### bcv <- coverage(b, mean.PriorDist = 0.2, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 50, nsim = 100)
    ### bcv <- coverage(b, A.or.r = 40, mean.PriorDist = 0.2, nsim = 100)

  ###############################################################
  # Poisson Regression Interactive Multi-level Modeling (PRIMM) #
  ###############################################################

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "poisson")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverageS, pcv$average.coverageRB, pcv$average.coverageS,
    ### pcv$minimum.coverageRB, pcv$raw.resultRB, pcv$raw.resultS

    ### pcv <- coverage(p, mean.PriorDist = 0.265, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 150, nsim = 100)
    ### pcv <- coverage(p, A.or.r = 150, mean.PriorDist = 0.265, nsim = 100)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coverage", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("coverage.plot")
### * coverage.plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: coverage.plot
### Title: Drawing the coverage plot
### Aliases: coverage.plot
### Keywords: methods

### ** Examples


  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)

  b <- gbp(z, n, model = "binomial")
  cov <- coverage(b, nsim = 10)  
  coverage.plot(cov)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("coverage.plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gbp")
### * gbp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gbp
### Title: Fitting Gaussian, Poisson, and Binomial Hierarchical Models
### Aliases: gbp
### Keywords: methods

### ** Examples


  # Loading datasets
  data(schools)
  y <- schools$y
  se <- schools$se

  # Arbitrary covariate for schools data
  x2 <- rep(c(-1, 0, 1, 2), 2)
  
  # baseball data where z is Hits and n is AtBats
  z <- c(18, 17, 16, 15, 14, 14, 13, 12, 11, 11, 10, 10, 10, 10, 10,  9,  8,  7)
  n <- c(45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45)

  # One covariate: 1 if a player is an outfielder and 0 otherwise
  x1 <- c(1,  1,  1,  1,  1,  0,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  0,  0)

  ################################################################
  # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
  ################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    g <- gbp(y, se, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A
    ### and a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    gcv <- coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  
 
    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A
    ### and regression coefficients, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    gcv <- coverage(g, A.or.r = 9, reg.coef = c(10, 1), nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8, model = "gaussian")
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A and
    ### 2nd level mean as true ones, not using estimated values as true ones.
    coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  ###############################################################
  # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
  ###############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r
    ### and a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r
    ### and regression coefficients, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    bcv <- coverage(b, A.or.r = 60, reg.coef = c(-1, 0), nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "binomial")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and
    ### 2nd level mean as true ones, not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  ##############################################################
  # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
  ##############################################################

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "poisson")
    p
    print(p, sort = FALSE)
    summary(p)
    plot(p)
    plot(p, sort = FALSE)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and
    ### 2nd level mean as true ones, not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gbp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hospital")
### * hospital

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hospital
### Title: Thirty-one Hospital Data
### Aliases: hospital
### Keywords: datasets

### ** Examples

  data(hospital)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hospital", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.gbp")
### * plot.gbp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.gbp
### Title: Drawing Shrinkage and Posterior Interval Plots
### Aliases: plot.gbp
### Keywords: methods

### ** Examples


  data(hospital)

  z <- hospital$d
  n <- hospital$n
  y <- hospital$y
  se <- hospital$se
  
  
  ###################################################################################
  # We do not have any covariates and do not know a mean of the prior distribution. #
  ###################################################################################

    ###############################################################
    # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
    ###############################################################

    g <- gbp(y, se, model = "gaussian")
    plot(g)
    plot(g, sort = FALSE)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "binomial")
    plot(b)
    plot(b, sort = FALSE)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, mean.PriorDist = 0.03, model = "poisson")
    plot(p)
    plot(p, sort = FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.gbp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.gbp")
### * print.gbp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.gbp
### Title: Displaying 'gbp' Class
### Aliases: print.gbp
### Keywords: methods

### ** Examples


  data(hospital)

  z <- hospital$d
  n <- hospital$n
  y <- hospital$y
  se <- hospital$se
  
  ###################################################################################
  # We do not have any covariates and do not know a mean of the prior distribution. #
  ###################################################################################

    ###############################################################
    # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
    ###############################################################

    g <- gbp(y, se, model = "gaussian")
    g
    print(g, sort = FALSE)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "binomial")
    b
    print(b, sort = FALSE)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, mean.PriorDist = 0.03, model = "poisson")
    p
    print(p, sort = FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.gbp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.summary.gbp")
### * print.summary.gbp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.summary.gbp
### Title: Displaying 'summary.gbp' Class
### Aliases: print.summary.gbp
### Keywords: methods

### ** Examples


  data(hospital)

  z <- hospital$d
  n <- hospital$n
  y <- hospital$y
  se <- hospital$se
  
  ###################################################################################
  # We do not have any covariates and do not know a mean of the prior distribution. #
  ###################################################################################

    ###############################################################
    # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
    ###############################################################

    g <- gbp(y, se, model = "gaussian")
    summary(g)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "binomial")
    summary(b)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, mean.PriorDist = 0.03, model = "poisson")
    summary(p)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.summary.gbp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("schools")
### * schools

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: schools
### Title: Eight Schools Data
### Aliases: schools
### Keywords: datasets

### ** Examples

  data(schools)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("schools", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.gbp")
### * summary.gbp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.gbp
### Title: Summarizing Estimation Result
### Aliases: summary.gbp
### Keywords: method

### ** Examples


  data(hospital)

  z <- hospital$d
  n <- hospital$n
  y <- hospital$y
  se <- hospital$se
  
  ###################################################################################
  # We do not have any covariates and do not know a mean of the prior distribution. #
  ###################################################################################

    ###############################################################
    # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
    ###############################################################

    g <- gbp(y, se, model = "gaussian")
    summary(g)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "binomial")
    summary(b)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, mean.PriorDist = 0.03, model = "poisson")
    summary(p)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.gbp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
