pkgname <- "Rgbp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Rgbp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Rgbp-package")
### * Rgbp-package

flush(stderr()); flush(stdout())

### Name: Rgbp
### Title: Bayesian Hierarchical Modeling and Frequentist Method Check
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

    g <- gbp(y, se)
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
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    gcv <- coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gr")
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

    ### when we want to simulate pseudo datasets based on different values of A,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    gcv <- coverage(g, A.or.r = 9, reg.coef = c(10, 1), covariates = x2, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8)
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
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  ###############################################################
  # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
  ###############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "br")
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
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "br")
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

    ### when we want to simulate pseudo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    bcv <- coverage(b, A.or.r = 60, reg.coef = c(-1, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")
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
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  ##############################################################
  # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
  ##############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    p <- gbp(z, n, model = "pr")
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
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, reg.coef = -5, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    p <- gbp(z, n, x1, model = "pr")
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

    ### when we want to simulate pseudo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    pcv <- coverage(p, A.or.r = 60, reg.coef = c(-2, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")
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
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  




cleanEx()
nameEx("baseball")
### * baseball

flush(stderr()); flush(stdout())

### Name: baseball
### Title: Baseball Data
### Aliases: baseball
### Keywords: datasets

### ** Examples


  data(baseball)
  z <- baseball$Hits
  n <- baseball$At.Bats
  x <- ifelse(baseball$Position == "fielder", 1, 0)


  ###############################################################################
  # We do have one covariates and do not know a mean of the prior distribution. #
  ###############################################################################

  ###############################################################
  # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
  ###############################################################

  b <- gbp(z, n, x, model = "br")
  b
  print(b, sort = FALSE)
  summary(b)
  plot(b)
  plot(b, sort = FALSE)

  ##############################################################
  # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
  ##############################################################

  p <- gbp(z, n, x, model = "pr")
  p
  print(p, sort = FALSE)
  summary(p)
  plot(p)
  plot(p, sort = FALSE)




cleanEx()
nameEx("coverage")
### * coverage

flush(stderr()); flush(stdout())

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

    g <- gbp(y, se)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A and 
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    gcv <- coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gr")
 
    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  
 
    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    gcv <- coverage(g, A.or.r = 9, reg.coef = c(10, 1), covariates = x2, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8)

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    gcv <- coverage(g, nsim = 10)  

    ### gcv$coverageRB, gcv$coverage10, gcv$average.coverageRB, gcv$average.coverage10,
    ### gcv$minimum.coverageRB, gcv$minimum.coverage10, gcv$raw.resultRB, gcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of A and
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  ################################################################
  # Binomial Regression Interactive Multi-level Modeling (BRIMM) #
  ################################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "br")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and 
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "br")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    bcv <- coverage(b, A.or.r = 60, reg.coef = c(-1, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    bcv <- coverage(b, nsim = 10)  

    ### bcv$coverageRB, bcv$coverage10, bcv$average.coverageRB, bcv$average.coverage10,
    ### bcv$minimum.coverageRB, bcv$minimum.coverage10, bcv$raw.resultRB, bcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  ###############################################################
  # Poisson Regression Interactive Multi-level Modeling (PRIMM) #
  ###############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    p <- gbp(z, n, model = "pr")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and 
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, reg.coef = -5, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    p <- gbp(z, n, x1, model = "pr")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    pcv <- coverage(p, A.or.r = 60, reg.coef = c(-2, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")

    ### when we want to simulate pseudo datasets considering the estimated values 
    ### as true ones.
    pcv <- coverage(p, nsim = 10)  

    ### pcv$coverageRB, pcv$coverage10, pcv$average.coverageRB, pcv$average.coverage10,
    ### pcv$minimum.coverageRB, pcv$minimum.coverage10, pcv$raw.resultRB, pcv$raw.result10.

    ### when we want to simulate pseudo datasets based on different values of r and
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  




cleanEx()
nameEx("gbp.default")
### * gbp.default

flush(stderr()); flush(stdout())

### Name: gbp
### Title: Fitting Bayesian Hierarchical Models
### Aliases: gbp gbp.default
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

    g <- gbp(y, se)
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
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    gcv <- coverage(g, A.or.r = 9, reg.coef = 10, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    g <- gbp(y, se, x2, model = "gr")
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

    ### when we want to simulate pseudo datasets based on different values of A,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    gcv <- coverage(g, A.or.r = 9, reg.coef = c(10, 1), covariates = x2, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    g <- gbp(y, se, mean.PriorDist = 8)
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
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    coverage(g, A.or.r = 9, mean.PriorDist = 5, nsim = 10)  

  ###############################################################
  # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
  ###############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    b <- gbp(z, n, model = "br")
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
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, reg.coef = -1, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    b <- gbp(z, n, x1, model = "br")
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

    ### when we want to simulate pseudo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    bcv <- coverage(b, A.or.r = 60, reg.coef = c(-1, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    b <- gbp(z, n, mean.PriorDist = 0.265, model = "br")
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
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    bcv <- coverage(b, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  

  ##############################################################
  # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
  ##############################################################

    ####################################################################################
    # If we do not have any covariate and do not know a mean of the prior distribution #
    ####################################################################################

    p <- gbp(z, n, model = "pr")
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
    ### of a regression coefficient (intercept), 
    ### not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, reg.coef = -5, nsim = 10)  

    ##################################################################################
    # If we have one covariate and do not know a mean of the prior distribution yet, #
    ##################################################################################

    p <- gbp(z, n, x1, model = "pr")
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

    ### when we want to simulate pseudo datasets based on different values of r,
    ### of regression coefficients, and of covariate, not using estimated values 
    ### as true ones. Two values of reg.coef are for beta0 and beta1
    pcv <- coverage(p, A.or.r = 60, reg.coef = c(-2, 0), covariates = x1, nsim = 10)  

    ################################################
    # If we know a mean of the prior distribution, #
    ################################################

    p <- gbp(z, n, mean.PriorDist = 0.265, model = "pr")
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
    ### of 2nd level mean as true ones, not using estimated values as true ones.
    pcv <- coverage(p, A.or.r = 60, mean.PriorDist = 0.3, nsim = 10)  




cleanEx()
nameEx("hospital")
### * hospital

flush(stderr()); flush(stdout())

### Name: hospital
### Title: Thirty-one Hospital Data
### Aliases: hospital
### Keywords: datasets

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

    g <- gbp(y, se)
    g
    print(g, sort = FALSE)
    summary(g)
    plot(g)
    plot(g, sort = FALSE)    

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "br")
    b
    print(b, sort = FALSE)
    summary(b)
    plot(b)
    plot(b, sort = FALSE)    

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, model = "pr")
    print(p, sort = FALSE)
    p
    summary(p)
    plot(p)
    plot(p, sort = FALSE)    




cleanEx()
nameEx("plot.gbp")
### * plot.gbp

flush(stderr()); flush(stdout())

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

    g <- gbp(y, se)
    plot(g)
    plot(g, sort = FALSE)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "br")
    plot(b)
    plot(b, sort = FALSE)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, model = "pr")
    plot(p)
    plot(p, sort = FALSE)




cleanEx()
nameEx("print.gbp")
### * print.gbp

flush(stderr()); flush(stdout())

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

    g <- gbp(y, se)
    g
    print(g, sort = FALSE)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "br")
    b
    print(b, sort = FALSE)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, model = "pr")
    p
    print(p, sort = FALSE)




cleanEx()
nameEx("print.summary.gbp")
### * print.summary.gbp

flush(stderr()); flush(stdout())

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

    g <- gbp(y, se)
    summary(g)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "br")
    summary(b)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, model = "pr")
    summary(p)




cleanEx()
nameEx("schools")
### * schools

flush(stderr()); flush(stdout())

### Name: schools
### Title: Eight Schools Data
### Aliases: schools
### Keywords: datasets

### ** Examples

  data(schools)
  y <- schools$y
  se <- schools$se

  ###################################################################################
  # We do not have any covariates and do not know a mean of the prior distribution. #
  ###################################################################################

  ###############################################################
  # Gaussian Regression Interactive Multilevel Modeling (GRIMM) #
  ###############################################################

  g <- gbp(y, se)
  g
  print(g, sort = TRUE)
  summary(g)
  plot(g)




cleanEx()
nameEx("summary.gbp")
### * summary.gbp

flush(stderr()); flush(stdout())

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

    g <- gbp(y, se)
    summary(g)

    ###############################################################
    # Binomial Regression Interactive Multilevel Modeling (BRIMM) #
    ###############################################################

    b <- gbp(z, n, model = "br")
    summary(b)

    ##############################################################
    # Poisson Regression Interactive Multilevel Modeling (PRIMM) #
    ##############################################################

    p <- gbp(z, n, model = "pr")
    summary(p)




### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
