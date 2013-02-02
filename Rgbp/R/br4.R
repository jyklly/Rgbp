BRInitialValueKn <- function(given) {
  # This function makes the initial values needed to run BRIMM.
  # "Kn" means the descriptive second level mean (mean of Beta distribution) is known.
  if (given$prior.mean == 0) {
    r.ini <- mean(given$sample.mean) * (1 - mean(given$sample.mean)) / var(given$sample.mean) -1    
  } else {
    r.ini <- given$prior.mean * (1 - given$prior.mean) / var(given$sample.mean) -1
  }
  list(r.ini = r.ini, a.ini = -log(r.ini))
}

BRInitialValueUn <- function(given) {
  # This function makes the initial values needed to run BRIMM.
  # "Un" means the descriptive second level mean (mean of Beta distribution) is unknown.
  y <- given$sample.mean
  z <- given$z
  n <- given$n
  x.ini <- given$x.ini
  if (identical(x.ini, NA) & !given$intercept) {
    print("In case we do not designate prior mean, at least one beta should be estimated (either intercept=T or at least one covariate is needed)")
    stop()
  } else if (identical(x.ini, NA) & given$intercept) {
    x <- matrix(1, length(y), 1)
    b.ini <- as.vector(log(mean(y) / (1 - mean(y))))
  } else if (!identical(x.ini, NA) & given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("cbind(z, n - z) ~ ", paste(xname, collapse = "+")))
    b.ini <- glm(formula, family = binomial)$coefficients
    x <- as.matrix(cbind(rep(1, length(y)), x.ini))
  } else if (!identical(x.ini, NA) & !given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("cbind(z, n - z) ~ ", paste(xname, collapse = "+"), "- 1"))
    b.ini <- glm(formula, family = binomial)$coefficients
    x <- x.ini
  }	
  p0.ini <- mean(exp(x %*% b.ini) / (1 + exp(x %*% b.ini)))
  r.ini <- p0.ini * (1 - p0.ini) / var(y) - 1
  list(x = x, b.ini = b.ini, a.ini = -log(r.ini))
}

BRAlphaEstKn <- function(given, ini) {
  # Alpha estimation of BRIMM when the second level mean is known
  z <- given$z
  n <- given$n
  k <- length(n)
  a.ini <- ini$a.ini
  p0 <- given$prior.mean
  q0 <- 1 - p0 

  BRDerivAlpha <- function(a) {
    # The first and second order derivatives of log likelihood with respect to alpha
    digamma.z.r.p <- digamma(z + exp(-a) * p0)
    digamma.r.p <- digamma(exp(-a) * p0)
    digamma.n.z.r.q <- digamma(n - z + exp(-a) * q0)
    digamma.r.q <- digamma(exp(-a) * q0)
    digamma.r <- digamma(exp(-a))
    digamma.n.r <- digamma(n + exp(-a))
    trigamma.z.r.p <- trigamma(z + exp(-a) * p0)
    trigamma.r.p <- trigamma(exp(-a) * p0)
    trigamma.n.z.r.q <- trigamma(n - z + exp(-a) * q0)
    trigamma.r.q <- trigamma(exp(-a) * q0)
    trigamma.r <- trigamma(exp(-a))
    trigamma.n.r <- trigamma(n + exp(-a))

    dig.part2 <- ((digamma.z.r.p - digamma.r.p) * p0 + (digamma.n.z.r.q - digamma.r.q) * q0
                  + digamma.r - digamma.n.r)
    trig.part3 <- ((trigamma.z.r.p - trigamma.r.p) * p0^2 + (trigamma.n.z.r.q - trigamma.r.q) * q0^2
                   + trigamma.r - trigamma.n.r)

    const1 <- dig.part2
    const3 <- trig.part3

    out <- c(1 - exp(-a) * sum(const1), exp(-a * 2) * sum(const3))
    out
  }

  dif <- 1
  eps <- 0.0001
  while (abs(dif) > eps) { 
    out1 <- BRDerivAlpha(a.ini)
    score <- out1[1]
    hessian <- out1[2]
    updated <- a.ini - score / hessian
    dif <- a.ini - updated
    a.ini <- updated
  }

  list(a.new = a.ini, a.var = - 1 / hessian)
}

BRAlphaBetaEstUn <- function(given, ini) {
  # Alpha and Beta estimation of BRIMM when the second level mean is unknown
  z <- given$z
  n <- given$n
  x <- ini$x
  k <- length(n)
  a.ini <- ini$a.ini
  b.ini <- ini$b.ini
  m <- ncol(x)

  BRDerivBeta <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to beta
    p <- exp(x %*% b) / (1 + exp(x %*% b))
    q <- 1 - p
    vec <- (digamma(z + exp(-a) * p) - digamma(exp(-a) * p)
            - digamma(n - z + exp(-a) * q) + digamma(exp(-a) * q)) * exp(-a) * p * q
    diag <- ((trigamma(z + exp(-a) * p) - trigamma(exp(-a) * p) 
              + trigamma(n - z + exp(-a) * q) - trigamma(exp(-a) * q)) * exp(-a) * p * q +
             (digamma(z + exp(-a) * p) - digamma(exp(-a) *p)
              - digamma(n - z + exp(-a) * q) + digamma(exp(-a) * q)) * (q - p)) * exp(-a) * p * q
    out <- cbind(t(x) %*% as.vector(vec), t(x) %*% diag(as.numeric(diag)) %*% x)
    out
  }
 
  BRDerivAlpha <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to alpha
    p <- exp(x %*% b) / (1 + exp(x %*% b))
    q <- 1 - p
    
    digamma.z.r.p <- digamma(z + exp(-a) * p)
    digamma.r.p <- digamma(exp(-a) * p)
    digamma.n.z.r.q <- digamma(n - z + exp(-a) * q)
    digamma.r.q <- digamma(exp(-a) * q)
    digamma.r <- digamma(exp(-a))
    digamma.n.r <- digamma(n + exp(-a))
    trigamma.z.r.p <- trigamma(z + exp(-a) * p)
    trigamma.r.p <- trigamma(exp(-a) * p)
    trigamma.n.z.r.q <- trigamma(n - z + exp(-a) * q)
    trigamma.r.q <- trigamma(exp(-a) * q)
    trigamma.r <- trigamma(exp(-a))
    trigamma.n.r <- trigamma(n + exp(-a))
    fourgamma.z.r.p <- psigamma(z + exp(-a) * p, deriv = 2)
    fourgamma.r.p <- psigamma(exp(-a) * p, deriv = 2)
    fourgamma.n.z.r.q <- psigamma(n - z + exp(-a) * q, deriv = 2)
    fourgamma.r.q <- psigamma(exp(-a) * q, deriv = 2)
    fifgamma.z.r.p <- psigamma(z + exp(-a) * p, deriv = 3)
    fifgamma.r.p <- psigamma(exp(-a) * p, deriv = 3)
    fifgamma.n.z.r.q <- psigamma(n - z + exp(-a) * q, deriv = 3)
    fifgamma.r.q <- psigamma(exp(-a) * q, deriv = 3)

    dig.part1 <- digamma.z.r.p - digamma.r.p - digamma.n.z.r.q + digamma.r.q
    dig.part2 <- ((digamma.z.r.p - digamma.r.p) * p + (digamma.n.z.r.q - digamma.r.q) * q 
                  + digamma.r - digamma.n.r)
    trig.part1 <- trigamma.z.r.p - trigamma.r.p + trigamma.n.z.r.q - trigamma.r.q
    trig.part2 <- (trigamma.z.r.p - trigamma.r.p) * p - (trigamma.n.z.r.q - trigamma.r.q) * q
    trig.part3 <- ((trigamma.z.r.p - trigamma.r.p) * p^2 + (trigamma.n.z.r.q - trigamma.r.q) * q^2
                   + trigamma.r - trigamma.n.r)
    fourg.part1 <- (fourgamma.z.r.p - fourgamma.r.p) * p + (fourgamma.n.z.r.q - fourgamma.r.q) * q
    fourg.part2 <- (fourgamma.z.r.p - fourgamma.r.p) * p^2 + (fourgamma.n.z.r.q - fourgamma.r.q) * q^2
    fifg.part1 <- (fifgamma.z.r.p - fifgamma.r.p) * p^2 + (fifgamma.n.z.r.q - fifgamma.r.q) * q^2

    const1 <- dig.part2
    const2 <- ((2 * trig.part1 + exp(-a) * fourg.part1) * exp(-a) * p^2 * q^2 
               + (dig.part1 + exp(-a) * trig.part2) * p * q * (q - p))
    const3 <- trig.part3    
    const4 <- ((2 * (trig.part1 + 2 * exp(-a) * fourg.part1) + exp(-a * 2) * fifg.part1) * p^2 * q^2
               + (2 * trig.part2 + exp(-a) * fourg.part2) * p * q * (q - p))
    sum.diag <- (trig.part1 * exp(-a) * p * q + dig.part1 * (q - p)) * exp(-a) * p * q
    
    if (m == 1) {
      out <- c(1 - exp(-a) * (sum(const1) - sum(const2) / sum(sum.diag) / 2), 
               exp(-a * 2) * (sum(const3) - (sum(const4) / sum(sum.diag) 
               - (sum(const2) / sum(sum.diag))^2) / 2))
      out      
    } else {
      out <- c(1 - exp(-a) * (sum(const1) - m / (2 * k) * sum(const2 / sum.diag)), 
               exp(-a * 2) * (sum(const3) - m / (2 * k) * sum(const4 / sum.diag - (const2 / sum.diag)^2)))
      out
    }
  }

  BRDerivAlphaBeta <- function(a, b) {
    # The cross derivative of log likelihood with respect to alpha and beta
    p <- exp(x %*% b) / (1 + exp(x %*% b))
    q <- 1 - p
    vec <- ((digamma(z + exp(-a) * p) - digamma(exp(-a) * p)
            - digamma(n - z + exp(-a) *q) + digamma(exp(-a) * q)) * p * q
           + ((trigamma(z + exp(-a) * p) - trigamma(exp(-a) * p)) * p
              - (trigamma(n - z + exp(-a) * q) - trigamma(exp(-a) * q)) * q) * exp(-a) * p * q)
    out <- -exp(-a) * t(x) %*% as.vector(vec)
    out
  }

  ini.value <- c(a.ini, b.ini)
  dif <- 1
  eps <- 0.0001
  while (max(abs(dif)) > eps) { 
    out1 <- BRDerivAlpha(ini.value[1], ini.value[2 : (m+1)])
    out2 <- BRDerivBeta(ini.value[1], ini.value[2 : (m+1)])
    out3 <- BRDerivAlphaBeta(ini.value[1], ini.value[2 : (m+1)])
    score <- c(out1[1], out2[, 1])
    hessian <- cbind(c(out1[2], out3), rbind(as.vector(out3), out2[, 2 : (m + 1)]))
    updated <- ini.value - solve(hessian) %*% score
    dif <- ini.value - updated
    ini.value <- updated
  }
  var.a.b <- -solve(hessian)
  list(a.new = ini.value[1], beta.new = ini.value[2 : (m + 1)], 
       a.var = var.a.b[1, 1], beta.var = var.a.b[2 : (m + 1), 2 : (m + 1)])
}

ShrinkageEst <- function(a.res, given) {	
  # This function calculates the shrinkage-related estimates
  a.new <- a.res$a.new
  a.var <- a.res$a.var
  B.hat <- exp(-a.new) / (given$n + exp(-a.new))  # shriankge = B
  inv.info <- 1 / a.var
  var.B.hat <- (B.hat * (1 - B.hat))^2 / ((B.hat * (1 - B.hat)) + inv.info)
  a1.beta <- inv.info / (1 - B.hat)
  a0.beta <- inv.info / B.hat
  central3.B <- 2 * (1 - 2 * B.hat) * B.hat * (1 - B.hat) / 
                (a1.beta + a0.beta + 1) / (a1.beta + a0.beta + 2)
  B.hat.low <- qbeta((1 - given$Alpha) / 2, a1.beta, a0.beta)
  B.hat.upp <- qbeta((1 + given$Alpha) / 2, a1.beta, a0.beta)
  list(B.hat = B.hat, inv.info = inv.info, var.B.hat = var.B.hat, central3.B = central3.B)
}

BRPosteriorEstKn <- function(B.res, given) {
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is known.
  B.hat <- B.res$B.hat
  var.B.hat <- B.res$var.B.hat
  p0 <- given$prior.mean
  central3.B <- B.res$central3.B
  y <- given$sample.mean
  n <- given$n
  p.hat <- y - B.hat * (y - p0)
  c1 <- y * (1 - y)
  c2 <- -3 * y^2 + 2 * (1 + p0) * y - p0
  c3 <- (y - p0) * (3 * y - 1 - p0)
  c4 <- (y - p0)^2
  d1 <- c1 + c3 * B.hat^2
  d2 <- c2 + 2 * c3 * B.hat - 3 * c4 * B.hat^2
  d3 <- c3 - 3 * c4 * B.hat
  d4 <- c4
  var.p.hat <- (d1 - d2 * B.hat - d3 * var.B.hat + d4 * central3.B) / n
  a1.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * p.hat
  a0.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * (1 - p.hat)
  p.hat.low <- qbeta((1 - given$Alpha) / 2, a1.beta.p, a0.beta.p)
  p.hat.upp <- qbeta((1 + given$Alpha) / 2, a1.beta.p, a0.beta.p)
  list(post.mean = p.hat, post.sd = sqrt(var.p.hat), 
       post.intv.low = p.hat.low, post.intv.upp = p.hat.upp, prior.mean = p0)
}

BRPosteriorEstUn <- function(B.res, a.res, ini, given){
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is unknown.
  x <- as.matrix(ini$x)
  n <- given$n
  beta.new <- a.res$beta.new
  beta.var <- a.res$beta.var
  B.hat <- B.res$B.hat
  y <- given$sample.mean
  xVx <- diag(x %*% beta.var %*% t(x))
  mu0 <- exp(x %*% beta.new + xVx / 2)
  b0 <- (1 + mu0) / (mu0 * (exp(xVx) - 1)) + 2
  b1 <- mu0 * (b0 - 1)

  k1p <- b1 / (b1 + b0)
  k2p <- as.vector(k1p * (1 - k1p) / (b1 + b0 + 1))
  k1b <- as.vector(B.hat)
  k2b <- as.vector(B.res$var.B.hat)
  k3b <- as.vector(B.res$central3.B)

  p0.hat <- k1p
  p.hat <- (1 - B.hat)*y + B.hat*p0.hat

  c1 <- as.vector(y - y^2 - (y - 3 * y^2) * B.hat^2 - 2 * y^2 * B.hat^3
                  - (B.hat^2 - 2 * B.hat^3) * k1p^2)
  c2 <- as.vector(-2 * y + 3 * y^2 + 2 * (y - 3 * y^2) * B.hat + 3 * y^2 * B.hat^2
                  - (-2 * B.hat + 3 * B.hat^2) * k1p^2)
  c3 <- as.vector(y - 3 * y^2 + 3 * y^2 * B.hat - (3 * B.hat - 1) * k1p^2)
  c4 <- as.vector(y^2 - k1p^2)
  c5 <- as.vector(4 * y * B.hat^3 - (4 * y - 1) * B.hat^2 + 2 * (B.hat^2 - 2 * B.hat^3) * k1p)
  c6 <- as.vector(2 * (4 * y - 1) * B.hat - 6 * y * B.hat^2 - 2 * y + 1 
                  + 2 * (-2 * B.hat + 3 * B.hat^2) * k1p)
  c7 <- as.vector(4 * y - 1 - 6 * y * B.hat + 2 * (3 * B.hat - 1) * k1p)
  c8 <- as.vector(-2 * y + 2 * k1p)
  c9 <- as.vector(B.hat^2 - 2 * B.hat^3)
  c10 <- as.vector(-2 * B.hat + 3 * B.hat^2)
  c11 <- as.vector(3 * B.hat - 1)

  var.p.hat <- ((c1 + c2 * k1b + c3 * k2b + c4 * k3b + c5 * k1p + c6 * k1p * k1b
                 + c7 * k1p * k2b + c8 * k1p * k3b + c9 * k2p + c10 * k1b * k2p
                 + c11 * k2b * k2p + k3b * k2p) / n 
                + k2b * k2p + k2b * k1p^2 - 2 * y * k2b * k1p + y^2 * k2b + k1b^2 * k2p)
  a1.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * p.hat
  a0.beta.p <- (p.hat * (1 - p.hat) / var.p.hat - 1) * (1 - p.hat)
  p.hat.low <- qbeta((1 - given$Alpha) / 2, a1.beta.p, a0.beta.p)
  p.hat.upp <- qbeta((1 + given$Alpha) / 2, a1.beta.p, a0.beta.p)
  list(post.mean = p.hat, post.sd = sqrt(var.p.hat),
       post.intv.low = p.hat.low, post.intv.upp = p.hat.upp, prior.mean = p0.hat)
}

BR <- function(z, n, X, prior.mean, intercept = T, Alpha = 0.9167){
  # The main function of PRIMM
  X <- if (missing(X)) {
         X <- NA
       } else {
         X <- X
       }
  prior.mean <- if (missing(prior.mean)) {
                  prior.mean <- NA
                } else {
                  prior.mean <- prior.mean
                }
  given <- list(z = z, n = n, sample.mean = z / n, x.ini = X,
                prior.mean = prior.mean, intercept = intercept, Alpha = Alpha)
  if (is.na(prior.mean)) {
    ini <- BRInitialValueUn(given)
  } else {
    ini <- BRInitialValueKn(given)
  }
  a.res <- if (is.na(prior.mean)) {
             BRAlphaBetaEstUn(given, ini)
           } else {
             BRAlphaEstKn(given, ini)
           }
  B.res <- ShrinkageEst(a.res, given)
  if (is.na(prior.mean)) {
    post.res <- BRPosteriorEstUn(B.res, a.res, ini, given)
  } else {
    post.res <- BRPosteriorEstKn(B.res,given)
  }

  output <- list(sample.mean = given$sample.mean, se = given$n, prior.mean = post.res$prior.mean,
                 shrinkage = B.res$B.hat, sd.shrinkage = sqrt(B.res$var.B.hat), 
                 post.mean = post.res$post.mean, post.sd = post.res$post.sd, 
                 post.intv.low = post.res$post.intv.low, post.intv.upp = post.res$post.intv.upp,
                 model="br", x = X, beta.new = a.res$beta.new, beta.hess = a.res$beta.hess,
                 intercept = intercept, a.new = a.res$a.new, a.var = 1 / B.res$inv.info)
  output
}


# test1
n<-rep(100,10)
z<-c(0,0,1,1,2,2,3,3,4,4)
x1<-X<-c(1,1,1,1,0,0,0,0,0,0)

system.time(b<-BR(z,n,prior.mean=0.02))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x1))
mean(b$shrinkage)

# test2
n<-rep(100,10)
z<-c(1,1,1,1,2,2,3,3,3,3)
x<-c(1,1,1,1,0,0,0,0,0,0)

system.time(b<-BR(z,n,prior.mean=0.02))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)

# test3
n<-rep(100,10)
z<-c(rep(1,5),rep(3,5))
x<-c(1,1,1,1,0,0,0,0,0,0)

system.time(b<-BR(z,n,prior.mean=0.02))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)

# test4
n<-rep(100,10)
z<-c(36, 37, 37,38, 42, 42, 42, 45, 45, 50)
x<-c(1,1,1,1,0,0,0,0,0,0)

system.time(b<-BR(z,n,prior.mean=0.414))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)

# test5
n<-rep(10000,10)
z<-c(0, 0, 0, 0, 1, 1, 1, 2, 2, 3)
x<-c(1,1,1,1,0,0,0,0,0,0)

system.time(b<-BR(z,n,prior.mean=0.0001))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)

# test6
n<-rep(1000,10)
z<-c(476, 482, 488, 494, 500,500, 506, 512, 518, 524)
x<-c(1,1,1,1,0,0,0,0,0,0)


system.time(b<-BR(z,n,prior.mean=0.5))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)

# test7
n<-c(10,10,10,10,10)
z<-c( 3, 4, 1, 2, 0)
x<-c( 1, 1, 0, 0, 0)

system.time(b<-BR(z,n,prior.mean=0.2))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)


# test8
n<-rep(10,10)
z<-c(2,2,2,2,3,3,3,4,5,5)
x<-c(1,1,1,1,0,0,0,0,0,0)

system.time(b<-BR(z,n,prior.mean=0.4))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)


# test9
z<-baseball[,4]
n<-baseball[,3]

system.time(b<-BR(z,n))
mean(b$shrinkage)



# test10
nba<-read.csv("demo_nba.csv",header=T)
nba$ThreePA<-round(nba$GamePlayed*nba$X3PA,0)
nba$ThreePM<-round(nba$GamePlayed*nba$X3PM,0)
nba$totalFGA<-nba$TwoPA+nba$ThreePA
nba$totalFGM<-nba$TwoPM+nba$ThreePM
nba$TAR<-round(nba$TwoPA/nba$totalFGA,3)
nba$Threeratio<-round(nba$ThreePA/nba$totalFGA,3)
name<-nba[c(1:12),1]
nba$logMIN<-round(log(nba$avgMIN),3)
nba<-nba[c(1:12),c(1,7,4,5,6,38,35,34,2,3,36)]
z<-nba$totalFGM
n<-nba$totalFGA
x<-X<-cbind(log(nba$avgMIN),nba$TAR)


system.time(b<-BR(z,n,prior.mean=0.4))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)


# test11
n<-c(18,24,28,37,42,48,56,60,68,73,75,79,91,99,104)
z<-c(4,1,3,1,2,1,2,3,6,0,8,3,3,9,7)
x<-c(0.251,-0.021,0.208,0.019,-0.169,0.164,0.296,0.199,0.209,0.093,0.002,0.064,0.105,0.073,0.209)

system.time(b<-BR(z,n,prior.mean=0.05))
mean(b$shrinkage)
system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)


# test12
system.time(b<-BR(test$z1,test$n,prior.mean=0.49))
mean(b$shrinkage)
system.time(b<-BR(test$z1,test$n))
mean(b$shrinkage)
system.time(b<-BR(test$z1,test$n,test[,4:6]))
mean(b$shrinkage)


system.time(b<-bp(test$z1,test$n,prior.mean=0.49,model="br"))
mean(b$shrinkage)
system.time(b<-bp(test$z1,test$n,model="br"))
mean(b$shrinkage)
system.time(b<-bp(test$z1,test$n,test[,4:6],model="br"))
mean(b$shrinkage)

# test13
bb<-read.csv("2005dataset.csv",header=T)
head(bb)
n1<-ab_1st<-bb$AB.4.+bb$AB.5.+bb$AB.6.
z1<-hit_1st<-bb$H.4.+bb$H.5.+bb$H.6.
n2<-ab_2nd<-bb$AB.7.+bb$AB.8.+bb$AB.9.10.
z2<-hit_2nd<-bb$H.7.+bb$H.8.+bb$H.9.10.
y1<-y_1st<-hit_1st/ab_1st
y2<-y_2nd<-hit_2nd/ab_2nd
x<-bb$isPitcher
data<-as.data.frame(cbind(hit_1st, ab_1st, y_1st, hit_2nd, ab_2nd, y_2nd, x))
data<-data[order(data$ab_1st,data$ab_2nd),]
first.half.0<-(data[,2]==0)
second.half.0<-(data[,5]==0)
excluded_without_justification<-(data[,2]>0 & data[,2]<=10)
fit.index0<-(data[,2]>0)
eval.index0<-(data[,2] >0 & data[,5]>0)
fit.index10<-(data[,2]>10)
eval.index10<-(data[,2] >10 & data[,5]>10)
data0<-data[fit.index0,]
index<-1:nrow(data0)
data0<-cbind(data0,index)
z0<-data0[,1]
n0<-data0[,2]
y0<-z0/n0
x0<-data0[,7]
data$eval<-fit.index10*eval.index10
data3<-data[fit.index10,]
data.eval<-data3[which(data3[,8]==1 ),]
z10<-data3[,1]
n10<-data3[,2]
y10<-z10/n10
x10<-data3[,7]

system.time(b<-BR(z0,n0))
mean(b$shrinkage)
system.time(b<-BR(z0,n0,x0))
mean(b$shrinkage)
system.time(b<-BR(z0,n0,cbind(x0,log(n0))))
mean(b$shrinkage)


# test14
system.time(b<-BR(z10,n10))
mean(b$shrinkage)
system.time(b<-BR(z10,n10,x10))
mean(b$shrinkage)
system.time(b<-BR(z10,n10,cbind(x10,log(n10))))
mean(b$shrinkage)


# test15
ALL<-c(147,134,134,122,118,117,115,114,111,111,109,100,100,97,95,94,88,87,86,84,84,83,81,79,78,77,74,69,66,62)
ALLPG<-c(1.62,1.46,1.46,1.33,1.3,1.3,1.25,1.24,1.22,1.25,1.18,1.09,1.09,1.07,1.03,1.03,0.98,0.97,0.93,0.95,0.91,0.9,0.88,0.9,0.85,0.84,0.8,0.78,0.73,0.67)
ngp<-ALL/ALLPG
exps<-ngp*40
MPF<-c(26,21,33,21,27,21,24,14,29,28,19,24,21,17,29,13,18,19,22,20,22,24,25,6,14,17,12,12,4,17)
dvs<-c(0,0,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1,0,1,1,1,1,0,1,1,0,0,1,0,1)
row.name<-c("FLA","ARI","ATL","WAS","KCR","TBR","SFG","CHW","COL","MIN","CIN","HOU","LAD","SEA","SDP","NYM","BOS","PIT","DET","BAL","TOR","OAK","MIL","NYY","TEX","CHC","STL","CLE","PHI","LAA")

system.time(b<-bp(MPF,exps,model="br"))
mean(b$shrinkage)

system.time(b<-BR(MPF,exps))
mean(b$shrinkage)



# test16
n<-c(33,50,70,63)
z<-c(23,40,62,40)
x1<-c(13.93,2.27)
x2<-c(9.19,1.5)
x3<-c(19.13,2.09)
x4<-c(18.97,1.94)
x<-rbind(x1,x2,x3,x4)

system.time(b<-bp(z,n,x,model="br"))
mean(b$shrinkage)

system.time(b<-BR(z,n,x))
mean(b$shrinkage)

# test17
n<-rep(2,27)
z<-c( 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1)
left<-rep(c(1,1,1,0,0,0,0,0,0),3)
right<-rep(c(0,0,0,1,1,1,0,0,0),3)
x<-cbind(left,right)

system.time(b<-BR(z,n))
mean(b$shrinkage)
system.time(b<-BR(z,n,x))
mean(b$shrinkage)

# test18
new.n<-rep(3,9)
new.z<-c(1,1,1,2,0,1,1,2,3)
new.left<-rep(c(1,0,0),3)
new.right<-rep(c(0,1,0),3)
new.x<-cbind(new.left,new.right)

system.time(b<-bp(new.z,new.n,model="br"))
mean(b$shrinkage)
system.time(b<-bp(new.z,new.n,new.x,model="br"))
mean(b$shrinkage)

system.time(b<-BR(new.z,new.n))
mean(b$shrinkage)
system.time(b<-BR(new.z,new.n,new.x))
mean(b$shrinkage)

# test19
hp<-read.csv("hospitals.csv",head=T)
z<-hp[,4]
n<-hp[,5]

system.time(b<-bp(z,n,model="br"))
mean(b$shrinkage)

system.time(b<-BR(z,n))
mean(b$shrinkage)

# test20
rat_dat <- read.table("rats_alex.txt", header=TRUE)
z <- rat_dat$y
n <- rat_dat$N

system.time(b<-bp(z,n,model="br"))
mean(b$shrinkage)

system.time(b<-BR(z,n))
mean(b$shrinkage)


# test21
ab<-c(150,358,605,581,573,506,591,553,428,345,563,542,588,421,497,428,507,526,556,625,415,233,490,353)
hit<-c(36,113,212,188,216,194,248,226,167,127,208,201,225,161,191,143,197,211,189,211,157,79,175,114)

system.time(b<-bp(hit,ab,model="br"))
mean(b$shrinkage)

system.time(b<-BR(hit,ab))
mean(b$shrinkage)
