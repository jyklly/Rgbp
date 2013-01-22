BRInitialValueKn <- function(given) {
  # This function makes the initial values needed to run BRIMM.
  # "Kn" means the descriptive second level mean (mean of Beta distribution) is known.
  r.ini <- given$prior.mean * (1 - given$prior.mean) / var(given$sample.mean) -1
  list(r.ini = r.ini, a.ini = -log(r.ini))
}

BRInitialValueUn <- function(given) {
  # This function makes the initial values needed to run BRIMM.
  # "Un" means the descriptive second level mean (mean of Beta distribution) is unknown.
  y <- given$sample.mean
  x.ini <- given$x.ini
  if (identical(x.ini, NA) & !given$intercept) {
    print("In case we do not designate prior mean, at least one beta should be estimated (either intercept=T or at least one covariate is needed)")
    stop()
  } else if (identical(x.ini, NA) & given$intercept) {
    x <- matrix(1, length(y), 1)
    b.ini <- as.vector(log(mean(y) / (1 - mean(y))))
  } else if (!identical(x.ini, NA) & given$intercept) {
    if (any(y == 0)) {
      y[which(y == 0)] <- 0.00001
    }
    if (any(y==1)) {
      y[which(y == 1)] <- 0.99999
    }		
    x <- as.matrix(cbind(rep(1, length(y)), x.ini))
    b.ini <- solve(t(x) %*% x) %*% t(x) %*% log(y/(1 - y))
    if (any(y==0.00001)) {
      y[which(y == 0.00001)] <- 0
    }
    if (any(y == 0.99999)) {
      y[which(y == 0.99999)] <- 1
    }
  } else if (!identical(x.ini, NA) & !given$intercept) {
    if (any(y==0)) {
      y[which(y == 0)] <- 0.00001
    }
    if (any(y==1)) {
      y[which(y == 1)] <- 0.99999
    }
    x <- x.ini
    b.ini <- solve(t(x) %*% x) %*% t(x) %*% log(y / (1 - y))
    if (any(y == 0.00001)) {
      y[which(y == 0.00001)] <- 0
    }
    if (any(y == 0.99999)) {
      y[which(y == 0.99999)] <- 1
    }
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

  list(a.new = a.ini, a.hess = hessian)
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

    out <- c(1 - exp(-a) * (sum(const1) - m / (2 * k) * sum(const2 / sum.diag)), 
             exp(-a * 2) * (sum(const3) - m / (2 * k) * sum(const4 / sum.diag - (const2 / sum.diag)^2)))
    out
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

  list(a.new = ini.value[1], beta.new = ini.value[2 : (m + 1)], 
       a.hess = hessian[1, 1], beta.hess = hessian[2 : (m + 1), 2 : (m + 1)])
}

ShrinkageEst <- function(a.res, given) {	
  # This function calculates the shrinkage-related estimates
  a.new <- a.res$a.new
  a.hess <- a.res$a.hess
  B.hat <- exp(-a.new) / (given$n + exp(-a.new))  # shriankge = B
  inv.info <- -a.hess
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
  beta.hess <- a.res$beta.hess
  B.hat <- B.res$B.hat
  y <- given$sample.mean
  xHx <- diag(x %*% solve(-beta.hess) %*% t(x))
  mu0 <- exp(x %*% beta.new + xHx / 2)
  b0 <- (1 + mu0) / (mu0 * (exp(xHx) - 1)) + 2
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