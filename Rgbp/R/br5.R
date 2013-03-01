BRInitialValueKn <- function(given) {
  # This function makes the initial values needed to run BRIMM.
  # "Kn" means the descriptive second level mean (mean of Beta distribution) is known.
  if (given$prior.mean == 0) {
    r.ini <- mean(given$sample.mean) * (1 - mean(given$sample.mean)) / var(given$sample.mean)
  } else {
    r.ini <- given$prior.mean * (1 - given$prior.mean) / var(given$sample.mean)
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
  r.ini <- p0.ini * (1 - p0.ini) / var(y) 
  list(x = x, b.ini = b.ini, a.ini = -log(r.ini))
}

BRAlphaEstKn2 <- function(given, ini) {
  # Alpha estimation of BRIMM when the second level mean is known
  z <- given$z
  n <- given$n
  k <- length(n)
  a.ini <- ini$a.ini
  p0 <- given$prior.mean
  q0 <- 1 - p0 

  BRDerivAlpha <- function(a) {
    # The first and second order derivatives of log likelihood with respect to alpha
    zap <- z + exp(-a) * p0 
    nzaq <- n - z + exp(-a) * q0
    ap <- exp(-a) * p0 
    aq <- exp(-a) * q0
    al <- exp(-a)

    const1 <- ((digamma(zap) - digamma(ap)) * p0 + (digamma(nzaq) - digamma(aq)) * q0
               + digamma(al) - digamma(n + al))
    const3 <- ((trigamma(zap) - trigamma(ap)) * p0^2 + (trigamma(nzaq) - trigamma(aq)) * q0^2
               + trigamma(al) - trigamma(n + al))

    out <- c(1 - exp(-a) * sum(const1), exp(-a) * (sum(const1) + exp(-a) * sum(const3)))
    out
  }

  dif <- 1
  eps <- 0.0001
  while (abs(dif) > eps) { 
    out1 <- BRDerivAlpha(a.ini)
    score <- out1[1]
    hessian <- out1[2]
    dif <- score / hessian
    a.ini <- a.ini - dif
  }
  list(a.new = a.ini, a.var = - 1 / hessian)
}


BRAlphaBetaEstUn2 <- function(given, ini) {
  # Alpha and Beta estimation of BRIMM when the second level mean is unknown
  z <- given$z
  n <- given$n
  x <- ini$x
  k <- length(n)
  a.ini <- ini$a.ini
  b.ini <- ini$b.ini
  m <- ncol(x)

  BRLogLikUn <- function(a, b) {
    # Log likelihood function of alpha and beta (regression coefficients) for BRIMM 
    # when the descriptive second level mean is unknown.

    p0.hat <- exp(x %*% as.matrix(b)) / (1 + exp(x %*% as.matrix(b)))
    t0 <- p0.hat * exp(-a)
    t1 <- (1 - p0.hat) * exp(-a)
    t2 <- exp(-a)
    if (any(c(t0, t1, t2, n - z + t1, t0 + z) <= 0)) {
      print("The components of lgamma should be positive")
      stop()
    } else {
      sum(lgamma(t0 + z) - lgamma(t0) + lgamma(n - z + t1) - lgamma(t1) + lgamma(t2) - lgamma(n + t2))
    }
  }

  BRDerivBeta <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to beta
    p <- exp(x %*% as.matrix(b)) / (1 + exp(x %*% as.matrix(b)))
    q <- 1 - p
    zap <- z + exp(-a) * p 
    nzaq <- n - z + exp(-a) * q
    ap <- exp(-a) * p 
    aq <- exp(-a) * q
    if (any(c(ap, aq)==0)) {
      vec[(ap == 0 | aq == 0)] <- 0
      tmp <- !(ap == 0 | aq == 0)
      vec[tmp] <- (digamma(zap[tmp]) - digamma(ap[tmp]) - digamma(nzaq[tmp]) 
                   + digamma(aq[tmp])) * exp(-a) * p[tmp] * q[tmp]
    } else {
      vec <- (digamma(zap) - digamma(ap) - digamma(nzaq) + digamma(aq)) * exp(-a) * p * q
    }
    if (any(c(ap, aq)==0)) {
      diag[(ap == 0 | aq == 0)] <- 0
      tmp <- !(ap == 0 | aq == 0)
      diag[tmp] <- ((trigamma(zap[tmp]) - trigamma(ap[tmp]) + trigamma(nzaq[tmp]) - trigamma(aq[tmp])) 
                    * exp(-a) * p[tmp] * q[tmp] 
                    + (digamma(zap[tmp]) - digamma(ap[tmp]) - digamma(nzaq[tmp]) + digamma(aq[tmp])) 
                    * (q[tmp] - p[tmp])) * exp(-a) * p[tmp] * q[tmp]
    } else {
      diag <- ((trigamma(z + exp(-a) * p) - trigamma(exp(-a) * p) 
                + trigamma(n - z + exp(-a) * q) - trigamma(exp(-a) * q)) * exp(-a) * p * q +
               (digamma(z + exp(-a) * p) - digamma(exp(-a) *p)
                - digamma(n - z + exp(-a) * q) + digamma(exp(-a) * q)) * (q - p)) * exp(-a) * p * q
    }
    out <- cbind(t(x) %*% as.vector(vec), t(x) %*% diag(as.numeric(diag)) %*% x)
    out
  } 

  BRDerivBeta2order <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to beta
    p <- exp(x %*% as.matrix(b)) / (1 + exp(x %*% as.matrix(b)))
    q <- 1 - p
    zap <- z + exp(-a) * p 
    nzaq <- n - z + exp(-a) * q
    ap <- exp(-a) * p 
    aq <- exp(-a) * q
    if (any(c(ap, aq)==0)) {
      diag[(ap == 0 | aq == 0)] <- 0
      tmp <- !(ap == 0 | aq == 0)
      diag[tmp] <- ((trigamma(zap[tmp]) - trigamma(ap[tmp]) + trigamma(nzaq[tmp]) - trigamma(aq[tmp])) 
                    * exp(-a) * p[tmp] * q[tmp] 
                    + (digamma(zap[tmp]) - digamma(ap[tmp]) - digamma(nzaq[tmp]) + digamma(aq[tmp])) 
                    * (q[tmp] - p[tmp])) * exp(-a) * p[tmp] * q[tmp]
    } else {
      diag <- ((trigamma(z + exp(-a) * p) - trigamma(exp(-a) * p) 
                + trigamma(n - z + exp(-a) * q) - trigamma(exp(-a) * q)) * exp(-a) * p * q +
               (digamma(z + exp(-a) * p) - digamma(exp(-a) *p)
                - digamma(n - z + exp(-a) * q) + digamma(exp(-a) * q)) * (q - p)) * exp(-a) * p * q
    }
    out <- t(x) %*% diag(as.numeric(diag)) %*% x
    out
  } 

  BetaHatSubAlpha <- function(a) {
    dif <- 1
    eps <- 0.0001
    n.iter <- 1
    while (any(abs(dif) > eps)) { 
      out <- BRDerivBeta(a, b.ini)
      score <- out[, 1]
      hessian <- out[, 2 : (m + 1)]
      dif <- solve(hessian) %*% score
      b.ini <- b.ini - dif
      n.iter <- n.iter + 1
      if (n.iter > 50) {
        stop()  
        print("Alpha estimate does not converge in Newton-Raphson")
      }
    }
    list(beta.new = b.ini, beta.hessian = BRDerivBeta2order(a, b.ini))
  }

  MarginalPostAlpha <- function(a) {
    b.sub.a <-  BetaHatSubAlpha(a)$beta.new
    a + BRLogLikUn(a, b.sub.a) - 0.5 * log(det(-BRDerivBeta2order(a, b.sub.a)))
  }
  a.temp <- optim(a.ini, MarginalPostAlpha, control = list(fnscale = -1), method= "L-BFGS-B",
                  hessian = T,  lower = -Inf, upper = Inf)
  a.new <- a.temp$par
  a.hess <- a.temp$hessian
  b.temp.result <- BetaHatSubAlpha(a.new)
  b.new <- b.temp.result$beta.new
  b.hessian <- b.temp.result$beta.hessian

  list(a.new = a.new, beta.new = b.new, 
       a.var = -1/a.hess, beta.var = -solve(b.hessian))
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

BR5 <- function(z, n, X, prior.mean, intercept = T, Alpha = 0.95){
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
             BRAlphaBetaEstUn2(given, ini)
           } else {
             BRAlphaEstKn2(given, ini)
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
                 model="br", x = X, beta.new = a.res$beta.new, beta.var = a.res$beta.var,
                 intercept = intercept, a.new = a.res$a.new, a.var = a.res$a.var)
  output
}
