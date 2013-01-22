PRInitialValueKn <- function(given) {
  # This function makes the initial values needed to run PRIMM.
  # "Kn" means the descriptive second level mean (mean of Beta distribution) is known.
  r.ini <- given$prior.mean / var(given$sample.mean)
  list(a.ini = -log(r.ini))
}

PRInitialValueUn <- function(given) {
  # This function makes the initial values needed to run PRIMM.
  # "Un" means the descriptive second level mean (mean of Beta distribution) is unknown.
  y <- given$sample.mean
  x.ini <- given$x.ini
  if (identical(x.ini, NA) & !given$intercept) {
    print("In case we do not designate prior mean, at least one beta should be estimated (either intercept=T or at least one covariate is needed)")
    stop()
  } else if (identical(x.ini, NA) & given$intercept) {
    x <- matrix(1, length(y), 1)
    b.ini <- mean(y)
  } else if (!identical(x.ini, NA) & given$intercept) {
    if (any(y == 0)) {
      y[which(y == 0)] <- 0.00001
    }
    x <- as.matrix(cbind(rep(1, length(y)), x.ini))
    b.ini <- solve(t(x) %*% x) %*% t(x) %*% log(y)
    if (any(y == 0.00001)) {
      y[which(y == 0.00001)] <- 0
    }
  } else if (!identical(x.ini,NA) & !given$intercept) {
    if (any(y == 0)) {
      y[which(y == 0)] <- 0.00001
    }
    x <- x.ini
    b.ini <- solve(t(x) %*% x) %*% t(x) %*% log(y)
    if (any(y == 0.00001)){
      y[which(y == 0.00001)] <- 0
    }
  }	
  mu0.ini <- mean(exp(x %*% b.ini))
  r.ini <- mu0.ini / var(y)
  list(x=x, b.ini=b.ini, a.ini= -log(r.ini))
}

PRLogLikKn <- function(a, given) {
  # Log likelihood function of alpha for PRIMM when the descriptive second level mean is known.
  z <- given$z
  n <- given$n
  prior.mean <- given$prior.mean
  sum(dnbinom(z, exp(-a) * prior.mean, prob = exp(-a) / (exp(-a) + n), log=T))
}

PRLogLikUn <- function(a, b, given, ini) {
  # Log likelihood function of alpha and beta (regression coefficients) for PRIMM 
  # when the descriptive second level mean is unknown.
  z <- given$z
  n <- given$n
  x <- ini$x
  sum(dnbinom(z, exp(-a + x %*% as.vector(b)), prob = exp(-a) / (exp(-a) + n), log=T))
}

PRAlphaEstKn <- function(given, ini) {
  # Alpha estimation of PRIMM when the descriptive second level mean is known
  z <- given$z
  n <- given$n
  mu0 <- given$prior.mean
  k <- length(n)
  a.ini <- ini$a.ini

  PRDerivAlpha <- function(a) {
    # The first and second order derivatives of log likelihood with respect to alpha
    digamma.z.r.mu0 <- digamma(z + exp(-a) * mu0)
    digamma.r.mu0 <- digamma(exp(-a) * mu0)
    trigamma.z.r.mu0 <- trigamma(z + exp(-a) * mu0)
    trigamma.r.mu0 <- trigamma(exp(-a) * mu0)

    dig.part1 <- digamma.z.r.mu0 - digamma.r.mu0 - log(1 + n * exp(a))
    trig.part1 <- trigamma.z.r.mu0 - trigamma.r.mu0

    const1 <- (dig.part1 + n / (exp(-a) + n)) * mu0 - z / (exp(-a) + n)
    const3 <- (trig.part1 * mu0 + n^2 * exp(a) / (exp(-a) + n)^2) * mu0 + z / (exp(-a) + n)^2

    out <- c(1 - exp(-a) * sum(const1), exp(-a * 2) * sum(const3))
    out
  }

  dif <- 1
  eps <- 0.0001
  while (abs(dif) > eps) { 
    out1 <- PRDerivAlpha(a.ini)
    score <- out1[1]
    hessian <- out1[2]
    updated <- a.ini - score / hessian
    dif <- a.ini - updated
    a.ini <- updated
  }
  list(a.new = a.ini, a.hess = hessian)
}

PRAlphaBetaEstUn <- function(given, ini) {
  # Alpha and Beta estimation of PRIMM when the descriptive second level mean is unknown
  z <- given$z
  n <- given$n
  x <- ini$x
  k <- length(n)
  a.ini <- ini$a.ini
  b.ini <- ini$b.ini
  m <- ncol(x)

  PRDerivBeta <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to beta
    mu0 <- exp(x %*% b)
    vec <- (digamma(z + exp(-a) * mu0) - digamma(exp(-a) * mu0) 
            - log(1 + n * exp(a))) * exp(-a) * mu0
    diag <- (trigamma(z + exp(-a) * mu0) - trigamma(exp(-a) * mu0)) * exp(-a * 2) * mu0^2 + vec
    out <- cbind(t(x) %*% as.vector(vec), t(x) %*% diag(as.numeric(diag)) %*% x)
    out
  }
 
  PRDerivAlpha <- function(a, b) {
    # The first and second order derivatives of log likelihood with respect to alpha
    mu0 <- exp(x %*% b)
    digamma.z.r.mu0 <- digamma(z + exp(-a) * mu0)
    digamma.r.mu0 <- digamma(exp(-a) * mu0)
    trigamma.z.r.mu0 <- trigamma(z + exp(-a) * mu0)
    trigamma.r.mu0 <- trigamma(exp(-a) * mu0)
    fourgamma.z.r.mu0 <- psigamma(z + exp(-a) * mu0, deriv = 2)
    fourgamma.r.mu0 <- psigamma(exp(-a) * mu0, deriv = 2)
    fifgamma.z.r.mu0 <- psigamma(z + exp(-a) * mu0, deriv = 3)
    fifgamma.r.mu0 <- psigamma(exp(-a) * mu0, deriv = 3)

    dig.part1 <- digamma.z.r.mu0 - digamma.r.mu0 - log(1 + n * exp(a))
    trig.part1 <- trigamma.z.r.mu0 - trigamma.r.mu0
    trig.part2 <- trig.part1 * mu0 + n * exp(a) / (exp(-a) + n)
    fourg.part1 <- (fourgamma.z.r.mu0 - fourgamma.r.mu0) * mu0
    fourg.part2 <- ((fourgamma.z.r.mu0 - fourgamma.r.mu0) * mu0^2 
                    - n * (2 * exp(-a) + n) * exp(2 * a) / (exp(-a) + n)^2)
    fifg.part1 <- (fifgamma.z.r.mu0 - fifgamma.r.mu0) * mu0^2 * exp(-a)

    sum.diag <- (trig.part1 * exp(-a) * mu0 + dig.part1) * exp(-a) * mu0
    const1 <- (dig.part1 + n / (exp(-a) + n)) * mu0 - z / (exp(-a) + n)
    const2 <- ((fourg.part1 * exp(-a) + 2 * trig.part1) * exp(-a) * mu0^2
               + (trig.part2 * exp(-a) + dig.part1) * mu0)
    const3 <- (trig.part1 * mu0 + n^2 * exp(a) / (exp(-a) + n)^2) * mu0 + z / (exp(-a) + n)^2
    const4 <- ((fifg.part1 + 4 * fourg.part1) * exp(-a) * mu0 + 2 * mu0 * trig.part1 
               + 2 * trig.part2 + exp(-a) * fourg.part2) * mu0
    out <- c(1 - exp(-a) * (sum(const1) - m / (2 * k) * sum(const2 / sum.diag)), 
             exp(-2 * a) * (sum(const3) - m / (2 * k) * sum(const4 / sum.diag - (const2 / sum.diag)^2)))
    out
  }

  PRDerivAlphaBeta <- function(a, b) {
    # The cross derivative of log likelihood with respect to alpha and beta
    mu0 <- exp(x %*% b)
    digamma.z.r.mu0 <- digamma(z + exp(-a) * mu0)
    digamma.r.mu0 <- digamma(exp(-a) * mu0)
    trigamma.z.r.mu0 <- trigamma(z + exp(-a) * mu0)
    trigamma.r.mu0 <- trigamma(exp(-a) * mu0)

    dig.part1 <- digamma.z.r.mu0 - digamma.r.mu0 - log(1 + n * exp(a))
    trig.part1 <- trigamma.z.r.mu0 - trigamma.r.mu0
    trig.part2 <- trig.part1 * mu0 + n * exp(a) / (exp(-a) + n)
    vec <- (dig.part1 + exp(-a) * trig.part2) * mu0
    out <- -exp(-a) * t(x) %*% as.vector(vec)
    out
  }

  ini.value <- c(a.ini, b.ini)
  dif <- 1
  eps <- 0.0001
  while (max(abs(dif)) > eps) { 
    out1 <- PRDerivAlpha(ini.value[1], ini.value[2 : (m+1)])
    out2 <- PRDerivBeta(ini.value[1], ini.value[2 : (m+1)])
    out3 <- PRDerivAlphaBeta(ini.value[1], ini.value[2 : (m+1)])
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
  B.hat.low <- qbeta((1 - given$Alpha) / 2, a1.beta, a0.beta)
  B.hat.upp <- qbeta((1 + given$Alpha) / 2, a1.beta, a0.beta)
  list(B.hat = B.hat, inv.info = inv.info, var.B.hat = var.B.hat)
}

PRPosteriorEstKn <- function(B.res, given) {
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is known.
  B.hat <- B.res$B.hat
  var.B.hat <- B.res$var.B.hat
  mu0 <- given$prior.mean
  n <- given$n
  y <- given$sample.mean
  lambda.hat <- y - B.hat * (y - mu0)
  var.lambda.hat <- (lambda.hat - B.hat * y + (var.B.hat + B.hat^2) * y - (var.B.hat + B.hat^2) * mu0) / n
                    + (y - mu0)^2 * var.B.hat
  u.gamma <- lambda.hat^2 / var.lambda.hat
  v.gamma <- lambda.hat / var.lambda.hat
  lambda.hat.low <- qgamma((1 - given$Alpha) / 2, shape = u.gamma, rate = v.gamma)
  lambda.hat.upp <- qgamma((1 + given$Alpha) / 2, shape = u.gamma, rate = v.gamma)
  list(post.mean = lambda.hat, post.sd = sqrt(var.lambda.hat), post.intv.low = lambda.hat.low,
       post.intv.upp = lambda.hat.upp, prior.mean = mu0)
}

PRPosteriorEstUn <- function(B.res, a.res, ini, given){
  # This function calculates posterior-related estimates
  # when the descriptive second level mean is unknown.
  x <- as.matrix(ini$x)
  n <- given$n
  y <- given$sample.mean
  beta.new <- a.res$beta.new
  B.hat <- B.res$B.hat
  var.B.hat <- B.res$var.B.hat
  V <- solve(-a.res$beta.hess)
  xSx <- diag(x %*% V %*% t(x))
  mu0.hat <- exp(x%*%beta.new + xSx/2)
  var.mu0 <- mu0.hat^2 * (exp(xSx) - 1)
  lambda.hat <- y - B.hat * (y - mu0.hat)
  var.lambda.hat <- ((lambda.hat - B.hat * y + (var.B.hat + B.hat^2) * y 
                      - (var.B.hat + B.hat^2) * mu0.hat) / n
                     + (var.B.hat + B.hat^2) * (y^2 - 2 * y * mu0.hat + var.mu0 + mu0.hat^2)
                     - (B.hat * (y - mu0.hat))^2)
  u.gamma <- lambda.hat^2 / var.lambda.hat
  v.gamma <- lambda.hat / var.lambda.hat
  lambda.hat.low <- qgamma((1 - given$Alpha) / 2, shape = u.gamma, rate = v.gamma)
  lambda.hat.upp <- qgamma((1 + given$Alpha) / 2, shape = u.gamma, rate = v.gamma)
  list(post.mean = lambda.hat, post.sd = sqrt(var.lambda.hat), 
       post.intv.low = lambda.hat.low, post.intv.upp = lambda.hat.upp, 
       prior.mean = mu0.hat)
}

# main function
PR <- function(z, n, X, prior.mean, intercept = T, Alpha = 0.9167) {
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

  given <- list(z = z, n = n, sample.mean = z/n, x.ini = X, prior.mean = prior.mean,
                intercept = intercept, Alpha = Alpha)

  if (is.na(prior.mean)) {
    ini <- PRInitialValueUn(given)
  }else{
    ini <- PRInitialValueKn(given)
  }

  a.res <- if (is.na(prior.mean)) {
             PRAlphaBetaEstUn(given, ini)
           } else {
             PRAlphaEstKn(given, ini)
           }

  B.res <- ShrinkageEst(a.res, given)

  if (is.na(prior.mean)) {
    post.res <- PRPosteriorEstUn(B.res, a.res, ini, given)
  }else{
    post.res <- PRPosteriorEstKn(B.res, given)
  }

  output <- list(sample.mean = given$sample.mean, se = given$n, prior.mean = post.res$prior.mean,
                 shrinkage = B.res$B.hat, sd.shrinkage = sqrt(B.res$var.B.hat), 
                 post.mean = post.res$post.mean, post.sd = post.res$post.sd, 
                 post.intv.low = post.res$post.intv.low, post.intv.upp = post.res$post.intv.upp,
                 model = "pr", x = X, beta.new = a.res$beta.new, beta.hess = a.res$beta.hess,
                 intercept = intercept, a.new = a.res$a.new, a.var = 1 / B.res$inv.info)
  output
}