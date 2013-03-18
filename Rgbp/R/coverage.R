
coverage <- function(gbp.object, A.or.r, reg.coef, covariates, mean.PriorDist, nsim = 10) {

  # Rao-Blackwellized criterion
  coverageRB <- matrix(NA, nrow = length(gbp.object$se), ncol = nsim)

  # 1-0 criterion that is 1 if interval includes true parameter, 0 if not
  coverage10 <- matrix(NA, nrow = length(gbp.object$se), ncol = nsim)

  if (missing(A.or.r)) {
    only.gbp.result <- TRUE
  } else {
    only.gbp.result <- FALSE
  }

  # if model=BRIMM	
  if (gbp.object$model == "br") {
    if (only.gbp.result) {
      # 1. initial values
      if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(gbp.object$beta.new)
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
      } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
        betas <- as.vector(gbp.object$beta.new)
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        X <- gbp.object$X
      } else if (!is.na(gbp.object$prior.mean)) {
        p0 <- gbp.object$prior.mean
      }
  
      n <- gbp.object$se
      r <- exp(-gbp.object$a.new)

      # 2. generate p matrix
      sim.p <- matrix(rbeta(length(n) * nsim, r * p0, r * (1 - p0)), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rbinom(nrow(sim.p) * nsim, n, sim.p), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, model = "br")
                 } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, X, model = "br")
                 } else if (!is.na(gbp.object$prior.mean)) {
                   gbp(sim.z[, i], n, mean.PriorDist = p0, model = "br")
                 }

          a1 <- r * p0 + sim.z[, i]
          a0 <- r * (1 - p0) + n - sim.z[, i]
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pbeta(upp, a1, a0) - pbeta(low, a1, a0)
          coverage10[, i] <- ifelse(low <= sim.p[, i] & sim.p[, i] <= upp, 1, 0)
        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {
      # 1. initial values
      if (missing(A.or.r)) {
        print("(A.or.r, reg.coef, covariates (if any) or 
               (A.or.r, mean.PriorDist) should be designated")
        stop()
      } else if (!missing(mean.PriorDist)) {
        p0 <- mean.PriorDist
      } else if (!missing(reg.coef) & missing(covariates)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(reg.coef)
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
      } else if (!missing(reg.coef) & !missing(covariates)) { 
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), covariates))
        betas <- as.vector(reg.coef)
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
      }
  
      n <- gbp.object$se
      r <- A.or.r
 
      # 2. generate p matrix
      sim.p <- matrix(rbeta(length(n) * nsim, r * p0, r * (1 - p0)), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rbinom(nrow(sim.p) * nsim, n, sim.p), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (!missing(reg.coef) & missing(covariates)) {
                   gbp(sim.z[, i], n, model = "br")
                 } else if (!missing(reg.coef) & !missing(covariates)) {
                   gbp(sim.z[, i], n, covariates, model = "br")
                 } else if (!missing(mean.PriorDist)) {
                   gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "br")
                 }
          a1 <- r * p0 + sim.z[, i]
          a0 <- r * (1 - p0) + n - sim.z[, i]
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pbeta(upp, a1, a0) - pbeta(low, a1, a0)
          coverage10[, i] <- ifelse(low <= sim.p[, i] & sim.p[, i] <= upp, 1, 0)
        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	
    }

  } else if(gbp.object$model == "pr") {

    if (only.gbp.result) {

      # 1. initial values
      if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(gbp.object$beta.new)
        lambda0 <- exp(temp.x %*% betas)
      } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
        betas <- as.vector(gbp.object$beta.new)
        lambda0 <- exp(temp.x %*% betas)
        X <- gbp.object$X
      } else if (!is.na(gbp.object$prior.mean)) {
        lambda0 <- gbp.object$prior.mean
      }
  
      n <- gbp.object$se
      r <- exp(-gbp.object$a.new)

      # 2. generate lambda matrix
      sim.lambda <- matrix(rgamma(length(n) * nsim, r * lambda0, r), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rpois(nrow(sim.lambda) * nsim, n * sim.lambda), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, model = "pr")
                 } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, X, model = "pr")
                 } else if (!is.na(gbp.object$prior.mean)) {
                   gbp(sim.z[, i], n, mean.PriorDist = lambda0, model = "pr")
                 }
          
          sh <- r * lambda0 + sim.z[, i]
          rt <- n + r
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pgamma(upp, sh, rt) - pgamma(low, sh, rt)
          coverage10[, i] <- ifelse(low <= sim.lambda[, i] & sim.lambda[, i] <= upp, 1, 0)
        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {
 
     # 1. initial values
      if (missing(A.or.r)) {
        print("(A.or.r, reg.coef, covariates (if any) or 
               (A.or.r, mean.PriorDist) should be designated")
        stop()
      } else if (missing(reg.coef)) {
        lambda0 <- mean.PriorDist
      } else if (missing(mean.PriorDist) & missing(covariates)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(reg.coef)
        lambda0 <- exp(temp.x %*% betas)
      } else if (missing(mean.PriorDist) & !missing(covariates)) { 
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), covariates))
        betas <- as.vector(reg.coef)
        lambda0 <- exp(temp.x %*% betas)
      }

      n <- gbp.object$se
      r <- A.or.r

      # 2. generate lambda matrix
      sim.lambda <- matrix(rgamma(length(n) * nsim, r * lambda0, r), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rpois(nrow(sim.lambda) * nsim, n * sim.lambda), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (missing(mean.PriorDist) & missing(covariates)) {
                   gbp(sim.z[, i], n, model = "pr")
                 } else if (missing(mean.PriorDist) & !missing(covariates)) {
                   gbp(sim.z[, i], n, covariates, model = "pr")
                 } else if (missing(reg.coef)) {
                   gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "pr")
                 }
          
          sh <- r * lambda0 + sim.z[, i]
          rt <- n + r
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pgamma(upp, sh, rt) - pgamma(low, sh, rt)
          coverage10[, i] <- ifelse(low <= sim.lambda[, i] & sim.lambda[, i] <= upp, 1, 0)
        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	
    }

  } else if (gbp.object$model == "gr") {	

     if (only.gbp.result) {
      # 1. initial values
      if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(gbp.object$beta.new)
        mu0 <- temp.x %*% betas
      } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
        betas <- as.vector(gbp.object$beta.new)
        mu0 <- temp.x %*% betas
        X <- gbp.object$X
      } else if (!is.na(gbp.object$prior.mean)) {
        mu0 <- gbp.object$prior.mean
      }
  
      se <- gbp.object$se
      A <- exp(gbp.object$a.new)

      # 2. generate mu matrix
      sim.mu <- matrix(rnorm(length(se) * nsim, mu0, sqrt(A)), nrow = length(se))

      # 3. generate y (data) matrix
      sim.y <- matrix(rnorm(nrow(sim.mu) * nsim, sim.mu, se), nrow = length(se))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                   gbp(sim.y[, i], se)
                 } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                   gbp(sim.y[, i], se, X)
                 } else if (!is.na(gbp.object$prior.mean)) {
                   gbp(sim.y[, i], se, mean.PriorDist = mu0)
                 }
          postmean <- mu0 * (se^2 / (se^2 + A)) + sim.y[, i] * (A / (se^2 + A))
          postsd <- sqrt(se^2 * (A / (se^2 + A)))
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pnorm(upp, postmean, postsd) - pnorm(low, postmean, postsd)
          coverage10[, i] <- ifelse(low <= sim.mu[, i] & sim.mu[, i] <= upp, 1, 0)
        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {
      # 1. initial values
      if (missing(A.or.r)) {
        print("(A.or.r, reg.coef, covariates (if any) or 
               (A.or.r, mean.PriorDist) should be designated")
        stop()
      } else if (!missing(mean.PriorDist)) {
        mu0 <- mean.PriorDist
      } else if (!missing(reg.coef) & missing(covariates)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(reg.coef)
        mu0 <- temp.x %*% betas
      } else if (!missing(reg.coef) & !missing(covariates)) { 
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), covariates))
        betas <- as.vector(reg.coef)
        mu0 <- temp.x %*% betas
      }
  
      se <- gbp.object$se
      A <- A.or.r
 
      # 2. generate mu matrix
      sim.mu <- matrix(rnorm(length(se) * nsim, mu0, sqrt(A)), nrow = length(se))

      # 3. generate y (data) matrix
      sim.y <- matrix(rnorm(nrow(sim.mu) * nsim, sim.mu, se), nrow = length(se))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (!missing(reg.coef) & missing(covariates)) {
                   gbp(sim.y[, i], se)
                 } else if (!missing(reg.coef) & !missing(covariates)) {
                   gbp(sim.y[, i], se, covariates)
                 } else if (!missing(mean.PriorDist)) {
                   gbp(sim.y[, i], se, mean.PriorDist = mean.PriorDist)
                 }
          postmean <- mu0 * (se^2 / (se^2 + A)) + sim.y[, i] * (A / (se^2 + A))
          postsd <- sqrt(se^2 * (A / (se^2 + A)))
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pnorm(upp, postmean, postsd) - pnorm(low, postmean, postsd)
          coverage10[, i] <- ifelse(low <= sim.mu[, i] & sim.mu[, i] <= upp, 1, 0)
        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	
    }
  }

  # average coverage probability
  result <- round(rowMeans(coverageRB, na.rm = TRUE), 3)
  avr.cov <- round(mean(result), 3)
  min.cov <- round(min(result), 3)

  result2 <- round(rowMeans(coverage10, na.rm = TRUE), 3)
  avr.cov2 <- round(mean(result2), 3)
  min.cov2 <- round(min(result2), 3)

  # plotting coverage graph
  par(xaxs = "r", yaxs = "r", mai = c(1, 0.6, 1, 0.3))
  plot(1 : length(gbp.object$se), result, ylim = c(0.7, 1), type = "l", col = 2, ylab = "",
       xlab = "Group_i, i = 1, ..., k", main = "Estimated Coverage Probability for Each Group",
       lwd = 3, lty = 1)
  abline(h = 0.95)
  points(1 : length(gbp.object$se), result2, type = "l", lty = 2, col = 4, lwd = 2)
  if (is.na(gbp.object$prior.mean)){
    if (gbp.object$model == "gr") {
      legend("bottomleft", c("Red Line: Rao-Blackwellized",
                             "Blue Dotted Line: (Unbiased)",
                             paste("A =", round(A, 2)), 
                             paste("beta", 0 : (length(betas) - 1), "=", round(betas, 3)), 
                             paste("AvgCoverage =", avr.cov, "(", avr.cov2, ")"), 
                             paste("MinCoverage =", min.cov, "(", min.cov2, ")")))
    } else {
      legend("bottomleft", c("Red Line: Rao-Blackwellized",
                             "Blue Dotted Line: (Unbiased)",
                             paste("r =", round(r, 2)), 
                             paste("beta", 0 : (length(betas) - 1), "=", round(betas, 3)), 
                             paste("AvgCoverage =", avr.cov, "(", avr.cov2, ")"),
                             paste("MinCoverage =", min.cov, "(", min.cov2, ")")))
    }
  } else {  # if prior mean is assigned
    if (gbp.object$model == "gr") {
      legend("bottomleft", c("Red Line: Rao-Blackwellized",
                             "Blue Dotted Line: (Unbiased)",
                             paste("A =", round(A, 2)), 
                             paste("AvgCoverage =", avr.cov, "(", avr.cov2, ")"), 
                             paste("MinCoverage =", min.cov, "(", min.cov2, ")")))
    } else {
      legend("bottomleft", c("Red Line: Rao-Blackwellized",
                             "Blue Dotted Line: (Unbiased)",
                             paste("r =", round(r, 2)), 
                             paste("AvgCoverage =", avr.cov, "(", avr.cov2, ")"),
                             paste("MinCoverage =", min.cov, "(", min.cov2, ")")))
    }
  }


  # print output
  output <- list(coverageRB = result, coverage10 = result2, 
                 average.coverageRB = avr.cov, average.coverage10 = avr.cov2, 
                 minimum.coverageRB = min.cov, minimum.coverage10 = min.cov2, 
                 raw.resultRB = coverageRB, raw.result10 = coverage10)
  return(output)
}
