
coverage <- function(gbp.object, A.or.r, reg.coef, mean.PriorDist, nsim = 100) {

  # Rao-Blackwellized criterion
  coverageRB <- matrix(NA, nrow = length(gbp.object$se), ncol = nsim)

  # 1-0 criterion that is 1 if interval includes true parameter, 0 if not
  coverageS <- matrix(NA, nrow = length(gbp.object$se), ncol = nsim)

  if (missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
    only.gbp.result <- TRUE
  } else {
    only.gbp.result <- FALSE
  }

######
  if (sum(is.na(gbp.object$weight)) == 1) {
    IS <- FALSE
  } else {
    IS <- TRUE
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
      r <- as.numeric(exp(-gbp.object$a.new))
      priormeanused <- p0

      # 2. generate p matrix
      sim.p <- matrix(rbeta(length(n) * nsim, r * p0, r * (1 - p0)), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rbinom(nrow(sim.p) * nsim, n, sim.p), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
######
          if(IS == 0) {
            out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, model = "binomial", Alpha = gbp.object$Alpha)
                   } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, X, model = "binomial", Alpha = gbp.object$Alpha)
                   } else if (!is.na(gbp.object$prior.mean)) {
                     gbp(sim.z[, i], n, mean.PriorDist = p0, model = "binomial", Alpha = gbp.object$Alpha)
                   }
          } else {
            out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, model = "binomial", Alpha = gbp.object$Alpha,
                         n.IS = length(gbp.object$weight))
                   } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, X, model = "binomial", Alpha = gbp.object$Alpha,
                         n.IS = length(gbp.object$weight))
                   } else if (!is.na(gbp.object$prior.mean)) {
                     gbp(sim.z[, i], n, mean.PriorDist = p0, model = "binomial", Alpha = gbp.object$Alpha, 
                         n.IS = length(gbp.object$weight))
                   }
          }


          a1 <- r * p0 + sim.z[, i]
          a0 <- r * (1 - p0) + n - sim.z[, i]
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pbeta(upp, a1, a0) - pbeta(low, a1, a0)
          coverageS[, i] <- ifelse(low <= sim.p[, i] & sim.p[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {

      # 1. initial values
      if (missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        p0 <- mean.PriorDist
        r <- as.numeric(exp(-gbp.object$a.new))
      } else if (missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because there is no covariate.")
          stop()
        } else if (identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(reg.coef)
        } else {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(reg.coef)
        }
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        r <- as.numeric(exp(-gbp.object$a.new))
      } else if (!missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          p0 <- gbp.object$prior.mean
        } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(gbp.object$beta.new)
          p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(gbp.object$beta.new)
          p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        }
        r <- A.or.r
      } else if (missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      } else if (!missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        p0 <- mean.PriorDist
        r <- A.or.r
      } else if (!missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because second-level mean is known in the gbp object to begin with.")
          stop()
        } else if (identical(gbp.object$prior.mean, NA)) {
          if (identical(gbp.object$X, NA)) {
            temp.x <- as.matrix(rep(1, length(gbp.object$se)))
            betas <- as.vector(reg.coef)
          } else {
            temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
            betas <- as.vector(reg.coef)
          }
          p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
          r <- A.or.r
        }
      } else if (!missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      }

      n <- gbp.object$se
      priormeanused <- p0
 
      # 2. generate p matrix
      sim.p <- matrix(rbeta(length(n) * nsim, r * p0, r * (1 - p0)), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rbinom(nrow(sim.p) * nsim, n, sim.p), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
######
          if(IS == 0) {
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "binomial", 
                         Alpha = gbp.object$Alpha)
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.z[, i], n, mean.PriorDist = gbp.object$prior.mean, model = "binomial", 
                           Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", Alpha = gbp.object$Alpha)
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", Alpha = gbp.object$Alpha)
                     }
                   }
          } else {
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "binomial", 
                         Alpha = gbp.object$Alpha, n.IS = length(gbp.object$weight))
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.z[, i], n, mean.PriorDist = gbp.object$prior.mean, model = "binomial", 
                           Alpha = gbp.object$Alpha, n.IS = length(gbp.object$weight))
                     } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", Alpha = gbp.object$Alpha,
                           n.IS = length(gbp.object$weight))
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", Alpha = gbp.object$Alpha,
                           n.IS = length(gbp.object$weight))
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", Alpha = gbp.object$Alpha,
                           n.IS = length(gbp.object$weight))
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", Alpha = gbp.object$Alpha,
                           n.IS = length(gbp.object$weight))
                     }
                   }
          }

          a1 <- r * p0 + sim.z[, i]
          a0 <- r * (1 - p0) + n - sim.z[, i]
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pbeta(upp, a1, a0) - pbeta(low, a1, a0)
          coverageS[, i] <- ifelse(low <= sim.p[, i] & sim.p[, i] <= upp, 1, 0)

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
      priormeanused <- lambda0

      # 2. generate lambda matrix
      sim.lambda <- matrix(rgamma(length(n) * nsim, r * lambda0, r), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rpois(nrow(sim.lambda) * nsim, n * sim.lambda), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, model = "poisson", Alpha = gbp.object$Alpha)
                 } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, X, model = "poisson", Alpha = gbp.object$Alpha)
                 } else if (!is.na(gbp.object$prior.mean)) {
                   gbp(sim.z[, i], n, mean.PriorDist = lambda0, model = "poisson", Alpha = gbp.object$Alpha)
                 }
          
          sh <- r * lambda0 + sim.z[, i]
          rt <- n + r
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pgamma(upp, sh, rt) - pgamma(low, sh, rt)
          coverageS[, i] <- ifelse(low <= sim.lambda[, i] & sim.lambda[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {

      # 1. initial values
      if (missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        lambda0 <- mean.PriorDist
        r <- exp(-gbp.object$a.new)
      } else if (missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because there is no covariate.")
          stop()
        } else if (identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(reg.coef)
        } else {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(reg.coef)
        }
        lambda0 <- exp(temp.x %*% betas)
        r <- exp(-gbp.object$a.new)
      } else if (!missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          lambda0 <- gbp.object$prior.mean
        } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(gbp.object$beta.new)
          lambda0 <- exp(temp.x %*% betas)
        } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(gbp.object$beta.new)
          lambda0 <- exp(temp.x %*% betas)
        }
        r <- A.or.r
      } else if (missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      } else if (!missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        lambda0 <- mean.PriorDist
        r <- A.or.r
      } else if (!missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because second-level mean is known in the gbp object to begin with.")
          stop()
        } else if (identical(gbp.object$prior.mean, NA)) {
          if (identical(gbp.object$X, NA)) {
            temp.x <- as.matrix(rep(1, length(gbp.object$se)))
            betas <- as.vector(reg.coef)
          } else {
            temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
            betas <- as.vector(reg.coef)
          }
          lambda0 <- exp(temp.x %*% betas)
          r <- A.or.r
        }
      } else if (!missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      }

      n <- gbp.object$se
      priormeanused <- lambda0
 
      # 2. generate lambda matrix
      sim.lambda <- matrix(rgamma(length(n) * nsim, r * lambda0, r), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rpois(nrow(sim.lambda) * nsim, n * sim.lambda), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "poisson", 
                         Alpha = gbp.object$Alpha)
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.z[, i], n, mean.PriorDist = gbp.object$prior.mean, model = "poisson", 
                           Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "poisson", Alpha = gbp.object$Alpha)
                     } else if (!identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "poisson", Alpha = gbp.object$Alpha)
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "poisson", Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "poisson", Alpha = gbp.object$Alpha)
                     }
                   }
          
          sh <- r * lambda0 + sim.z[, i]
          rt <- n + r
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pgamma(upp, sh, rt) - pgamma(low, sh, rt)
          coverageS[, i] <- ifelse(low <= sim.lambda[, i] & sim.lambda[, i] <= upp, 1, 0)

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
      priormeanused <- mu0

      # 2. generate mu matrix
      sim.mu <- matrix(rnorm(length(se) * nsim, mu0, sqrt(A)), nrow = length(se))

      # 3. generate y (data) matrix
      sim.y <- matrix(rnorm(nrow(sim.mu) * nsim, sim.mu, se), nrow = length(se))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                   gbp(sim.y[, i], se, Alpha = gbp.object$Alpha)
                 } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                   gbp(sim.y[, i], se, X, Alpha = gbp.object$Alpha)
                 } else if (!is.na(gbp.object$prior.mean)) {
                   gbp(sim.y[, i], se, mean.PriorDist = mu0, Alpha = gbp.object$Alpha)
                 }
          postmean <- mu0 * (se^2 / (se^2 + A)) + sim.y[, i] * (A / (se^2 + A))
          postsd <- sqrt(se^2 * (A / (se^2 + A)))
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pnorm(upp, postmean, postsd) - pnorm(low, postmean, postsd)
          coverageS[, i] <- ifelse(low <= sim.mu[, i] & sim.mu[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {

      # 1. initial values
      if (missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        mu0 <- mean.PriorDist
        A <- exp(gbp.object$a.new)
      } else if (missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because there is no covariate.")
          stop()
        } else if (identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(reg.coef)
        } else {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(reg.coef)
        }
        mu0 <- temp.x %*% betas
        A <- exp(gbp.object$a.new)
      } else if (!missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          mu0 <- gbp.object$prior.mean
        } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(gbp.object$beta.new)
          mu0 <- temp.x %*% betas
        } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(gbp.object$beta.new)
          mu0 <- temp.x %*% betas
        }
        A <- A.or.r
      } else if (missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      } else if (!missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        mu0 <- mean.PriorDist
        A <- A.or.r
      } else if (!missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because second-level mean is known in the gbp object to begin with.")
          stop()
        } else if (identical(gbp.object$prior.mean, NA)) {
          if (identical(gbp.object$X, NA)) {
            temp.x <- as.matrix(rep(1, length(gbp.object$se)))
            betas <- as.vector(reg.coef)
          } else {
            temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
            betas <- as.vector(reg.coef)
          }
          mu0 <- temp.x %*% betas
          A <- A.or.r
        }
      } else if (!missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      }

      se <- gbp.object$se
      priormeanused <- mu0
 
      # 2. generate mu matrix
      sim.mu <- matrix(rnorm(length(se) * nsim, mu0, sqrt(A)), nrow = length(se))

      # 3. generate y (data) matrix
      sim.y <- matrix(rnorm(nrow(sim.mu) * nsim, sim.mu, se), nrow = length(se))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.y[, i], se, mean.PriorDist = mean.PriorDist, model = "gaussian", 
                         Alpha = gbp.object$Alpha)
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.y[, i], se, mean.PriorDist = gbp.object$prior.mean, model = "gaussian", 
                           Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, model = "gaussian", Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, gbp.object$X, model = "gaussian", Alpha = gbp.object$Alpha)
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, model = "gaussian", Alpha = gbp.object$Alpha)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, gbp.object$X, model = "gaussian", Alpha = gbp.object$Alpha)
                     }
                   }

          postmean <- mu0 * (se^2 / (se^2 + A)) + sim.y[, i] * (A / (se^2 + A))
          postsd <- sqrt(se^2 * (A / (se^2 + A)))
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pnorm(upp, postmean, postsd) - pnorm(low, postmean, postsd)
          coverageS[, i] <- ifelse(low <= sim.mu[, i] & sim.mu[, i] <= upp, 1, 0)

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
  se.cov <- round(sqrt(apply(coverageRB, 1, var) / nsim), 4)
  result2 <- round(rowMeans(coverageS, na.rm = TRUE), 3)
  avr.cov2 <- round(mean(result2), 3)
  se.cov2 <- round(sqrt(apply(coverageS, 1, var) / nsim), 4)
  effective.n <- nsim - sum(is.na(coverageS[1, ]))

  # plotting coverage graph
  par(xaxs = "r", yaxs = "r", mai = c(1, 0.6, 1, 0.3))
  n.units <- length(gbp.object$se)
  plot(1 : length(gbp.object$se), result, ylim = c(0.6, 1), type = "l", col = 2, ylab = "",
       xlab = paste("Unit_j", ", ", "j = 1, ...,", n.units), 
       main = "Estimated Coverage Probability for Each Unit",
       lwd = 3, lty = 1)
  abline(h = gbp.object$Alpha)

  if (is.na(gbp.object$prior.mean) & missing(mean.PriorDist)) {
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Gaussian"), 
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True A =", round(A, 2)), 
                             paste("Given True beta", 0 : (length(betas) - 1), "=", round(betas, 3)), 
                             paste("Overall Coverage =", avr.cov)))
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial", "Poisson")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True r =", round(r, 2)), 
                             paste("Given True beta", 0 : (length(betas) - 1), "=", round(betas, 3)), 
                             paste("Overall Coverage =", avr.cov)))
    }

  } else if (is.na(gbp.object$prior.mean) & !missing(mean.PriorDist)) {
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Gaussian"),
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True A =", round(A, 2)), 
                             paste("Known Prior Mean: ", round(priormeanused, 2)), 
                             paste("Overall Coverage =", avr.cov)))
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial", "Poisson")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True r =", round(r, 2)), 
                             paste("Known Prior Mean: ", round(priormeanused, 2)), 
                             paste("Overall Coverage =", avr.cov)))
    }

  } else if (!is.na(gbp.object$prior.mean) & !missing(mean.PriorDist)) {  # if prior mean is assigned
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Gaussian"),
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True A =", round(A, 2)), 
                             paste("Known Prior Mean: ", round(priormeanused, 2)), 
                             paste("Overall Coverage =", avr.cov)))
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial", "Poisson")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True r =", round(r, 2)), 
                             paste("Known Prior Mean: ", round(priormeanused, 2)), 
                             paste("Overall Coverage =", avr.cov)))
    }
  } else if (!is.na(gbp.object$prior.mean) & missing(mean.PriorDist)) {  # if prior mean is assigned
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Gaussian"),
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True A =", round(A, 2)), 
                             paste("Known Prior Mean: ", round(priormeanused, 2)), 
                             paste("Overall Coverage =", avr.cov)))
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial", "Poisson")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red Line: RB coverage estimates",
                             paste("# of Simulations per Unit: ", effective.n),
                             paste("Given True r =", round(r, 2)), 
                             paste("Known Prior Mean: ", round(priormeanused, 2)), 
                             paste("Overall Coverage =", avr.cov)))
    }
  }

  # print output
  output <- list(coverageRB = result, coverageS = result2, 
                 average.coverageRB = avr.cov, average.coverageS = avr.cov2, 
                 se.coverageRB = se.cov, se.coverageS = se.cov2, 
                 raw.resultRB = coverageRB, raw.resultS = coverageS)
  return(output)
}
