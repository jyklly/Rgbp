gbp <- function(x, w, covariates, mean.PriorDist, model, intercept, Alpha) UseMethod("gbp")

gbp.default <- function(x, w, covariates, mean.PriorDist, model = "gr", intercept = TRUE, Alpha = 0.95) {
  res <- switch(model, 
           gr = gr(x, w, X = covariates, mu = mean.PriorDist, Alpha = Alpha, intercept = intercept), 
           br = br(x, w, X = covariates, prior.mean = mean.PriorDist, intercept = intercept, Alpha = Alpha), 
           pr = pr(x, w, X = covariates, prior.mean = mean.PriorDist, intercept = intercept, Alpha = Alpha))
  
  class(res) <- "gbp"	
  res
}

print.gbp <- function(x, sort = TRUE, ...) {
  
  if (any(is.na(x$prior.mean)) & !identical(x$X, NA)) {

	cova <- as.matrix(x$X)
	colnames(cova) <- paste(1 : dim(cova)[2])
    
    if (x$model == "gr") {
      temp <- data.frame(obs.mean = x$sample.mean, se = x$se, x = cova, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(obs.mean = x$sample.mean, n = x$se, x = cova, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }

  } else if (any(is.na(x$prior.mean)) & identical(x$X, NA)) { 
  # if there are neither prior.mean and X assigned

    if (x$model == "gr") {
      temp <- data.frame(obs.mean = x$sample.mean, se = x$se, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(obs.mean = x$sample.mean, n = x$se, 
                         prior.mean = x$prior.mean.hat, shrinkage = x$shrinkage, 
                         intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }


  } else if (!any(is.na(x$prior.mean))) {
    if (x$model == "gr") {
      temp <- data.frame(obs.mean = x$sample.mean, se = x$se, 
                         prior.mean = x$prior.mean, shrinkage = x$shrinkage, 
                         intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(obs.mean = x$sample.mean, n = x$se, 
                         prior.mean = x$prior.mean, shrinkage = x$shrinkage, 
                         intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }
  }
  
  if (sort == TRUE) {
    temp <- temp[order(temp[, 2]), ]
  }

  temp.mean <- colMeans(temp)
  temp <- data.frame(rbind(temp, temp.mean), row.names = c(rownames(temp), "colMeans"))

  cat("Summary for whole observations: \n")
  cat("\n")
  print(round(temp, 3))
}

summary.gbp <- function(object, ...) {

  if (any(is.na(object$prior.mean)) & !identical(object$X, NA)) {

	cova <- as.matrix(object$X)
	colnames(cova) <- paste(1 : dim(cova)[2])
    
    if (object$model == "gr") {
      temp <- data.frame(obs.mean = object$sample.mean, se = object$se, x = cova, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         intv.low = object$post.intv.low, post.mean = object$post.mean, 
                         intv.upp = object$post.intv.upp, post.sd = object$post.sd)
    } else {
      temp <- data.frame(obs.mean = object$sample.mean, n = object$se, x = cova, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         intv.low = object$post.intv.low, post.mean = object$post.mean, 
                         intv.upp = object$post.intv.upp, post.sd = object$post.sd)
    }

  } else if (any(is.na(object$prior.mean)) & identical(object$X, NA)) { 
  # if there are neither prior.mean and X assigned

    if (object$model == "gr") {
      temp <- data.frame(obs.mean = object$sample.mean, se = object$se, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         intv.low = object$post.intv.low, post.mean = object$post.mean, 
                         intv.upp = object$post.intv.upp, post.sd = object$post.sd)
    } else {
      temp <- data.frame(obs.mean = object$sample.mean, n = object$se, 
                         prior.mean = object$prior.mean.hat, shrinkage = object$shrinkage, 
                         intv.low = object$post.intv.low, post.mean = object$post.mean, 
                         intv.upp = object$post.intv.upp, post.sd = object$post.sd)
    }


  } else if (!any(is.na(object$prior.mean))) {
    if (object$model == "gr") {
      temp <- data.frame(obs.mean = object$sample.mean, se = object$se, 
                         prior.mean = object$prior.mean, shrinkage = object$shrinkage, 
                         intv.low = object$post.intv.low, post.mean = object$post.mean, 
                         intv.upp = object$post.intv.upp, post.sd = object$post.sd)
    } else {
      temp <- data.frame(obs.mean = object$sample.mean, n = object$se, 
                         prior.mean = object$prior.mean, shrinkage = object$shrinkage,
                         intv.low = object$post.intv.low, post.mean = object$post.mean, 
                         intv.upp = object$post.intv.upp, post.sd = object$post.sd)
    }
  }



  if (min(object$se) == max(object$se)) {  # if se or n are all the same

    temp2 <- temp[order(temp$obs.mean), ]

    if (length(object$se) %% 2) {  # if number of groups is odd
      summary.table <- temp2[c(1, (dim(temp2)[1] + 1) / 2, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (number.of.medians == 1) {
        row.names(summary.table) <- c("Unit w/ min(obs.mean)", 
                                      "Unit w/ median(obs.mean)", 
                                      "Unit w/ max(obs.mean)")
      } else {  # if there are more than one median
        row.names(summary.table) <- c("Unit w/ min(obs.mean)", 
                                     paste("Unit w/ median(obs.mean)", 1 : number.of.medians, sep = ""),
                                     "Unit w/ max(obs.mean)")
      }
    } else {  # if number of groups is even
      summary.table <- temp2[c(1, dim(temp2)[1] / 2, dim(temp2)[1] / 2 + 1, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (number.of.medians == 1) {
        row.names(summary.table) <- c("Unit w/ min(obs.mean)", 
                                      "Unit w/ median(obs.mean)", 
                                      "Unit w/ max(obs.mean)")
      } else {  # if there are more than one median
        row.names(summary.table) <- c("Unit w/ min(obs.mean)", 
                                     paste("Unit w/ median(obs.mean)", 1 : number.of.medians, sep = ""),
                                     "Unit w/ max(obs.mean)")
      }
    }
  } else { # if n or se are different from each group

    temp2 <- temp[order(temp[, 2], temp$obs.mean), ]

    if (length(object$se) %% 2) {  # if number of groups is odd
      summary.table <- temp2[c(1, (dim(temp2)[1] + 1) / 2, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (object$model == "gr") {
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Unit w/ min(se)", "Unit w/ median(se)", "Unit w/ max(se)")
        } else {
          row.names(summary.table) <- c("Unit w/ min(se)", 
                                       paste("Unit w/ median(se)", 1 : number.of.medians, sep = ""),
                                       "Unit w/ max(se)")
        }
      } else {  # if model is not "gr"
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Unit w/ min(n)", "Unit w/ median(n)", "Unit w/ max(n)")
        } else {
          row.names(summary.table) <- c("Unit w/ min(n)", 
                                       paste("Unit w/ median(n)", 1 : number.of.medians, sep = ""),
                                       "Unit w/ max(n)")
        }
      }

    } else {   # if number of groups is even
      summary.table <- temp2[c(1, dim(temp2)[1] / 2, dim(temp2)[1] / 2 + 1, dim(temp2)[1]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      if (object$model == "gr") {
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Unit w/ min(se)", "Unit w/ median(se)", "Unit w/ max(se)")
        } else {
          row.names(summary.table) <- c("Unit w/ min(se)", 
                                       paste("Unit w/ median(se)", 1 : number.of.medians, sep = ""),
                                       "Unit w/ max(se)")
        }
      } else {  # if model is not "gr"
        if (number.of.medians == 1) {
          row.names(summary.table) <- c("Unit w/ min(n)", "Unit w/ median(n)", "Unit w/ max(n)")
        } else {
          row.names(summary.table) <- c("Unit w/ min(n)", 
                                       paste("Unit w/ median(n)", 1 : number.of.medians, sep = ""),
                                       "Unit w/ max(n)")
        }
      }
    }
  }

  temp.mean <- colMeans(temp)
  
  if (object$model == "gr") {
    tamp.mean[2] <- round(tamp.mean[2], 1)
    summary.table[, 2] <- round(summary.table[, 2], 1)
  } else {
    tamp.mean[2] <- round(tamp.mean[2])
    summary.table[, 2] <- round(summary.table[, 2])
  }

  summary.table <- data.frame(rbind(summary.table, temp.mean), 
                              row.names = c(rownames(summary.table), "Overall Mean"))

  post.mode.alpha <- object$a.new
  post.sd.alpha <- sqrt(object$a.var)
  result2 <- data.frame(post.mode.alpha, post.sd.alpha)

  if (any(is.na(object$prior.mean))) {
    estimate <- as.vector(object$beta.new)
    if (object$intercept == TRUE) {
      names(estimate) <- paste("beta", 0 : (length(estimate) - 1), sep = "")
    } else {
      names(estimate) <- paste("beta", 1 : (length(estimate)), sep = "")
    }
    se <- as.vector(sqrt(diag(object$beta.var)))
    z.val <- estimate / se
    p.val <- 2 * pnorm(-abs(z.val))
    beta.result <- data.frame(estimate, se, z.val, p.val)
    res <- list(main = summary.table, sec.var = result2, reg = beta.result)
  } else {
    res <- list(main = summary.table, sec.var = result2, reg = NA)
  }
  
  class(res) <- "summary.gbp"
  res
}


print.summary.gbp <- function(x, ...) {
  if (identical(x$reg, NA)) {
    cat("Main summary:\n")
    cat("\n")
    print(round(x$main, 3))
    cat("\n")
    cat("\n")
    cat("Second-level Variance Component Estimation Summary:\n")
    cat("alpha = log(A) for Gaussian and log(1/r) for Binomial and Poisson data:\n")
    cat("\n")
    print(round(x$sec.var, 3))
  } else {
    cat("Main summary:\n")
    cat("\n")
    print(round(x$main, 3))
    cat("\n")
    cat("\n")
    cat("Second-level Variance Component Estimation Summary:\n")
    cat("alpha = log(A) for Gaussian and log(1/r) for Binomial and Poisson data:\n")
    cat("\n")
    print(round(x$sec.var, 3))
    cat("\n")
    cat("\n")
    cat("Regression Summary:\n")
    cat("\n")
    print(round(x$reg, 3))
  }
}

plot.gbp <- function(x, sort = TRUE, ...) {

  y <- x$sample.mean
  se <- x$se
  if (any(is.na(x$prior.mean))) {
    pr.m <- x$prior.mean.hat
  } else {
    pr.m <- x$prior.mean
  }
  po.m <- x$post.mean
  po.sd <- x$post.sd
  po.low <- x$post.intv.low
  po.upp <- x$post.intv.upp
  xlabel <- "Indices (Groups) by the order of data input"

  if (sort == TRUE) {
    temp.data <- as.data.frame(cbind(y, se, pr.m, po.m, po.sd, po.low, po.upp))
    temp.data <- temp.data[order(temp.data$se), ]
    y <- temp.data$y
    se <- temp.data$se
    pr.m <- temp.data[, 3]
    po.m <- temp.data[, 4]
    po.sd <- temp.data[, 5]
    po.low <- temp.data[, 6]
    po.upp <- temp.data[, 7]
    if (x$model == "gr") {
      xlabel <- "Indices (Groups) sorted by the order of se"
    } else {
      xlabel <- "Indices (Groups) sorted by the order of n"
    }
  }

  cx <- (mean(log(se + 2)) + min(log(se + 2))) / 2
  index <- 1 : length(se)
  ylim.low <- ifelse(min(po.low, y) >= 0, 0.8 * min(po.low, y), 1.2 * min(po.low, y))
  ylim.upp <- ifelse(max(po.upp, y) >= 0, 1.2 * max(po.upp, y), 0.8 * max(po.upp, y))

  par(mfrow = c(1, 2), xaxs = "r", yaxs = "r", mai = c(1, 0.6, 1, 0.3), las = 1, ps = 13,oma=c(0,5,0,0))
  
  sqrtV <- se
  sdlens <- sqrtV / max(sqrtV)
  postlens <- po.sd / max(po.sd)
  xmin <- min(c(y, po.m, pr.m))
  xmax <- max(c(y, po.m, pr.m))
  sunflowerplot(rep(4, length(y)) ~ y, ylim = c(-1, 5), xlim = c(xmin - abs(xmin) * 0.1, 
                xmax + abs(xmax) * 0.1), yaxt = "n", col.lab = "white", main = "Shrinkage Plot", pch=1,cex=1)
  
  if (length(unique(pr.m)) == 1) {
    abline(v = pr.m, col = 4)
  } else {
    points(pr.m, rep(0,length(pr.m)), col = 4, pch = "-", cex = 2)
  }
  
  sunflowerplot(rep(0, length(y)) ~ po.m, add = TRUE,col="red",cex=1,pch=16)
  abline(h = 4)
  abline(h = 0)
 # axis(2, c(0, 4), c(expression(hat(theta)), expression(bar(y))), cex.axis = 1.1)
  sapply(1 : length(y), function(i) {
    lines(c(y[i], po.m[i]), c(4, 0))
    lines(c(y[i], y[i] + sdlens[i] * sd(y) * 0.4), c(4, 4 + sdlens[i]), col = "darkviolet")
    ##posterior variance lines
    lines(c(po.m[i] - postlens[i] * sd(y) * 0.4, po.m[i]), c(0 - postlens[i], 0), col = "darkgreen")
    xcord <- ((4 * po.m[i] / (y[i] - po.m[i]) - 4 * po.m / (y - po.m)) / 
              (4 / (y[i] - po.m[i]) - 4 / (y - po.m)))
    ycord <- 4 / (y - po.m) * xcord - 4 / (y - po.m) * po.m
    coords <- subset(cbind(xcord, ycord), ycord > 0 & ycord < 4)
    points(coords, pch=0)
  })

  
  plot(index, po.m, ylim = c(ylim.low, ylim.upp), xlab = xlabel, ylab = expression(theta),
       main = "100*Alpha% Intervals", cex = log(se + 2) / cx, col = 2, pch = 19)
  sapply(1 : length(y), function(j) {
    lines(rep(index[j], 2), c(po.low[j], po.upp[j]), lwd = 0.5)
  })
  points(index, po.low, cex = 1.5, pch = "-")
  points(index, po.upp, cex = 1.5, pch = "-")
  points(index, y, cex = log(se + 2) / cx)
  if (length(unique(pr.m)) == 1) {
    abline(h = pr.m, col = 4)
  } else {
    points(index, pr.m, col = 4, pch = "-", cex = 2)
  }

  se.or.n <- switch(x$model, "gr" = "standard error", "br" = "n", "pr" = "n")
  par(new=TRUE,mfrow=c(1,1))
  par(xpd=NA)
  plot(1, type="n", axes=F, xlab="", ylab="")
  legend("topleft", pch = c(19, 1, NA, NA, NA,0), col = c(2, 1, 4,"darkviolet", "darkgreen",1), 
         lwd = c(NA, NA, 2, 2, 2), 
         c("posterior mean", "sample mean", "prior mean", se.or.n, "posterior sd", "crossover"),
         seg.len = 0.5, bty = "n", inset = c(-0.175, 0))
 # legend("topleft", pch = c(19, 1, NA,NA), col = c(2, 1, 4,"darkviolet"), lwd = c(NA, NA, 2,2), c("posterior mean", "sample mean", "prior mean","se"), seg.len = 0.5, bty = "n")
  par(xpd = FALSE)


}
