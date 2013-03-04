gbp <- function(x, ...) UseMethod("gbp")

gbp.default <- function(x, y, covariates, mu, model = "gr", intercept = T, Alpha = 0.95, ...) {
  res <- switch(model, 
                gr = gr(x, y, X = covariates, mu, Alpha), 
                br = br(x, y, X = covariates, prior.mean = mu, intercept = intercept, Alpha = Alpha), 
                pr = pr(x, y, X = covariates, prior.mean = mu, intercept = intercept, Alpha = Alpha))
  
  class(res) <- "gbp"	
  res
}

print.gbp <- function(x, ...) {
  
  if (is.na(x$prior.mean) & !identical(x$X, NA)) {

	cova <- as.matrix(x$X)
	colnames(cova) <- paste(1 : dim(cova)[2])
    
    if (x$model == "gr") {
      temp <- data.frame(sample.mean = x$sample.mean, se = x$se, x = cova, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(sample.mean = x$sample.mean, n = x$se, x = cova, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }

  } else if (is.na(x$prior.mean) & identical(x$X, NA)) { # if there are neither prior.mean and X assigned

    if (x$model == "gr") {
      temp <- data.frame(sample.mean = x$sample.mean, se = x$se, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(sample.mean = x$sample.mean, n = x$se, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }


  } else if (!is.na(x$prior.mean)) {
    if (x$model == "gr") {
      temp <- data.frame(sample.mean = x$sample.mean, se = x$se, prior.mean = x$prior.mean, 
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(sample.mean = x$sample.mean, n = x$se, prior.mean = x$prior.mean, 
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }
  }

  cat("Summary for whole observations:\n")
  cat("\n")
  print(round(temp,3))
}

summary.gbp <- function(object, ...) {

  x <- object
  if (is.na(x$prior.mean) & !identical(x$X, NA)) {

	cova <- as.matrix(x$X)
	colnames(cova) <- paste(1 : dim(cova)[2])
    
    if (x$model == "gr") {
      temp <- data.frame(sample.mean = x$sample.mean, se = x$se, x = cova, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(sample.mean = x$sample.mean, n = x$se, x = cova, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }

  } else if (is.na(x$prior.mean) & identical(x$X, NA)) { # if there are neither prior.mean and X assigned

    if (x$model == "gr") {
      temp <- data.frame(sample.mean = x$sample.mean, se = x$se, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(sample.mean = x$sample.mean, n = x$se, prior.mean = x$prior.mean.hat,
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }


  } else if (!is.na(x$prior.mean)) {
    if (x$model == "gr") {
      temp <- data.frame(sample.mean = x$sample.mean, se = x$se, prior.mean = x$prior.mean, 
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    } else {
      temp <- data.frame(sample.mean = x$sample.mean, n = x$se, prior.mean = x$prior.mean, 
                         shrinkage = x$shrinkage, sd.shrinkage = x$sd.shrinkage, 
                         post.intv.low = x$post.intv.low, post.mean = x$post.mean, 
                         post.intv.upp = x$post.intv.upp, post.sd = x$post.sd)
    }
  }
  if (min(x$se) == max(x$se)) {
    if (length(x$se) %% 2) {
      summary.index <- (x$sample.mean == min(x$sample.mean) | 
                        x$sample.mean == max(x$sample.mean) | 
                        x$sample.mean == sort(x$sample.mean)[(length(x$se) + 1) / 2])
      summary.table <- temp[summary.index, ]
      summary.table <- summary.table[order(summary.table$sample.mean), ]
      number.of.medians <- dim(summary.table)[1] - 2
      row.names(summary.table) <- c("Group w/ min(sample.mean)", 
                                    paste("Group w/ median(sample.mean)", 1 : number.of.medians, sep = ""),
                                    "Group w/ max(sample.mean)")
    } else { 
      summary.index <- (x$sample.mean == min(x$sample.mean) | 
                        x$sample.mean == max(x$sample.mean) | 
                        x$sample.mean == sort(x$sample.mean)[length(x$se) / 2] |
                        x$sample.mean == sort(x$sample.mean)[length(x$se) / 2 + 1])
      summary.table <- temp[summary.index, ]
      summary.table <- summary.table[order(summary.table$sample.mean), ]
      number.of.medians <- dim(summary.table)[1] - 2
      row.names(summary.table) <- c("Group w/ min(sample.mean)", 
                                    paste("Group w/ median(sample.mean)", 1 : number.of.medians, sep = ""),
                                    "Group w/ max(sample.mean)")
    }

  } else { # if n is different from each group
    if (length(x$se) %% 2) {
      summary.index <- (x$se == min(x$se) | x$se == max(x$se) | x$se == sort(x$se)[(length(x$se) + 1) / 2])
      summary.table <- temp[summary.index, ]
      summary.table <- summary.table[order(summary.table[, 2]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      row.names(summary.table) <- c("Group w/ min(n)", 
                                    paste("Group w/ median(n)", 1 : number.of.medians, sep = ""),
                                    "Group w/ max(n)")
    } else { 
      summary.index <- (x$se == min(x$se) | x$se == max(x$se) | x$se == sort(x$se)[length(x$se) / 2] |
                        x$se == sort(x$se)[length(x$se) / 2 + 1])
      summary.table <- temp[summary.index, ]
      summary.table <- summary.table[order(summary.table[, 2]), ]
      number.of.medians <- dim(summary.table)[1] - 2
      row.names(summary.table) <- c("Group w/ min(n)", 
                                    paste("Group w/ median(n)", 1 : number.of.medians, sep = ""),
                                    "Group w/ max(n)")
    }
  }

  alpha.hat <- x$a.new
  alpha.hat.sd <- sqrt(x$a.var)
  result2 <- data.frame(alpha.hat, alpha.hat.sd)

  if (is.na(x$prior.mean)) {
    estimate <- as.vector(x$beta.new)
    if (x$intercept == T) {
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
    cat("\n")
    print(round(x$sec.var, 3))
  } else {
    cat("Main summary:\n")
    cat("\n")
    print(round(x$main, 3))
    cat("\n")
    cat("\n")
    cat("Second-level Variance Component Estimation Summary:\n")
    cat("\n")
    print(round(x$sec.var, 3))
    cat("\n")
    cat("\n")
    cat("Regression Summary:\n")
    cat("\n")
    print(round(x$reg, 3))
  }
}

plot.gbp <- function(x,...) {

  y <- x$sample.mean
  se <- x$se
  if (is.na(x$prior.mean)) {
    pr.m <- x$prior.mean.hat
  } else {
    pr.m <- x$prior.mean
  }
  po.m <- x$post.mean
  po.sd <- x$post.sd
  po.low <- x$post.intv.low
  po.upp <- x$post.intv.upp
  cx <- (mean(log(se + 2)) + min(log(se + 2))) / 2
  index <- 1 : length(se)
  ylim.low <- ifelse(min(po.low, y) >= 0, 0.8 * min(po.low, y), 1.2 * min(po.low, y))
  ylim.upp <- ifelse(max(po.upp, y) >= 0, 1.2 * max(po.upp, y), 0.8 * max(po.upp, y))

  par(mfrow = c(1, 2), xaxs = "r", yaxs = "r", mai = c(1, 0.6, 1, 0.3), las = 1, ps = 13)
  xl <- c("Indices (Groups) by the order of data input")
  
  sqrtV <- se
  sdlens <- sqrtV / max(sqrtV)
  postlens <- po.sd / max(sqrtV)
  xmin <- min(c(y, po.m, pr.m))
  xmax <- max(c(y, po.m, pr.m))
  sunflowerplot(rep(4, length(y)) ~ y, ylim = c(-1, 5), xlim = c(xmin - abs(xmin) * 0.1, 
                xmax + abs(xmax) * 0.1), yaxt = "n", col.lab = "white", main = "Shrinkage Plot")
  if (min(pr.m) == max(pr.m)) {
    if(length(unique(pr.m)) == 1)
      points(pr.m[1], 0, col = "darkviolet", pch = 2, cex = 4)
    legend("bottomright", c(ifelse(x$model == "gr", "se", "n"), "prior mean"), col = c(4, "darkviolet"),
           lty = 1, lwd = c(2, NA), pch = c(NA, 2), seg.len = 0.5, bty = "n")
  } else {
    if(length(unique(pr.m)) == 1)
      points(pr.m[1], 0, col = "darkviolet", pch = 2, cex = 4)
    legend("bottomright", c(ifelse(x$model == "gr", "se", "n")), col = 4,
           lty = 1, lwd = 2, seg.len = 0.5, bty = "n")
  }
  sunflowerplot(rep(0, length(y)) ~ po.m, add = T)
  abline(h = 4)
  abline(h = 0)
  axis(2, c(0, 4), c(expression(hat(theta)), expression(bar(y))), cex.axis = 1.1)
  sapply(1 : length(y), function(i) {
    lines(c(y[i], po.m[i]), c(4, 0))
    lines(c(y[i], y[i] + sdlens[i] * sd(y) * 0.4), c(4, 4 + sdlens[i]), col = 4)
    ##posterior variance lines
    lines(c(po.m[i] - postlens[i] * sd(y) * 0.4, po.m[i]), c(0 - postlens[i], 0), col = 4)
    xcord <- ((4 * po.m[i] / (y[i] - po.m[i]) - 4 * po.m / (y - po.m)) / 
              (4 / (y[i] - po.m[i]) - 4 / (y - po.m)))
    ycord <- 4 / (y - po.m) * xcord - 4 / (y - po.m) * po.m
    coords <- subset(cbind(xcord, ycord), ycord > 0 & ycord < 4)
    points(coords, col = 2)
  })
  
  plot(index, po.m, ylim = c(ylim.low, ylim.upp), xlab = xl, ylab = expression(theta),
       main = "100*Alpha% Intervals for Posterior Mean", cex = log(se + 2) / cx, col = 2, pch = 19)
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
  legend("topright", pch = c(19, 1, NA), col = c(2, 1, 4), lwd = c(NA, NA, 2), 
         c("posterior mean", "sample mean", "prior mean"), seg.len = 0.5, bty = "n")
}
