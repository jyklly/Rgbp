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
  z <- given$z
  n <- given$n
  x.ini <- given$x.ini
  if (identical(x.ini, NA) & !given$intercept) {
    print("In case we do not designate prior mean, at least one beta should be estimated (either intercept=T or at least one covariate is needed)")
    stop()
  } else if (identical(x.ini, NA) & given$intercept) {
    x <- matrix(1, length(y), 1)
    b.ini <- mean(y)
  } else if (!identical(x.ini, NA) & given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("z ~ ", paste(xname, collapse = "+")))
    b.ini <- glm(formula, offset = log(n), family = poisson)$coefficients
    x <- as.matrix(cbind(rep(1, length(y)), x.ini))
  } else if (!identical(x.ini,NA) & !given$intercept) {
    x.ini <- as.matrix(x.ini)
    xname <- paste("x.ini[, ", 1 : ncol(x.ini), "]", sep = "")
    formula <- as.formula(paste("z ~ ", paste(xname, collapse = "+"), "- 1"))
    b.ini <- glm(formula, offset = log(n), family = poisson)$coefficients
    x <- x.ini
  }	
  mu0.ini <- mean(exp(x %*% b.ini))
  r.ini <- mu0.ini / var(y)
  list(x=x, b.ini=b.ini, a.ini= -log(r.ini))
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
#
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
  dif <- rep(1, m + 1)
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


# test1
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
PR(z,n,x)
PR(z,n)
PR(z,n,prior.mean=0.4)

# test2
n<-c(10,10,10,10,10)
z<-c( 3, 4, 1, 2, 0)
x<-c( 1, 1, 0, 0, 0)
PR(z,n,x)  # not working
PR(z,n)  # not working
PR(z,n,prior.mean=0.2)

# test3
n<-c(18,24,28,37,42,48,56,60,68,73,75,79,91,99,104)
z<-c(4,1,3,1,2,1,2,3,6,0,8,3,3,9,7)
x<-c(0.251,-0.021,0.208,0.019,-0.169,0.164,0.296,0.199,0.209,0.093,0.002,0.064,0.105,0.073,0.209)
PR(z,n,x)
PR(z,n)  # not working 
PR(z,n,prior.mean=0.05)

# test4
PR(test$z1,test$n,test[,4:6])
PR(test$z1,test$n)
PR(test$z1,test$n,prior.mean=0.49)
PR(test$z3,test$n,test[,4:6])
PR(test$z3,test$n)  # not working
PR(test$z3,test$n,prior.mean=0.5)

# test5
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
system.time(PR(z0,n0,x0))
PR(z0,n0,cbind(x0,log(n0)))
system.time(PR(z0,n0))  # not working
system.time(PR(z10,n10,x10))
PR(z10,n10,cbind(x10,log(n10)))
PR(z10,n10)

# test6
n<-rep(10,10)
z<-c(2,2,2,2,3,3,3,4,5,5)
x<-c(1,1,1,1,0,0,0,0,0,0)
PR(z,n,x)
PR(z,n)
PR(z,n,prior.mean=0.4)

# test7
ALL<-c(147,134,134,122,118,117,115,114,111,111,109,100,100,97,95,94,88,87,86,84,84,83,81,79,78,77,74,69,66,62)
ALLPG<-c(1.62,1.46,1.46,1.33,1.3,1.3,1.25,1.24,1.22,1.25,1.18,1.09,1.09,1.07,1.03,1.03,0.98,0.97,0.93,0.95,0.91,0.9,0.88,0.9,0.85,0.84,0.8,0.78,0.73,0.67)
ngp<-ALL/ALLPG
exps<-ngp*40
MPF<-c(26,21,33,21,27,21,24,14,29,28,19,24,21,17,29,13,18,19,22,20,22,24,25,6,14,17,12,12,4,17)
dvs<-c(0,0,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1,0,1,1,1,1,0,1,1,0,0,1,0,1)
row.name<-c("FLA","ARI","ATL","WAS","KCR","TBR","SFG","CHW","COL","MIN","CIN","HOU","LAD","SEA","SDP","NYM","BOS","PIT","DET","BAL","TOR","OAK","MIL","NYY","TEX","CHC","STL","CLE","PHI","LAA")
PR(MPF,exps)  # not working

# test8
n<-rep(100,10)
z<-c(0,0,1,1,2,2,3,3,4,4)
x1<-c(0,0,0,0,1,1,1,1,1,1)
x2<-rnorm(10)
PR(z,n)  # not working
PR(z,n,prior.mean=0.02)
PR(z,n,x1)
PR(z,n,x2)

# test9
n<-rep(100,10)
z<-c(1,1,1,1,2,2,3,3,3,3)
x<-c(1,1,1,1,0,0,0,0,0,0)
PR(z,n)  # not working
PR(z,n,prior.mean=0.02)
PR(z,n,x)

# test10
n<-rep(100,10)
z<-c(rep(1,5),rep(3,5))
x<-c(1,1,1,1,0,0,0,0,0,0)
PR(z,n)  # not working
PR(z,n,prior.mean=0.02)
PR(z,n,x)

# test11
n<-rep(100,10)
z<-c(36, 37, 37,38, 42, 42, 42, 45, 45, 50)
x<-c(1,1,1,1,0,0,0,0,0,0)
PR(z,n)
PR(z,n,prior.mean=0.4)
PR(z,n,x)

# test12
n<-rep(10000,10)
z<-c(0, 0, 0, 1, 1, 1, 1, 2, 2, 3)
x<-c(1,1,1,1,0,0,0,0,0,0)
PR(z,n)  # not working
PR(z,n,prior.mean=0.0001)
PR(z,n,x)

# test13
n<-rep(1000,10)
z<-c(476, 482, 488, 494, 500,500, 506, 512, 518, 524)
x<-c(1,1,1,1,0,0,0,0,0,0)
PR(z,n)
PR(z,n,prior.mean=0.5)
PR(z,n,x)

# test14
z<-baseball[,4]
n<-baseball[,3]
PR(z,n)
PR(z,n,prior.mean=0.5)  # not working

# test15
n<-c(33,50,70,63)
z<-c(23,40,62,40)
x1<-c(13.93,2.27)
x2<-c(9.19,1.5)
x3<-c(19.13,2.09)
x4<-c(18.97,1.94)
x<-rbind(x1,x2,x3,x4)
PR(z,n,x)  # not working
PR(z,n)  # not working

# test16
n<-rep(2,27)
z<-c( 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1)
left<-rep(c(1,1,1,0,0,0,0,0,0),3)
right<-rep(c(0,0,0,1,1,1,0,0,0),3)
x<-cbind(left,right)
PR(z,n,x)
PR(z,n)

# test17
new.n<-rep(3,9)
new.z<-c(1,1,1,2,0,1,1,2,3)
new.left<-rep(c(1,0,0),3)
new.right<-rep(c(0,1,0),3)
new.x<-cbind(new.left,new.right)
PR(new.z,new.n,new.x)
PR(new.z,new.n)  # not working

# test18
hp<-read.csv("hospitals.csv",head=T)
z<-hp[,4]
n<-hp[,5]
PR(z,n)  # not working

# test19
rat_dat <- read.table("rats_alex.txt", header=TRUE)
z <- rat_dat$y
n <- rat_dat$N
PR(z,n)  # not working

# test20
ab<-c(150,358,605,581,573,506,591,553,428,345,563,542,588,421,497,428,507,526,556,625,415,233,490,353)
hit<-c(36,113,212,188,216,194,248,226,167,127,208,201,225,161,191,143,197,211,189,211,157,79,175,114)
PR(hit,ab)  # not working
