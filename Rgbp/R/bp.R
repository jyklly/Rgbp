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
  lambda0.ini <- mean(exp(x %*% b.ini))
  r.ini <- lambda0.ini / var(y)
  list(x=x, b.ini=b.ini, a.ini= -log(r.ini))
}

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

PRLogLikKn <- function(a, given) {
  # Log likelihood function of alpha for PRIMM when the descriptive second level mean is known.
  z <- given$z
  n <- given$n
  prior.mean <- given$prior.mean
  sum(dnbinom(z, exp(-a) * prior.mean, prob = exp(-a) / (exp(-a) + n), log=T))
}

PRLogLikUn <- function(a, b, given, ini) {
  # Log likelihood function of alpha and beta (regression coefficients) for PRIMM when the descriptive second level mean is unknown.
  z <- given$z
  n <- given$n
  x <- ini$x
  sum(dnbinom(z, exp(-a + x %*% as.vector(b)), prob = exp(-a) / (exp(-a) + n), log=T))
}

BRLogLikKn <- function(a, given) {
  # Log likelihood function of alpha for BRIMM when the descriptive second level mean is known.
  n <- given$n
  z <- given$z
  prior.mean <- given$prior.mean
  t0 <- prior.mean * exp(-a)
  t1 <- (1 - prior.mean) * exp(-a)
  t2 <- exp(-a)
  if (any(c(t0, t1, t2, n-z+t1, t0+z)<=0)) {
    print("The components of lgamma should be positive")
    stop()
  } else {
    sum(lgamma(t0 + z) - lgamma(t0) + lgamma(n - z + t1) - lgamma(t1) + lgamma(t2) - lgamma(n + t2))
  }
}

BRLogLikUn <- function(a, b, given, ini) {
  # Log likelihood function of alpha and beta (regression coefficients) for BRIMM when the descriptive second level mean is unknown.
  n <- given$n
  z <- given$z
  x <- ini$x
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

# pr hessian matrix with respect to beta, regression coefficients
pr.deriv2b<-function(a,b,given,ini){
  z<-given$z
  n<-given$n
  x<-ini$x
  temp<-(trigamma(z+exp(-a+x%*%b))*exp(-a+x%*%b)+digamma(z+exp(-a+x%*%b))-trigamma(exp(-a+x%*%b))*exp(-a+x%*%b)-digamma(exp(-a+x%*%b))-log(1+n*exp(a)))*exp(x%*%b)
  exp(-a)*t(x)%*%diag(as.numeric(temp))%*%x
}






# alpha estimation 2
BRAlphaBetaEstUn <- function(given, ini) {
  # Alpha and Beta estimation of BRIMM when the second level mean is unknown
  z <- given$z
  n <- given$n
  x <- ini$x
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

    trig.part1 <- trigamma.z.r.p - trigamma.r.p + trigamma.n.z.r.q - trigamma.r.q
    trig.part2 <- (trigamma.z.r.p - trigamma.r.p) * p - (trigamma.n.z.r.q - trigamma.r.q) * q
    trig.part3 <- ((trigamma.z.r.p - trigamma.r.p) * p^2 + (trigamma.n.z.r.q - trigamma.r.q) * q^2
                   + trigamma.r - trigamma.n.r)
    dig.part1 <- digamma.z.r.p - digamma.r.p - digamma.n.z.r.q + digamma.r.q
    dig.part2 <- ((digamma.z.r.p - digamma.r.p) * p + (digamma.n.z.r.q - digamma.r.q) * q 
                  + digamma.r - digamma.n.r)
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

    out <- c(1 - exp(-a) * (sum(const1) - m / 2 * sum(const2) / sum.diag), 
             exp(-a * 2) * (sum(const3) - m / 2 * (sum(const4) / sum.diag - (sum(const2) / sum.diag)^2)))
    out
  }

  BRDerivAlphaBeta <- function(a, b) {
    # The cross derivative of log likelihood with respect to alpha and beta
    p <- exp(x %*% b) / (1 + exp(x %*% b))
    q <- 1 - p
    vec <- (digamma(z + exp(-a) * p) - digamma(exp(-a) * p)
            - digamma(n - z + exp(-a) *q) + digamma(exp(-a) * q)) * p * q
           + ((trigamma(z + exp(-a) * p) - trigamma(exp(-a) * p)) * p
              - (trigamma(n - z + exp(-a) * q) - trigamma(exp(-a) * q)) * q) * exp(-a) * p * q
    exp(-a) * t(x) %*% as.vector(vec)
  }

  ini.value <- c(a.ini, b.ini)
  dif <- 1
  eps <- 0.0001
  n.iter <- 0
  while (max(abs(dif)) > eps) { 
    out1 <- BRDerivAlpha(ini.value[1], ini.value[2 : (m+1)])
    out2 <- BRDerivBeta(ini.value[1], ini.value[2 : (m+1)])
    out3 <- BRDerivAlphaBeta(ini.value[1], ini.value[2 : (m+1)])
    hessian <- cbind(c(out1[2], out3), rbind(as.vector(out3), out2[, 2 : (m + 1)]))
    score <- c(out1[1], out2[, 1])
    updated <- ini.value - solve(hessian) %*% score
    dif <- ini.value - updated
    ini.value <- updated
    n.iter <- n.iter + 1
  }
  list(a.hat = ini.value[1], b.hat = ini.value[2 : (m+1)], hessian = hessian)
}

BRAlphaBetaEstUn(given, ini)












# alpha estimation (when prior.mean is known)
alpha.est.prior.kn<-function(given,ini){
  # change log likelihood depending on the model: br or pr
  llik<-function(a){
    switch(given$model,br=BRLoglikKn(a,given),pr=PRLogLikKn(a,given))
  }
  opti <- optim(ini$a.ini,function(a){a+llik(a)},control=list(fnscale=-1),method="L-BFGS-B",hessian=T,lower=-Inf,upper=Inf)
  list(a.new=opti$par,a.hess=opti$hessian,beta.new=NA,beta.hess=NA)
}


alpha.est.prior.un<-function(given,ini){

  # switch log likelihood
  loglik<-function(a,b){switch(given$model,br=BRLogLikUn(a,b,given,ini),pr=PRLogLikUn(a,b,given,ini))}

  # switch hessian with respect to beta, regression coefficients (function of alpha)
  hess_a<-function(a,b){switch(given$model,br=br.deriv2b(a,b,given,ini),pr=pr.deriv2b(a,b,given,ini))}

  # objective function of alpha
  objfun_a<-function(a){
  	btemp_a<-optim(ini$b.ini,function(b){loglik(a,b)},control=list(fnscale=-1),method="L-BFGS-B",hessian=T,lower=-Inf,upper=Inf)$par
	a+loglik(a,btemp_a)-0.5*log(det(-hess_a(a,btemp_a)))
  }

  # estimation of alpha
  a.new<-optim(ini$a.ini,objfun_a,control=list(fnscale=-1),method="L-BFGS-B",hessian=T,lower=-Inf,upper=Inf)$par
  a.hess<-optim(ini$a.ini,objfun_a,control=list(fnscale=-1),method="L-BFGS-B",hessian=T,lower=-Inf,upper=Inf)$hessian

  # objective function of beta
  objfun_b<-function(b){
	loglik(a.new,b)
  }

  # estimation of beta
  if(ncol(as.matrix(ini$x))>1){
    beta.new<-optim(ini$b.ini,objfun_b,control=list(fnscale=-1),hessian=T)$par
    beta.hess<-optim(ini$b.ini,objfun_b,control=list(fnscale=-1),hessian=T)$hessian
  }else{
    beta.new<-optim(ini$b.ini,objfun_b,control=list(fnscale=-1),method="L-BFGS-B",hessian=T,lower=-Inf,upper=Inf)$par
    beta.hess<-optim(ini$b.ini,objfun_b,control=list(fnscale=-1),method="L-BFGS-B",hessian=T,lower=-Inf,upper=Inf)$hessian
  }
  list(a.new=a.new,a.hess=a.hess,beta.new=beta.new,beta.hess=beta.hess)
}

# shrinkage estimation
shrink.est<-function(a.res,given){	
  a.new<-a.res$a.new
  a.hess<-a.res$a.hess
  B.hat<-exp(-a.new)/(given$n+exp(-a.new))
  inv.info<- -a.hess
  var.B.hat<-(B.hat*(1-B.hat))^2/((B.hat*(1-B.hat))+inv.info)
  a1.beta<-inv.info/(1-B.hat)
  a0.beta<-inv.info/B.hat
  central3.B<-2*(1-2*B.hat)*B.hat*(1-B.hat)/(a1.beta+a0.beta+1)/(a1.beta+a0.beta+2) # for brimm
  B.hat.low<-qbeta((1-given$Alpha)/2,a1.beta,a0.beta)
  B.hat.upp<-qbeta((1+given$Alpha)/2,a1.beta,a0.beta)
  list(B.hat=B.hat,inv.info=inv.info,var.B.hat=var.B.hat,central3.B=central3.B)
}

# brimm posterior estimation (when prior.mean is known)
br.post.est.prior.kn<-function(B.res,given){
  B.hat<-B.res$B.hat
  var.B.hat<-B.res$var.B.hat
  p0<-given$prior.mean
  central3.B<-B.res$central3.B
  y<-given$sample.mean; n<-given$n
  p.hat<- y-B.hat*(y-p0)
  c1<- y*(1-y)
  c2<- -3*y^2+2*(1+p0)*y-p0
  c3<-(y-p0)*(3*y-1-p0)
  c4<-(y-p0)^2
  d1<-c1+c3*B.hat^2
  d2<-c2+2*c3*B.hat-3*c4*B.hat^2
  d3<-c3-3*c4*B.hat
  d4<-c4
  var.p.hat<-1/n*(d1-d2*B.hat-d3*var.B.hat+d4*central3.B)
  sd.p.hat<-sqrt(var.p.hat)
  a1.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*p.hat
  a0.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*(1-p.hat)
  p.hat.low<- qbeta((1-given$Alpha)/2,a1.beta.p,a0.beta.p)
  p.hat.upp<- qbeta((1+given$Alpha)/2,a1.beta.p,a0.beta.p)
  list(post.mean=p.hat,post.sd=sd.p.hat,post.intv.low=p.hat.low,post.intv.upp=p.hat.upp,prior.mean=p0)
}

# brimm posterior estimation (when prior.mean is unknown)
br.post.est.prior.un<-function(B.res,a.res,ini,given){
  x<-as.matrix(ini$x)
  n<-given$n
  beta.new<-a.res$beta.new
  beta.hess<-a.res$beta.hess
  B.hat<-B.res$B.hat
  y<-given$sample.mean
  xHx <- diag(x%*%solve(-beta.hess)%*%t(x))
  mu0 <- exp(x %*% beta.new + xHx/2)
  b0<- (1+mu0)/(mu0*(exp(xHx)-1))+2
  b1<- mu0*(b0-1)
  k1p <- b1/(b1+b0)
  # cumulants	
  k2p<-as.vector(k1p*(1-k1p)/(b1+b0+1))
  k1b<-as.vector(B.hat)
  k2b<-as.vector(B.res$var.B.hat)
  k3b<-as.vector(B.res$central3.B)
  p0.hat<-k1p
  p.hat<-(1-B.hat)*y+B.hat*p0.hat	# posterior p
  # coefficients		
  c1<-as.vector(y-y^2-(y-3*y^2)*B.hat^2-2*y^2*B.hat^3-(B.hat^2-2*B.hat^3)*k1p^2)
  c2<-as.vector(-2*y+3*y^2+2*(y-3*y^2)*B.hat+3*y^2*B.hat^2-(-2*B.hat+3*B.hat^2)*k1p^2)
  c3<-as.vector(y-3*y^2+3*y^2*B.hat-(3*B.hat-1)*k1p^2)
  c4<-as.vector(y^2-k1p^2)
  c5<-as.vector(4*y*B.hat^3-(4*y-1)*B.hat^2+2*(B.hat^2-2*B.hat^3)*k1p)
  c6<-as.vector(2*(4*y-1)*B.hat-6*y*B.hat^2-2*y+1+2*(-2*B.hat+3*B.hat^2)*k1p)
  c7<-as.vector(4*y-1-6*y*B.hat+2*(3*B.hat-1)*k1p)
  c8<-as.vector(-2*y+2*k1p)
  c9<-as.vector(B.hat^2-2*B.hat^3)
  c10<-as.vector(-2*B.hat+3*B.hat^2)
  c11<-as.vector(3*B.hat-1)
  # posterior variance
  var.p.hat<-1/n*(c1+c2*k1b+c3*k2b+c4*k3b+c5*k1p+c6*k1p*k1b+c7*k1p*k2b+c8*k1p*k3b+c9*k2p+c10*k1b*k2p+c11*k2b*k2p+k3b*k2p)+k2b*k2p+k2b*k1p^2-2*y*k2b*k1p+y^2*k2b+k1b^2*k2p;
  sd.p.hat<-sqrt(var.p.hat)
  a1.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*p.hat
  a0.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*(1-p.hat)
  p.hat.low<- qbeta((1-given$Alpha)/2,a1.beta.p,a0.beta.p)
  p.hat.upp<- qbeta((1+given$Alpha)/2,a1.beta.p,a0.beta.p)
  list(post.mean=p.hat,post.sd=sd.p.hat,post.intv.low=p.hat.low,post.intv.upp=p.hat.upp,prior.mean=p0.hat)
}

# primm1 posterior estimation (when prior.mean is known)
pr.post.est.prior.kn<-function(B.res,given){
  B.hat<-B.res$B.hat
  var.B.hat<-B.res$var.B.hat
  lambda0<-given$prior.mean
  n<-given$n
  y<-given$sample.mean;lambda.hat<- y-B.hat*(y-lambda0)
  var.lambda.hat<-1/n*(lambda.hat-B.hat*y+(var.B.hat+B.hat^2)*y-(var.B.hat+B.hat^2)*lambda0)+(y-lambda0)^2*var.B.hat
  sd.lambda.hat<-sqrt(var.lambda.hat)
  u.gamma<-lambda.hat^2/var.lambda.hat
  v.gamma<-lambda.hat/var.lambda.hat
  lambda.hat.low<-qgamma((1-given$Alpha)/2,shape=u.gamma,rate=v.gamma)
  lambda.hat.upp<-qgamma((1+given$Alpha)/2,shape=u.gamma,rate=v.gamma)
  list(post.mean=lambda.hat, post.sd=sd.lambda.hat, post.intv.low=lambda.hat.low, post.intv.upp=lambda.hat.upp, prior.mean=lambda0)
}


# primm1 posterior estimation (when prior.mean is unknown)
pr.post.est.prior.un<-function(B.res,a.res,ini,given){
  x<-as.matrix(ini$x)
  n<-given$n
  y<-given$sample.mean
  beta.new<-a.res$beta.new
  B.hat<-B.res$B.hat
  var.B.hat<-B.res$var.B.hat
  V<-solve(-a.res$beta.hess)
  u<-sapply(1:length(n),function(t){x[t,]%*%V%*%x[t,]})
  lambda0.hat<-exp(x%*%beta.new+1/2*u)
  lambda.hat<- y-B.hat*(y-lambda0.hat)
  vr<-exp(2*x%*%beta.new+u)*(exp(u)-1)
  var.lambda.hat<-1/n*(lambda.hat-B.hat*y+(var.B.hat+B.hat^2)*y-(var.B.hat+B.hat^2)*lambda0.hat)+(var.B.hat+B.hat^2)*(y^2-2*y*lambda0.hat+vr+lambda0.hat^2)-(B.hat*(y-lambda0.hat))^2
  sd.lambda.hat<-sqrt(var.lambda.hat)
  u.gamma<-lambda.hat^2/var.lambda.hat
  v.gamma<-lambda.hat/var.lambda.hat
  lambda.hat.low<-qgamma((1-given$Alpha)/2,shape=u.gamma,rate=v.gamma)
  lambda.hat.upp<-qgamma((1+given$Alpha)/2,shape=u.gamma,rate=v.gamma)
  list(post.mean=lambda.hat, post.sd=sd.lambda.hat, post.intv.low=lambda.hat.low, post.intv.upp=lambda.hat.upp, prior.mean=lambda0.hat)
}

# main function
bp <- function(z, n, X, prior.mean, model="br", intercept=T, Alpha=0.95){

  X<-if(missing(X)){X<-NA}else{X<-X}
  prior.mean<-if(missing(prior.mean)){prior.mean<-NA}else{prior.mean<-prior.mean}

  given<-list(z=z,n=n,sample.mean=z/n,x.ini=X,prior.mean=prior.mean,model=model,intercept=intercept,Alpha=Alpha)

  # initial values
  if(is.na(prior.mean)){
    ini<-switch(model,br=BRInitialValueUn(given),pr=PRInitialValueUn(given))
  }else{
    ini<-switch(model,br=BRInitialValueKn(given),pr=PRInitialValueKn(given))
  }

  # alpha (and beta if prior.mean is unknown) estimation including hessian
  a.res<-if(is.na(prior.mean)){alpha.est.prior.un(given,ini)}else{alpha.est.prior.kn(given,ini)}

  # shrinkage estimation
  B.res<-shrink.est(a.res,given)

  # posterior estimation
  if(is.na(prior.mean)){
    post.res<-switch(model,br=br.post.est.prior.un(B.res,a.res,ini,given),pr=pr.post.est.prior.un(B.res,a.res,ini,given))
  }else{
    post.res<-switch(model,br=br.post.est.prior.kn(B.res,given),pr=pr.post.est.prior.kn(B.res,given))
  }

  output<-list(sample.mean=given$sample.mean,se=given$n,prior.mean=post.res$prior.mean, shrinkage=B.res$B.hat, sd.shrinkage=sqrt(B.res$var.B.hat), post.mean=post.res$post.mean, post.sd=post.res$post.sd, post.intv.low=post.res$post.intv.low, post.intv.upp=post.res$post.intv.upp, model=model, x=X, beta.new=a.res$beta.new, beta.hess=a.res$beta.hess, intercept=intercept, a.new=a.res$a.new, a.var=1/B.res$inv.info)
  output
}

n<-rep(10,10)
z<-rbinom(10,10,0.5)
X<-rnorm(10)
bp(z,n,x[,2])
BRDerivBeta(ini$a.ini, ini$b.ini)