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

  list(a.ini = a.ini, r.ini = exp(-a.ini))
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
              - digamma(n - z + exp(-a) * q) + digamma(exp(-a) * q)) * (p - q)) * exp(-a) * p * q
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
               + (dig.part1 + exp(-a) * trig.part2) * p * q * (p - q))
    const3 <- trig.part3    
    const4 <- ((2 * (trig.part1 + 2 * exp(-a) * fourg.part1) + exp(-a * 2) * fifg.part1) * p^2 * q^2
               + (2 * trig.part2 + exp(-a) * fourg.part2) * p * q * (p - q))
    sum.diag <- (trig.part1 * exp(-a) * p * q + dig.part1 * (p - q)) * exp(-a) * p * q
    
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
    hessian.inverse <- solve(hessian)
    updated <- ini.value - hessian.inverse %*% score
    dif <- ini.value - updated
    ini.value <- updated
  }
  list(x = x, b.ini = ini.value[2 : (m + 1)], a.ini = ini.value[1])
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
  # Log likelihood function of alpha and beta (regression coefficients) for BRIMM 
  # when the descriptive second level mean is unknown.
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

# br hessian matrix with respect to beta, regression coefficients
br.deriv2b<-function(a,b,given,ini){
  z<-given$z
  n<-given$n
  x<-ini$x
  p.<-exp(x%*%b)/(1+exp(x%*%b))
  q.<-1-p.
  temp<-((trigamma(z+exp(-a)*p.)-trigamma(exp(-a)*p.)+trigamma(n-z+exp(-a)*q.)-trigamma(exp(-a)*q.))*exp(-a)*p.*q.+(digamma(z+exp(-a)*p.)-digamma(exp(-a)*p.)-digamma(n-z+exp(-a)*q.)+digamma(exp(-a)*q.))*(q.-p.))*p.*q.*exp(-a)
  t(x)%*%diag(as.numeric(temp))%*%x
}

# alpha estimation (when prior.mean is known)
alpha.est.prior.kn<-function(given,ini){

  llik<-function(a){
    BRLogLikKn(a,given)
  }
  opti <- optim(ini$a.ini,function(a){a+llik(a)},control=list(fnscale=-1),method="L-BFGS-B",hessian=T,lower=-Inf,upper=Inf)
  list(a.new=opti$par,a.hess=opti$hessian,beta.new=NA,beta.hess=NA)
}

alpha.est.prior.un<-function(given,ini){

  # switch log likelihood
  loglik<-function(a,b){BRLogLikUn(a,b,given,ini)}

  # switch hessian with respect to beta, regression coefficients (function of alpha)
  hess_a<-function(a,b){br.deriv2b(a,b,given,ini)}

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
  xVx <- diag(x%*%solve(-beta.hess)%*%t(x))
  mu0 <- exp(x %*% beta.new + xVx/2)
  b0<- (1+mu0)/(mu0*(exp(xVx)-1))+2
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

# main function
BR <- function(z, n, X, prior.mean, intercept=T, Alpha=0.95){

  X<-if(missing(X)){X<-NA}else{X<-X}
  prior.mean<-if(missing(prior.mean)){prior.mean<-NA}else{prior.mean<-prior.mean}

  given<-list(z=z,n=n,sample.mean=z/n,x.ini=X,prior.mean=prior.mean,intercept=intercept,Alpha=Alpha)

  # initial values
  if(is.na(prior.mean)){
    ini<-BRInitialValueUn(given)
    ini<-BRAlphaBetaEstUn(given,ini)    
  }else{
    ini<-BRInitialValueKn(given)
    ini<-BRAlphaEstKn(given,ini)
  }

  # alpha (and beta if prior.mean is unknown) estimation including hessian
  a.res<-if(is.na(prior.mean)){alpha.est.prior.un(given,ini)}else{alpha.est.prior.kn(given,ini)}

  # shrinkage estimation
  B.res<-shrink.est(a.res,given)

  # posterior estimation
  if(is.na(prior.mean)){
    post.res<-br.post.est.prior.un(B.res,a.res,ini,given)
  }else{
    post.res<-br.post.est.prior.kn(B.res,given)
  }
  a.var=1/B.res$inv.info
  output<-list(sample.mean=given$sample.mean,se=given$n,prior.mean=post.res$prior.mean, shrinkage=B.res$B.hat, sd.shrinkage=sqrt(B.res$var.B.hat), post.mean=post.res$post.mean, post.sd=post.res$post.sd, post.intv.low=post.res$post.intv.low, post.intv.upp=post.res$post.intv.upp, model="br", x=X, beta.new=a.res$beta.new, beta.hess=a.res$beta.hess, intercept=intercept, a.new=a.res$a.new, a.var=1/B.res$inv.info)
  output
}
