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
  }else{
    ini<-BRInitialValueKn(given)
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
