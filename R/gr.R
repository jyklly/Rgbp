#Note this R code fits the hierarchical normal-normal model using ADM
#Author: Joseph Kelly

gr<-function(y,se,X,mu,CI=0.95,intercept=T,eps=0.0001){
	
  ##define some values that will be used often
  muknown <- !missing(mu)
  type <- 1
  k <- length(y)
  V <- se^2
  if(missing(X)){
    X.ini<-NA
    X <- as.matrix(rep(1,k))
  }else{
    X <- X.ini <- as.matrix(X)
    if(intercept)
      X <- cbind(rep(1,k),X)
  }
	
  ##optimize depending on if mu is specified. 
  if(muknown){
    r<-0
    est <- optim(mean(log(V)),function(x){gr.ll.muknown(exp(x),y,V,mu,type)},lower=0,upper=Inf,control=list(fnscale=-1),method="L-BFGS-B",hessian=T)
    est.var <- solve(-1*as.matrix(est$hessian))
    Ahat <- exp(est$par)
    Betahat<-NA
    hess<-NA
  }else{
    r <- dim(X)[2]
    ## chang<-1
    ## count<-0
    ## mu <- rep(mean(y),length(y))
    ## while(chang > eps){
    ##   est<-optim(log(mean(V)),allk,hessian=TRUE)
    ##   Ahat<-exp(est$par)
    ##   Bhat<-V/(V+Ahat)
    ##   DVAhat<-diag(V+Ahat)
    ##   Betahat<-solve(t(X)%*%solve(DVAhat)%*%X)%*%t(X)%*%solve(DVAhat)%*%y
    ##   munew<-X%*%Betahat
    ##   chang<-max(abs((mu-munew)/mu))
    ##   mu<-munew
    ##   count<-count+1
    ## }
    est <- optim(c(mean(log(V)),rep(0,r)), function(x){gr.ll.muunknown(exp(x[1]),y,V,X,type)},control=list(fnscale=-1),method="L-BFGS-B",hessian=T)
    est.var <- -1/est$hessian
    Ahat <- exp(est$par[1])
  }

  ##calculate estimates
  ninfo <- -1*as.matrix(est$hessian)[1,1]
  Avar<-est.var[1,1]*Ahat^2
  Bhat<-V/(V+Ahat)
  a0<-ninfo/Bhat
  a1<-ninfo/(1-Bhat)
  v<-Bhat*(1-Bhat)/(a1+a0+1)
  seB<-sqrt(v)
  skewB<-2*(a0-a1)*sqrt(1+a0+a1)/(sqrt(a0*a1)*(2+a0+a1))
  LCLB <- qbeta((1-CI)/2,a1,a0)
  UCLB <- qbeta(1 - (1-CI)/2,a1,a0)
  
  ##if mu is estimated use different formula for posterior sd and estimate mu
  if(muknown){
    shat<-sqrt(V*(1-Bhat)+v*(y-mu)^2)
  }else{
    DVAhat<-diag(V+Ahat)
    DVAhati <- diag(1/(V+Ahat))
    Betahat <- solve(t(X)%*%DVAhati%*%X)%*%t(X)%*%DVAhati%*%y
    BetahatSE <- solve(t(X)%*%DVAhati%*%X)
    hess <- -1/BetahatSE^2
    mu<-X%*%Betahat
    ##calculate posterior variance
    E <- eigen(DVAhat)
    W <- E$values
    Q <- E$vectors
    Z <- Q%*%diag(1/sqrt(W))%*%t(Q)
    Phat<-Z%*%X%*%solve(t(X)%*%solve(DVAhat)%*%X)%*%t(X)%*%Z
    p<-diag(Phat)
    shat<-sqrt((1-(1-p)*Bhat)*V+v*(y-mu)^2)
  }
  
  ##calculate some other moments
  thetahat<-(1-Bhat)*y+Bhat*mu
  mu3<-skewB*seB^3
  skew<-((mu-y)^3*mu3-3*V*(mu-y)*v)/shat^3
  sqrtV<-sqrt(V)
  hmB <- length(Bhat)/sum(1/Bhat)
  m4B <- (6*(a0^3 - a0^2*(2*a1-1) + a1^2*(a1+1)-2*a0*a1*(a1+2))/(a0*a1*(a0+a1+2)*(a0+a1+3)) + 3)*v^2
  m4 <- 3*V^2*beta(a1+2,a0)/beta(a1,a0) + (mu-y)^4*m4B +6*(beta(a1+3,a0)/beta(a1,a0)*V*(y-mu)^2 - 2*beta(a1+2,a0)/beta(a1,a0)*(y-mu)^2*(1-a0/(a1+a0))*V + V*(1-a0/(a0+a1))^3*(y-mu)^2)
  kurt <- m4/shat^2
	
  ## calculate CIs using skewed normal
  ## TODO: use apply not loop
  snparam <- matrix(nrow=length(thetahat),ncol=3)
  skewedmat <- c()
  for(i in 1:length(thetahat)){    
    snparam[i,] <- gr.cp.to.dp(c(thetahat[i],shat[i],sign(skew[i])*min(abs(skew[i]),0.94)))
    row <- c(qsn((1-CI)/2,snparam[i,1],snparam[i,2],snparam[i,3]),thetahat[i],qsn(1-(1-CI)/2,snparam[i,1],snparam[i,2],snparam[i,3]))
    skewedmat <- rbind(skewedmat,row)
  }
  skewedmat <- as.matrix(skewedmat)
  
  ## return output
  ## TODO: discuss with Tak and sync output
  output<- list(sample.mean=y,se=se,prior.mean=mu,shrinkage=Bhat,se.shrinkage=seB,post.intv.low=skewedmat[,1],post.mean=thetahat,post.intv.upp=skewedmat[,3],post.sd=shat,model="gr",x=X.ini,beta.new=Betahat,beta.hess=hess,intercept=T,a.new=log(Ahat),a.var=est.var[1,1])
  return(output)
}

##log adjusted posterior for known mu
gr.ll.muknown<-function(A,y,V,mu,type){
  lla<-type*log(A)+sum(dnorm(x=y,mean=mu,sd=sqrt(V+A),log=TRUE))
  return(lla)
}

##log adjusted posterior for known mu
gr.ll.muunknown<-function(A,y,V,X,type){
  DVA<-diag(V+A)
  DVAi <- diag(1/(V+A))
  BetaA <- solve(t(X)%*%DVAi%*%X)%*%t(X)%*%DVAi%*%y
  lla<-type*log(A) - 1/2*sum(log(V+A)) - 1/2*log(det(t(X)%*%DVAi%*%X)) - 1/2*t(y-X%*%BetaA)%*%DVAi%*%(y-X%*%BetaA)
  return(lla)
}

##this function is a slightly modified version of cp.to.dp in the sn package
##Author: Adelchi Azzalini
##Webpage: http://azzalini.stat.unipd.it/SN/
gr.cp.to.dp <- function (param) { 
  b <- sqrt(2/pi)
  m <- length(param) - 2
  gamma1 <- param[m + 2]
  if (abs(gamma1) > 0.995271746431) 
    gamma1 <- 0.9952  ##this line was altered from original
  A <- sign(gamma1) * (abs(2 * gamma1/(4 - pi)))^(1/3)
  delta <- A/(b * sqrt(1 + A^2))
  lambda <- delta/sqrt(1 - delta^2)
  E.Z <- b * delta
  sd.Z <- sqrt(1 - E.Z^2)
  location <- param[1:m]
  location[1] <- param[1] - param[m + 1] * E.Z/sd.Z
  scale <- param[m + 1]/sd.Z
  dp <- c(location, scale, lambda)
  names(dp)[(m + 1):(m + 2)] <- c("scale", "shape")
  if (m == 1) 
    names(dp)[1] <- "location"
  dp
}
