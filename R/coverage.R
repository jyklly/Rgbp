
coverage<-function(x, nsim=100, detail=F,...){

  # assign empty space for the result to be input
  coverage<-matrix(NA, nrow=length(x$se), ncol=nsim)
  coverage.l<-matrix(NA, nrow=length(x$se), ncol=nsim)
  coverage.u<-matrix(NA, nrow=length(x$se), ncol=nsim)

  # 10 means 1-0 criterion that is 1 if interval icludes true parameter, 0 if not
  coverage10<-matrix(NA, nrow=length(x$se), ncol=nsim)
  coverage10.l<-matrix(NA, nrow=length(x$se), ncol=nsim)
  coverage10.u<-matrix(NA, nrow=length(x$se), ncol=nsim)

  # covariate matrix
  if(identical(x$x,NA)){
    pre.x<-as.matrix(rep(1,length(x$se)))
  }else{
    pre.x<-as.matrix(cbind(rep(1,length(x$se)),x$x))
  }

  # if model=BRIMM	
  if(x$model=="br"){

    # 1. initial values
    betas<-as.vector(x$beta.new)
    pre.p<-exp(pre.x%*%betas)/(1+exp(pre.x%*%betas))
    n<-x$se
    X<-x$x
    r<-exp(-x$a.new)
    r.l<-as.numeric(exp(-(x$a.new+qnorm(1-0.0833/2)*sqrt(x$a.var))))
    r.u<-as.numeric(exp(-(x$a.new-qnorm(1-0.0833/2)*sqrt(x$a.var))))

    # 2. generate p matrix
    sim.p<-matrix(rbeta(length(pre.p)*nsim, r*pre.p, r*(1-pre.p)), nrow=length(n))
    sim.p.l<-matrix(rbeta(length(pre.p)*nsim, r.l*pre.p, r.l*(1-pre.p)), nrow=length(n))
    sim.p.u<-matrix(rbeta(length(pre.p)*nsim, r.u*pre.p, r.u*(1-pre.p)), nrow=length(n))

    # 3. generate z (data) matrix
    sim.z<-matrix(rbinom(nrow(sim.p)*nsim, n, sim.p), nrow=length(n))
    sim.z.l<-matrix(rbinom(nrow(sim.p.l)*nsim, n, sim.p.l), nrow=length(n))
    sim.z.u<-matrix(rbinom(nrow(sim.p.u)*nsim, n, sim.p.u), nrow=length(n))

    # 4. simulation
    for (i in 1:nsim){
      tryCatch({
        out<-if(identical(X,NA)){gbp(sim.z[,i],n,model="br")}else{gbp(sim.z[,i],n,X,model="br")}
        out.l<-if(identical(X,NA)){gbp(sim.z.l[,i],n,model="br")}else{gbp(sim.z.l[,i],n,X,model="br")}
        out.u<-if(identical(X,NA)){gbp(sim.z.u[,i],n,model="br")}else{gbp(sim.z.u[,i],n,X,model="br")}
        a1<-r*pre.p+sim.z[,i]
        a1.l<-r.l*pre.p+sim.z.l[,i]
        a1.u<-r.u*pre.p+sim.z.u[,i]
        a0<-r*(1-pre.p)+n-sim.z[,i]
        a0.l<-r.l*(1-pre.p)+n-sim.z.l[,i]
        a0.u<-r.u*(1-pre.p)+n-sim.z.u[,i]
        low<-out$post.intv.low
        low.l<-out.l$post.intv.low
        low.u<-out.u$post.intv.low
        upp<-out$post.intv.upp
        upp.l<-out.l$post.intv.upp
        upp.u<-out.u$post.intv.upp
        coverage[,i]<-pbeta(upp,a1,a0)-pbeta(low,a1,a0)
        coverage.l[,i]<-pbeta(upp.l,a1.l,a0.l)-pbeta(low.l,a1.l,a0.l)
        coverage.u[,i]<-pbeta(upp.u,a1.u,a0.u)-pbeta(low.u,a1.u,a0.u)
        coverage10[,i]<-ifelse(low<=sim.p[,i] & sim.p[,i]<=upp,1,0)
        coverage10.l[,i]<-ifelse(low.l<=sim.p.l[,i] & sim.p.l[,i]<=upp.l,1,0)
        coverage10.u[,i]<-ifelse(low.u<=sim.p.u[,i] & sim.p.u[,i]<=upp.u,1,0)
      },error=function(x){print(c(i,"error"))},warning=function(x){print(c(i,"warning"))})
    }
		
  # if model=PRIMM
  }else if(x$model=="pr"){
	
    # 1. initial values
    betas<-as.vector(x$beta.new)
    pre.lambda<-exp(pre.x%*%betas)
    n<-x$se
    X<-x$x
    r<-exp(-x$a.new)
    r.l<-as.numeric(exp(-(x$a.new+qnorm(1-0.0833/2)*sqrt(x$a.var))))
    r.u<-as.numeric(exp(-(x$a.new-qnorm(1-0.0833/2)*sqrt(x$a.var))))

    # 2. generate lambda matrix
    sim.lambda<-matrix(rgamma(length(pre.lambda)*nsim, pre.lambda*r, r), nrow=length(n))
    sim.lambda.l<-matrix(rgamma(length(pre.lambda)*nsim, pre.lambda*r.l, r.l), nrow=length(n))
    sim.lambda.u<-matrix(rgamma(length(pre.lambda)*nsim, pre.lambda*r.u, r.u), nrow=length(n))

    # 3. generate z (data) matrix
    sim.z<-matrix(rpois(nrow(sim.lambda)*nsim, n*sim.lambda), nrow=length(n))
    sim.z.l<-matrix(rpois(nrow(sim.lambda.l)*nsim, n*sim.lambda.l), nrow=length(n))
    sim.z.u<-matrix(rpois(nrow(sim.lambda.u)*nsim, n*sim.lambda.u), nrow=length(n))

    # 4. simulation
    for (i in 1:nsim){
      tryCatch({
        out<-if(identical(X,NA)){gbp(sim.z[,i],n,model="pr")}else{gbp(sim.z[,i],n,X,model="pr")}
        out.l<-if(identical(X,NA)){gbp(sim.z.l[,i],n,model="pr")}else{gbp(sim.z.l[,i],n,X,model="pr")}
        out.u<-if(identical(X,NA)){gbp(sim.z.u[,i],n,model="pr")}else{gbp(sim.z.u[,i],n,X,model="pr")}
        sh<-(r*pre.lambda+sim.z[,i])
        sh.l<-(r.l*pre.lambda+sim.z.l[,i])
        sh.u<-(r.u*pre.lambda+sim.z.u[,i])
        rt<-(n+r)
        rt.l<-(n+r.l)
        rt.u<-(n+r.u)
        low<-out$post.intv.low
        low.l<-out.l$post.intv.low
        low.u<-out.u$post.intv.low
        upp<-out$post.intv.upp
        upp.l<-out.l$post.intv.upp
        upp.u<-out.u$post.intv.upp
        coverage[,i]<-pgamma(upp,sh,rt)-pgamma(low,sh,rt)
        coverage.l[,i]<-pgamma(upp.l,sh.l,rt.l)-pgamma(low.l,sh.l,rt.l)
        coverage.u[,i]<-pgamma(upp.u,sh.u,rt.u)-pgamma(low.u,sh.u,rt.u)
        coverage10[,i]<-ifelse(low<=sim.lambda[,i] & sim.lambda[,i]<=upp,1,0)
        coverage10.l[,i]<-ifelse(low.l<=sim.lambda.l[,i] & sim.lambda.l[,i]<=upp.l,1,0)
        coverage10.u[,i]<-ifelse(low.u<=sim.lambda.u[,i] & sim.lambda.u[,i]<=upp.u,1,0)
      },error=function(x){print(c(i,"error"))},warning=function(x){print(c(i,"warning"))})
    }

  # if model=GRIMM
  }else if(x$model=="gr"){

    # 1. initial values
    betas<-x$beta.new
    pre.mu<-as.matrix(pre.x)%*%as.vector(betas)
    se<-x$se
    X<-x$x
    A<-exp(x$a.new)
    A.l<-exp((x$a.new-qnorm(1-0.0833/2)*sqrt(x$a.var)))
    A.u<-exp((x$a.new+qnorm(1-0.0833/2)*sqrt(x$a.var)))

    # 2. generate mu matrix
    sim.mu<-matrix(rnorm(length(pre.mu)*nsim, pre.mu, sqrt(A)), nrow=length(se))
    sim.mu.l<-matrix(rnorm(length(pre.mu)*nsim, pre.mu, sqrt(A.l)), nrow=length(se))
    sim.mu.u<-matrix(rnorm(length(pre.mu)*nsim, pre.mu, sqrt(A.u)), nrow=length(se))

    # 3. generate y (data) matrix
    sim.y<-matrix(rnorm(nrow(sim.mu)*nsim, sim.mu, se), nrow=length(se))
    sim.y.l<-matrix(rnorm(nrow(sim.mu.l)*nsim, sim.mu.l, se), nrow=length(se))
    sim.y.u<-matrix(rnorm(nrow(sim.mu.u)*nsim, sim.mu.u, se), nrow=length(se))

    # 4. simulation
    for (i in 1:nsim){
      tryCatch({
        out<-if(identical(x$x,NA)){gbp(sim.y[,i], se)}else{gbp(sim.y[,i], se, X)}
        out.l<-if(identical(x$x,NA)){gbp(sim.y.l[,i], se)}else{gbp(sim.y.l[,i], se, X)}
        out.u<-if(identical(x$x,NA)){gbp(sim.y.u[,i], se)}else{gbp(sim.y.u[,i], se, X)}
        low<-out$post.intv.low
        low.l<-out.l$post.intv.low
        low.u<-out.u$post.intv.low
        upp<-out$post.intv.upp
        upp.l<-out.l$post.intv.upp
        upp.u<-out.u$post.intv.upp
        postmean<-pre.mu*(se^2/(se^2+A))+sim.y[,i]*(A/(se^2+A))
        postmean.l<-pre.mu*(se^2/(se^2+A.l))+sim.y.l[,i]*(A.l/(se^2+A.l))
        postmean.u<-pre.mu*(se^2/(se^2+A.u))+sim.y.u[,i]*(A.u/(se^2+A.u))
        postsd<-sqrt(se^2*(A/(se^2+A)))
        postsd.l<-sqrt(se^2*(A.l/(se^2+A.l)))
        postsd.u<-sqrt(se^2*(A.u/(se^2+A.u)))
        coverage[,i]<-pnorm(upp,postmean,postsd)-pnorm(low,postmean,postsd)
        coverage.l[,i]<-pnorm(upp.l,postmean.l,postsd.l)-pnorm(low.l,postmean.l,postsd.l)
        coverage.u[,i]<-pnorm(upp.u,postmean.u,postsd.u)-pnorm(low.u,postmean.u,postsd.u)
        coverage10[,i]<-ifelse(low<=sim.mu[,i] & sim.mu[,i]<=upp,1,0)
        coverage10.l[,i]<-ifelse(low.l<=sim.mu.l[,i] & sim.mu.l[,i]<=upp.l,1,0)
        coverage10.u[,i]<-ifelse(low.u<=sim.mu.u[,i] & sim.mu.u[,i]<=upp.u,1,0)
      },error=function(x){print(c(i,"error"))},warning=function(x){print(c(i,"warning"))})
    }
  }

  # average coverage probability
  result<-round(apply(cbind(coverage,coverage.l,coverage.u),1,mean,na.rm=T),3)
  avr.cov<-round(apply(as.matrix(result),2,mean),3)
  min.cov<-round(apply(as.matrix(result),2,min),3)

  result2<-round(apply(cbind(coverage10,coverage10.l,coverage10.u),1,mean,na.rm=T),3)
  avr.cov2<-round(apply(as.matrix(result2),2,mean),3)
  min.cov2<-round(apply(as.matrix(result2),2,min),3)

  # plotting coverage graph
  par(xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3));
  plot(1:length(x$se),result,ylim=c(0.7,1),type="l",col=2,ylab="",xlab="Group_i, i=1,...,k",main="Coverage for Each Group",lwd=3,lty=1)
  abline(h=0.95)
  points(1:length(x$se),result2,type="l",lty=2,col=4,lwd=2)
  if(x$model=="gr"){
    legend("bottomleft",c(paste("A =",paste( round(A.l,2),",", round(A,2),",", round(A.u,2) )) ,paste("beta",0:(length(betas)-1),"=",round(betas,3)),paste("average.coverage =",avr.cov),paste("minimum.coverage =",min.cov),"Rao-Blackwellized","Unbiased"),lty=c(1,rep(1,length(betas)),1,1,1,2),col=c(1,rep(1,length(betas)),1,1,2,4),lwd=c(0,rep(0,length(betas)),0,0,2,2))
  }else{
    legend("bottomleft",c(paste("r =",paste(round(r.l,2),",",round(r,2),",",round(r.u,2))),paste("beta",0:(length(betas)-1),"=",round(betas,3)),paste("average.coverage =",avr.cov),paste("minimum.coverage =",min.cov),"Rao-Blackwellized","Unbiased"),lty=c(1,rep(1,length(betas)),1,1,1,2),col=c(1,rep(1,length(betas)),1,1,2,4),lwd=c(0,rep(0,length(betas)),0,0,2,2))
  }

  # print output
  if(identical(detail,T)){
    output<-list(coverage=result,coverage.1or0=result2,average.coverage=avr.cov,average.coverage.1or0=avr.cov2,minimum.coverage=min.cov,minimum.coverage.1or0=min.cov2,raw.result=cbind(coverage,coverage.l,coverage.u),raw.result.1or0=cbind(coverage10,coverage10.l,coverage10.u))
  }else{
    output<-list(coverage=result,coverage.1or0=result2,average.coverage=avr.cov,average.coverage.1or0=avr.cov2,minimum.coverage=min.cov,minimum.coverage.1or0=min.cov2)
  }
  return(output)
}
