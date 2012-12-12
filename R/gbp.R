gbp<-function(x, ...) UseMethod("gbp")

gbp.default<-function(x,y,X,mu,model="gr",Alpha=0.95,intercept=T, ...){

  res<-switch(model, 
              gr=gr(x,y,X,mu,Alpha), 
              br=bp(x,y,X,prior.mean=mu,model="br",Alpha=Alpha,intercept=intercept), 
              pr=bp(x,y,X,prior.mean=mu,model="pr",Alpha=Alpha,intercept=intercept) )
  
  class(res)<-"gbp"	
  res
}


print.gbp<-function(x,...){
  
  if(!identical(x$x,NA)){
	cova<-as.matrix(x$x)
	colnames(cova)<-paste(1:dim(cova)[2])
    
    if(x$model=="gr"){
      temp<-data.frame(sample.mean=x$sample.mean, se=x$se, x=cova, prior.mean=x$prior.mean, shrinkage=x$shrinkage, sd.shrinkage=x$sd.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }else{
      temp<-data.frame(sample.mean=x$sample.mean, n=x$se, x=cova, prior.mean=x$prior.mean, shrinkage=x$shrinkage, sd.shrinkage=x$sd.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }
  }else{
    if(x$model=="gr"){
      temp<-data.frame(sample.mean=x$sample.mean, se=x$se, prior.mean=x$prior.mean, shrinkage=x$shrinkage, sd.shrinkage=x$sd.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }else{
      temp<-data.frame(sample.mean=x$sample.mean, n=x$se, prior.mean=x$prior.mean, shrinkage=x$shrinkage, sd.shrinkage=x$sd.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }
  }

  cat("Summary for whole observations:\n")
  cat("\n")
  print(round(temp,3))
  
}


summary.gbp<-function(object,...){

  if(!identical(object$x,NA)){

	cova<-as.matrix(object$x)
	colnames(cova)<-paste(1:dim(cova)[2])

    if(object$model=="gr"){
      temp<-data.frame(sample.mean=object$sample.mean, se=object$se, x=cova, prior.mean=object$prior.mean, shrinkage=object$shrinkage, sd.shrinkage=object$sd.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }else{
      temp<-data.frame(sample.mean=object$sample.mean, n=object$se, x=cova, prior.mean=object$prior.mean, shrinkage=object$shrinkage, sd.shrinkage=object$sd.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }
  }else{
    if(object$model=="gr"){
      temp<-data.frame(sample.mean=object$sample.mean, se=object$se, prior.mean=object$prior.mean, shrinkage=object$shrinkage, sd.shrinkage=object$sd.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }else{
      temp<-data.frame(sample.mean=object$sample.mean, n=object$se, prior.mean=object$prior.mean, shrinkage=object$shrinkage, sd.shrinkage=object$sd.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }
  }
  
  temp.summary<-apply(temp,2,summary)
  
  a.new<-object$a.new
  a.var<-object$a.var
  if(object$model=="gr"){
    A.hat<-exp(a.new)
    sd.A<-sqrt(a.var*A.hat^2)
    low.A.hat<-exp(a.new-qnorm(1-0.0833/2)*sqrt(a.var))
    upp.A.hat<-exp(a.new+qnorm(1-0.0833/2)*sqrt(a.var))
    result2<-data.frame(A.hat,sd.A,low.A.hat,upp.A.hat)
  }else{
    r.hat<-exp(-a.new)
    sd.r<-sqrt(a.var*exp(-2*a.new))
    low.r.hat<-exp(-(a.new+qnorm(1-0.0833/2)*sqrt(a.var)))
    upp.r.hat<-exp(-(a.new-qnorm(1-0.0833/2)*sqrt(a.var)))
    result2<-data.frame(r.hat,sd.r,low.r.hat,upp.r.hat)
  }
  
  if(!identical(object$beta.new,NA)){
    estimate<-as.vector(object$beta.new)
    if(object$intercept==T){
      names(estimate)<-paste("beta", 0:(length(estimate)-1), sep = "")
    }else{
      names(estimate)<-paste("beta", 1:(length(estimate)), sep = "")
    }
    se<-as.vector(sqrt(diag(solve(-object$beta.hess))))
    z.val<-estimate/se
    p.val<-2*pnorm(-abs(z.val))
    beta.result<-data.frame(estimate,se,z.val,p.val)
  }else{
    beta.result<-NA
  }
  
  res<-list(main=temp.summary,sec.var=result2,reg=beta.result)
  class(res)<-"summary.gbp"
  res
}


print.summary.gbp<-function(x,...){
  if(identical(x$reg,NA)){
    cat("Main summary:\n")
    cat("\n")
    print(round(x$main,3))
    cat("\n")
    cat("\n")
    cat("Second-level Variance Component Estimation Summary:\n")
    cat("\n")
    print(round(x$sec.var,3))
  }else{
    cat("Main summary:\n")
    cat("\n")
    print(round(x$main,3))
    cat("\n")
    cat("\n")
    cat("Second-level Variance Component Estimation Summary:\n")
    cat("\n")
    print(round(x$sec.var,3))
    cat("\n")
    cat("\n")
    cat("Regression Summary:\n")
    cat("\n")
    print(round(x$reg,3))
  }
}


plot.gbp<-function(x,...){
  
  y<-x$sample.mean
  se<-x$se
  pr.m<-x$prior.mean
  po.m<-x$post.mean
  po.sd<-x$post.sd
  po.low<-x$post.intv.low
  po.upp<-x$post.intv.upp
  cx<-(mean(log(se+2))+min(log(se+2)))/2
  index<-1:length(se)
  ylim.low<-ifelse(min(po.low,y)>=0, 0.8*min(po.low,y), 1.2*min(po.low,y))
  ylim.upp<-ifelse(max(po.upp,y)>=0, 1.2*max(po.upp,y), 0.8*max(po.upp,y))

  par(mfrow=c(1,2),xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3),las=1,ps=13)
  xl<-c("Indices (Groups) by the order of data input")
  
  sqrtV <- se
  sdlens <- sqrtV/max(sqrtV)
  postlens <- po.sd/max(sqrtV)
  xmin <- min(c(y,po.m,pr.m))
  xmax <- max(c(y,po.m,pr.m))
  sunflowerplot(rep(4,length(y))~y,ylim=c(-1,5),xlim=c(xmin-abs(xmin)*0.1,xmax+abs(xmax)*0.1), yaxt="n", col.lab="white", main="Shrinkage Plot")
  if(length(unique(pr.m))==1)
    points(pr.m[1],0,col="darkviolet",pch=2,cex=4)
  legend("bottomright",ifelse(x$model=="gr","se","n"),col="blue",lty=1,seg.len=0.5,lwd=2)
  sunflowerplot(rep(0,length(y))~po.m,add=T)
  abline(h=4)
  abline(h=0)
  axis(2,c(0,4),c(expression(hat(theta)),expression(bar(y))), cex.axis=1.1)
  sapply(1:length(y), function(i){
    lines(c(y[i],po.m[i]),c(4,0))
    lines(c(y[i],y[i]+sdlens[i]*sd(y)*0.4),c(4,4+sdlens[i]),col="blue")
    ##posterior variance lines
    lines(c(po.m[i]-postlens[i]*sd(y)*0.4,po.m[i]),c(0-postlens[i],0),col="blue")
    xcord <- (4*po.m[i]/(y[i]-po.m[i])-4*po.m/(y-po.m))/(4/(y[i]-po.m[i])-4/(y-po.m))
    ycord <- 4/(y-po.m)*xcord - 4/(y-po.m)*po.m
    coords <- subset(cbind(xcord,ycord),ycord>0 & ycord<4)
    points(coords,col="red")
  })
  
  plot(index,po.m,ylim=c(ylim.low,ylim.upp),xlab=xl,ylab=expression(theta),main="100*Alpha% Intervals for Posterior Mean",cex=log(se+2)/cx,col="red",pch=19)
  sapply(1:length(y),function(j){lines(rep(index[j],2),c(po.low[j],po.upp[j]),lwd=0.5)})
  points(index,po.low,cex=1.5,pch="-")
  points(index,po.upp,cex=1.5,pch="-")
  points(index,y,cex=log(se+2)/cx)
  if(length(unique(pr.m))==1){
    abline(h=pr.m,col=4)
  }else{
    points(index,pr.m,col=4,pch="-",cex=2)
  }
  legend("topright",pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior mean","sample mean","prior mean"),seg.len=0.5)
  
}


coverage.gbp<-function(x, nsim=100, detail=F,...){

  # assign empty space for the result to be input
  coverage<-matrix(NA, nr=length(x$se), nc=nsim)
  coverage.l<-matrix(NA, nr=length(x$se), nc=nsim)
  coverage.u<-matrix(NA, nr=length(x$se), nc=nsim)

  # 10 means 1-0 criterion that is 1 if interval icludes true parameter, 0 if not
  coverage10<-matrix(NA, nr=length(x$se), nc=nsim)
  coverage10.l<-matrix(NA, nr=length(x$se), nc=nsim)
  coverage10.u<-matrix(NA, nr=length(x$se), nc=nsim)

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
    sim.p<-matrix(rbeta(length(pre.p)*nsim, r*pre.p, r*(1-pre.p)), nr=length(n))
    sim.p.l<-matrix(rbeta(length(pre.p)*nsim, r.l*pre.p, r.l*(1-pre.p)), nr=length(n))
    sim.p.u<-matrix(rbeta(length(pre.p)*nsim, r.u*pre.p, r.u*(1-pre.p)), nr=length(n))

    # 3. generate z (data) matrix
    sim.z<-matrix(rbinom(nrow(sim.p)*nsim, n, sim.p), nr=length(n))
    sim.z.l<-matrix(rbinom(nrow(sim.p.l)*nsim, n, sim.p.l), nr=length(n))
    sim.z.u<-matrix(rbinom(nrow(sim.p.u)*nsim, n, sim.p.u), nr=length(n))

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
    sim.lambda<-matrix(rgamma(length(pre.lambda)*nsim, pre.lambda*r, r), nr=length(n))
    sim.lambda.l<-matrix(rgamma(length(pre.lambda)*nsim, pre.lambda*r.l, r.l), nr=length(n))
    sim.lambda.u<-matrix(rgamma(length(pre.lambda)*nsim, pre.lambda*r.u, r.u), nr=length(n))

    # 3. generate z (data) matrix
    sim.z<-matrix(rpois(nrow(sim.lambda)*nsim, n*sim.lambda), nr=length(n))
    sim.z.l<-matrix(rpois(nrow(sim.lambda.l)*nsim, n*sim.lambda.l), nr=length(n))
    sim.z.u<-matrix(rpois(nrow(sim.lambda.u)*nsim, n*sim.lambda.u), nr=length(n))

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
    sim.mu<-matrix(rnorm(length(pre.mu)*nsim, pre.mu, sqrt(A)), nr=length(se))
    sim.mu.l<-matrix(rnorm(length(pre.mu)*nsim, pre.mu, sqrt(A.l)), nr=length(se))
    sim.mu.u<-matrix(rnorm(length(pre.mu)*nsim, pre.mu, sqrt(A.u)), nr=length(se))

    # 3. generate y (data) matrix
    sim.y<-matrix(rnorm(nrow(sim.mu)*nsim, sim.mu, se), nr=length(se))
    sim.y.l<-matrix(rnorm(nrow(sim.mu.l)*nsim, sim.mu.l, se), nr=length(se))
    sim.y.u<-matrix(rnorm(nrow(sim.mu.u)*nsim, sim.mu.u, se), nr=length(se))

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
