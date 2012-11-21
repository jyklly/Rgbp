gbp<-function(x, ...) UseMethod("gbp")

gbp.default<-function(x,y,X,mu,model="gr",CI=0.95,intercept=T,eps=0.0001, ...){

  res<-switch(model, 
              gr=gr(x,y,X,mu,CI), 
              br=bp(x,y,X,prior.mean=mu,model="br",CI=CI,intercept=intercept,eps=eps), 
              pr=bp(x,y,X,prior.mean=mu,model="pr",CI=CI,intercept=intercept,eps=eps) )
  
  class(res)<-"gbp"	
  res
}


print.gbp<-function(x,...){
  
  if(!identical(x$x,NA)){
    if(x$model=="gr"){
      temp<-data.frame(sample.mean=x$sample.mean, se=x$se, x=x$x, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }else{
      temp<-data.frame(sample.mean=x$sample.mean, n=x$se, x=x$x, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }
  }else{
    if(x$model=="gr"){
      temp<-data.frame(sample.mean=x$sample.mean, se=x$se, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }else{
      temp<-data.frame(sample.mean=x$sample.mean, n=x$se, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.sd=x$post.sd)
    }
  }

  temp.mean<-colMeans(temp)
  res<-data.frame(rbind(temp,temp.mean),row.names=c(rownames(temp),"= Mean ="))
  
  cat("Summary for whole observations:\n")
  cat("\n")
  print(round(res,3))
  
}


summary.gbp<-function(object,...){

  if(!identical(object$x,NA)){
    if(object$model=="gr"){
      temp<-data.frame(sample.mean=object$sample.mean, se=object$se, x=object$x, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }else{
      temp<-data.frame(sample.mean=object$sample.mean, n=object$se, x=object$x, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }
  }else{
    if(object$model=="gr"){
      temp<-data.frame(sample.mean=object$sample.mean, se=object$se, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }else{
      temp<-data.frame(sample.mean=object$sample.mean, n=object$se, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.sd=object$post.sd)
    }
  }
  
  temp.summary<-apply(temp,2,summary)
  
  a.new<-object$a.new
  a.var<-object$a.var
  if(object$model=="gr"){
    A.hat<-exp(a.new)
    se.A<-sqrt(a.var*A.hat^2)
    low.A.hat<-exp(a.new-1.96*sqrt(a.var))
    upp.A.hat<-exp(a.new+1.96*sqrt(a.var))
    result2<-data.frame(A.hat,se.A,low.A.hat,upp.A.hat)
  }else{
    r.hat<-exp(-a.new)
    se.r<-sqrt(a.var*exp(-2*a.new))
    low.r.hat<-exp(-(a.new+1.96*sqrt(a.var)))
    upp.r.hat<-exp(-(a.new-1.96*sqrt(a.var)))
    result2<-data.frame(r.hat,se.r,low.r.hat,upp.r.hat)
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
  po.se<-x$post.sd
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
  postlens <- po.se/max(sqrtV)
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
  
  plot(index,po.m,ylim=c(ylim.low,ylim.upp),xlab=xl,ylab=expression(theta),main="100*CI% Intervals for Posterior Mean",cex=log(se+2)/cx,col="red",pch=19)
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
