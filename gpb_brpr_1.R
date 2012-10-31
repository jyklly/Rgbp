######### BRIMM (Hyung Suk Tak and Carl N. Morris, June 2012). This code is only for Stat324. Any feedback on it will be greatly appreciated.

# we define many functions outside the main function pb
ini.r.a<-function(model,q,w){
	r.ini<-switch(model,br=q*(1-q)/var(w)-1,pr1=q/var(w))
	a.ini<- -log(r.ini)
	list(r.ini=r.ini,a.ini=a.ini)
}

# brimm log likelihood function of alpha (when prior.mean is known)
br.llik.prior.known<-function(a,prior.mean){
	t0<-prior.mean*exp(-a)
	t1<-(1-prior.mean)*exp(-a)
	t2<-exp(-a)
	if(any(c(t0,t1,t2,n-z+t1,t0+z)<=0)){
		print("The components of lgamma should be positive");stop()
	}else{
		sum(lgamma(t0+z)-lgamma(t0)+lgamma(n-z+t1)-lgamma(t1)+lgamma(t2)-lgamma(n+t2))
	}
}

# primm1 log likelihood function of alpha (when prior.mean is known)
pr1.llik.prior.known<-function(a,prior.mean){
	sum(dnbinom(z,exp(-a)*prior.mean,exp(-a)/(exp(-a)+n),log=T))
}

# alpha estimation (when prior.mean is known)
alpha.est.prior.known<-function(dif,eps,a.ini,model,prior.mean,n.iter){

	# change log likelihood depending on the model: br or pr1
	llik<-function(a){
		switch(model,br=br.llik.prior.known(a,prior.mean),pr1=pr1.llik.prior.known(a,prior.mean))
	}

	while(abs(dif)>eps){
		opti <- optim(a.ini,function(a){a+llik(a)},control=list(fnscale=-1),method="BFGS",hessian=T)
		a.new <- opti$par
		a.hessian <- opti$hessian
		dif<-a.new-a.ini
		a.ini<-a.new
		n.iter<-n.iter+1
	}
	list(a.new=a.new,a.hess=a.hessian,n.iter=n.iter)
}

# shrinkage estimation (when prior.mean is known)
shrink.est<-function(x,model,grp.name,sample.mean,result.order){
	B.hat<-exp(-x$a.new)/(n+exp(-x$a.new))
	inv.info<-(-1)*x$a.hess
	var.B.hat<-(B.hat*(1-B.hat))^2/((B.hat*(1-B.hat))+inv.info)
	se.B.hat<-sqrt(var.B.hat)
	a1.beta<-inv.info/(1-B.hat)
	a0.beta<-inv.info/B.hat
	central3.B<-ifelse(model=="br",2*(1-2*B.hat)*B.hat*(1-B.hat)/(a1.beta+a0.beta+1)/(a1.beta+a0.beta+2),NA)
	B.hat.low<-qbeta(0.025,a1.beta,a0.beta)
	B.hat.upp<-qbeta(0.975,a1.beta,a0.beta)
	if(!is.na(grp.name)){
		B.result<-data.frame(n,sample.mean,B.hat,se.B.hat,B.hat.low,B.hat.upp,a1.beta,a0.beta,row.names=grp.name)
	}else{
		B.result<-data.frame(n,sample.mean,B.hat,se.B.hat,B.hat.low,B.hat.upp,a1.beta,a0.beta)
	}
	if(result.order==T){
		B.result<-B.result[order(B.result$n),]
	}
	B.result.mean<-colMeans(B.result)
	Bf<-data.frame(rbind(B.result,B.result.mean),row.names=c(rownames(B.result),"= Mean ="))
	Bf[,1]<-round(Bf[,1],0)
	Bf[,2:8]<-round(Bf[,2:8],3)
	list(B.table=Bf,B.res=B.result,var.B.hat=var.B.hat,central3.B=central3.B)
}

# r estimation (when prior.mean is known)
r.est<-function(x,r.alpha){
	r.hat<-round(exp(-x$a.new),0)
	se.r.hat<-sqrt(solve(-x$a.hess)*exp(-2*x$a.new))
	alpha.hat<-x$a.new
	se.alpha.hat<-sqrt(solve(-x$a.hess))
	low.r.hat<-round(exp(-(alpha.hat-qnorm(r.alpha/2)*se.alpha.hat)),0)
	upp.r.hat<-round(exp(-(alpha.hat+qnorm(r.alpha/2)*se.alpha.hat)),0)
	r.result<-round(data.frame(r.hat,se.r.hat,low.r.hat,upp.r.hat,alpha.hat,se.alpha.hat),3)
	list(r.table=r.result)
}

# brimm posterior estimation (when prior.mean is known)
br.post.est<-function(x,prior.mean,y,grp.name,result.order){
	B.hat<-x$B.res[,3]
	var.B.hat<-x$var.B.hat
	p0<-prior.mean
	central3.B<-x$central3.B
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
	se.p.hat<-sqrt(var.p.hat)
	a1.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*p.hat
	a0.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*(1-p.hat)
	p.hat.low<- qbeta(0.025,a1.beta.p,a0.beta.p)
	p.hat.upp<- qbeta(0.975,a1.beta.p,a0.beta.p)
	if(!is.na(grp.name)){
		p.result<-data.frame(n,sample.mean=y,p0,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp,row.names=grp.name)
	}else{
		p.result<-data.frame(n,sample.mean=y,p0,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp)
	}
	if(result.order==T){
		p.result<-p.result[order(p.result$n),]
	}
	p.result.mean<-colMeans(p.result)
	pf<-data.frame(rbind(p.result,p.result.mean),row.names=c(rownames(p.result),"= Mean ="))
	pf[,1]<-round(pf[,1],0)
	pf[,2:ncol(pf)]<-round(pf[,2:ncol(pf)],3)
	list(post.table=pf,post.res=p.result)
}

# primm1 posterior estimation (when prior.mean is known)
pr1.post.est<-function(x,prior.mean,y,grp.name,result.order){
	B.hat<-x$B.res[,3]
	var.B.hat<-x$var.B.hat
	lambda0<-prior.mean

	lambda.hat<- y-B.hat*(y-lambda0)
	var.lambda.hat<-1/n*(lambda.hat-B.hat*y+(var.B.hat+B.hat^2)*y-(var.B.hat+B.hat^2)*lambda0)+(y-lambda0)^2*var.B.hat
	se.lambda.hat<-sqrt(var.lambda.hat)
	u.gamma<-lambda.hat^2/var.lambda.hat
	v.gamma<-lambda.hat/var.lambda.hat
	lambda.hat.low<-qgamma(0.025,shape=u.gamma,rate=v.gamma)
	lambda.hat.upp<-qgamma(0.975,shape=u.gamma,rate=v.gamma)
	if(!is.na(grp.name)){
		lambda.result<-data.frame(n,sample.mean=y,lambda0,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp,row.names=grp.name)
	}else{
		lambda.result<-data.frame(n,sample.mean=y,lambda0,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp)
	}
	if(result.order==T){
		lambda.result<-lambda.result[order(lambda.result$n),]
	}
	lambda.result.mean<-colMeans(lambda.result)
	lambdaf<-data.frame(rbind(lambda.result,lambda.result.mean),row.names=c(rownames(lambda.result),"= Mean ="))
	lambdaf[,1]<-round(lambdaf[,1],0)
	lambdaf[,2:ncol(lambdaf)]<-round(lambdaf[,2:ncol(lambdaf)],3)
	list(post.table=lambdaf,post.res=lambda.result)
}

# graph function
post.graph<-function(x,graph.legend.position,graph.order,prior.mean){
	k<-x$post.res
	legend.position<-switch(graph.legend.position,"topleft","topright","bottomleft","bottomright")
	temp<-data.frame(n=k[,1],y=k[,2],p.hat=k[,5],p.hat.low=k[,7],p.hat.upp=k[,8],p0.hat=k[,3])
	xl<-c("Indices (Groups) with input order")
	if(graph.order==T){
		temp<-temp[order(temp$n),]
		xl<-c("Indices (Groups) ordered by increasing n")
	}
	par(xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3))
	index<-1:length(n)
	cx<-(mean(log(temp$n+2))+min(log(temp$n+2)))/2
	plot(index,temp$p.hat,ylim=c(0.8*min(c(temp$p.hat.low,temp$y)),1.2*max(c(temp$p.hat.upp,temp$y))),xlab=xl,ylab="Posterior p",main="95% Probability Interval for p",cex=log(temp$n+2)/cx,col="red",pch=19)
	for (j in 1:length(n)){
		lines(rep(index[j],2),c(temp$p.hat.low[j],temp$p.hat.upp[j]),lwd=0.5)
	}
	points(index,temp$p.hat.low,cex=1.5,pch="-")
	points(index,temp$p.hat.upp,cex=1.5,pch="-")
	points(index,temp$y,cex=log(temp$n+2)/cx)
	if(!is.na(prior.mean)){
		abline(h=prior.mean,col=4)
	}else if(is.na(prior.mean) & ncol(x)>1){
		points(index,temp$p0.hat,col=4,pch="-",cex=2)
	}else{
		points(index,temp$p0.hat,type="l",col=4)
	}
	legend(legend.position,pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior p","sample mean","prior mean"),seg.len=0.5)
}

# main function
pb<-function(z,n,x=NA,prior.mean=NA,model="br",intercept=T,result.order=F,graph=T,graph.order=F,graph.legend.position=2,grp.name=NA,eps=0.0001,r.alpha=0.0833){

	y<-sample.mean<-as.vector(z/n)
	dif<-1
	n.iter<-0

	# initial values of r (depending on the model) and alpha
	a.ini<-ini.r.a(model,prior.mean,y)$a.ini

	# alpha estimation including hessian
	a.res<-alpha.est.prior.known(dif,eps,a.ini,model,prior.mean,n.iter)

	# shrinkage estimation
	B.res<-shrink.est(a.res,model,grp.name,sample.mean,result.order)
	
	# r estimation
	r.res<-r.est(a.res,r.alpha)

	# posterior estimation
	post.res<-switch(model,br=br.post.est(B.res,prior.mean,y,grp.name,result.order),pr1=pr1.post.est(B.res,prior.mean,y,grp.name,result.order))

	# graph
	if(graph==T){
		post.graph(post.res,graph.legend.position,graph.order,prior.mean)
	}

	options(scipen=3)
	output<-list(post.est=post.res$post.table,shrinkage=B.res$B.table,r.hat=r.res$r.table,n.iterations=a.res$n.iter)
	return(output)
	options(digits=7)
}