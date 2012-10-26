######### PRIMM1 (Hyung Suk Tak and Carl N. Morris, June 2012). This code is only for Stat324. Any feedback on it will be greatly appreciated.

primm1<-function(z,n,x=NA,prior.mean=NA,intercept=T,result.order=F,graph=T,graph.order=F,graph.legend.position=2,grp.name=NA,eps=0.0001,r.alpha=0.0833){
	y<-sample.mean<-z/n
	niter<-1
######### If we know lambda0 (prior.mean), we do not need to estimate it
	if(!is.na(prior.mean)){
# initial values
		lambda0<-prior.mean
		r<-lambda0/var(y)
		a<- -log(r)
# log likelihood function of alpha
		llik<-function(a){
			sum(dnbinom(z,exp(-a)*lambda0,exp(-a)/(exp(-a)+n),log=T))
		}
######### Estimation of alpha	
		while(TRUE){
			adm.a<-function(a){
				a+llik(a)
			}
			opti <- optim(a,adm.a,control=list(fnscale=-1),method="BFGS",hessian=T)
			a.new <- opti$par
			a.hessian <- opti$hessian
			if(any(abs(a.new-a)>eps)){
				a<-a.new;niter<-niter+1
			}else{
				break;
			}
		}
# output: shrinkage (ADM)
		B.hat<-exp(-a.new)/(n+exp(-a.new))
		inv.info<-(-1)*a.hessian
		var.B.hat<-(B.hat*(1-B.hat))^2/((B.hat*(1-B.hat))+inv.info)
		se.B.hat<-sqrt(var.B.hat)
		a1.beta<-inv.info/(1-B.hat)
		a0.beta<-inv.info/B.hat
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
# output: r.hat (with (1-r.alpha)*100% interval) and alpha.hat
		r.hat<-round(exp(-a.new),0)
		se.r.hat<-round(sqrt(solve(inv.info)*exp(-2*a.new)),3)
		alpha.hat<-round(a.new,3)
		se.alpha.hat<-round(sqrt(solve(inv.info)),3)
		low.r.hat<-round(exp(-(alpha.hat-qnorm(r.alpha/2)*se.alpha.hat)),0)
		upp.r.hat<-round(exp(-(alpha.hat+qnorm(r.alpha/2)*se.alpha.hat)),0)
		r.result<-data.frame(r.hat,se.r.hat,low.r.hat,upp.r.hat,alpha.hat,se.alpha.hat)

# output: posterior mean and variance
# posterior lambda
		lambda.hat<-(y-B.hat*(y-lambda0))
		var.lambda.hat<-1/n*(lambda.hat-B.hat*y+(var.B.hat+B.hat^2)*y-(var.B.hat+B.hat^2)*lambda0)+(y-lambda0)^2*var.B.hat
		se.lambda.hat<-sqrt(var.lambda.hat)
		u.gamma<-lambda.hat^2/var.lambda.hat
		v.gamma<-lambda.hat/var.lambda.hat
		lambda.hat.low<-qgamma(0.025,shape=u.gamma,rate=v.gamma)
		lambda.hat.upp<-qgamma(0.975,shape=u.gamma,rate=v.gamma)
		if(!is.na(grp.name)){
			lambda.result<-data.frame(n,sample.mean,lambda0,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp,row.names=grp.name)
		}else{
			lambda.result<-data.frame(n,sample.mean,lambda0,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp)
		}
		if(result.order==T){
			lambda.result<-lambda.result[order(lambda.result$n),]
		}
		lambda.result.mean<-colMeans(lambda.result)
		lambdaf<-data.frame(rbind(lambda.result,lambda.result.mean),row.names=c(rownames(lambda.result),"= Mean ="))
		lambdaf[,1]<-round(lambdaf[,1],0)
		lambdaf[,2:ncol(lambdaf)]<-round(lambdaf[,2:ncol(lambdaf)],3)
		lambda0.hat<-lambda0
######### Print!
		options(scipen=3)
		output<-list(posterior.lambda=lambdaf,shrinkage=Bf,r.hat=r.result,n.iterations=niter)
	}else{
######### If we do not know lambda0, we need to estimate it by regression
# covariate matrix
		x.ini<-x
		if(identical(x,NA) & intercept==F){
			print("In case we do not designate prior mean, at least one beta should be estimated (either intercept or covariate is needed)");stop()
		}else if(identical(x,NA) & intercept==T){
			x<-matrix(1, length(n), 1)
			colnames(x)="b0"
			b<-as.vector(log(mean(y)))
		}else if(!identical(x,NA) & intercept==T){
			if(any(y==0)){y[which(y==0)]<-0.00001}
			x<-as.matrix(cbind(rep(1,length(z)),x))
			b<-solve(t(x)%*%x)%*%t(x)%*%log(y)
			if(any(y==0.00001)){y[which(y==0.00001)]<-0}
		}else if(!identical(x,NA) & intercept==F){
			if(any(y==0)){y[which(y==0)]<-0.00001}
			b<-solve(t(x)%*%x)%*%t(x)%*%log(y)
			if(any(y==0.00001)){y[which(y==0.00001)]<-0};
		}	
# initial values
		m<-ncol(x)		
		lambda0.ini<-mean(exp(x%*%b))
		r<-lambda0.ini/var(y)
		a<- -log(r)
		beta.hessian<- -solve(t(x)%*%x)
# log Likelihood Function
		llik<-function(a,b){
			sum(dnbinom(z,exp(-a+x%*%b),exp(-a)/(exp(-a)+n),log=T))
		}
# estimating beta and r
		while(TRUE){
# laplace approximation
			adm.a<-function(a){
				a+llik(a,b)-0.5*log(det(-beta.hessian))
			}
			opti <- optim(a,adm.a,control=list(fnscale=-1),method="BFGS",hessian=T)
			a.new <- opti$par
			a.hessian <- opti$hessian
# laplace approximation
			nr.b<-function(b){
				llik(a.new,b)-0.5*log(det(-a.hessian))
			}
			if(ncol(x)>1){
				opti2 <- optim(b,nr.b,control=list(fnscale=-1),hessian=T)
				beta.new <- opti2$par
				beta.hessian <- opti2$hessian
			}else{
				opti2 <- optim(b,nr.b,control=list(fnscale=-1),method="BFGS",hessian=T)
				beta.new <- opti2$par
				beta.hessian <- opti2$hessian
			}
			if(any(c(abs(a.new-a),abs(beta.new-b))>eps)){
				b<-beta.new;a<-a.new;niter<-niter+1
			}else{
				break;
			}
		}
# output: regression coefficients
		estimate<-as.vector(beta.new)
		if(intercept==T){
			names(estimate)<-paste("beta", 0:(length(beta.new)-1), sep = "")
     	}else{
			names(estimate)<-paste("beta", 1:(length(beta.new)), sep = "")
		}
		se<-sqrt(diag(solve(-beta.hessian)))
		z.val<-estimate/se
		p.val<-ifelse(z.val<=0,2*pnorm(z.val),2*pnorm(-z.val))
		beta.result<-round(data.frame(estimate,se,z.val,p.val),3)	
# output: shrinkage
		B.hat<-exp(-a.new)/(n+exp(-a.new));
		inv.info<-(-1)*a.hessian;
		var.B.hat<-(B.hat*(1-B.hat))^2/((B.hat*(1-B.hat))+inv.info);
		se.B.hat<-sqrt(var.B.hat);
		a1.beta<-inv.info/(1-B.hat);
		a0.beta<-inv.info/B.hat;
		B.hat.low<-qbeta(0.025,a1.beta,a0.beta);
		B.hat.upp<-qbeta(0.975,a1.beta,a0.beta);
		if(!is.na(grp.name)){
			B.result<-data.frame(n,y,B.hat,se.B.hat,B.hat.low,B.hat.upp,a1.beta,a0.beta,row.names=grp.name)
		}else{
			B.result<-data.frame(n,y,B.hat,se.B.hat,B.hat.low,B.hat.upp,a1.beta,a0.beta)
		}
		if(result.order==T){
			B.result<-B.result[order(B.result$n),]
		}
		B.result.mean<-colMeans(B.result);
		Bf<-data.frame(rbind(B.result,B.result.mean),row.names=c(rownames(B.result),"= Mean ="))
		Bf[,1]<-round(Bf[,1],0)
		Bf[,2:8]<-round(Bf[,2:8],3)
# output: r.hat
		r.hat<-round(exp(-a.new),0)
		se.r.hat<-round(sqrt(solve(inv.info)*exp(-2*a.new)),3)
		alpha.hat<-round(a.new,3)
		se.alpha.hat<-round(sqrt(solve(inv.info)),3)
		low.r.hat<-round(exp(-(alpha.hat-qnorm((1-r.alpha)/2)*se.alpha.hat)),0)
		upp.r.hat<-round(exp(-(alpha.hat+qnorm((1-r.alpha)/2)*se.alpha.hat)),0)
		r.result<-data.frame(r.hat,se.r.hat,low.r.hat,upp.r.hat,alpha.hat,se.alpha.hat)
# output: posterior lambda
		V<-solve(-beta.hessian)
		u<-sapply(1:length(n),function(t){x[t,]%*%V%*%x[t,]})
		lambda.hat<-(y-B.hat*(y-exp(x%*%beta.new+1/2*u)))
		lambda0.hat<-exp(x%*%beta.new+1/2*u)
		vr<-exp(2*x%*%beta.new+u)*(exp(u)-1)
		var.lambda.hat<-1/n*(lambda.hat-B.hat*y+(var.B.hat+B.hat^2)*y-(var.B.hat+B.hat^2)*lambda0.hat)+(var.B.hat+B.hat^2)*(y^2-2*y*lambda0.hat+vr+lambda0.hat^2)-(B.hat*(y-lambda0.hat))^2
		se.lambda.hat<-sqrt(var.lambda.hat)
		u.gamma<-lambda.hat^2/var.lambda.hat
		v.gamma<-lambda.hat/var.lambda.hat
		lambda.hat.low<-qgamma(0.025,shape=u.gamma,rate=v.gamma)
		lambda.hat.upp<-qgamma(0.975,shape=u.gamma,rate=v.gamma)
		if(identical(x.ini,NA)){
			if(!is.na(grp.name)){
				lambda.result<-data.frame(n,sample.mean,lambda0.hat,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp,row.names=grp.name)
			}else{
				lambda.result<-data.frame(n,sample.mean,lambda0.hat,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp)
			}
		}else{
			if(!is.na(grp.name)){
				lambda.result<-data.frame(n,x.ini,sample.mean,lambda0.hat,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp,row.names=grp.name)
			}else{
				lambda.result<-data.frame(n,x.ini,sample.mean,lambda0.hat,B.hat,lambda.hat,se.lambda.hat,lambda.hat.low,lambda.hat.upp)
			}
		}
		if(result.order==T){
			lambda.result<-lambda.result[order(lambda.result$n),]
		}
		lambda.result.mean<-colMeans(lambda.result)
		lambdaf<-data.frame(rbind(lambda.result,lambda.result.mean),row.names=c(rownames(lambda.result),"= Mean ="))
		lambdaf[,1]<-round(lambdaf[,1],0)
		lambdaf[,2:ncol(lambdaf)]<-round(lambdaf[,2:ncol(lambdaf)],3)
######### Print!
		options(scipen=3)
		output<-list(posterior.lambda=lambdaf,shrinkage=Bf,reg.coeff=beta.result,r.hat=r.result,n.iterations=niter)
	}
######### Plotting Probability Interval Graph ordered by increasing post.lambda
	if(graph==T){
		legend.position<-switch(graph.legend.position,"topleft","topright","bottomleft","bottomright")
		temp<-data.frame(n,y,lambda.hat,lambda.hat.low,lambda.hat.upp,lambda0.hat)
		xl<-c("Indices (Groups)")
		if(graph.order==T){
			temp<-temp[order(temp$n),]
			xl<-c("Indices (Groups) ordered by increasing n")
		}
		par(xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3));
		index<-1:length(n)
		cx<-(mean(log(temp$n+2))+min(log(temp$n+2)))/2
		plot(index,temp$lambda.hat,ylim=c(0.8*min(c(lambda.hat.low,y)),1.2*max(c(lambda.hat.upp,y))),xlab=xl,ylab="Posterior Lambda",main="95% Probability Interval for Lambda",cex=log(temp$n+2)/cx,col="red",pch=19)
		for (j in 1:length(n)){
			lines(rep(index[j],2),c(temp$lambda.hat.low[j],temp$lambda.hat.upp[j]),lwd=0.5)
		}
		points(index,temp$lambda.hat.low,cex=1.5,pch="-")
		points(index,temp$lambda.hat.upp,cex=1.5,pch="-")
		points(index,temp$y,cex=log(temp$n+2)/cx)
		if(!is.na(prior.mean)){
			abline(h=lambda0.hat,col=4)
		}else if(is.na(prior.mean) & ncol(x)>1){
			points(index,temp$lambda0.hat,col=4,pch="-",cex=2)
		}else{
			points(index,temp$lambda0.hat,type="l",col=4)
		}
		legend(legend.position,pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior lambda","sample mean","prior mean"),seg.len=0.5)
	}
	return(output)
	options(digits=7)
}

# 28 June: skewness and kurtosis are removed
# 1  July: se.r.hat is added
# 5  July: missing(x) -> identical(x,NA)
# 15 July: alpha.hat, se.alpha.hat, and (1-r.alpha)*100% interval of r.hat are added
# 30 Aug : r.alpha=0.9163 -> r.alpha=0.0837 -> r.alpha=0.0833, and corrected the case p0 known (r.alpha)
# 2  Sep : y->sample.mean in the output
# 2  Sep : The length of blue line in the legend of the graph changed
# 22 Sep : y->sample.mean in the output
# 25 Oct : some of identical -> is.na