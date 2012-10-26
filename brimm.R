######### BRIMM (Hyung Suk Tak and Carl N. Morris, June 2012). This code is only for Stat324. Any feedback on it will be greatly appreciated.

brimm<-function(z,n,x=NA,prior.mean=NA,intercept=T,result.order=F,graph=T,graph.order=F,graph.legend.position=2,grp.name=NA,eps=0.0001,r.alpha=0.0833){

	y<-sample.mean<-z/n
	niter<-1

# if we know p0, we need not estimate it
	if(!is.na(prior.mean)){
# initial values
		p0<-prior.mean		
		r<-p0*(1-p0)/var(y)-1
		a<- -log(r)
# log likelihood function of alpha
		llik<-function(a){
			t0<-p0*exp(-a)
			t1<-(1-p0)*exp(-a)
			t2<-exp(-a)
			if(any(c(t0,t1,t2,n-z+t1,t0+z)<=0)){
				print("The components of lgamma should be positive");stop()
			}else{
				sum(lgamma(t0+z)-lgamma(t0)+lgamma(n-z+t1)-lgamma(t1)+lgamma(t2)-lgamma(n+t2))
			}
		}
# estimation of alpha
		while(TRUE){
			adm.a<-function(a){
				a+llik(a);
			}
			opti <- optim(a,adm.a,control=list(fnscale=-1),method="BFGS",hessian=T)
			a.new <- opti$par
			a.hessian <- opti$hessian
			if(any(abs(a.new-a)>eps)){
				a<-a.new;niter<-niter+1
			}else{
				break
			}
		}
# output: shrinkage (ADM)
		B.hat<-exp(-a.new)/(n+exp(-a.new))
		inv.info<-(-1)*a.hessian
		var.B.hat<-(B.hat*(1-B.hat))^2/((B.hat*(1-B.hat))+inv.info)
		se.B.hat<-sqrt(var.B.hat)
		a1.beta<-inv.info/(1-B.hat)
		a0.beta<-inv.info/B.hat
		central3.B<-2*(1-2*B.hat)*B.hat*(1-B.hat)/(a1.beta+a0.beta+1)/(a1.beta+a0.beta+2)
		B.hat.low<- qbeta(0.025,a1.beta,a0.beta)
		B.hat.upp<- qbeta(0.975,a1.beta,a0.beta)
		if(!is.na(grp.name)){
			B.result<-data.frame(n,sample.mean,B.hat,se.B.hat,B.hat.low,B.hat.upp,a1.beta,a0.beta,row.names=grp.name)
		}else{
			B.result<-data.frame(n,sample.mean,B.hat,se.B.hat,B.hat.low,B.hat.upp,a1.beta,a0.beta)
		}
		if(result.order==T){
			B.result<-B.result[order(B.result$n),]
		}
		B.result.mean<-colMeans(B.result);
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

# output: Posterior mean and variance
# posterior p
		p.hat<-(y-B.hat*(y-p0))
# for posterior variance
		c1<-y*(1-y)
		c2<-(-3*y^2+2*(1+p0)*y-p0)
		c3<-(y-p0)*(3*y-1-p0)
		c4<-(y-p0)^2
		d1<-c1+c3*B.hat^2
		d2<-c2+2*c3*B.hat-3*c4*B.hat^2
		d3<-c3-3*c4*B.hat
		d4<-c4;
# posterior variance
		var.p.hat<-1/n*(d1-d2*B.hat-d3*var.B.hat+d4*central3.B)
		se.p.hat<-sqrt(var.p.hat)
		a1.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*p.hat
		a0.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*(1-p.hat)
		p.hat.low<- qbeta(0.025,a1.beta.p,a0.beta.p)
		p.hat.upp<- qbeta(0.975,a1.beta.p,a0.beta.p)
		p0.hat<-p0
		if(!is.na(grp.name)){
			p.result<-data.frame(n,sample.mean,p0,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp,row.names=grp.name)
		}else{
			p.result<-data.frame(n,sample.mean,p0,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp)
		}
		if(result.order==T){
			p.result<-p.result[order(p.result$n),]
		}
		p.result.mean<-colMeans(p.result)
		pf<-data.frame(rbind(p.result,p.result.mean),row.names=c(rownames(p.result),"= Mean ="))
		pf[,1]<-round(pf[,1],0)
		pf[,2:ncol(pf)]<-round(pf[,2:ncol(pf)],3)

######### Print!
		options(scipen=3)
		output<-list(posterior.p=pf,shrinkage=Bf,r.hat=r.result,n.iterations=niter)
	}
######### If we do not know p0, we need to estimate it
	if(is.na(prior.mean)){
# covariate matrix
		x.ini<-x
		if(identical(x,NA) & intercept==F){
			print("In case we do not designate prior mean, at least one beta should be estimated (either intercept or covariate is needed)");stop()
		}else if(identical(x,NA) & intercept==T){
			x<-matrix(1, length(n), 1)
			colnames(x)="b0"
			b<-as.vector(log(mean(y)/(1-mean(y))))
		}else if(!identical(x,NA) & intercept==T){
			if(any(y==0)){y[which(y==0)]<-0.00001}
			if(any(y==1)){y[which(y==1)]<-0.99999}		
			x<-as.matrix(cbind(rep(1,length(z)),x))
			b<-solve(t(x)%*%x)%*%t(x)%*%log(y/(1-y))
			if(any(y==0.00001)){y[which(y==0.00001)]<-0}
			if(any(y==0.99999)){y[which(y==0.99999)]<-1}
		}else if(!identical(x,NA) & intercept==F){
			if(any(y==0)){y[which(y==0)]<-0.00001}
			if(any(y==1)){y[which(y==1)]<-0.99999}
			b<-solve(t(x)%*%x)%*%t(x)%*%log(y/(1-y))
			if(any(y==0.00001)){y[which(y==0.00001)]<-0}
			if(any(y==0.99999)){y[which(y==0.99999)]<-1}
		}
# initial values
		m<-ncol(x)		
		p0.ini<-mean(exp(x%*%b)/(1+exp(x%*%b)))
		r<-p0.ini*(1-p0.ini)/var(y)-1
		a<- -log(r)
		beta.hessian<- -solve(t(x)%*%x)
# log likelihood function of alpha and beta
	    llik<-function(a,b){
	    	p1<-exp(x%*%b)/(1+exp(x%*%b))
			t0<-p1*exp(-a)
			t1<-(1-p1)*exp(-a)
			t2<-exp(-a)
			if(any(c(t0,t1,t2,n-z+t1,t0+z)<=0)){
				print("The components of lgamma should be positive");stop()
			}else{
				sum(lgamma(t0+z)-lgamma(t0)+lgamma(n-z+t1)-lgamma(t1)+lgamma(t2)-lgamma(n+t2))
			}
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
				break
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
		B.hat<-exp(-a.new)/(n+exp(-a.new))
		inv.info<-(-1)*a.hessian
		var.B.hat<-(B.hat*(1-B.hat))^2/((B.hat*(1-B.hat))+inv.info)
		se.B.hat<-sqrt(var.B.hat)
		a1.beta<-inv.info/(1-B.hat)
		a0.beta<-inv.info/B.hat
		central3.B<-2*(1-2*B.hat)*B.hat*(1-B.hat)/(a1.beta+a0.beta+1)/(a1.beta+a0.beta+2)
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
		B.result.mean<-as.vector(colMeans(B.result))
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
# output: posterior p
# cumulants		
		b1<-sapply(1:length(n),function(i){(1+exp(x[i,]%*%beta.new))/(x[i,]%*%solve(-beta.hessian)%*%x[i,])+0.5})
		b0<-sapply(1:length(n),function(i){(1+exp(x[i,]%*%beta.new))/(x[i,]%*%solve(-beta.hessian)%*%x[i,])/exp(x[i,]%*%beta.new)+0.5})
		k1p<-as.vector(b1/(b1+b0))
		p0.hat<-k1p
		k2p<-as.vector(k1p*(1-k1p)/(b1+b0+1))
		k1b<-as.vector(B.hat)
		k2b<-as.vector(var.B.hat)
		k3b<-as.vector(central3.B)
# posterior p
		p.hat<-(1-B.hat)*y+B.hat*p0.hat
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
		se.p.hat<-sqrt(var.p.hat)
		a1.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*p.hat
		a0.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*(1-p.hat)
		p.hat.low<- qbeta(0.025,a1.beta.p,a0.beta.p)
		p.hat.upp<- qbeta(0.975,a1.beta.p,a0.beta.p)

		if(identical(x.ini,NA)){
			if(!is.na(grp.name)){
				p.result<-data.frame(n,sample.mean,p0.hat,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp,row.names=grp.name)
			}else{
				p.result<-data.frame(n,sample.mean,p0.hat,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp)
			}
		}else{
			if(!is.na(grp.name)){
				p.result<-data.frame(n,x.ini,sample.mean,p0.hat,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp,row.names=grp.name)
			}else{
				p.result<-data.frame(n,x.ini,sample.mean,p0.hat,B.hat,p.hat,se.p.hat,p.hat.low,p.hat.upp)
			}
		}		
		if(result.order==T){
			p.result<-p.result[order(p.result$n),]
		}
		p.result.mean<-colMeans(p.result)
		pf<-data.frame(rbind(p.result,p.result.mean),row.names=c(rownames(p.result),"= Mean ="))
		pf[,1]<-round(pf[,1],0)
		pf[,2:ncol(pf)]<-round(pf[,2:ncol(pf)],3)		
######### Print!
		options(scipen=3)
		output<-list(posterior.p=pf,shrinkage=Bf,reg.coeff=beta.result,r.hat=r.result,n.iterations=niter)
	}
######### Plotting Probability Interval Graph
	if(graph==T){
		legend.position<-switch(graph.legend.position,"topleft","topright","bottomleft","bottomright")
		temp<-data.frame(n,y,p.hat,p.hat.low,p.hat.upp,p0.hat)
		xl<-c("Indices (Groups)")
		if(graph.order==T){
			temp<-temp[order(temp$n),]
			xl<-c("Indices (Groups) ordered by increasing n")
		}
		par(xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3))
		index<-1:length(n)
		cx<-(mean(log(temp$n+2))+min(log(temp$n+2)))/2
		plot(index,temp$p.hat,ylim=c(0.8*min(c(p.hat.low,y)),1.2*max(c(p.hat.upp,y))),xlab=xl,ylab="Posterior p",main="95% Probability Interval for p",cex=log(temp$n+2)/cx,col="red",pch=19)
		for (j in 1:length(n)){
			lines(rep(index[j],2),c(temp$p.hat.low[j],temp$p.hat.upp[j]),lwd=0.5)
		}
		points(index,temp$p.hat.low,cex=1.5,pch="-")
		points(index,temp$p.hat.upp,cex=1.5,pch="-")
		points(index,temp$y,cex=log(temp$n+2)/cx)
		if(!is.na(prior.mean)){
			abline(h=p0,col=4)
		}else if(is.na(prior.mean) & ncol(x)>1){
			points(index,temp$p0.hat,col=4,pch="-",cex=2)
		}else{
			points(index,temp$p0.hat,type="l",col=4)
		}
		legend(legend.position,pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior p","sample mean","prior mean"),seg.len=0.5)
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