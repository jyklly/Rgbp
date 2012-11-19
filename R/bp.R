######### BRIMM (Hyung Suk Tak and Carl N. Morris, Nov 2012). This code is only for Stat324. Any feedback on it will be greatly appreciated.

# primm1 initial values (when prior.mean is known)
pr.ini.kn<-function(given){
	r.ini<-given$prior.mean/var(given$sample.mean); a.ini<- -log(r.ini)
	list(r.ini=r.ini,a.ini=a.ini)
}

# primm1 initial values (when prior.mean is unknown)
pr.ini.un<-function(given){
	y<-given$sample.mean; x.ini<-given$x.ini
	if(identical(x.ini,NA) & !given$intercept){
		print("In case we do not designate prior mean, at least one beta should be estimated (either intercept=T or at least one covariate is needed)");stop()
	}else if(identical(x.ini,NA) & given$intercept){
		x<-matrix(1, length(y), 1);	b.ini<-mean(y)
	}else if(!identical(x.ini,NA) & given$intercept){
		if(any(y==0)){y[which(y==0)]<-0.00001}
		x<-as.matrix(cbind(rep(1,length(y)),x.ini))
		b.ini<-solve(t(x)%*%x)%*%t(x)%*%log(y)
		if(any(y==0.00001)){y[which(y==0.00001)]<-0}
	}else if(!identical(x.ini,NA) & !given$intercept){
		if(any(y==0)){y[which(y==0)]<-0.00001}
		x<-x.ini
		b.ini<-solve(t(x)%*%x)%*%t(x)%*%log(y)
		if(any(y==0.00001)){y[which(y==0.00001)]<-0}
	}	
	lambda0.ini<-mean(exp(x%*%b.ini)); r.ini<-lambda0.ini/var(y)
	list(x=x,b.ini=b.ini,lambda0.ini=lambda0.ini,r.ini=r.ini,a.ini=-log(r.ini),beta.hess.ini=-solve(t(x)%*%x))
}


# brimm initial values (when prior.mean is known)
br.ini.kn<-function(given){
	r.ini<-given$prior.mean*(1-given$prior.mean)/var(given$sample.mean)-1
	a.ini<- -log(r.ini)
	list(r.ini=r.ini,a.ini=a.ini)
}

# brimm initial values (when prior.mean is unknown)
br.ini.un<-function(given){
	y<-given$sample.mean; x.ini<-given$x.ini
	if(identical(x.ini,NA) & !given$intercept){
		print("In case we do not designate prior mean, at least one beta should be estimated (either intercept=T or at least one covariate is needed)");stop()
	}else if(identical(x.ini,NA) & given$intercept){
		x<-matrix(1, length(y), 1)
		b.ini<-as.vector(log(mean(y)/(1-mean(y))))
	}else if(!identical(x.ini,NA) & given$intercept){
		if(any(y==0)){y[which(y==0)]<-0.00001}
		if(any(y==1)){y[which(y==1)]<-0.99999}		
		x<-as.matrix(cbind(rep(1,length(y)),x.ini))
		b.ini<-solve(t(x)%*%x)%*%t(x)%*%log(y/(1-y))
		if(any(y==0.00001)){y[which(y==0.00001)]<-0}
		if(any(y==0.99999)){y[which(y==0.99999)]<-1}
	}else if(!identical(x.ini,NA) & !given$intercept){
		if(any(y==0)){y[which(y==0)]<-0.00001}
		if(any(y==1)){y[which(y==1)]<-0.99999}
		x<-x.ini
		b.ini<-solve(t(x)%*%x)%*%t(x)%*%log(y/(1-y))
		if(any(y==0.00001)){y[which(y==0.00001)]<-0}
		if(any(y==0.99999)){y[which(y==0.99999)]<-1}
	}	
	p0.ini<-mean(exp(x%*%b.ini)/(1+exp(x%*%b.ini)))
	r.ini<-p0.ini*(1-p0.ini)/var(y)-1
	list(x=x,b.ini=b.ini,p0.ini=p0.ini,r.ini=r.ini,a.ini=-log(r.ini),beta.hess.ini=-solve(t(x)%*%x))
}

# brimm log likelihood function of alpha (when prior.mean is known)
br.llik.prior.kn<-function(a,given){
	n<-given$n; z<-given$z; prior.mean<-given$prior.mean
	t0<-prior.mean*exp(-a); t1<-(1-prior.mean)*exp(-a); t2<-exp(-a)
	if(any(c(t0,t1,t2,n-z+t1,t0+z)<=0)){
		print("The components of lgamma should be positive");stop()
	}else{
		sum(lgamma(t0+z)-lgamma(t0)+lgamma(n-z+t1)-lgamma(t1)+lgamma(t2)-lgamma(n+t2))
	}
}

# brimm log likelihood function of alpha and beta (when prior.mean is unknown)
br.llik.prior.un<-function(a,b,given,ini){
	n<-given$n; z<-given$z; x<-ini$x; 
	p0.hat<-exp(x%*%as.matrix(b))/(1+exp(x%*%as.matrix(b)))
	t0<-p0.hat*exp(-a); t1<-(1-p0.hat)*exp(-a); t2<-exp(-a)
	if(any(c(t0,t1,t2,n-z+t1,t0+z)<=0)){
		print("The components of lgamma should be positive");stop()
	}else{
		sum(lgamma(t0+z)-lgamma(t0)+lgamma(n-z+t1)-lgamma(t1)+lgamma(t2)-lgamma(n+t2))
	}
}

# primm1 log likelihood function of alpha (when prior.mean is known)
pr.llik.prior.kn<-function(a,given){
	z<-given$z; n<-given$n; prior.mean<-given$prior.mean
	sum(dnbinom(z,exp(-a)*prior.mean,exp(-a)/(exp(-a)+n),log=T))
}

# primm1 log likelihood function of alpha (when prior.mean is unknown)
pr.llik.prior.un<-function(a,b,given,ini){
	z<-given$z; n<-given$n; x<-ini$x
	sum(dnbinom(z,exp(-a+x%*%as.matrix(b)),exp(-a)/(exp(-a)+n),log=T))
}

# alpha estimation (when prior.mean is known)
alpha.est.prior.kn<-function(given,ini){
	# change log likelihood depending on the model: br or pr
	llik<-function(a){
		switch(given$model,br=br.llik.prior.kn(a,given),pr=pr.llik.prior.kn(a,given))
	}
	dif<-1; n.iter<-0; a.ini<-ini$a.ini; eps<-given$eps
	while(dif>eps){
		opti <- optim(a.ini,function(a){a+llik(a)},control=list(fnscale=-1),method="BFGS",hessian=T)
		a.new <- opti$par
		a.hess <- opti$hessian
		dif<-abs(a.new-a.ini)
		a.ini<-a.new
		n.iter<-n.iter+1
	}
	list(a.new=a.new,a.hess=a.hess,n.iter=n.iter,beta.new=NA,beta.hess=NA)
}

# alpha estimation (when prior.mean is unknown)
alpha.est.prior.un<-function(given,ini){
	# change log likelihood depending on the model: br or pr
	llik<-function(a,b){
		switch(given$model,br=br.llik.prior.un(a,b,given,ini),pr=pr.llik.prior.un(a,b,given,ini))
	}
	dif<-1; n.iter<-0; a.ini<-ini$a.ini; b.ini<-ini$b.ini; beta.hess<-ini$beta.hess.ini; eps<-given$eps
	while(dif>eps){
		adm.a<-function(a){
			a+llik(a,b.ini)-0.5*log(det(-beta.hess))
		} # laplace approximation
		opti <- optim(a.ini,adm.a,control=list(fnscale=-1),method="BFGS",hessian=T)
		a.new <- opti$par; a.hess <- opti$hessian
		marginal.b<-function(b){
			llik(a.new,b)-0.5*log(det(-a.hess))
		} # laplace approximation
		if(ncol(as.matrix(ini$x))>1){
			opti2 <- optim(b.ini,marginal.b,control=list(fnscale=-1),hessian=T)
			beta.new <- opti2$par; beta.hess <- opti2$hessian
		}else{
			opti2 <- optim(b.ini,marginal.b,control=list(fnscale=-1),method="BFGS",hessian=T)
			beta.new <- opti2$par; beta.hess <- opti2$hessian
		}
		dif<-max(abs(a.new-a.ini),abs(beta.new-b.ini))
		b.ini<-beta.new; a.ini<-a.new; n.iter<-n.iter+1
	}
	list(a.new=a.new,a.hess=a.hess,beta.new=beta.new,beta.hess=beta.hess,n.iter=n.iter)
}

# shrinkage estimation
shrink.est<-function(a.res,given){	
	a.new<-a.res$a.new
	a.hess<-a.res$a.hess
	B.hat<-exp(-a.new)/(given$n+exp(-a.new))
	inv.info<- -a.hess
	var.B.hat<-(B.hat*(1-B.hat))^2/((B.hat*(1-B.hat))+inv.info)
	se.B.hat<-sqrt(var.B.hat)
	a1.beta<-inv.info/(1-B.hat)
	a0.beta<-inv.info/B.hat
	central3.B<-2*(1-2*B.hat)*B.hat*(1-B.hat)/(a1.beta+a0.beta+1)/(a1.beta+a0.beta+2) # for brimm
	B.hat.low<-qbeta((1-given$CI)/2,a1.beta,a0.beta)
	B.hat.upp<-qbeta((1+given$CI)/2,a1.beta,a0.beta)
	list(B.hat=B.hat,inv.info=inv.info,var.B.hat=var.B.hat,se.B.hat=se.B.hat,central3.B=central3.B)
}

# brimm posterior estimation (when prior.mean is known)
br.post.est.prior.kn<-function(B.res,given){
	B.hat<-B.res$B.hat; var.B.hat<-B.res$var.B.hat; p0<-given$prior.mean; central3.B<-B.res$central3.B
	y<-given$sample.mean; n<-given$n
	p.hat<- y-B.hat*(y-p0)	# posterior p
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
	p.hat.low<- qbeta((1-given$CI)/2,a1.beta.p,a0.beta.p)
	p.hat.upp<- qbeta((1+given$CI)/2,a1.beta.p,a0.beta.p)
	list(post.mean=p.hat,post.se=se.p.hat,post.intv.low=p.hat.low,post.intv.upp=p.hat.upp,prior.mean=p0)
}

# brimm posterior estimation (when prior.mean is unknown)
br.post.est.prior.un<-function(B.res,a.res,ini,given){
	x<-as.matrix(ini$x); n<-given$n; beta.new<-a.res$beta.new; beta.hess<-a.res$beta.hess; B.hat<-B.res$B.hat;
	y<-given$sample.mean
	b1<-sapply(1:length(n),function(i){(1+exp(x[i,]%*%beta.new))/(x[i,]%*%solve(-beta.hess)%*%x[i,])+0.5})
	b0<-sapply(1:length(n),function(i){(1+exp(x[i,]%*%beta.new))/(x[i,]%*%solve(-beta.hess)%*%x[i,])/exp(x[i,]%*%beta.new)+0.5})
	# cumulants	
	k1p<-as.vector(b1/(b1+b0))
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
	se.p.hat<-sqrt(var.p.hat)
	a1.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*p.hat
	a0.beta.p<- (p.hat*(1-p.hat)/var.p.hat-1)*(1-p.hat)
	p.hat.low<- qbeta((1-given$CI)/2,a1.beta.p,a0.beta.p)
	p.hat.upp<- qbeta((1+given$CI)/2,a1.beta.p,a0.beta.p)
	list(post.mean=p.hat,post.se=se.p.hat,post.intv.low=p.hat.low,post.intv.upp=p.hat.upp,prior.mean=p0.hat)
}

# primm1 posterior estimation (when prior.mean is known)
pr.post.est.prior.kn<-function(B.res,given){
	B.hat<-B.res$B.hat; var.B.hat<-B.res$var.B.hat; lambda0<-given$prior.mean; n<-given$n;
	y<-given$sample.mean;lambda.hat<- y-B.hat*(y-lambda0)
	var.lambda.hat<-1/n*(lambda.hat-B.hat*y+(var.B.hat+B.hat^2)*y-(var.B.hat+B.hat^2)*lambda0)+(y-lambda0)^2*var.B.hat
	se.lambda.hat<-sqrt(var.lambda.hat)
	u.gamma<-lambda.hat^2/var.lambda.hat
	v.gamma<-lambda.hat/var.lambda.hat
	lambda.hat.low<-qgamma((1-given$CI)/2,shape=u.gamma,rate=v.gamma)
	lambda.hat.upp<-qgamma((1+given$CI)/2,shape=u.gamma,rate=v.gamma)
	list(post.mean=lambda.hat, post.se=se.lambda.hat, post.intv.low=lambda.hat.low, post.intv.upp=lambda.hat.upp, prior.mean=lambda0)
}


# primm1 posterior estimation (when prior.mean is unknown)
pr.post.est.prior.un<-function(B.res,a.res,ini,given){
	x<-as.matrix(ini$x); n<-given$n; y<-given$sample.mean; beta.new<-a.res$beta.new; B.hat<-B.res$B.hat;
	var.B.hat<-B.res$var.B.hat
	V<-solve(-a.res$beta.hess)
	u<-sapply(1:length(n),function(t){x[t,]%*%V%*%x[t,]})
	lambda0.hat<-exp(x%*%beta.new+1/2*u)
	lambda.hat<- y-B.hat*(y-lambda0.hat)
	vr<-exp(2*x%*%beta.new+u)*(exp(u)-1)
	var.lambda.hat<-1/n*(lambda.hat-B.hat*y+(var.B.hat+B.hat^2)*y-(var.B.hat+B.hat^2)*lambda0.hat)+(var.B.hat+B.hat^2)*(y^2-2*y*lambda0.hat+vr+lambda0.hat^2)-(B.hat*(y-lambda0.hat))^2
	se.lambda.hat<-sqrt(var.lambda.hat)
	u.gamma<-lambda.hat^2/var.lambda.hat
	v.gamma<-lambda.hat/var.lambda.hat
	lambda.hat.low<-qgamma((1-given$CI)/2,shape=u.gamma,rate=v.gamma)
	lambda.hat.upp<-qgamma((1+given$CI)/2,shape=u.gamma,rate=v.gamma)
	list(post.mean=lambda.hat, post.se=se.lambda.hat, post.intv.low=lambda.hat.low, post.intv.upp=lambda.hat.upp, prior.mean=lambda0.hat)
}

# main function
bp<-function(z,n,X,prior.mean,model="br",intercept=T,eps=0.0001,CI=0.95){

	X<-if(missing(X)){X<-NA}else{X<-X}
	prior.mean<-if(missing(prior.mean)){prior.mean<-NA}else{prior.mean<-prior.mean}

	given<-list(z=z,n=n,sample.mean=z/n,x.ini=X,prior.mean=prior.mean,model=model,intercept=intercept,eps=eps,CI=CI)

	# initial values
	if(is.na(prior.mean)){
		ini<-switch(model,br=br.ini.un(given),pr=pr.ini.un(given))
	}else{
		ini<-switch(model,br=br.ini.kn(given),pr=pr.ini.kn(given))
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

	output<-list(sample.mean=given$sample.mean,se=given$n,prior.mean=post.res$prior.mean, shrinkage=B.res$B.hat, se.shrinkage=B.res$se.B.hat, post.mean=post.res$post.mean, post.se=post.res$post.se, post.intv.low=post.res$post.intv.low, post.intv.upp=post.res$post.intv.upp, model=model, x=X, beta.new=a.res$beta.new, beta.hess=a.res$beta.hess, intercept=intercept, a.new=a.res$a.new, a.var=1/B.res$inv.info)
	output
}
