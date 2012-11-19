summary.gbp<-function(object,...){

	if(!identical(object$x,NA)){
		if(object$model=="gr"){
			temp<-data.frame(sample.mean=object$sample.mean, se=object$se, x=object$x, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.se=object$post.se)
		}else{
			temp<-data.frame(sample.mean=object$sample.mean, n=object$se, x=object$x, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.se=object$post.se)
		}
	}else{
		if(object$model=="gr"){
			temp<-data.frame(sample.mean=object$sample.mean, se=object$se, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.se=object$post.se)
		}else{
			temp<-data.frame(sample.mean=object$sample.mean, n=object$se, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.se=object$post.se)
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
		p.val<-ifelse(z.val<=0,2*pnorm(z.val),2*pnorm(-z.val))
		beta.result<-data.frame(estimate,se,z.val,p.val)
		print(list(main.summary=round(temp.summary,3),second.level.variance.summary=round(result2,3),regression.summary=round(beta.result,3)))
	}else{
		print(list(main.summary=round(temp.summary,3),second.level.variance.summary=round(result2,3)))
	}
}