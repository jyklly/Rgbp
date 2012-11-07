print.gbp<-function(object,...){

######## if model is brimm

	if(object$model=="br"){
		if(is.na(object$prior.mean)){
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean.hat=object$p0.hat, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$p.hat, se.post.mean=object$se.p.hat, intv.low=object$p.hat.low, intv.upp=object$p.hat.upp)
		}else{
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean=object$prior.mean, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$p.hat, se.post.mean=object$se.p.hat, intv.low=object$p.hat.low, intv.upp=object$p.hat.upp)	}

######## if model is primm

	}else if (object$model=="pr"){
		if(is.na(object$prior.mean)){
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean.hat=object$lambda0.hat, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$lambda.hat, se.post.mean=object$se.lambda.hat, intv.low=object$lambda.hat.low, intv.upp=object$lambda.hat.upp)
		}else{
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean=object$prior.mean, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$lambda.hat, se.post.mean=object$se.lambda.hat, intv.low=object$lambda.hat.low, intv.upp=object$lambda.hat.upp)	
		}

######## if model is grimm (it is autometically grimm)

	}else if(object$model){
		if(is.na(object$mu)){
			temp<-data.frame(y=object$y, se=object$se, mu=object$mu, Bhat=object$Bhat, seB=object$seB, LCL=object$LCL, theta=object$theta, UCL=object$UCL, shat=object$shat)	
		}else{
			temp<-data.frame(y=object$y, se=object$se, mu=object$mu, Bhat=object$Bhat, seB=object$seB, LCL=object$LCL, theta=object$theta, UCL=object$UCL, shat=object$shat)	
		}
	}
	temp.mean<-colMeans(temp)
	res<-data.frame(rbind(temp,temp.mean),row.names=c(rownames(temp),"= Mean ="))
	print(round(res,3))
}