print.gbp<-function(object,...){

	if(object$model=="br"){
		if(is.na(object$prior.mean)){
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean.hat=object$p0.hat, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$p.hat, se.post.mean=object$se.p.hat, intv.low=object$p.hat.low, intv.upp=object$p.hat.upp)
		}else{
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean.hat=object$prior.mean, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$p.hat, se.post.mean=object$se.p.hat, intv.low=object$p.hat.low, intv.upp=object$p.hat.upp)	}

	}else if (object$model=="pr"){
		if(is.na(object$prior.mean)){
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean.hat=object$lambda0.hat, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$lambda.hat, se.post.mean=object$se.lambda.hat, intv.low=object$lambda.hat.low, intv.upp=object$lambda.hat.upp)
		}else{
			temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean.hat=object$prior.mean, shrinkage=object$B.hat, se.shrinkage=object$se.B.hat, post.mean=object$lambda.hat, se.post.mean=object$se.lambda.hat, intv.low=object$lambda.hat.low, intv.upp=object$lambda.hat.upp)	
		}
	}
	temp.mean<-colMeans(temp)
	res<-data.frame(rbind(temp,temp.mean),row.names=c(rownames(temp),"= Mean ="))
	print(round(res,3))
}