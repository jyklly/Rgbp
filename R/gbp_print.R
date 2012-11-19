print.gbp<-function(x,...){
	if(!identical(x$x,NA)){
		if(x$model=="gr"){
			temp<-data.frame(sample.mean=x$sample.mean, se=x$se, x=x$x, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.se=x$post.se)
		}else{
			temp<-data.frame(sample.mean=x$sample.mean, n=x$se, x=x$x, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.se=x$post.se)
		}
	}else{
		if(x$model=="gr"){
			temp<-data.frame(sample.mean=x$sample.mean, se=x$se, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.se=x$post.se)
		}else{
			temp<-data.frame(sample.mean=x$sample.mean, n=x$se, prior.mean=x$prior.mean, shrinkage=x$shrinkage, se.shrinkage=x$se.shrinkage, post.intv.low=x$post.intv.low, post.mean=x$post.mean, post.intv.upp=x$post.intv.upp, post.se=x$post.se)
		}
	}

	temp.mean<-colMeans(temp)
	res<-data.frame(rbind(temp,temp.mean),row.names=c(rownames(temp),"= Mean ="))
	print(round(res,3))

}