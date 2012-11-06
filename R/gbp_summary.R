summary.gbp<-function(object,...){

	if(is.na(object$prior.mean)){
		temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean.hat=object[[36]], shrinkage=object[[21]], se.shrinkage=object[[24]], post.mean=object[[30]], se.post.mean=object[[31]], intv.low=object[[32]], intv.upp=object[[33]])
	}else{
		temp<-data.frame(n=object$n, sample.mean=object$sample.mean, prior.mean=object[[5]], shrinkage=object[[21]], se.shrinkage=object[[24]], post.mean=object[[30]], se.post.mean=object[[31]], intv.low=object[[32]], intv.upp=object[[33]])
	}

	temp.mean<-colMeans(temp)
	res<-data.frame(rbind(temp,temp.mean),row.names=c(rownames(temp),"= Mean ="))
	print(round(res,3))

}