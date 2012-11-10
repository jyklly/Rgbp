print.gbp<-function(object,...){

	temp<-data.frame(sample.mean=object$sample.mean, se=object$se, prior.mean=object$prior.mean, shrinkage=object$shrinkage, se.shrinkage=object$se.shrinkage, post.intv.low=object$post.intv.low, post.mean=object$post.mean, post.intv.upp=object$post.intv.upp, post.se=object$post.se)
	temp.mean<-colMeans(temp)
	res<-data.frame(rbind(temp,temp.mean),row.names=c(rownames(temp),"= Mean ="))
	print(round(res,3))

}