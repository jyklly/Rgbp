plot.gbp<-function(object,legend.pos=2,gp.ord=F){

	y<-object$sample.mean
	se<-object$se
	pr.m<-object$prior.mean
	po.m<-object$post.mean
	po.se<-object$post.se
	po.low<-object$post.intv.low
	po.upp<-object$post.intv.upp
	cx<-(mean(log(se+2))+min(log(se+2)))/2
	index<-1:length(se)
	ylim.low<-ifelse(min(po.low,y)>=0, 0.8*min(po.low,y), 1.2*min(po.low,y))
	ylim.upp<-ifelse(max(po.upp,y)>=0, 1.2*max(po.upp,y), 0.8*max(po.upp,y))

	par(mfrow=c(1,2),xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3),las=1,ps=13)
	leg.pos<-switch(legend.pos,"topleft","topright","bottomleft","bottomright")
	xl<-c("Indices (Groups) by the order of data input")

	sqrtV <- se
	sdlens <- sqrtV/max(sqrtV)
	postlens <- po.se/max(sqrtV)
	sunflowerplot(rep(4,length(y))~y,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1), yaxt="n", col.lab="white", main="Shrinkage Plot")
	if(length(unique(pr.m))==1)
		points(pr.m[1],0,col="darkviolet",pch=2,cex=4)
	legend(leg.pos,c("se"),col="blue",lty=1,seg.len=0.5,lwd=2)
	sunflowerplot(rep(0,length(y))~po.m,add=T)
	abline(h=4)
	abline(h=0)
	axis(2,c(0,4),c(expression(hat(theta)),expression(bar(y))), cex.axis=1.1)
	for(i in 1:length(y)){
		lines(c(y[i],po.m[i]),c(4,0),xlim=c(min(y)-abs(min(y))*0,max(y)+abs(max(y))*0))
		lines(c(y[i],y[i]+sdlens[i]*sd(y)*0.4),c(4,4+sdlens[i]),col="blue")
	##posterior variance lines
	lines(c(po.m[i]-postlens[i]*sd(y)*0.4,po.m[i]),c(0-postlens[i],0),col="blue")
	xcord <- (4*po.m[i]/(y[i]-po.m[i])-4*po.m/(y-po.m))/(4/(y[i]-po.m[i])-4/(y-po.m))
	ycord <- 4/(y-po.m)*xcord - 4/(y-po.m)*po.m
	coords <- subset(cbind(xcord,ycord),ycord>0 & ycord<4)
	points(coords,col="red")
	}

	plot(index,po.m,ylim=c(ylim.low,ylim.upp),xlab=xl,ylab=expression(theta),main="100(1-CI)% Intervals for Posterior Mean",cex=log(se+2)/cx,col="red",pch=19)
	sapply(1:length(y),function(j){lines(rep(index[j],2),c(po.low[j],po.upp[j]),lwd=0.5)})
	points(index,po.low,cex=1.5,pch="-")
	points(index,po.upp,cex=1.5,pch="-")
	points(index,y,cex=log(se+2)/cx)
	if(length(unique(pr.m))==1){
		abline(h=pr.m,col=4)
	}else{
		points(index,pr.m,col=4,pch="-",cex=2)
	}
	legend(leg.pos,pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior mean","sample mean","prior mean"),seg.len=0.5)

}