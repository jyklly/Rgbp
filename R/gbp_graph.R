plot.gbp<-function(object,legend.pos=2,gp.ord=F,many.grp=F){

	leg.pos<-switch(legend.pos,"topleft","topright","bottomleft","bottomright")
	prior.mean<-object$prior.mean
	if(is.na(object$prior.mean)){
		temp<-data.frame(n=object$n,y=object$sample.mean,p.hat=object[[25]],p.hat.low=object[[27]],p.hat.upp=object[[28]],p0.hat=object[[31]])
	}else{
		temp<-data.frame(n=object$n,y=object$sample.mean,p.hat=object[[19]],p.hat.low=object[[21]],p.hat.upp=object[[22]],p0.hat=object[[26]])
	}
	xl<-c("Indices (Groups) by the order of data input")
	if(gp.ord==T){
		temp<-temp[order(temp$n),]
		xl<-c("Indices (Groups) ordered by increasing n")
	}
	par(xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3))
	index<-1:length(temp$n)
	cx<-(mean(log(temp$n+2))+min(log(temp$n+2)))/2
	if(many.grp){
		plot(index,temp$p.hat,ylim=c(0.8*min(c(temp$p.hat.low,temp$y)),1.2*max(c(temp$p.hat.upp,temp$y))),xlab=xl,ylab="Posterior Mean",main="95% Probability Intervals for Posterior Mean",type="l",col="red",lwd=1.5)
		sapply(1:length(temp$n),function(j){lines(rep(index[j],2),c(temp$p.hat.low[j],temp$p.hat.upp[j]),lwd=0.5)})
		points(index,temp$p.hat.low,type="l")
		points(index,temp$p.hat.upp,type="l")
		points(index,temp$y,type="l",lty=2)
		if(!is.na(prior.mean)){
			abline(h=prior.mean,col=4)
		}else if(!identical(object$x.ini,NA)){
			points(index,temp$p0.hat,col=4,type="l")
		}else{
			points(index,temp$p0.hat,type="l",col=4)
		}
		legend(leg.pos,pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior mean","sample mean","prior mean"),seg.len=0.5)
	}else{
		plot(index,temp$p.hat,ylim=c(0.8*min(c(temp$p.hat.low,temp$y)),1.2*max(c(temp$p.hat.upp,temp$y))),xlab=xl,ylab="Posterior Mean",main="95% Probability Intervals for Posterior Mean",cex=log(temp$n+2)/cx,col="red",pch=19)
		sapply(1:length(temp$n),function(j){lines(rep(index[j],2),c(temp$p.hat.low[j],temp$p.hat.upp[j]),lwd=0.5)})
		points(index,temp$p.hat.low,cex=1.5,pch="-")
		points(index,temp$p.hat.upp,cex=1.5,pch="-")
		points(index,temp$y,cex=log(temp$n+2)/cx)
		if(!is.na(prior.mean)){
			abline(h=prior.mean,col=4)
		}else if(!identical(object$x.ini,NA)){
			points(index,temp$p0.hat,col=4,pch="-",cex=2)
		}else{
			points(index,temp$p0.hat,type="l",col=4)
		}
		legend(leg.pos,pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior mean","sample mean","prior mean"),seg.len=0.5)
	}
}