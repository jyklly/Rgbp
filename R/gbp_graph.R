plot.gbp<-function(object,legend.pos=2,gp.ord=F){
	par(mfrow=c(1,2),xaxs = "r", yaxs = "r", mai=c(1,0.6,1,0.3),las=1,ps=13)
	leg.pos<-switch(legend.pos,"topleft","topright","bottomleft","bottomright")
	xl<-c("Indices (Groups) by the order of data input")

	if(object$model=="gr"){
		if(is.na(object$mu)){
			# not available now
		}else{
			y<-object$y
			se<-object$se
			mu<-object$mu
			thetahat<-object$theta
			shat<-object$shat
			LCL<-object$LCL
			UCL<-object$UCL
			index<-1:length(y)

			fs <- 1.3
			sqrtV <- se
			sdlens <- sqrtV/max(sqrtV)
			postlens <- shat/max(sqrtV)
			sunflowerplot(rep(4,length(y))~y,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n", col.lab="white",cex.axis=fs,main="Shrinkage Plot")
			if(length(unique(mu))==1)
				points(mu[1],0,col="darkviolet",pch=2,cex=4)
			legend(leg.pos,c("SD"),col="blue",lty=1)
			sunflowerplot(rep(0,length(y))~thetahat,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n",col.lab="white",add=T)
			abline(h=4)
			abline(h=0)
			axis(2,c(0,4),c(expression(hat(theta)),"y"), cex.axis=fs*1.1)
			for(i in 1:length(y)){
				lines(c(y[i],thetahat[i]),c(4,0),xlim=c(min(y)-abs(min(y))*0,max(y)+abs(max(y))*0))
				lines(c(y[i],y[i]+sdlens[i]*sd(y)*0.4),c(4,4+sdlens[i]),col="blue")
				##posterior variance lines
				lines(c(thetahat[i]-postlens[i]*sd(y)*0.4,thetahat[i]),c(0-postlens[i],0),col="blue")
				xcord <- (4*thetahat[i]/(y[i]-thetahat[i]) - 4*thetahat/(y-thetahat))/(4/(y[i]-thetahat[i]) - 4/(y-thetahat))
				ycord <- 4/(y-thetahat)*xcord - 4/(y-thetahat)*thetahat
				coords <- subset(cbind(xcord,ycord),ycord>0 & ycord<4)
				points(coords,col="red")
			}

			plot(index,thetahat,ylim=c(0.8*min(c(LCL,y)),1.2*max(c(UCL,y))),xlab=xl,ylab=expression(mu),main=expression(paste("100(1-CI)% Posterior Intervals for",~mu)),cex=sqrt(se/length(y)),col="red",pch=19,cex.axis=fs,cex.main=fs)
j<-1
			sapply(1:length(y),function(j){lines(rep(index[j],2),c(LCL[j],UCL[j]),lwd=0.5)})
			points(index,LCL,cex=1.5,pch="-")
			points(index,UCL,cex=1.5,pch="-")
			points(index,y,cex=sqrt(se/length(y)))
			if(!is.na(mu)){
				abline(h=mu,col=4)
			}else if(is.na(mu) & ncol(x)>1){
				# not available now
			}else{
				# not available now
			}
			legend(leg.pos,pch=c(19,1,NA),col=c(2,1,4),lwd=c(NA,NA,2),c("posterior mean","sample mean","prior mean"),seg.len=0.5)
		}
	}else if(object$model=="br"){
		if(is.na(object$prior.mean)){
			# not available now
		}else{
			y<-object$sample.mean
			n<-object$n
			prior.mean<-object$prior.mean
			post.mean<-object$p.hat
			post.sd<-object$se.p.hat
			intv.low<-object$p.hat.low
			intv.upp<-object$p.hat.upp
			index<-1:length(z)
			temp<-data.frame(n,y,p.hat=post.mean,p.hat.low=intv.low,p.hat.upp=intv.upp)
			cx<-(mean(log(n+2))+min(log(n+2)))/2
			fs <- 1.3
			sqrtV <- sqrt(1/n)
			sdlens <- sqrtV/max(sqrtV)
			postlens <- post.sd/max(sqrtV)
			sunflowerplot(rep(4,length(y))~y,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n", col.lab="white",cex.axis=fs,main="Shrinkage Plot")
			if(length(unique(prior.mean))==1)
				points(prior.mean[1],0,col="darkviolet",pch=2,cex=4)
			legend(leg.pos,c("SD"),col="blue",lty=1)
			sunflowerplot(rep(0,length(y))~post.mean,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n",col.lab="white",add=T)
			abline(h=4)
			abline(h=0)
			axis(2,c(0,4),c("p","y"), cex.axis=fs*1.1)
			for(i in 1:length(y)){
				lines(c(y[i],post.mean[i]),c(4,0),xlim=c(min(y)-abs(min(y))*0,max(y)+abs(max(y))*0))
				lines(c(y[i],y[i]+sdlens[i]*sd(y)*0.4),c(4,4+sdlens[i]),col="blue")
				##posterior variance lines
				lines(c(post.mean[i]-postlens[i]*sd(y)*0.4,post.mean[i]),c(0-postlens[i],0),col="blue")
				xcord <- (4*post.mean[i]/(y[i]-post.mean[i]) - 4*post.mean/(y-post.mean))/(4/(y[i]-post.mean[i]) - 4/(y-post.mean))
				ycord <- 4/(y-post.mean)*xcord - 4/(y-post.mean)*post.mean
				coords <- subset(cbind(xcord,ycord),ycord>0 & ycord<4)
				points(coords,col="red")
			}
			xl<-c("Indices (Groups) by the order of data input")
			if(gp.ord==T){
				temp<-temp[order(temp$n),]
				xl<-c("Indices (Groups) ordered by increasing n")
			}
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
	}else{
		if(is.na(object$prior.mean)){
			# not available now
		}else{
			y<-object$sample.mean
			n<-object$n
			prior.mean<-object$prior.mean
			post.mean<-object$lambda.hat
			post.sd<-object$se.lambda.hat
			intv.low<-object$lambda.hat.low
			intv.upp<-object$lambda.hat.upp
			index<-1:length(z)
			temp<-data.frame(n,y,p.hat=post.mean,p.hat.low=intv.low,p.hat.upp=intv.upp)
			cx<-(mean(log(n+2))+min(log(n+2)))/2
			fs <- 1.3
			sqrtV <- sqrt(1/n)
			sdlens <- sqrtV/max(sqrtV)
			postlens <- post.sd/max(sqrtV)
			sunflowerplot(rep(4,length(y))~y,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n", col.lab="white",cex.axis=fs,main="Shrinkage Plot")
			if(length(unique(prior.mean))==1)
				points(prior.mean[1],0,col="darkviolet",pch=2,cex=4)
			legend(leg.pos,c("SD"),col="blue",lty=1)
			sunflowerplot(rep(0,length(y))~post.mean,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n",col.lab="white",add=T)
			abline(h=4)
			abline(h=0)
			axis(2,c(0,4),c("p","y"), cex.axis=fs*1.1)
			for(i in 1:length(y)){
				lines(c(y[i],post.mean[i]),c(4,0),xlim=c(min(y)-abs(min(y))*0,max(y)+abs(max(y))*0))
				lines(c(y[i],y[i]+sdlens[i]*sd(y)*0.4),c(4,4+sdlens[i]),col="blue")
				##posterior variance lines
				lines(c(post.mean[i]-postlens[i]*sd(y)*0.4,post.mean[i]),c(0-postlens[i],0),col="blue")
				xcord <- (4*post.mean[i]/(y[i]-post.mean[i]) - 4*post.mean/(y-post.mean))/(4/(y[i]-post.mean[i]) - 4/(y-post.mean))
				ycord <- 4/(y-post.mean)*xcord - 4/(y-post.mean)*post.mean
				coords <- subset(cbind(xcord,ycord),ycord>0 & ycord<4)
				points(coords,col="red")
			}
			xl<-c("Indices (Groups) by the order of data input")
			if(gp.ord==T){
				temp<-temp[order(temp$n),]
				xl<-c("Indices (Groups) ordered by increasing n")
			}
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
}
