
gbp.shrinkage.graph.R <- function(y,se,mu,thetahat,shat){
  fs <- 1.3
  sqrtV <- se
  sdlens <- sqrtV/max(sqrtV)
  postlens <- shat/max(sqrtV)
  sunflowerplot(rep(4,length(y))~y,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n", col.lab="white",cex.axis=fs)
  if(length(unique(mu))==1)
    points(mu[1],0,col="darkviolet",pch=2,cex=4)
  legend("topright",c("SD"),col="blue",lty=1)
  par(new=TRUE)
  sunflowerplot(rep(0,length(y))~thetahat,ylim=c(-1,5),xlim=c(min(y)-abs(min(y))*0.1,max(y)+abs(max(y))*0.1),yaxt="n",cex.axis=fs,col.lab="white")
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
}
