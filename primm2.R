# All Broken Bats
ALL<-c(147,134,134,122,118,117,115,114,111,111,109,100,100,97,95,94,88,87,86,84,84,83,81,79,78,77,74,69,66,62)

# All Broken Bats Per Game
ALLPG<-c(1.62,1.46,1.46,1.33,1.3,1.3,1.25,1.24,1.22,1.25,1.18,1.09,1.09,1.07,1.03,1.03,0.98,0.97,0.93,0.95,0.91,0.9,0.88,0.9,0.85,0.84,0.8,0.78,0.73,0.67)

# Number of Games played and Exposure
ngp<-ALL/ALLPG
exps<-ngp*40

# MPF broken bats and Sample Mean
MPF<-c(26,21,33,21,27,21,24,14,29,28,19,24,21,17,29,13,18,19,22,20,22,24,25,6,14,17,12,12,4,17)
p_hat<-MPF/exps

# One covariate: Division, 1 for American league, 0 for National league
dvs<-c(0,0,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1,0,1,1,1,1,0,1,1,0,0,1,0,1)

# For convenience
z<-MPF
n<-exps
k<-length(z)
x<-as.matrix(cbind(rep(1, k), dvs))
p<-ncol(x)
z.mean <- z/n
beta0<-as.vector(glm(z.mean~dvs,family=poisson)$coefficients)
mu0<-mean(exp(x%*%beta))
r0<-mu0/var(z.mean)

m<-max(z)


trick<-function(r,b,i){
	sum(log(r*exp(x[which(z>=i),]%*%b)+i-1))		
}

trick1<-function(r,b,i){
	sum((r*exp(x[i,]%*%b)+z[i])*log(1+n[i]/r))		
}

loglik<-function(r,b){
	trick(r,b,1)+trick(r,b,2)+trick(r,b,3)+trick(r,b,4)+trick(r,b,5)+trick(r,b,6)+trick(r,b,7)+trick(r,b,8)+trick(r,b,9)+trick(r,b,10)+trick(r,b,11)+trick(r,b,12)+trick(r,b,13)+trick(r,b,14)+trick(r,b,15)+trick(r,b,16)+trick(r,b,17)+trick(r,b,18)+trick(r,b,19)+trick(r,b,20)+trick(r,b,21)+trick(r,b,22)+trick(r,b,23)+trick(r,b,24)+trick(r,b,25)+trick(r,b,26)+trick(r,b,27)+trick(r,b,28)+trick(r,b,29)+trick(r,b,30)+trick(r,b,31)+trick(r,b,32)+trick(r,b,33)-log(r)*(sum(z))-(trick1(r,b,1)+trick1(r,b,2)+trick1(r,b,3)+trick1(r,b,4)+trick1(r,b,5)+trick1(r,b,6)+trick1(r,b,7)+trick1(r,b,8)+trick1(r,b,9)+trick1(r,b,10)+trick1(r,b,11)+trick1(r,b,12)+trick1(r,b,13)+trick1(r,b,14)+trick1(r,b,15)+trick1(r,b,16)+trick1(r,b,17)+trick1(r,b,18)+trick1(r,b,19)+trick1(r,b,20)+trick1(r,b,21)+trick1(r,b,22)+trick1(r,b,23)+trick1(r,b,24)+trick1(r,b,25)+trick1(r,b,26)+trick1(r,b,27)+trick1(r,b,28)+trick1(r,b,29)+trick1(r,b,30))
}	
	
# Start making function!
primm2<-function(z,n,x,n0=NA,r0=NA,intercept=T){
	
	k <- length(z)
	r <- ncol(x)
	x <- as.matrix(cbind(rep(1, k), x))
	z.mean <- z/n
	pre.poi.reg<-glm(z.mean~x,family=poisson)
	pre.poi.reg.coeff<-as.vector(glm(z.mean~x,family=poisson)$coefficients)
	loglik<-function(r,beta=pre.poi.reg.coeff){
		m<-max(z)
		for (i in 1:m){
			sum(log(r*exp(x[which(z>=i),]%*%beta)+i-1))
		}
	}	
}

primm2<-function(z,n,x,n0=NA,r0=NA,beta0=NA,intercept=T){
