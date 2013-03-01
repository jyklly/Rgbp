library(gbp)
data(schools)
attach(schools)




tr <- function(M){
  sum(diag(M))
}

V = se^2
X = matrix(rnorm(length(y)*2),ncol=2)
source("gr.R")
a <- gr(y,se,X)

derval = function(alpha,y,V,X){
  
  A = exp(alpha)
  wv = 1/(V+A)
  Wm = diag(wv)
  Sigmam = solve(t(X)%*%Wm%*%X)
  Betahat <- solve(t(X)%*%Wm%*%X)%*%t(X)%*%Wm%*%y

  l2 = log(A) - 1/2*sum(log(V+A)) + 1/2*log(det(Sigmam)) - 1/2*sum(wv*(y-X%*%Betahat)^2)

  dlralphaBEVAL = 1 - A/2*sum(wv) + A/2*tr(Sigmam%*%t(X)%*%Wm^2%*%X) + A/2*sum(wv^2*(y-X%*%Betahat)^2)
  dlralphaBEVAL2 = -A/2*sum(wv^2*V) + A/2*tr(Sigmam%*%t(X)%*%Wm^2%*%X) + A/2*sum(wv^2*(y-X%*%Betahat)^2) +
    A^2/2*tr((Sigmam%*%t(X)%*%Wm^2%*%X)%*%(Sigmam%*%t(X)%*%Wm^2%*%X)) - A^2*tr(Sigmam%*%t(X)%*%Wm^3%*%X) - A^2*sum(wv^3*(y-X%*%Betahat)^2
)
  dbetahatA = Sigmam%*%t(X)%*%Wm^2%*%(X%*%Betahat - y)
  dbetahatA2 = 2*Sigmam%*%t(X)%*%Wm^2%*%X%*%dbetahatA - 2*Sigmam%*%t(X)%*%Wm^3%*%(X%*%Betahat-y)
  dl2alpha = dlralphaBEVAL + A*sum(wv*(y-X%*%Betahat)*X%*%dbetahatA)
  dl2alpha2 = dlralphaBEVAL2 + A*sum(wv^2*V*(y-X%*%Betahat)*X%*%dbetahatA) - A^2*sum(wv*(X%*%dbetahatA)^2) +
    A^2*sum(wv*(y-X%*%Betahat)*X%*%dbetahatA2) - A^2*sum(wv^2*(y-X%*%Betahat)*X%*%dbetahatA)
  return(list(alpha=alpha,l2=l2,dl2alpha=dl2alpha,dl2alpha2=dl2alpha2))
}

NR = function(vals,tol = 0.01,maxits = 100, weight = 1){

  newsor = vals$alpha
  change = 1
  its = 0
  while(change > tol & its < maxits){
    new = vals$alpha - vals$dl2alpha/vals$dl2alpha2
    newsor = (1-weight)*newsor + new*weight
    change = abs((newsor - vals$alpha)/vals$alpha)
    vals = derval(newsor,y,V,X)
    its = its + 1
    print(newsor)
    print(as.matrix(c(vals$dl2alpha,vals$dl2alpha2)))
  }
  return(list(est = vals$alpha, der2 = vals$dl2alpha2))
}

test = NR(derval(log(5),y,V,X),maxits=5,weight=0.1,tol=0.0001)
#golden section

alphavec = seq(-5,15,length.out=1000)
res = lapply(alphavec,derval)
l2 <- as.numeric(getattr(res,"l2"))
dl2alpha <- as.numeric(getattr(res,"dl2alpha"))
par(mfrow = c(2,1))
plot(l2~alphavec)
plot(dl2alpha ~ alphavec)

