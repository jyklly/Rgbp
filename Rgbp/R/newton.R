data(schools)
attach(schools)

tr <- function(M){
  sum(diag(M))
}

V = se^2
X = matrix(rnorm(length(y)*2),ncol=2)

derval = function(alpha){
  
  A = exp(alpha)
  wv = 1/(V+A)
  Wm = diag(wv)
  Sigmam = solve(t(X)%*%Wm%*%X)
  Betahat <- solve(t(X)%*%Wm%*%X)%*%t(X)%*%Wm%*%y

  l2 = log(A) - 1/2*sum(log(V+A)) + 1/2*log(det(Sigmam)) - 1/2*sum(wv*(y-X%*%Betahat))

  dlralphaBEVAL = 1 - A/2*sum(wv) + A/2*tr(Sigmam%*%t(X)%*%Wm^2%*%X) + A/2*sum(wv^2*(y-X%*%Betahat)^2)
  dlralphaBEVAL2 = -A/2*sum(wv^2*V) + A/2*tr(Sigmam%*%t(X)%*%Wm^2%*%X) + A/2*sum(wv^2*(y-X%*%Betahat)^2) +
    A^2/2*tr((Sigmam%*%t(X)%*%Wm^2%*%X)^2) - A^2*tr(Sigmam%*%t(X)%*%Wm^3%*%X) - A^2*sum(wv^3*(y-X%*%Betahat)^2)

  dbetahatA = Sigmam%*%t(X)%*%Wm^2%*%(X%*%Betahat - y)
  dbetahatA2 = 2*Sigmam%*%t(X)%*%Wm^2%*%X%*%dbetahatA - 2*Sigmam%*%t(X)%*%Wm^3%*%(X%*%Betahat-y)

  dl2alpha = dlralphaBEVAL + A*sum(wv*(y-X%*%Betahat)*X%*%dbetahatA)
  dl2alpha2 = dlralphaBEVAL2 + A*sum(wv^2*V*(y-X%*%Betahat)*X%*%dbetahatA) - A^2*sum(wv*(X%*%dbetahatA)^2) +
    A^2*sum(wv*(y-X%*%Betahat)*X%*%dbetahatA2) - A^2*sum(wv^2*(y-X%*%Betahat)*X%*%dbetahatA)
  return(list(alpha=alpha,l2=l2,dl2alpha=dl2alpha,dl2alpha2=dl2alpha2))
}

NR = function(vals,tol = 0.000001,maxits = 10){

  change = 1
  its = 0
  while(change > tol & its < maxits){
    new = vals$alpha - vals$l2/vals$dl2alpha
    print(new)
    change = abs((new - vals$alpha)/vals$alpha)
    vals = derval(new)
    its = its + 1
  }

  return(list(est = vals$alpha, der2 = vals$dl2alpha2))
}

test = NR(derval(-2))

