gbp<-function(arg1,arg2,mu,X,model="gr",CI=0.95,intercept=T,eps=0.0001){

  if(missing(X))
    x <- NA
  if(missing(mu))
    prior.mean <- NA
  
  print(x)
  res<-switch(model, 
              gr=gr(arg1,arg2,X,mu,CI), 
              br=bp(arg1,arg2,prior.mean,x,model="br",CI,intercept=intercept,eps=eps), 
              pr=bp(arg1,arg2,prior.mean,x,model="pr",CI,intercept=intercept,eps=eps) )
  
  class(res)<-"gbp"	
  res
}
