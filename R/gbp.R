gbp <- function(x, ...) UseMethod("gbp")

gbp.default<-function(arg1,arg2,mu,X,model="gr",CI=0.95,intercept=T,eps=0.0001){

  res<-switch(model, 
              gr=gr(arg1,arg2,X,mu,CI), 
              br=bp(arg1,arg2,X,prior.mean=mu,model="br",CI=CI,intercept=intercept,eps=eps), 
              pr=bp(arg1,arg2,X,prior.mean=mu,model="pr",CI=CI,intercept=intercept,eps=eps) )
  
  class(res)<-"gbp"	
  res

}
