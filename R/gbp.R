gbp<-function(arg1,arg2,mu=NA,X=NA,model="gr",CI=0.95,intercept=T,eps=0.0001){
  res<-switch(model, 
              gr=gr(arg1,arg2,X,mu,CI), 
              br=bp(arg1,arg2,prior.mean=mu,x=X,model="br",CI,intercept=intercept,eps=eps), 
              pr=bp(arg1,arg2,prior.mean=mu,x=X,model="pr",CI,intercept=intercept,eps=eps) )
  
  class(res)<-"gbp"	
  res
}
