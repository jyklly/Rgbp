gbp<-function(arg1,arg2,mu=NA,X=NA,model="gr",CI=0.95,intercept=T){

	res<-switch(model, 
			gr=gr(arg1,arg2,X,mu,CI), 
			br=bp(arg1,arg2,prior.mean=mu,x=X,model="br",CI,intercept=intercept), 
			pr=bp(arg1,arg2,prior.mean=mu,x=X,model="pr",CI,intercept=intercept) )

	class(res)<-"gbp"	
	res
}

g<-gbp(y,se,mu=8.75,model="gr")
b<-gbp(z,n,mu=0.265,model="br")
p<-gbp(z,n,mu=0.265,model="pr")
