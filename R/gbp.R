gbp<-function(arg1,arg2,arg3,arg4,model="gr",intv.alpha,intercept=T){

	switch(model, 
		gr=adm(arg1,arg2,arg3,arg4,intv.alpha,intercept=T), 
		br=bp(arg1,arg2,arg3,arg4,model="br",intv.alpha,intercept=T), 
		pr=bp(arg1,arg2,arg3,arg4,model="pr",intv.alpha,intercept=T))
}
