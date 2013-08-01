nits <- 1000
results <- lapply(1:nits, function(s){load(file=paste("output/mcmcout",s,".RData", sep = ""));return(out)})
ind <- apply(do.call("rbind",lapply(results,function(x){x$coverage.ind})),2,mean)
rb <- apply(do.call("rbind",lapply(results,function(x){x$coverage.RB})),2,mean)
