library(Rgbp)
data(schools)
Rprof("gr.out")
lapply(1:1000,function(x){gbp(schools$y, schools$se)})
Rprof(NULL)
summaryRprof("gr.out")


Rprof("qsn.out")
lapply(1:10000,function(x){qsn(0.75, 1, 2, 1,engine = "T.Owen")})
Rprof(NULL)
summaryRprof("qsn.out")

