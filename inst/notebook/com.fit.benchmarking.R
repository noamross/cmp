library("microbenchmark")
library("compoisson")
library("dplyr")

set.seed(0)

tdata = as.matrix(table(rcom(10000, 1.6, 0.3))) %>% cbind(as.integer(rownames(.)), . )
compoissonOld::com.fit(tdata)
com.fit(tdata)


Rprof("old.prof")
replicate(10, a <- compoissonOld::com.fit(tdata))
Rprof("new.prof")
replicate(1000, a <- com.fit(tdata))
Rprof(NULL)

microbenchmark(com.fit(tdata), compoissonOld::com.fit(tdata), times=10)

library(noamtools)
proftable("old.prof")
proftable("new.prof")
