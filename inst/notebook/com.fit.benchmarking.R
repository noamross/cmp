delibrary("microbenchmark")
library("compoisson")
library("dplyr")

set.seed(0)

tdata = as.matrix(table(rcom(10000, 1.6, 2))) %>% cbind(as.integer(rownames(.)), . )
compoisson::com.fit(tdata)
com.fit(tdata)


Rprof("old.prof")
replicate(10, a <- compoisson::com.fit(tdata))
Rprof("new.prof")
replicate(10000, a <- com.fit(tdata))
Rprof(NULL)

microbenchmark(com.fit(tdata), compoisson::com.fit(tdata), times=10)

library(noamtools)
proftable("old.prof")
proftable("new.prof")

ttdata = rcom(100, 1.6, 0.3)
com_loglik(tdata, 1.6, 0.3)
pois_loglik(tdata, 1.6)
nb_loglik(tdata, 1.6, 1)
pois_fit(tdata)
nb_fit(tdata)
com_fit(tdata)
