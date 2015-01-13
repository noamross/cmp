devtools::load_all()
data("insurance")
Rprof("comfit.prof")
x = replicate(500, com.fit(Lemaire))
Rprof(NULL)
noamtools::proftable("comfit.prof", 30)

#---

library(microbenchmark)
logcomp = microbenchmark(com.log.factorial(1:100), lfactorial(1:100))
library(ggplot2)
autoplot(logcomp)

#---

# TODO: Replace com.log.factorial with lfactorial