devtools::load_all()
data("insurance")
Rprof("comfit.prof")
x = replicate(1500, com.fit(Lemaire))
Rprof(NULL)
noamtools::proftable("comfit.prof", 30)

#---

data("insurance")
library(microbenchmark)
com_fit_compare = microbenchmark(compoisson::com.fit(Lemaire), compoisson2::com.fit(Lemaire), times=1000L)
summary(com_fit_compare)
library(ggplot2)

model = compoisson::com.fit(Lemaire)
model2 = compoisson2::com.fit(Lemaire)
compoisson::com.mean(model$lambda, model$nu)
compoisson2::com.mean(model$lambda, model$nu)
#---


# TODO: Replace com.log.factorial with lfactorial