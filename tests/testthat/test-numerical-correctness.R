context('Numerical Correctness')

test_that("Distribution Functions Match", {
  testvals = expand.grid(lambda=c(0.1, 1), nu=10^(-3:3))
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_equal(
        dcom(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        compoissonOld::dcom(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])))))
  }
  
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_equal(
        pcom(10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        sum(compoissonOld::dcom(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i]))))))
  }
  
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_equal(
        com.log.density(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        compoissonOld::com.log.density(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])))))
  }
  
})

