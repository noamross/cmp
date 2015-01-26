context('Numerical Correctness')

test_that("Distribution Functions Match", {
  testvals = expand.grid(lambda=seq(0.1, 1.9, by=0.2), nu=seq(0.1, 1.9, by=0.2))
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_equal(
        com_compute_log_z(lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        compoissonOld::com.compute.log.z(lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])))))
  }
  
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

test_that("pcom is bounded", {
  testvals = expand.grid(lambda=seq(0.1, 1.9, by=0.2), nu=seq(0.1, 1.9, by=0.2))
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_less_than(pcom(1000, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 1.000001))) }})

