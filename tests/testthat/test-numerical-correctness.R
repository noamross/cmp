context('Numerical Correctness')

test_that("Distribution Functions Match", {
  testvals = expand.grid(lambda=seq(0.1, 1.9, by=0.2), nu=seq(0.1, 1.9, by=0.2))
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_equal(
        compute_log_z(lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        compoisson::com.compute.log.z(lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])))))
  }
  
  for(i in 1:nrow(testvals)) {
    
    eval(bquote(
      expect_equal(
        dcmp(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        compoisson::dcom(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])))))
  }
  
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_equal(
        pcmp(10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        sum(compoisson::dcmp(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i]))))))
  }
  
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_equal(
        com.log.density(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 
        compoisson::com.log.density(0:10, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])))))
  }
  
})

test_that("pcmp is bounded", {
  testvals = expand.grid(lambda=seq(0.1, 1.9, by=0.2), nu=seq(0.1, 1.9, by=0.2))
  for(i in 1:nrow(testvals)) {
    eval(bquote(
      expect_less_than(pcmp(1000, lambda = .(testvals$lambda[i]), nu = .(testvals$nu[i])), 1.000001))) }})

