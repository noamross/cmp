context('Backwards Compatability')

test_that("CMPControl compatible", {
  library("CMPControl")
  graphics.off()
  expect_equal_to_reference(ControlCharts(nonconformities,"Sample Number", "Number of nonconformities"), 'CMPControlChart.rds')
  })

