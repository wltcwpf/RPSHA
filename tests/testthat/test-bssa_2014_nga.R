test_that("test bssa_2014_nga vs NGA-West2 spreadsheet", {
  expect_equal(bssa_2014_nga(M = 5, period = -1, Rjb = 85, Fault_Type = 3, region = 1,
                             z1 = -999, Vs30 = 350, coeffs = as.matrix(bssa_2014_coeffs))$med, 0.269344818)

  expect_equal(round(bssa_2014_nga(M = 5, period = 0, Rjb = 85, Fault_Type = 3, region = 1,
                             z1 = -999, Vs30 = 350, coeffs = as.matrix(bssa_2014_coeffs))$sigm, digits = 4), 0.7022)

  expect_equal(round(bssa_2014_nga(M = 6.89, period = 0, Rjb = 251, Fault_Type = 1, region = 1,
                             z1 = 1.5, Vs30 = 487, coeffs = as.matrix(bssa_2014_coeffs))$med, digits = 6), round(0.005315, digits = 6))

})
