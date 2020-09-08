library(hyperSpec)

test_that("simple integration works", {
  expect_equal(integrate_trapezoid(c(0, 1), c(0, 1)), 0.5)
  expect_equal(integrate_rectangular(c(0, 1), c(0, 1)), 0.5)
})


test_that("integration of a single wavelength returns 0", {
  expect_true(all(integrate_wl(chondro, test = 1600)[[]][, "test"] == 0))
})
