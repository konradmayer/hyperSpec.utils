library(hyperSpec)
testarray <- spcmap2array(chondro)


test_that("the returned object is of class array", {
  expect_is(testarray, "array")
})

test_that("the array has the correct attributes attached", {
  expect_named(attributes(testarray), c("dim", "wavelength"))
})

test_that("dimensions of the array equal to those of the hyperSpec object", {
  expect_equal(dim(testarray), c(x = 35, y = 25, nwl = 300))
})
