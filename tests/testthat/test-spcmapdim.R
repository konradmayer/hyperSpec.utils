library(hyperSpec)



reference_dims <- c(
  x = length(unique(chondro$x)), y = length(unique(chondro$y)),
  nwl = length(chondro@wavelength)
)

test_that("reference_dims returns correct dimensions", {
  expect_equal(spcmap_dim(chondro), reference_dims)
})

test_that("reference_dims returns an error when no spatial dimensions are present", {
  expect_error(spcmap_dim(paracetamol))
})
