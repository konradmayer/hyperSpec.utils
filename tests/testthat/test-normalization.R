library(hyperSpec)

test_that("the range after minmax_normalization is 0:1", {
  expect_true(
    all(
      apply(
        hyperSpec::apply(
          minmax_normalization(chondro), 1, range
        )$spc, 1,
        function(x) all.equal(x, c(0, 1), check.attributes = FALSE)
      )
    )
  )
})

test_that("the sum of intensities after area normalization is 1", {
  expect_true(
    all(
      apply(
        hyperSpec::apply(area_normalization(chondro, "sum"), 1, sum)$spc, 1,
        function(x) all.equal(x, 1, check.attributes = FALSE)
      )
    )
  )
})


test_that("the mean after SNV normalization is 0", {
  expect_true(
    all(
      apply(
        hyperSpec::apply(snv_normalization(chondro), 1, mean)$spc, 1,
        function(x) all.equal(x, 0, check.attributes = FALSE)
      )
    )
  )
})


test_that("the sd after SNV normalization is 1", {
  expect_true(
    all(
      apply(
        hyperSpec::apply(snv_normalization(chondro), 1, sd)$spc, 1,
        function(x) all.equal(x, 1, check.attributes = FALSE)
      )
    )
  )
})


test_that("the intensity for the single wavelength used for band normalization equals 1", {
  expect_true(
    all(
      apply(
        band_normalization(chondro, 856)[, , 856]$spc, 1,
        function(x) all.equal(x, 1, check.attributes = FALSE)
      )
    )
  )
})


test_that("the mean intensity for the wavelength range used for band normalization equals 1", {
  expect_true(
    all(
      apply(
        hyperSpec::apply(band_normalization(chondro, 800 ~ 900)[, , 800 ~ 900], 1, mean)$spc, 1,
        function(x) all.equal(x, 1, check.attributes = FALSE)
      )
    )
  )
})
