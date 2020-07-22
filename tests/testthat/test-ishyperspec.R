library(hyperSpec)

test_that("is_hyperSpec returns TRUE for a hyperSpec object with x and y variable", {
  expect_true(
    is.hyperSpec(chondro)
  )
})

test_that("is_hyperSpecMap returns TRUE for a hyperSpec object with x and y variables", {
  expect_true(
    is.hyperSpecMap(chondro)
  )
})

test_that("is_hyperSpec returns TRUE for a hyperSpec object", {
  expect_true(
    is.hyperSpec(paracetamol)
  )
})

test_that("is_hyperSpecMap returns TRUE for a hyperSpec object if x and y variables are missing", {
  expect_false(
    is.hyperSpecMap(paracetamol)
  )
})
