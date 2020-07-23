library(hyperSpec)
chondro_spike <- chondro
chondro_spike[[20, , 800]] <- 5000
chondro_spike[[50, , 1400]] <- 4000
chondro_spike[[50, , 1200]] <- 5000

despiked <- crr(chondro_spike)
is_spike <- unlist(lapply(despiked$crr, length)) > 0

test_that("cosmic ray removal finds the spectra having manually introduced spikes", {
  expect_equal(which(is_spike), c(20, 50))
})


test_that("found spikes are at the correct wavenumber", {
  expect_equal(round(unlist(despiked$crr[is_spike]), digits = -1), c(800, 1200, 1400))
})
