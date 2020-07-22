library(tidyverse)
library(hyperSpec)



# normalization -----------------------------------------------------------

carb <- read.jdx("/media/konrad/raman2/6 Team/14 Konrad/datenaustausch/mixture_carb_r/carb.dx")

plot(minmax_normalization(carb))
apply(minmax_normalization(carb), 1, range)$spc # minmax 0 to 1

apply(area_normalization(carb, "sum"), 1, sum)$spc # should be 1

apply(snv_normalization(carb), 1, mean)$spc # should be 0
apply(snv_normalization(carb), 1, sd)$spc # should be 1

plot(band_normalization(carb, 856))
abline(v = 856)
abline(h = 1)

band_normalization(carb, 856)[, , 856]$spc # band 856 should be 1

plot(band_normalization(carb, 2900 ~ 2970))
abline(v = c(2900, 2970))
apply(band_normalization(carb, 2900 ~ 2970)[, , 2900 ~ 2970], 1, mean)$spc # mean over wavelengths used for standardization should be 1


plot(vector_normalization(carb))
apply(vector_normalization(carb), 1, function(x) sum(x^2))$spc # sum of squares should be 1


# is_hyperSpec ------------------------------------------------------------

spcmap <- read.txt.Witec.Graph("/media/konrad/Elements/arbeit/image-registration/test/graph-asci/G120_oct_10mV 0,01s_001_Spec.Data 1 (Header).txt",
  type = "map", encoding = "latin1"
)

is.hyperSpec(spcmap) # should be TRUE
is.hyperSpecMap(spcmap) # should be TRUE

is.hyperSpec(carb) # should be TRUE
is.hyperSpecMap(carb) # should be FALSE


# spcmap2array ------------------------------------------------------------

testarray <- spcmap2array(spcmap)
testarray %>% dim() # x, y should multiply to 8100 pixels; 1024 wavelengths
testarray %>% class() # should be array
attributes(testarray) %>% str() # should contain dim and wavelengths
