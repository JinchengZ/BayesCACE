rm(list = ls())
## R package BayeCACE installation 
library(devtools)
Sys.setenv("TAR" = "internal")
devtools::install_github("JinchengZ/BayesCACE")

Sys.getenv("PATH")
devtools::install_local("BayesCACE_1.0.tar.gz", dependencies = NA, upgrade = "never")
# remove.packages("BayesCACE")
# BayesCACE:::plot.cacebayes

library("BayesCACE")
data("epidural_c", package = "BayesCACE")
epidural_c
data("epidural_ic", package = "BayesCACE")
head(epidural_ic)

plot.noncomp(data = epidural_c, overall = TRUE)

set.seed(123)
out.study <- cace.study(data = epidural_c, conv.diag = TRUE, mcmc.samples = 
                          +   TRUE, two.step = TRUE)
out.study$CACE
out.study$conv.out[[1]]
out.study$meta

set.seed(123)
out.meta.c <- cace.meta.c(data = epidural_c, conv.diag = TRUE, mcmc.samples= TRUE, study.specific = TRUE)
out.meta.c$smry
out.meta.c$DIC

set.seed(123)
out.meta.ic <- cace.meta.ic(data = epidural_ic, conv.diag = TRUE, 
                            mcmc.samples = TRUE, study.specific = TRUE)

plot.cacebayes(obj = out.meta.ic)

plot.forest(data = epidural_ic, obj = out.meta.ic)

plot.forest(data = epidural_c, obj = out.study, obj2 = out.meta.c)

