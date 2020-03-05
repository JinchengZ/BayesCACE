## R package BayeCACE installation 
library(devtools)
devtools::install_github("JinchengZ/BayesCACE")
setwd("C:/Users/jzhou02/Documents/Jincheng files/BayesCACE")
Sys.setenv("TAR" = "internal")
Sys.getenv("PATH")
devtools::install_local("BayesCACE_0.1.tar.gz", dependencies = NA, upgrade = "never")
# remove.packages("BayesCACE")
# BayesCACE::plot.cacebayes

library("BayesCACE")
data("epidural_c", package = "BayesCACE")
epidural_c
data("epidural_ic", package = "BayesCACE")
head(epidural_ic)

noncomp <- plot.noncomp(data = epidural_c, overall = TRUE)

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

plotsets <- plot.cacebayes(obj = out.meta.ic)

forest_ic <- plot.forest(data = epidural_ic, obj = out.meta.ic)

forest_c <- plot.forest(data = epidural_c, obj = out.study, obj2 = out.meta.c)

