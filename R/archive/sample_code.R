library("BayesCACE")
data("epidural_c", package = "BayesCACE")
epidural_c
data("epidural_ic", package = "BayesCACE")
head(epidural_ic)


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


smry<-signif(summary(out.meta.ic$samples), digits = 2)


smry<-summary(jags.out$samples)
smry<-cbind(smry$statistics[,c("Mean","SD")],smry$quantiles[,c("2.5%","50%","97.5%")], 
            smry$statistics[,c("Naive SE","Time-series SE")])
smry<-signif(smry,digits=3)
smry<-round(smry,digits=3)
out$smry <- smry