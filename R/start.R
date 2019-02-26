rm(list=ls())
library(MASS)
library(mvtnorm)
library(coda)
library(rjags)

source('fun_coda_dic.R')
obs <- read.table("../data/epidural.txt",header=T)
attach(obs)

Ntol <- n000+n001+n010+n011+n100+n101+n110+n111
is.vector(Ntol)
N0 <- n000+n001+n010+n011
N1 <- n100+n101+n110+n111
is.vector(N0)
R <- cbind(n000,n001,n010,n011, n100,n101,n110,n111)
n <- length(Ntol)
r <- (n100+n101+n110+n111)/Ntol
data <- list(N0=N0, N1=N1, R=R, n=n)