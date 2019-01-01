rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

setwd('U:\\GIT_models\\github_SBM_gamma')
sourceCpp('rcpp_func.cpp')
source('gibbs functions.R')
source('SBM_main.R')
dat=data.matrix(read.csv('fake data.csv',as.is=T))

ngibbs=2000
ngroup.loc=10
ngroup.spp=10
res=SBM(dat=dat,ngroup.loc=ngroup.loc,ngroup.spp=ngroup.spp,ngibbs=ngibbs,burnin=ngibbs/2)
