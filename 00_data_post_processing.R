## DON'T FORGET TO SET WORKING DIRECTORY TO SOURCE FILE LOCATION

rm(list=ls(all=TRUE))

library(future.apply)
library(LaplacesDemon)
library(matrixcalc)
library(moments)
library(mvnfast)
library(mvtnorm)
library(MCS)
library(NMOF)
library(plyr)
library(Rfast)
library(Rglpk)
library(truncnorm)
library(tidyquant)
library(xtable)
library(xts)
library(zoo)



load('data/realized_measures.Rdata')

my.date = as.Date(index(rets))
assets  = colnames(rets)
dm      = ncol(rets)

# standardize and estimate static marginals

stand   = stand.full = matrix(NA,ncol=dm,nrow=dim(rets)[1])

for(t in 1:dim(rets)[1]){
  stand[t,] = rets[t,]/sqrt(diag(RCov[,,t]))
}

stand_means = apply(stand, 2, mean)
stand_sds   = apply(stand, 2, sd)
marginals   = list(stand_means,stand_sds)
names(marginals) = c('mean','sd')

save(marginals,file='temp/marginals.Rdata')

for(i in 1:dm){
  stand.full[,i] = (stand[,i]-stand_means[i])/stand_sds[i]
}

##

Sigma = alply(RCor,3)
date  = as.Date(index(rets))
rets  = coredata(rets)
stand = stand.full
udata = pnorm(stand.full)

save(Sigma,rets,RVs,stand,udata,assets,RCov,RCor,date,file='data/FXdata.Rdata')

