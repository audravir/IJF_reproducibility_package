rm(list=ls(all=TRUE))

load("temp/individual_llh.RData")
library(truncnorm)
library(tidyquant)

##------
## Density combination
##------

ws_gew    = lL_gew = matrix(NA,ncol=K,nrow=post.sample)
ws_jore1  = lL_jore1 = matrix(NA,ncol=K,nrow=post.sample)
ws_jore10 = lL_jore10 = matrix(NA,ncol=K,nrow=post.sample)
lL_equal  = lL_DN = matrix(NA,ncol=K,nrow=post.sample)

##----------------------------
## Calculate DelNegro weights
##----------------------------

# load a Particle MCMC function for DelNegro weights
source('individual_models_and_functions/FUN_pmcmc_delnegro.R')

# Function arguments:
# first_llik,second_llik,Npart,prior,M,propsd

resDN = PMCMC_delNegro(lL_caw,lL_tdcc,1000,c(0,2),10000,propsd = 0.5)
ws_DN = t(pnorm(resDN$weights_xs))

# match the dimension of DelNegro weights with other outputs
if (dim(ws_DN)[1] > post.sample){
  dnind = round(seq(1,dim(ws_DN)[1],length=post.sample))
  ws_DN = ws_DN[dnind,]
}

# Gather DelNegro weight results
resDNxs            = list()
resDNxs$beta       = resDN$beta[dnind]
resDNxs$weights_xs = resDN$weights_xs[,dnind]
resDNxs$likelihood = resDN$likelihood[,dnind]
resDNxs$acc        = resDN$acc[dnind]

##----------------------------
## Calculate the rest of the weights
##----------------------------

lhf = exp(lL_caw)
llf = exp(lL_tdcc)

for(m in 1:post.sample){
  for (t in 1:K){
    # weights geweke
    a1p = lhf[m,1:t]
    a2p = llf[m,1:t]
    opw = function(x){
      -sum(log(x*a1p+(1-x)*a2p))
    }
    ws_gew[m,t] = optim(0.5, opw, gr = NULL,
                        method = c("L-BFGS-B"),
                        lower = 0, upper = 1, hessian = FALSE)$par
    lL_gew[m,t] = log(ws_gew[m,t]*lhf[m,t]+(1-ws_gew[m,t])*llf[m,t])
    
    # DN
    
    lL_DN[m,t] = log(ws_DN[m,t]*lhf[m,t]+(1-ws_DN[m,t])*llf[m,t])
    
    # jore 1 past
    ws_jore1[m,t] = a1p[t]/(a1p[t]+a2p[t])
    lL_jore1[m,t] = log(ws_jore1[m,t]*lhf[m,t]+(1-ws_jore1[m,t])*llf[m,t])
    
    # jore 10 past
    pst = 10
    ws_jore10[m,t] = exp(sum(log(a1p[max(t-pst+1,1):t])))/
      (exp(sum(log(a1p[max(t-pst+1,1):t])))+exp(sum(log(a2p[max(t-pst+1,1):t]))))
    lL_jore10[m,t] = log(ws_jore10[m,t]*lhf[m,t]+(1-ws_jore10[m,t])*llf[m,t])
    
    # equally-weighted
    lL_equal[m,t] = log(0.5*lhf[m,t]+0.5*llf[m,t])
  }
}


# Calculate market volatility
mkvol = apply(tail(apply(RVs^2,2,scale),K),1,mean)

# -----------------------------
# LPTS
# -----------------------------

# Find the quantiles Q=50,25,10,5

cuts = seq(0.6,1,length=5000)
perc = rep(NA,length(cuts))
for(i in 1:length(cuts)){
  perc[i] = mean(apply(tail(udata,K)<cuts[i],1,all))
}

ind50 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.5))],1,all))
ind25 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.25))],1,all))
ind10 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.10))],1,all))
ind05 = which(apply(tail(udata,K)<cuts[which.min(abs(perc - 0.05))],1,all))

# -----------------------------
# Correlation between weights and avrg volatility
# -----------------------------

datayahoo = getSymbols('^VIX', from = date[1],
                       to = tail(date,1)+1,warnings = FALSE,
                       auto.assign = FALSE,periodicity = "daily",
                       return.class = 'xts')

tmp      = datayahoo$VIX.Adjusted
ind.date = match(date,index(tmp))
vix      = tail(tmp[ind.date],K)

ws = rep(1/dm,dm)

avrg2 =  rep(NA,K)
for(t in 1:K){
  avrg2[t] = t(ws)%*%RCov[,,nn+t]%*%(ws)
}

marketvols = cbind(mkvol,avrg2,coredata(vix))

##---------------
## at mcmc
##---------------

corr.mat = array(NA,dim=c(post.sample,7,7))

for(m in 1:post.sample){
  tmp.matrix = cbind(marketvols,ws_gew[m,],ws_jore1[m,],ws_DN[m,],lL_caw[m,]-lL_tdcc[m,])
  corr.mat[m,,] = cor(tmp.matrix,use = "pairwise.complete.obs")
}

corrs_res = apply(corr.mat,c(2,3),median)
colnames(corrs_res) = c('avrg RV','Mkt:eql','VIX','Geweke','Jore1','DelNegro','logLik diff')
rownames(corrs_res) = c('avrg RV','Mkt:eql','VIX','Geweke','Jore1','DelNegro','logLik diff')
corrs_res[upper.tri(corrs_res)] = NA

save.image(file="temp/results_FX_pooling.Rdata")


