rm(list=ls(all=TRUE))
load('data/FXdata.Rdata')

data     = stand # Prediction etc is performed on STANDARDIZED returns
dm       = dim(data)[2] 
start    = 1

nn.all   = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2020")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}

nn       = end.date
Sig      = Sigma
K        = nn.all-end.date
post.sample = 1000

##------
## Static
#  did not estimate separately
#  estimator is here, with IWishart prior for the R matrix
#  the posterior predictive is multivariate t
##------

lL_static = matrix(NA,ncol=K,nrow=1)

for(t in (nn+1):(nn+K)){
  stdf = 10+(t-1)-dm+1
  scm  = ((10-dm-1)*diag(dm)+t(data[start:(t-1),])%*%data[start:(t-1),])/stdf 
  lL_static[(t-nn)] <- mvtnorm::dmvt(data[t,],delta=rep(0,dm),scm,df = stdf,log=TRUE)
}

sum(lL_static) #the LPS for K=797 

##------
## vector dcc Gaussian Copula
##------

load('temp/results_vectordcc.Rdata')
M   = dim(res$resdcc)[1] # size of MCMC
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q       = R = array(NA,c(dm, dm, nn+K))
Q[,,1] <- cor(data[start:nn,])
a      <- res$resdcc[ind,1:dm]
b      <- res$resdcc[ind,(dm+1):(dm*2)]
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
Sbar   = cov(data[start:nn,])

Vpred  <- list()
lL_dcc = matrix(NA,ncol=K,nrow=post.sample)

for(m in 1:post.sample){
  A      = Outer(a[m,],a[m,])
  B      = Outer(b[m,],b[m,])

  for(t in 2:(nn+K)){
    Q[,,t]   <- (Oiota-A-B)*Sbar+A*Outer(data[t-1,],data[t-1,])+B*Q[,,(t-1)]
    t.ma  = Q[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    R[,,t] = t.R
    if(t>nn){
      lL_dcc[m,(t-nn)] <- mvnfast::dmvn(data[t,], rep(0,dm), t.R, log=TRUE)
    }
  }
}

##-----------------------
## vector dcc t Copula
##-----------------------

load('temp/results_vectordcc_tcop.Rdata')

M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth

Q  = R = array(NA,c(dm, dm, nn+K))
a  <- res$r[ind,2:(dm+1)]
b  <- res$r[ind,(dm+2):(2*dm+1)]
nu <- res$r[ind,1]

lL_tdcc = matrix(NA,ncol=K,nrow=post.sample)
iota    = rep(1,dm)
Oiota   = Outer(iota,iota)
INL.N   = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  tdata <- qt(udata,nu[m])
  Sbar  <- cova(tdata[start:nn,])
  A     <- Outer(a[m,],a[m,])
  B     <- Outer(b[m,],b[m,])
  B0    <- (Oiota-A-B)*Sbar
  INL   <- dt(tdata,df=nu[m],log=TRUE)
  Q[,,1] <- cora(tdata[start:nn,])
  
  for(t in 2:(nn+K)){
    Q[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Q[,,(t-1)]
    t.ma  = Q[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    R[,,t] = t.R
    if(t>nn){
      lL_tdcc[m,(t-nn)] = mvnfast::dmvt(tdata[t,], rep(0,dm),t.R, nu[m], log=TRUE)+sum(INL.N[t,])-sum(INL[t,])
    }
  }
}

#-----------------------
# dcc-HEAVY-scalar t Copula
#-----------------------

load('temp/results_heavy_t.Rdata')
M     = dim(res$r)[1]
ind   = round(seq(1,M,length=post.sample)) #thin every xth
lL_ht = matrix(NA,ncol=K,nrow=post.sample)

nu = res$r[ind,1]
a  = res$r[ind,2]
b  = res$r[ind,3]

Rpred = res$Rpred[ind]
Pbar  = Reduce('+',Sigma[1:nn])/nn
R     = array(NA,c(dm, dm, nn+K))
INL.N = dnorm(data,log=TRUE)

for(m in 1:post.sample){
  
  tdata  = qt(udata,nu[m])
  Rbar   = cor(tdata[1:nn,])
  R[,,1] = Rbar
  INL.T  = dt(tdata,df=nu[m],log=TRUE)
  
  for(t in 2:(nn+K)){
    R[,,t] = (1-b[m])*Rbar-a[m]*Pbar+a[m]*Sigma[[t-1]]+b[m]*R[,,t-1]
    if(t>nn){
      lL_ht[m,(t-nn)]= mvnfast::dmvt(tdata[t,], rep(0,dm),R[,,t],nu[m],log=TRUE)+sum(INL.N[t,])-sum(INL.T[t,])
    }
  }
}

##------
## XM
##------

Sbar   = Reduce('+',Sig[1:nn])/nn
iota   = rep(1,dm)

load('temp/results_xm.Rdata')
M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_xm = matrix(NA,ncol=K,nrow=post.sample)

lag = res$r[ind,1]
nu  = res$r[ind,2]
b1  = res$r[ind,3:(dm+2)]
b2  = res$r[ind,(dm+3):(2*dm+2)]

for(m in 1:post.sample){
  Vpred = vector(mode = "list", length = K)
  B1  = Outer(b1[m,],b1[m,])
  B2  = Outer(b2[m,],b2[m,])
  B0  = (Oiota-B1-B2)*Sbar
  
  for(t0 in 1:K){
    Vpred[[t0]] = B0+B1*Sig[[nn+t0-1]]+B2*(Reduce('+',Sig[(nn+t0-lag[m]):(nn+t0-1)])/lag[m])
    mtdf = nu[m]-dm
    lL_xm[m,t0] = mvnfast::dmvt(data[nn+t0,], rep(0,dm), (mtdf-1)/(mtdf+1)*Vpred[[t0]], df = mtdf+1, log=TRUE)
  }
}

##------
## CAW-iw
##------

Sbar   = Reduce('+',Sig[1:nn])/nn
iota   = rep(1,dm)

load('temp/results_caw_iw.Rdata')
M   = dim(res$r)[1]
ind = round(seq(1,M,length=post.sample)) #thin every xth
lL_caw = matrix(NA,ncol=K,nrow=post.sample)

nu  = res$r[ind,1]
b1  = res$r[ind,2:(dm+1)]
b2  = res$r[ind,(dm+2):(2*dm+1)]
Vlast = res$Vpred[ind]
R = vector(mode = "list", length = K+nn)
R[[1]] = Sbar

for(m in 1:post.sample){
  B1  = Outer(b1[m,],b1[m,])
  B2  = Outer(b2[m,],b2[m,])
  B0  = (Oiota-B1-B2)*Sbar
  mtdf = nu[m]-dm
  
  for(t in 2:(K+nn)){
    R[[t]] = B0+B1*Sig[[t-1]]+B2*R[[t-1]]
    if(t>nn){
      lL_caw[m,t-nn] = mvnfast::dmvt(data[t,], rep(0,dm), (mtdf-1)/(mtdf+1)*R[[t]], df = mtdf+1, log=TRUE)
    }
  }
}

##----------------------
## save here
##----------------------
rm(res,Sig,R,Q)

save.image("temp/individual_llh.RData")



