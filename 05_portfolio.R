rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')
load('temp/results_FX_pooling.Rdata')
load('temp/marginals.Rdata')
load('temp/RV_forc.Rdata')

# marginals$rvs=tail(RVs,K)
marginals$rvs=RV_forc$'1sa'

load('temp/results_vectordcc_tcop.Rdata')
reslf = res
load('temp/results_caw_iw.Rdata')
reshf = res
load('temp/results_heavy_t.Rdata')
resH  = res

rm(res)

# How many returns are simulated for numerical portfolio optimization 

MonteCarloSize=10000

# 1. jore's1
# 2. geweke's
# 3. equally w.
# 4. caw
# 5. dcct
# 6. dcc-heavy-t

models  = c('Jore1','Geweke','Equal','CAW','DCC-t','DCC-HEAVY-t')
ws_gmvr = ws_CVAR05r = ws_CVAR10r = array(NA,c(length(models),K,dm))

# for Low Frequency (rank-1 DCC-t) model
Q     = R = array(NA,c(dm, dm, nn+K))
a     = apply(reslf$r[ind,2:(dm+1)],2,median)
b     = apply(reslf$r[ind,(dm+2):(2*dm+1)],2,median)
nu    = median(reslf$r[ind,1])
iota  = rep(1,dm)
Oiota = Outer(iota,iota)
tdata = qt(udata,nu)
Sbar  = cova(tdata[1:nn,])
A     = Outer(a,a)
B     = Outer(b,b)
B0lf  = (Oiota-A-B)*Sbar
INL   = dt(tdata,df=nu,log=TRUE)
Q[,,1] = cora(tdata[1:nn,])

# for DCC-HEAVY-t model
Rh     = array(NA,c(dm, dm, nn+K))
nuh    = median(resH$r[ind,1])
ah     = median(resH$r[ind,2])
bh     = median(resH$r[ind,3])
Pbar   = Reduce('+',Sigma[1:nn])/nn
tdatah = qt(udata,nuh)
Rbar   = cora(tdatah[1:nn,])
Rh[,,1] = Rbar

# for High Frequency (CAW) model
nux   = mean(reshf$r[ind,1])
b1    = apply(reshf$r[ind,2:(dm+1)],2,median)
b2    = apply(reshf$r[ind,(dm+2):(2*dm+1)],2,median)
Rcaw  = array(NA,c(dm, dm, nn+K))
Rcaw[,,1] = Pbar
B1    = Outer(b1,b1)
B2    = Outer(b2,b2)
B0    = (Oiota-B1-B2)*Pbar
mtdf  = nux-dm

# some other initials
US    = matrix(runif(3*(nn+K)),ncol=nn+K,nrow=3)
A.tmp = B.tmp = C.tmp = matrix(NA,nrow=MonteCarloSize,ncol=dm)
cov1  = cov2 = cov3 = cov4 = cov5 = cov6 = array(NA,c(dm,dm,K))

class(A.tmp) <- "numeric"
class(B.tmp) <- "numeric"
class(C.tmp) <- "numeric"

# median weights
wj1 = apply(ws_jore1,2,median)
wgw = apply(ws_gew,2,median)

for(t in 2:(nn+K)){
  
  t0 = Sys.time()

  # for DCC-t model
  Q[,,t] = B0lf+A*Outer(tdata[t-1,],tdata[t-1,])+B*Q[,,(t-1)]
  t.ma   = Q[,,t]
  t.dv   = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R    = Outer(t.dv,t.dv)*t.ma
  R[,,t] = t.R
  
  # for DCC-HEAVY-t model
  Rh[,,t] = (1-bh)*Rbar-ah*Pbar+ah*Sigma[[t-1]]+bh*Rh[,,t-1]
  
  # for CAW model
  Rcaw[,,t] = B0+B1*Sigma[[t-1]]+B2*Rcaw[,,t-1]
  
  if(t>nn){

    mvnfast::rmvt(MonteCarloSize,rep(0,dm),(mtdf-1)/(mtdf+1)*Rcaw[,,t],df=mtdf+1,A=A.tmp)
    sample_uxm = pt(A.tmp,df = mtdf+1)

    mvnfast::rmvt(MonteCarloSize,rep(0,dm),R[,,t],df=nu,A=B.tmp)
    sample_udcct = pt(B.tmp,df = nu)

    mvnfast::rmvt(MonteCarloSize,rep(0,dm),Rh[,,t],df=nuh,A=C.tmp)
    sample_udccth = pt(C.tmp,df = nuh)

    sample_standretxm = qnorm(sample_uxm)
    sample_standretdcct = qnorm(sample_udcct)
    sample_standretdccth = qnorm(sample_udccth)

    sample_retsxm = ((t(sample_standretxm)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]
    sample_retsdcct = ((t(sample_standretdcct)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]
    sample_retsdccth = ((t(sample_standretdccth)*marginals$sd)+marginals$mean)*marginals$rvs[t-nn,]

    ##------------
    ## Portfolio opt.
    ##------------
    
    # CAW
    tmp  = t(sample_retsxm)
    SS   = cova(tmp)
    cov4[,,t-nn]=SS
    
    res <- minvar(SS, wmin = 0, wmax = 1)
    ws_gmvr[4,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10r[4,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05r[4,t-nn,]=c(res)

    # DCC-t
    tmp  = t(sample_retsdcct)
    SS   = cova(tmp)
    cov5[,,t-nn]=SS

    res <- minvar(SS, wmin = 0, wmax = 1)
    ws_gmvr[5,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10r[5,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05r[5,t-nn,]=c(res)

    # DCC-HEAVY-t
    tmp  = t(sample_retsdccth)
    SS   = cova(tmp)
    cov6[,,t-nn]=SS

    res <- minvar(SS, wmin = 0, wmax = 1)
    ws_gmvr[6,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10r[6,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05r[6,t-nn,]=c(res)

    # Jore's
    if(wj1[t-nn]>US[1,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

    tmp  = t(selected)
    SS   = cova(tmp)
    cov1[,,t-nn]=SS
    
    res <- minvar(SS, wmin = 0, wmax = 1)
    ws_gmvr[1,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10r[1,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05r[1,t-nn,]=c(res)
    
    # Geweke's
    if(wgw[t-nn]>US[2,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

    tmp  = t(selected)
    SS   = cova(tmp)
    cov2[,,t-nn]=SS

    res <- minvar(SS, wmin = 0, wmax = 1)
    ws_gmvr[2,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10r[2,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05r[2,t-nn,]=c(res)
 
    # Equally weighted
    if(0.5>US[3,t-nn]) selected = sample_retsxm else selected = sample_retsdcct

    tmp  = t(selected)
    SS   = cova(tmp)
    cov3[,,t-nn]=SS
    
    res <- minvar(SS, wmin = 0, wmax = 1)
    ws_gmvr[3,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.10, wmin = 0, wmax = 1)
    ws_CVAR10r[3,t-nn,]=c(res)
    
    res <- minCVaR(tmp, 0.05, wmin = 0, wmax = 1)
    ws_CVAR05r[3,t-nn,]=c(res)
    
    if(!t%%100){
      print(paste('Portfolio optimization: ', round(t/(nn+K)*100,2),"%",sep=""))
      print(Sys.time()-t0)
    }
  }
}
print("Done!")

rm(reshf,reslf,resH,Sigma,ws_DN,ws_gew,ws_jore1,Q,RCor,RCov,A.tmp,B.tmp,C.tmp,R,Rh,Rcaw,
   sample_retsdcct,sample_retsdccth,sample_retsxm,sample_standretdcct,sample_standretdccth,
   sample_standretxm,sample_udcct,sample_udccth,sample_uxm,selected,tmp)

save.image(file = 'temp/FX_portfolio_at_median.Rdata')

