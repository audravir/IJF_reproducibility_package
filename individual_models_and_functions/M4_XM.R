plan(multisession, workers = 4)

nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2020")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}

T0  = end.date  ## for estimation
K = nn-end.date # for oos evaluation

data = Sigma[1:T0]

propsdb  = 0.001
propsdnu = 0.01

t0   = Sys.time()
t1   = Sys.time()

diwish = function(Sig,nu,S){dinvwishart(Sig, nu, S, log=TRUE)}

Sig    = data
dm     = dim(Sig[[1]])[1]
TT     = length(Sig)
bi     = min(M,50000)
TIMING = rep(NA,M+bi)
resc   = matrix(NA,nrow=M+bi,ncol=dm*2+2)
Vpred  = vector(mode = "list", length = M)
nu     = 25 
lag    = 30
b1     = rep(0.46,dm) 
b2     = rep(0.76,dm) 
Sbar   = Reduce('+',Sig)/TT
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
B1     = Outer(b1,b1) 
B2     = Outer(b2,b2)
B0     = (Oiota-B1-B2)*Sbar
llo    = lln = rep(0,TT)
LLH    = rep(NA,M+bi)
accB1  = accnu = accl = accB2 = rep(0,bi+M)
V      = Vn = list()
G1     = G2 = G2n = c(list(matrix(0,nrow=dm,ncol=dm)),Sig[-TT])
V[[1]] = Vn[[1]] = (B0+(b1%*%t(b1))*G1[[1]]+(b2%*%t(b2))*G2[[1]])

for(t in 2:TT){
  G2[[t]] = Reduce('+',Sig[max(1,t-lag):(t-1)])/min(c(t-1,lag))
  V[[t]]  = (B0+(b1%*%t(b1))*G1[[t]]+(b2%*%t(b2))*G2[[t]])
}

diwish.t = function(x,y){LaplacesDemon::dinvwishart(x, nu, (nu-dm-1)*y, log=TRUE)}
llo <- future_mapply(diwish.t,Sig,V)

for(m in 1:(bi+M)){
  t2   = Sys.time()

  ##-----
  ## l
  ##-----
  repeat{
    pos  = rbinom(1,1,0.5)
    lagnew = lag+(-1)^(pos)
    if(lagnew>1) break
  }
  
  for(t in 2:TT){
    G2n[[t]] = Reduce('+',Sig[max(1,t-lagnew):(t-1)])/min(c(t-1,lagnew))
    V[[t]]   = (B0+B1*G1[[t]]+B2*G2n[[t]])
    # lln[t]   = diwish(Sig[[t]],nu,(nu-dm-1)*V[[t]])
  }
  
  diwish.t <- function(x,y){LaplacesDemon::dinvwishart(x, nu, (nu-dm-1)*y, log=TRUE)}
  lln      <- future_mapply(diwish.t,Sig,V)
  
  if((sum(lln)-sum(llo))>log(runif(1))){
    llo     = lln
    lag     = lagnew
    G2      = G2n
    accl[m] = 1
  }

  ##-----
  ## bs split randomly
  ##-----
  block1 <- sample(c(TRUE,FALSE),size=dm*2,replace = TRUE)
  block2 <- (!block1)
  
  # block 1
  repeat{
    b.prop = rnorm(dm*2,c(b1,b2),sd=propsdb)
    bn     = b.prop*block1+c(b1,b2)*block2
    
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar
    
    cond1 = b1n[1]>0
    cond2 = b2n[1]>0
    cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1
    cond4 = sum(abs(B1+B2)<1)==dm^2
    if(cond1 && cond2 && cond4) {
      break
    }
  }
  
  if(cond3){
    for(t in 1:TT){
      Vn[[t]]  = (B0+B1*G1[[t]]+B2*G2[[t]])
    }
      lln <- future_mapply(diwish.t,Sig,Vn)
  } else {lln=-Inf}
  
  if((sum(lln)-sum(llo)+
      sum(dnorm(b1n,0,sqrt(10),log=TRUE))-sum(dnorm(b1,0,sqrt(10),log=TRUE))+
      sum(dnorm(b2n,0,sqrt(10),log=TRUE))-sum(dnorm(b2,0,sqrt(10),log=TRUE)))>log(runif(1))){
    b1   = b1n
    b2   = b2n
    accB1[m] = 1
    llo  = lln
    V    = Vn
  }

  # block 2
  repeat{
    b.prop = rnorm(dm*2,c(b1,b2),sd=propsdb)
    bn     = b.prop*block2+c(b1,b2)*block1
    
    b1n = bn[1:dm]
    b2n = bn[(dm+1):(2*dm)]
    B1  = Outer(b1n,b1n)
    B2  = Outer(b2n,b2n)
    B0  = (Oiota-B1-B2)*Sbar
    
    cond1 = b1n[1]>0
    cond2 = b2n[1]>0
    cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1
    cond4 = sum(abs(B1+B2)<1)==dm^2
    if(cond1 && cond2 && cond4) {
      break
    }
  }
  
  if(cond3){
    for(t in 1:TT){
      Vn[[t]]  = (B0+B1*G1[[t]]+B2*G2[[t]])
    }
      lln <- future_mapply(diwish.t,Sig,Vn)
  } else {lln=-Inf}
  
  if((sum(lln)-sum(llo)+
      sum(dnorm(b1n,0,sqrt(10),log=TRUE))-sum(dnorm(b1,0,sqrt(10),log=TRUE))+
      sum(dnorm(b2n,0,sqrt(10),log=TRUE))-sum(dnorm(b2,0,sqrt(10),log=TRUE)))>log(runif(1))){
    b1   = b1n
    b2   = b2n
    accB2[m] = 1
    llo  = lln
    V    = Vn
  }

  B1  = Outer(b1,b1)
  B2  = Outer(b2,b2)
  B0  = (Oiota-B1-B2)*Sbar
  
  ##-----
  ## nu
  ##-----

  nun  = truncnorm::rtruncnorm(1,a = dm+1,mean = nu,sd = propsdnu) 
  
  diwish.t <- function(x,y){LaplacesDemon::dinvwishart(x, nun, (nun-dm-1)*y, log=TRUE)}
  lln      <- future_mapply(diwish.t,Sig,V)
  
  if((sum(lln)-sum(llo)+
      dexp(nun,0.01,log=TRUE)-dexp(nu,0.01,log=TRUE)+
      log(truncnorm::dtruncnorm(nun,a = dm+1,b=Inf,mean = nu,sd = propsdnu))-
      log(truncnorm::dtruncnorm(nu,a = dm+1,b=Inf,mean = nun,sd = propsdnu)))>log(runif(1))){
    llo   = lln
    accnu[m] = 1
    nu    = nun
  }
  
  ##-----
  ## Collect results and predict
  ##-----
  
  resc[m,] = c(lag,nu,b1,b2)
  LLH[m]   = sum(llo)
  
  if (m>bi){
    Vpred[[m-bi]] = B0+B1*Sig[[TT]]+B2*Reduce('+',Sig[(TT+1-lag):TT])/lag
  }
  
  if(!m%%100){
    print(paste('Model 4: ',round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}
print("Done!")

plan(sequential)

post.size = 5000
ind       = round(seq(1,M,length=post.size))
r         = resc[(bi+1):(bi+M),]


res = list(Vpred[ind],r[ind,],
           accl[(bi+1):(bi+M)][ind],
           accnu[(bi+1):(bi+M)][ind],
           accB1[(bi+1):(bi+M)][ind],
           accB2[(bi+1):(bi+M)][ind],LLH[(bi+1):(bi+M)][ind])
names(res) = c('Vpred','r','accl','accnu','accB1','accB2','LLH')
save(res,file='temp/results_xm.Rdata')


