nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2020")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}

T0 = end.date  # for estimation
K  = nn-end.date # for oos evaluation

data = stand[1:T0,]

# M    = 50000 is defined in mother file

propsd   = 0.0001
propsdnu = 0.5

###

t0   = Sys.time()
t1   = Sys.time()

TT   = dim(data)[1]
dm   = dim(data)[2]
bi   = min(M,50000)
TIMING = rep(NA,M+bi)

udata = udata[1:TT,]

Qold   = array(NA,c(dm, dm, TT))
aold   <- rep(0.12,dm)
bold   <- rep(0.99,dm)
nuold  <- 10
tdata  <- qt(udata,nuold)
llold  <- rep(0,TT)
LLH    <- rep(NA,M+bi)
Qold[,,1] = cora(tdata)
Rold   <- Qold

resdcc <- matrix(NA,ncol=dm*2+1,nrow=M+bi)
iota   = rep(1,dm)
Oiota  = Outer(iota,iota)
Sbar   = cov(tdata)
A      = Outer(aold,aold)
B      = Outer(bold,bold)
B0     = (Oiota-A-B)*Sbar

accdcc1 = accdcc2 = rep(0,bi+M)
accnu  <- rep(0,bi+M)
Vpred  = Qpred = vector(mode = "list", length = M)

for(t in 2:TT){
  Qold[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Qold[,,(t-1)]
  t.ma  = Qold[,,t]
  t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
  t.R   = Outer(t.dv,t.dv)*t.ma
  Rold[,,t] = t.R
  inlik = sum(dt(tdata[t,],df=nuold,log=TRUE))
  llold[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm),t.R, nuold, log=TRUE)-inlik
}

for(m in 1:(M+bi)){
  t2 = Sys.time()
  
  ##-----
  ## bs split randomly
  ##-----
  
  block1 <- sample(c(TRUE,FALSE),size=dm*2,replace = TRUE)
  block2 <- (!block1)
  
  ##-----
  ## bs
  ##-----

  # block 1
  repeat{
    b.prop = rnorm(dm*2,c(aold,bold),sd=propsd)
    bn     = b.prop*block1+c(aold,bold)*block2
    
    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    A    = Outer(anew,anew)
    B    = Outer(bnew,bnew)
    B0   = (Oiota-A-B)*Sbar
    
    cond1 = anew[1]>0
    cond2 = bnew[1]>0
    # cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 
    cond4 = (sum(abs(A+B)<1)==dm^2)
    if(cond1 && cond2 && cond4) {
      break
    }
  }
  
  IPD   <- rep(1,TT)
  llnew <- rep(0,TT)
  Qnew  = Qold
  Rnew  = Rold
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Qnew[,,(t-1)]
    t.ma  = Qnew[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    Rnew[,,t] = t.R
    IPD[t] =  prod(eigen(t.R,symmetric = TRUE,only.values = TRUE)$values>0)==1
  }
  
  if(sum(IPD)==TT){
    for(t in 2:TT){
      inlik    <- sum(dt(tdata[t,],df=nuold,log=TRUE))
      llnew[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), Rnew[,,t], df = nuold, log=TRUE)-
        inlik
    } 
  } else {llnew=-Inf}
  
  if((sum(llnew)-sum(llold)+
      sum(dnorm(anew,0,sqrt(10),log=TRUE))-sum(dnorm(aold,0,sqrt(10),log=TRUE))+
      sum(dnorm(bnew,0,sqrt(10),log=TRUE))-sum(dnorm(bold,0,sqrt(10),log=TRUE)))>log(runif(1)))
  {
    llold  = llnew
    aold   = anew
    bold   = bnew
    Rold   = Rnew
    Qold   = Qnew
    accdcc1[m] = 1
  }
  
  
  # block 2
  repeat{
    b.prop = rnorm(dm*2,c(aold,bold),sd=propsd)
    bn     = b.prop*block2+c(aold,bold)*block1
    
    anew = bn[1:dm]
    bnew = bn[(dm+1):(2*dm)]
    A    = Outer(anew,anew)
    B    = Outer(bnew,bnew)
    B0   = (Oiota-A-B)*Sbar
    
    cond1 = anew[1]>0
    cond2 = bnew[1]>0
    # cond3 = prod(eigen(B0,symmetric = TRUE,only.values = TRUE)$values>0)==1 
    cond4 = (sum(abs(A+B)<1)==dm^2)
    if(cond1 && cond2 && cond4) {
      break
    }
  }
  
  IPD   <- rep(1,TT)
  llnew <- rep(0,TT)
  Qnew  = Qold
  Rnew  = Rold
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Qnew[,,(t-1)]
    t.ma  = Qnew[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    Rnew[,,t] = t.R
    IPD[t] =  prod(eigen(t.R,symmetric = TRUE,only.values = TRUE)$values>0)==1
  }
  
  if(sum(IPD)==TT){
    for(t in 2:TT){
      inlik    <- sum(dt(tdata[t,],df=nuold,log=TRUE))
      llnew[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm), Rnew[,,t], df = nuold, log=TRUE)-
        inlik
    } 
  } else {llnew=-Inf}
  
  if((sum(llnew)-sum(llold)+
      sum(dnorm(anew,0,sqrt(10),log=TRUE))-sum(dnorm(aold,0,sqrt(10),log=TRUE))+
      sum(dnorm(bnew,0,sqrt(10),log=TRUE))-sum(dnorm(bold,0,sqrt(10),log=TRUE)))>log(runif(1)))
  {
    llold  = llnew
    aold   = anew
    bold   = bnew
    Rold   = Rnew
    Qold   = Qnew
    accdcc2[m] = 1
  }
  
  ##-----
  ## nu
  ##-----

  nunew  = truncnorm::rtruncnorm(1,a = 5,mean = nuold,sd = propsdnu)
  
  tdata  = qt(udata,nunew)
  Sbar   = cova(tdata)
  
  llnew <- rep(0,TT)
  Qnew  = Qold
  Qnew[,,1] = cora(tdata)
  Rnew  = Qnew
  A     = Outer(aold,aold)
  B     = Outer(bold,bold)
  B0    = (Oiota-A-B)*Sbar
  
  for(t in 2:TT){
    Qnew[,,t] <- B0+A*Outer(tdata[t-1,],tdata[t-1,])+B*Qnew[,,(t-1)]
    t.ma  = Qnew[,,t]
    t.dv  = t.ma[ col(t.ma)==row(t.ma) ]^{-1/2} 
    t.R   = Outer(t.dv,t.dv)*t.ma
    Rnew[,,t] = t.R
    inlik = sum(dt(tdata[t,],df=nunew,log=TRUE))
    llnew[t]  <- mvnfast::dmvt(tdata[t,], rep(0,dm), t.R, nunew, log=TRUE)-inlik
  }
  
  if((sum(llnew)-sum(llold)+
      dexp(nunew,0.01,log=TRUE)-dexp(nuold,0.01,log=TRUE)+
      log(truncnorm::dtruncnorm(nunew,a = 5,b=Inf,mean = nuold,sd = propsdnu))-
      log(truncnorm::dtruncnorm(nuold,a = 5,b=Inf,mean = nunew,sd = propsdnu)))>log(runif(1)))
  {
    llold  = llnew
    nuold  = nunew
    Qold   = Qnew
    Rold   = Rnew
    accnu[m] = 1
  }
  
  tdata  = qt(udata,nuold)
  Sbar   = cova(tdata)
  
  ##-----
  ## Collect results and prediction
  ##-----
  
  resdcc[m,] <- c(nuold,aold,bold) 
  LLH[m]     <- sum(llold)
  
  if(m>bi){
    A      = Outer(aold,aold)
    B      = Outer(bold,bold)
    B0     = (Oiota-A-B)*Sbar
    Qpred[[m-bi]] <- B0+A*Outer(tdata[TT,],tdata[TT,])+B*Qold[,,TT]
    Vpred[[m-bi]] <- diag(diag(Qpred[[m-bi]])^{-1/2})%*%Qpred[[m-bi]]%*%diag(diag(Qpred[[m-bi]])^{-1/2})
  }

  if(!m%%100){
    print(paste('Model 2: ',round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}  
print("Done!")


nu = resdcc[,1]
b1 = resdcc[,2:(dm+1)]
b2 = resdcc[,(dm+2):(2*dm+1)]

post.size = 5000
ind       = round(seq(1,M,length=post.size))

r   = resdcc[(bi+1):(bi+M),]
res = list(Vpred[ind],r[ind,],
           accdcc1[(bi+1):(bi+M)][ind],accdcc2[(bi+1):(bi+M)][ind],accnu[(bi+1):(bi+M)][ind],
           LLH[(bi+1):(bi+M)][ind])
names(res) = c('Vpred','r','accdcc1','accdcc2','accnu','LLH')
save(res,file=paste('temp/results_vectordcc_tcop.Rdata',sep=''))

