nn       = length(date)
end.date = which(zoo::as.yearmon(date)=="ene 2020")[1]-1
if(is.na(date[end.date])){end.date = which(zoo::as.yearmon(date)=="jan 2020")[1]-1}

T0  = end.date  #
K = nn-end.date

data  = stand[1:T0,]
udata = udata[1:T0,]
Sig   = Sigma[1:T0]

propsd   = 0.002
propsdnu = 0.5

t0   = Sys.time()
t1   = Sys.time()

TT   = dim(data)[1]
dm   = dim(data)[2]

R    = array(NA,c(dm, dm, TT))
pold   <- c(0.09,0.8)
nuold  <- 50
tdata  <- qt(udata,nuold)
Rbar   <- cor(tdata)
llold  <- rep(0,TT)
bi     = min(M,50000)
TIMING = rep(NA,M+bi)
LLH    <- rep(NA,M+bi)
R[,,1] <- Rbar
Rnew   = R
Pbar   = Reduce('+',Sig)/TT
Rtilde = (1-pold[2])*Rbar-pold[1]*Pbar
IPD    = rep(TRUE,TT)
is.positive.definite(Rtilde)

res <- matrix(NA,ncol=3,nrow=M+bi)
accdcc   = accnu = rep(0,bi+M)
Rpred    = vector(mode = "list", length = M)


for(t in 2:TT){
  R[,,t]   <- Rtilde+pold[1]*Sig[[t-1]]+pold[2]*R[,,t-1]
  llold[t] <- mvnfast::dmvt(tdata[t,], rep(0,dm),R[,,t],df=nuold,log=TRUE)-
    sum(dt(tdata[t,],df=nuold,log = TRUE))
}

for(m in 1:(M+bi)){
  t2 = Sys.time()
  
  repeat{
    pnew  = rnorm(2,pold,sd=propsd)
    if(all(pnew>0)) break
  }
  
  llnew = rep(0,TT)
  for(t in 2:TT){
    Rnew[,,t]   <- (1-pnew[2])*Rbar-pnew[1]*Pbar+pnew[1]*Sig[[t-1]]+pnew[2]*Rnew[,,t-1]
    IPD[t] = prod(eigen(Rnew[,,t],symmetric = TRUE,only.values = TRUE)$values>0)==1 
  }
  
  if (all(IPD==TRUE)){
    for(t in 2:TT){
      llnew[t] = mvnfast::dmvt(tdata[t,], rep(0,dm),Rnew[,,t],df=nuold,log=TRUE)-
        sum(dt(tdata[t,],df=nuold,log = TRUE))
    }
  } else {llnew = -Inf}
  
  if((sum(llnew)-sum(llold)+
      dbeta(pnew[1],3,10,log=TRUE)-dbeta(pold[1],3,10,log=TRUE)+
      dbeta(pnew[2],10,3,log=TRUE)-dbeta(pold[2],10,3,log=TRUE))>log(runif(1)))
  {
    llold  = llnew
    pold   = pnew
    R      = Rnew
    accdcc[m] = 1
  }
  
  ## nu
   nunew = truncnorm::rtruncnorm(1,5,b = Inf,mean = nuold,sd = propsdnu)
   tdata <- qt(udata,nunew)
   Rbar  <- cor(tdata)
   llnew = rep(0,TT)
   
   for(t in 2:TT){
     Rnew[,,t]   <- (1-pold[2])*Rbar-pold[1]*Pbar+pold[1]*Sig[[t-1]]+pold[2]*Rnew[,,t-1]
     IPD[t] = prod(eigen(Rnew[,,t],symmetric = TRUE,only.values = TRUE)$values>0)==1 
   }
   
   if (all(IPD==TRUE)){
     for(t in 2:TT){
       llnew[t] = mvnfast::dmvt(tdata[t,], rep(0,dm),Rnew[,,t],df=nunew,log=TRUE)-
         sum(dt(tdata[t,],df=nunew,log = TRUE))
     }
   } else {llnew = -Inf}
  
   if((sum(llnew)-sum(llold)+
       dexp(nunew,0.01,log=TRUE)-dexp(nuold,0.01,log=TRUE)+
       log(truncnorm::dtruncnorm(nunew,5,Inf,mean = nuold,sd = propsdnu))-
       log(truncnorm::dtruncnorm(nuold,5,Inf,mean = nunew,sd = propsdnu)))>log(runif(1)))
   {
     llold  = llnew
     nuold  = nunew
     R      = Rnew
     accnu[m] = 1
   }
   
   tdata  = qt(udata,nuold)
   Rbar   <- cor(tdata)
   
  ##-----
  ## Collect results and prediction
  ##-----
  
  res[m,] <- c(nuold,pold)
  LLH[m]  <- sum(llold)
  
  if(m>bi){
    Rpred[[m-bi]] <- (1-pold[2])*Rbar-pold[1]*Pbar+pold[1]*Sig[[TT]]+pold[2]*R[,,TT]
  }
  
  if(!m%%100){
    print(paste('Model 3: ',round(m/(M+bi)*100,2),"%",sep=""))
    print(Sys.time()-t1)
    print(Sys.time()-t0)
    t1   = Sys.time()
  }
  TIMING[m] = Sys.time()-t2
}
print("Done!")


post.size = 5000
ind       = round(seq(1,M,length=post.size))
r         = res[(bi+1):(bi+M),]

res = list(Rpred[ind],r[ind,],
           accdcc[(bi+1):(bi+M)][ind],accnu[(bi+1):(bi+M)][ind],
           LLH[(bi+1):(bi+M)][ind])
names(res) = c('Rpred','r','accdcc','accnu','LLH')

save(res,file='temp/results_heavy_t.Rdata')

