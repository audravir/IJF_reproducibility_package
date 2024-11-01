PMCMC_delNegro = function(first_llik,second_llik,Npart,prior,M,propsd){
  t0   = Sys.time()
  t1   = Sys.time()
  
  y1 = apply(exp(first_llik),2,median)
  y2 = apply(exp(second_llik),2,median)
  K = length(y1)
  priorm = prior[1]
  priorsd = prior[2]
  
  bi = M
  
  ##-------
  ## Particle MCMC
  ##-------
  
  betas = rep(NA,M+bi)
  xssp  = wssp = matrix(NA,ncol=M+bi,nrow=K)
  acc   = rep(0,M+bi)
  lams  = pnorm(rnorm(K,0,1))
  llo   = sum(log(lams*y1+(1-lams)*y2))
  beta  = truncnorm::rtruncnorm(1,a=-1,b=1,0.1,sqrt(0.1))
  tau2  = 1-beta^2
  xso = xsn = wsso = rep(0,K)
  
  for(m in 1:(M+bi)){
    betan  = truncnorm::rtruncnorm(1,a=-1,b=1,beta,propsd)
    tau2n  = 1-betan^2
    
    xn    = rnorm(Npart,0,1)
    wssn   = rep(0,K)
    
    for(t in 1:K){
      # blind propagate
      xs  = rnorm(Npart,betan*xn,sqrt(tau2n))
      lams = pnorm(xs)
      # weights
      ws  = lams*y1[t]+(1-lams)*y2[t]
      #if(sum(ws)==0){ws=1/Npart}
      ind = sample(1:Npart,Npart,replace = T,prob=ws/sum(ws))
      xn = xs[ind]
      wssn[t] = mean(ws[ind])
      xsn[t] = mean(xn)
    }
    
    lln = sum(log(wssn))
    if((lln-llo+log(dtruncnorm(betan,-1,1,priorm,priorsd))-
        log(dtruncnorm(beta,-1,1,priorm,priorsd)))>log(runif(1))){
      llo  = lln
      beta = betan
      tau2 = tau2n
      acc[m]  = 1
      xso = xsn
      wsso = wssn
    }
    if(!m%%100){
      print(paste('DelNegro weights: ',round(m/(M+bi)*100,2),"%",sep=""))
      print(Sys.time()-t1)
      print(Sys.time()-t0)
      t1   = Sys.time()
    }
    
    betas[m]  = beta
    xssp[,m]  = xso
    wssp[,m]  = wsso
  }
  print("Done!")
  resdn = list(betas[(bi+1):(bi+M)],xssp[,(bi+1):(bi+M)],
               wssp[,(bi+1):(bi+M)],acc[(bi+1):(bi+M)])
  names(resdn) = c('beta','weights_xs','likelihood','acc')
  return(resdn)
}