rm(list=ls(all=TRUE))

load('temp/FX_portfolio_at_median.Rdata')
load('data/FXdata.Rdata')

# load risk-free data as for Sharpe ratio

DGS10      = read.csv("data/DGS10.csv")
DGS10$date = as.Date(DGS10$DATE, "%Y-%m-%d")
DGS10$rf   = as.numeric(DGS10$DGS10)
DGS10$rf   = na.locf(DGS10$rf)
ind        = DGS10$date %in% tail(date,K)
rf         = DGS10$rf[ind]

# define expected shortfall function
esfun=function(x,p){
  es=mean(x[which(x<quantile(x,p))])
  return(es)
}

all.ws        = list(ws_gmvr,ws_CVAR05r,ws_CVAR10r)
names(all.ws) = c('gmvr','CVAR05r','CVAR10r')
all.prets     = all.co = all.to = all.portsd = list()
pret          = co = to = portsd = array(NA,dim=c(length(models),K))

for(j in 1:length(all.ws)){
  for(i in 1:length(models)){
    for(t in 1:K){
      pret[i,t]   = sum(all.ws[[j]][i,t,]*rets[nn+t,])
      portsd[i,t] = sqrt(t(all.ws[[j]][i,t,])%*%RCov[,,nn+t]%*%all.ws[[j]][i,t,])
      co[i,t]     = (sum((all.ws[[j]][i,t,])^2))^(1/2)
      if(t<K)  to[i,t]   = sum(abs(all.ws[[j]][i,t+1,]-all.ws[[j]][i,t,]*((1+rets[nn+t,])/(1+pret[i,t]))))
    }
  }
  all.prets[[j]] = pret
  all.co[[j]]    = co
  all.to[[j]]    = to
  all.portsd[[j]]    = portsd
}

names(all.prets) = names(all.co) = names(all.to) = names(all.portsd) = names(all.ws)

test.results = list()

# GMV
loss.sd <- sqrt(t((all.prets[[1]] - apply(all.prets[[1]],1,mean))^2))
colnames(loss.sd) = models
SSM.sd <- MCSprocedure(Loss = loss.sd, alpha = 0.05, B = 10000, statistic = "Tmax")
test.results[[1]] = SSM.sd@show[,6]

# CVAR05
loss.sd <- sqrt(t((all.prets[[2]] - apply(all.prets[[2]],1,mean))^2))
colnames(loss.sd) = models
SSM.sd <- MCSprocedure(Loss = loss.sd, alpha = 0.05, B = 10000, statistic = "Tmax")
test.results[[2]] = SSM.sd@show[,6]

# CVAR10
loss.sd <- sqrt(t((all.prets[[3]] - apply(all.prets[[3]],1,mean))^2))
colnames(loss.sd) = models
SSM.sd <- MCSprocedure(Loss = loss.sd, alpha = 0.05, B = 10000, statistic = "Tmax")
test.results[[3]] = SSM.sd@show[,6]

# Group all results
# ncol = number of models
# nrow = return + sd + (p-value) +SR+CVAR+CVAR+TO+CO (8)*3 (portfolios)

ALL.res = matrix(NA,nrow=24,ncol=length(models))

# 1-8 GMV portfolio
# 9-16 CVAR05
# 17-24 CVAR10

# 1 return 2 sd 3 (p-value) 4 SR 5 CVAR05 6 CVAR10 7 TO 8 CO 

for(i in 1:3){
  ALL.res[i*8-7,] = lapply(all.prets, function(x) apply(x,1,mean)*252)[[i]]
  ALL.res[i*8-6,] = lapply(all.prets, function(x) apply(x,1,sd)*sqrt(252))[[i]]
  if(length(test.results[[i]])<length(models)){
    y <- test.results[[i]][models]
    names(y) = models
    ALL.res[i*8-5,] = y
  }
  else {
    ALL.res[i*8-5,] = test.results[[i]]
  }
  ALL.res[i*8-4,] = lapply(all.prets, function(x) apply(sweep(x*252, 2, rf, '-'),1,mean)/(apply(x,1,sd)*sqrt(252)))[[i]]
  ALL.res[i*8-3,] = lapply(all.prets, function(x) apply(x*252,1,esfun,0.05))[[i]]
  ALL.res[i*8-2,] = lapply(all.prets, function(x) apply(x*252,1,esfun,0.10))[[i]]
  ALL.res[i*8-1,] = lapply(all.co, function(x) apply(x,1,mean))[[i]]
  ALL.res[i*8,] = lapply(all.to, function(x) apply(x,1,mean,na.rm=TRUE))[[i]]
}

colnames(ALL.res) = models
rownames(ALL.res) = rep(c('avrg','sd','(p-value)','SR','rCVaR05','rCVaR10','CO','TO'),3)

round(ALL.res,3)

save(ALL.res,file='temp/read_portfolio.RData')



