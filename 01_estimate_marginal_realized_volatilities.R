rm(list=ls(all=TRUE))

load('data/FXdata.Rdata')

nn  = dim(RVs)[1]
dm  = dim(RVs)[2]
oos = 797

# log-HAR(1,5,22)
RVsforc1_m3 = matrix(NA,ncol=dm,nrow=oos)
for(i in 1:dm){
  rvt   = log(RVs[,i])
  rvt1  = c(NA,rvt[-length(rvt)])
  rvt5  = rollmean(rvt1,5,fill = c(NA,NA,NA),align = 'right')
  rvt22 = rollmean(rvt1,22,fill = c(NA,NA,NA),align = 'right')
  dtf   = data.frame(rvt,rvt1,rvt5,rvt22)
  lrvm4 = lm(rvt~rvt1+rvt5+rvt22,data = dtf[1:(nn-oos),])
  RVsforc1_m3[,i]  = exp(predict(lrvm4, newdata = dtf[(nn-oos+1):nn,]))
}

RV_forc        = list(tail(date,oos),RVsforc1_m3)
names(RV_forc) = c('Date','1sa')
save(RV_forc,file='temp/RV_forc.Rdata')





