##-----------------------
## TABLE 1
##-----------------------
rm(list=ls(all=TRUE))

load("temp/individual_llh.RData")

lpbf = c(sum(lL_static[,]),
         sum(apply(lL_dcc[,],2,median)),
         sum(apply(lL_tdcc[,],2,median)),
         sum(apply(lL_xm[,],2,median)),
         sum(apply(lL_caw[,],2,median)),
         sum(apply(lL_ht,2,median)))
names(lpbf) = c('Static','DCC','DCC-t','AIW','CAW','DCC-HEAVY-t')
lpbfdf = as.data.frame(t(lpbf))

print(xtable(lpbfdf,align= 'ccccccc',
             caption = '1-step-ahead $LPS$ 
             for all individual models: Static, DCC, DCC-t, AIW, CAW,
and DCC-HEAVY model with $t$ copula for 2020/01/02 - 2023/01/31 out-of-sample period
(K = 797 observations).',
             label = 'table:lps_FX', digits = 2),
      file='tables_and_figures/lps_FX.tex',
      include.rownames = FALSE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE)

##-----------------------
## FIGURES 2, 3 and TABLE 2
##-----------------------

rm(list=ls(all=TRUE))

load(file="temp/results_FX_pooling.Rdata")

move.axis = 100
atx       =  seq(date[(nn+1)], date[(nn+K)], by=25)

# FIGURE 2

pdf('tables_and_figures/weights_FX.pdf',height=8,width=14)
par(mfrow=c(2,1), mar=c(4, 3, 1, 1) + 0.1)
plot(date[(nn+1):(nn+K)],mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
par(new = TRUE)
plot(date[(nn+1):(nn+K)],apply(ws_jore1,2,median),ylim=c(0,1),
     type='l',ylab='',xlab='',xaxt="n",lwd=3,xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]),col='pink')
lines(date[(nn+1):(nn+K)],apply(ws_jore10,2,median),col='pink4',lwd=2,lty=2)
lines(date[(nn+1):(nn+K)],apply(ws_DN,2,median),lty=1,lwd=2,col='blue')
lines(date[(nn+1):(nn+K)],apply(ws_gew,2,median),col=1,lwd=2,lty=1)
abline(h=0.5)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=1,col=c('pink','pink4','blue',1),
       lty=c(1,2,1,1),lwd=c(2,2,2,2),
       legend=c('Jore1','Jore10','DelNegro','Geweke'))
plot(tail(date,K),mkvol,type='l',axes = FALSE,
     col='gray90',lwd=3,ylab='',xlab='',
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
par(new = TRUE)
plot(tail(date,K),cumsum(apply(lL_caw,2,median))-cumsum(apply(lL_caw,2,median)),
     type='l',ylim=c(-500,500),ylab='',xlab='',xaxt="n",
     xlim=c(date[(nn+1)]-move.axis,date[(nn+K)]))
lines(tail(date,K),cumsum(apply(lL_tdcc,2,median))-cumsum(apply(lL_caw,2,median)),lwd=2,
      col='coral3')
lines(tail(date,K),cumsum(apply(lL_gew,2,median))-cumsum(apply(lL_caw,2,median)),
      col=1,lwd=2)
lines(tail(date,K),cumsum(apply(lL_jore1,2,median))-cumsum(apply(lL_caw,2,median)),
      col='pink',lwd=2)
lines(tail(date,K),cumsum(apply(lL_jore10,2,median))-cumsum(apply(lL_caw,2,median)),
      col='pink4',lwd=2,lty=2)
lines(tail(date,K),cumsum(apply(lL_DN,2,median))-cumsum(apply(lL_caw,2,median)),lwd=2,
      col='blue')
lines(tail(date,K),cumsum(apply((lL_equal),2,median))-cumsum(apply(lL_caw,2,median)),
      col='coral',lwd=2,lty=2)
lines(tail(date,K),cumsum(apply(lL_ht,2,median))-cumsum(apply(lL_caw,2,median)),
      col='violet',lwd=2)
axis(1, at=atx, labels=format(atx, "%Y/%m"),las=2,cex.axis=0.75)
legend(x=date[(nn+1)]-move.axis-35,y=500,col=c('coral3',1,'pink','pink4','blue','coral','violet'),
       lty=c(1,1,1,2,1,2,1),lwd=c(2,2,2,2,2,2,2),
       legend=c('DCC-t','Geweke','Jore1','Jore10','DelNegro','Equal','DCC-HEAVY-t'))
dev.off()


# FIGURE 3

pdf('tables_and_figures/lpts_FX.pdf',height=8,width=12)
par(mfrow=c(2,2))
plot(density(apply(lL_caw, 1,mean)),xlim=c(-4,-2.9),ylab='',xlab='',main='LPS',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc, 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht, 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew, 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1, 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal, 1,mean)),lwd=2,col='coral',lty=2)
legend(x=-4,y=110,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))
plot(density(apply(lL_caw[,ind25], 1,mean)),xlim=c(-0.65,0.25),ylab='',xlab='',main='LPTS(25%)',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc[,ind25], 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht[,ind25], 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew[,ind25], 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1[,ind25], 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal[,ind25], 1,mean)),lwd=2,col='coral',lty=2)
legend(x=-0.65,y=90,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))
plot(density(apply(lL_caw[,ind10], 1,mean)),xlim=c(0.2,1),ylab='',xlab='',main='LPTS(10%)',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc[,ind10], 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht[,ind10], 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew[,ind10], 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1[,ind10], 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal[,ind10], 1,mean)),lwd=2,col='coral',lty=2)
legend(x=0.2,y=90,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))
plot(density(apply(lL_caw[,ind05], 1,mean)),xlim=c(0.4,1.2),ylab='',xlab='',main='LPTS(5%)',
     lwd=2,lty=1,col='royalblue')
lines(density(apply(lL_tdcc[,ind05], 1,mean)),lwd=2,col='coral3')
lines(density(apply(lL_ht[,ind05], 1,mean)),lwd=2,col='violet')
lines(density(apply(lL_gew[,ind05], 1,mean)),lwd=2,col='black')
lines(density(apply(lL_jore1[,ind05], 1,mean)),lwd=2,col='pink')
lines(density(apply(lL_equal[,ind05], 1,mean)),lwd=2,col='coral',lty=2)
legend(x=0.4,y=90,lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,1,2),
       col=c('coral3','royalblue','violet','black','pink','coral'),
       legend=c('DCC-t','CAW','DCC-HEAVY-t','Geweke','Jore1','Equal'))
dev.off()

## TABLE 2

print(xtable::xtable(corrs_res,
                     caption = "Posterior medians of sample correlations 
             between three proxies for the market volatility and
             the preference for high-frequency model for
             2020/01/02 to 2023/01/31 out-of-sample period ($K=797$ observations).
             The proxies for the market volatility are: average standardized 
             realized volatility (avrg RV), 
             equally weighted market portfolio realized volatility (Mkt:eql) 
             and VIX index. The preference for the high-frequency model is measured as a 
             high-frequency 
             component weight in various  pooling schemes 
             as well as the difference between the daily log likelihood (logLik diff) 
             between the CAW and DCC-t models.",
                     label = 'table:corrs_FX', digits = 3,
                     align = "lccc|cccc"),
      file='tables_and_figures/corrs_FX.tex',
      include.rownames = TRUE,latex.environments = "center" ,
      caption.placement = "top",
      include.colnames= TRUE,
      rotate.colnames = FALSE,scalebox = 0.8,
      hline.after = c(-1, 0, 3,nrow(corrs_res)))


##-----------------------
## TABLE 3
##-----------------------

rm(list=ls(all=TRUE))

load('temp/read_portfolio.RData')

tableLines <- print(xtable(ALL.res,digits = 3,caption="Portfolio allocation results based on 1-step-ahead
              predictions for 2020/01/02 to 2023/01/31 out-of-sample period ($K$ = 797 observations).
              The three portfolios are:
             Global Minimum Variance (GMV), and minimum Conditional Value at Risk
             for lower 5 and 10 percentiles (CVaR05 and CVaR10), all with short-sale constraints.
             The table reports average portfolio return (avrg),
             overall standard deviation in \\% (sd),
             the p-value corresponding to the model confidence set of
             Hansen, Lunde and Nason (2011) using a significance level of 5\\%,
             adjusted Sharpe ratio (SR), realized Conditional Value at Risk
             for lower 5 and 10 percentiles  (rCVaR05 and rCVaR10),
             turnover (TO), and concentration (CO), all quantities annualized, for the
             pooled models (Geweke's, Jore's and
             equally weighted),
            two best individual models (CAW  and DCC-t) and a competitor model (DCC-HEAVY-t).
            In gray is highlighted the best performing portfolio and the second best is underlined.",
                           align = "lccc|ccc",label='table:gmvfull_FX_new2'),
                    scalebox=0.8,sanitize.text.function=function(x){x},
                    hline.after = c(-1, 0, 8,16,nrow(ALL.res)),
                    caption.placement = "top")
writeLines (tableLines, con = "tables_and_figures/gmvfull_FX_new2.tex")




