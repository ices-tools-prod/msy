
directory<- c("C:\\home\\CADIGAN\\stocks\\SR\\nonparam\\help_example\\")
setwd(directory)

path.to.file <- directory

source('Rfunc.txt')
require("numDeriv")
require("scam")
require("xtable")
require("nlstools")

########################################################;
## read data;
########################################################;


sscale = 1
s.name = c('SSB (t)') 
rscale = 1
r.name = c('  Recruit (millions)') 
p.name = c('Recruit per spawner  ')

## random sample;
sr.data = read.table(paste(directory,'sample.dat',sep=''),col.name=c('s','r'))  
## population - the truth;
sr.data.pop = read.table(paste(directory,'pop.dat',sep=''),col.name=c('s','r'))

n = length(sr.data$s)                 
sr.data$year=1:n
s <- sr.data$s
r <- sr.data$r
year = sr.data$year

start.alpha.bh <- max(r)
ms <- mean(s)
mr<- mean(r)
n = length(r)

x = s/sscale
y = r/rscale
ys = y/x
log.ys = log(ys)

wt.ext = rep(1,n) 

srdat = data.frame(logrec=log(y),logssb=log(x),ssb=x,weights=wt.ext)

### Beverton-Holt  #####

init.Rmax = 600
init.S50 = 250
init.alpha= init.Rmax
init.beta = init.S50
ls.slope = sum(x*y)/sum(x*x) 
Rmax = max(y)
Pmax = 2*max(y/x)

BH.fit <- nls(logrec ~ log(alpha) + logssb - log(beta + ssb),data=srdat,
           start = list(alpha = init.alpha,beta = init.beta),
           algorithm="port",lower=c(0.1,1),upper=c(Rmax,1000))

summary(BH.fit)
BH.parm = coef(BH.fit)

ssb.max=max(x)*3
ssb.pred.points = seq(0.1,ssb.max,length=500)
BH.pred.rec = exp(predict(BH.fit,list(ssb=ssb.pred.points,logssb=log(ssb.pred.points))))  
BH.pred.prod =  BH.pred.rec/ssb.pred.points

BH.raw.res=resid(BH.fit)
BH.sig2 = sum(BH.raw.res**2)/(sum(wt.ext)-2)
BH.res = BH.raw.res/sqrt(BH.sig2)

BH.boo <- nlsBoot(BH.fit)

plot(BH.boo)
file.name <- paste(path.to.file,'BHparms_boo', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)
dev.off()

### Hockey-Stick ##########

init.alpha= 0.5*init.Rmax/init.S50
init.delta = init.S50
gam2by4 = 0.1/4
delta.max = max(x) 
delta.min = min(x)

#rec = log(alpha) + log(ssb + sqrt(delta**2 + gam2by4) - sqrt((ssb-delta)**2 + gam2by4))
#plot(ssb,rec)

HS.fit <- nls(logrec ~ log(alpha) + log(ssb + sqrt(delta**2 + gam2by4) - sqrt((ssb-delta)**2 + gam2by4)),
           data=srdat,
           start = list(alpha = init.alpha,delta = init.delta),
           algorithm="port",lower=c(0,delta.min),upper=c(10,delta.max))

summary(HS.fit)
HS.parm = coef(HS.fit)

HS.pred.rec = exp(predict(HS.fit,list(ssb=ssb.pred.points)))  
HS.pred.prod =  HS.pred.rec/ssb.pred.points

HS.raw.res=resid(HS.fit)
HS.sig2 = sum(HS.raw.res**2)/(sum(wt.ext)-2)
HS.res = HS.raw.res/sqrt(HS.sig2)

HS.boo <- nlsBoot(HS.fit)

plot(HS.boo)
file.name <- paste(path.to.file,'HSparms_boo', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)
dev.off()  

### Ricker ##########

init.alpha= init.Rmax/init.S50
init.beta = init.alpha/(2.72*init.Rmax)
srdat$one=1

RK.fit <- nls(logrec ~ log(alpha) + log(beta) + one + logssb - beta*ssb,
           data=srdat,
           start = list(alpha = init.alpha,beta = init.beta),
           algorithm="port",lower=c(1e-8,1e-8),upper=c(Rmax,0.5))

summary(RK.fit)
RK.parm = coef(RK.fit)

RK.pred.rec = exp(predict(RK.fit,list(ssb=ssb.pred.points,logssb=log(ssb.pred.points),
one = rep(1,length(ssb.pred.points)))))  
RK.pred.prod =  RK.pred.rec/ssb.pred.points

RK.raw.res=resid(RK.fit)
RK.sig2 = sum(RK.raw.res**2)/(sum(wt.ext)-2)
RK.res = RK.raw.res/sqrt(RK.sig2)

RK.boo <- nlsBoot(RK.fit)

plot(RK.boo)
file.name <- paste(path.to.file,'RKparms_boo', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)
dev.off()

####  shape constrained Spline fits #########

rng = diff(range(x))/n
ssb.pred.points1 = c(seq(0.1,min(x),by=rng/4),seq(min(x),ssb.max,by=rng))
    
np = length(ssb.pred.points1)
dat = data.frame(stock.size=c(x,ssb.pred.points1),recruit=c(y,rep(1,np)),
                 wt=c(wt.ext,rep(0,np)))
dat$log.recruit=log(dat$recruit)
dat$offset = log(dat$stock.size)

dat.nz = subset(dat,wt==1)

## Nonparametric compensatory mortality SR model;


m = length(dat$stock.size)

nknots=50

gcv.fit = optim(-7,mygcv,upper=5,method="L-BFGS-B")

b.cm <- scam(log.recruit ~ s(stock.size,k=nknots,bs="mpd") + offset(offset),
        family=gaussian(link="identity"),data=dat,optimizer="nlm",
        weights=dat$wt,sp=exp(gcv.fit$par))

nboot=1000
boot.b.cm = exp(bootscam(b.cm,exp(gcv.fit$par),nboot))
boot.b.cmq = apply(boot.b.cm,1,quantile,probs=c(0.025,0.975))

my.df = 3
fit3 = optim(-7,mydf,upper=5,method="L-BFGS-B")
sp3 = exp(fit3$par)
b.cm3 <- scam(log.recruit ~ s(stock.size,k=nknots,bs="mpd",m=2) + offset(offset),
        family=gaussian(link="identity"),data=dat,optimizer="nlm",
        weights=dat$wt,sp=sp3)
        
boot.b.cm3 = exp(bootscam(b.cm3,sp3,nboot))
boot.b.cm3q = apply(boot.b.cm3,1,quantile,probs=c(0.025,0.975))


## scam does not adjust for zero-weighted data when computing df's;
se.scale = sqrt((np+n-sum(b.cm$edf))/(n-sum(b.cm$edf)))
b.cm$fitted.se = se.scale*sqrt(diag(b.cm$X %*% b.cm$Vp.t %*% t(b.cm$X)))
b.cm$pred.prod = exp(b.cm$fitted.value)/dat$stock.size

se.scale = sqrt((np+n-sum(b.cm3$edf))/(n-sum(b.cm3$edf)))
b.cm3$fitted.se = se.scale*sqrt(diag(b.cm3$X %*% b.cm3$Vp.t %*% t(b.cm3$X)))
b.cm3$pred.prod = exp(b.cm3$fitted.value)/dat$stock.size

 

##### Stock-Recruit Plot ##### 

ind = order(dat$stock.size)
  
b.cm.pred = data.frame(stock.size = dat$stock.size[ind],recruit = exp(b.cm$fitted.values)[ind],
  LCI.rec.asy = exp(b.cm$fitted.values - qnorm(0.975)*b.cm$fitted.se)[ind],
  UCI.rec.asy = exp(b.cm$fitted.values + qnorm(0.975)*b.cm$fitted.se)[ind],    
  LCI.rec = boot.b.cmq[1,ind],
  UCI.rec = boot.b.cmq[2,ind],
  prod = b.cm$pred.prod[ind])
b.cm.pred$LCI.prod = b.cm.pred$LCI.rec/b.cm.pred$stock.size  
b.cm.pred$UCI.prod = b.cm.pred$UCI.rec/b.cm.pred$stock.size
   

xlim = range(b.cm.pred$stock.size);
ylim = c(0,max(b.cm.pred$UCI.rec,HS.pred.rec,BH.pred.rec,RK.pred.rec,y))

file.name <- paste(path.to.file,'fit',c,'.png', sep = "") 
png(file=file.name,width=7.5,height=7.5,units='cm',res=600) 
par(mfrow=c(2,1),oma=c(3,1,0.5,0),mar=c(0.5,3,0,1),las=0,cex=0.85)

    
plot(b.cm.pred$stock.size,b.cm.pred$recruit,type='l',lwd=1.5,
  ylab='',xlab='',las=1,xaxt='n',ylim=ylim,xlim=xlim,log='x')  
points(x,y,cex=0.5)  
lines(b.cm.pred$stock.size,b.cm.pred$LCI.rec,type='l',lwd=1,lty=2)
lines(ssb.pred.points,HS.pred.rec,lty=1,col='green',lwd=1.5)
lines(ssb.pred.points,BH.pred.rec,lty=1,col='blue',lwd=1.5)  
lines(ssb.pred.points,RK.pred.rec,lty=1,col='red',lwd=1.5) 
lines(sr.data.pop$s,sr.data.pop$r,lty=1,col='grey',lwd=1.5) 

rect(7800,0,11450,700,col='white',border = NA)
p1 = which.min(abs(b.cm.pred$stock.size-7900))
text(9400,b.cm.pred$recruit[p1]-10,'Scam',cex=0.5)   
p1 = which.min(abs(ssb.pred.points-7900))
text(9400,HS.pred.rec[p1]-10,'HS',cex=0.5)  
text(9400,BH.pred.rec[p1]+10,'BH',cex=0.5)
text(9400,RK.pred.rec[p1],'RK',cex=0.5)   
p1 = which.min(abs(sr.data.pop$s-7900))
text(9400,sr.data.pop$r[p1]+10,'True',cex=0.5) 

box(lty=1)
 
mtext(side=2,line=3,outer=FALSE,r.name,cex=0.85)

ylim = c(0,max(b.cm.pred$UCI.prod,HS.pred.prod,BH.pred.prod,RK.pred.prod,y/x)) 
       
plot(b.cm.pred$stock.size,b.cm.pred$prod,type='l',lwd=1.5,
  ylab='',xlab='',las=1,ylim=ylim,xlim=xlim,log='x')        
points(x,y/x,cex=0.5)    
lines(b.cm.pred$stock.size,b.cm.pred$LCI.prod,type='l',lwd=1,lty=2)
lines(b.cm.pred$stock.size,b.cm.pred$UCI.prod,type='l',lwd=1,lty=2)
lines(ssb.pred.points,HS.pred.prod,lty=1,col='green',lwd=1.5)
lines(ssb.pred.points,BH.pred.prod,lty=1,col='blue',lwd=1.5) 
lines(ssb.pred.points,RK.pred.prod,lty=1,col='red',lwd=1.5)  
sr.data.pop$prod = sr.data.pop$r/sr.data.pop$s   
lines(sr.data.pop$s,sr.data.pop$prod,lty=1,col='grey',lwd=1.5) 


rect(41.2,0,60,3.5,col='white',border = NA)

p1 = which.min(abs(b.cm.pred$stock.size-42))
text(50,b.cm.pred$prod[p1],'Scam',cex=0.5)   
p1 = which.min(abs(ssb.pred.points-42))
text(50,HS.pred.prod[p1],'HS',cex=0.5)  
text(50,BH.pred.prod[p1]-0.15,'BH',cex=0.5)
text(50,RK.pred.prod[p1],'RK',cex=0.5)   
p1 = which.min(abs(sr.data.pop$s-42))
text(50,sr.data.pop$prod[p1]-0.15,'True',cex=0.5) 

box(lty=1)  
mtext(side=2,line=3,outer=FALSE,p.name,cex=0.85)     
mtext(side=1,line=2.2,outer=FALSE,s.name,cex=0.85)


dev.off()

#### compare smoothers with different df's #####

ind = order(dat$stock.size)

b.cm.comp = data.frame(stock.size = dat$stock.size[ind],
  recruit = exp(b.cm$fitted.values)[ind],     
  LCI.rec = boot.b.cmq[1,ind],
  UCI.rec = boot.b.cmq[2,ind],
  recruit3 = exp(b.cm3$fitted.values)[ind],     
  LCI.rec3 = boot.b.cm3q[1,ind],
  UCI.rec3 = boot.b.cm3q[2,ind])

ylim = c(0,max(b.cm.comp$UCI.rec,b.cm.comp$UCI.rec3,y))
  
par(mar=c(4,4,0.5,1),las=0)

 
plot(b.cm.comp$stock.size,b.cm.comp$recruit,type='l',lwd=2,
  ylab='',xlab='',las=1,ylim=ylim)
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec,type='l',lwd=1,lty=2)
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec,type='l',lwd=1,lty=2)
    
lines(b.cm.comp$stock.size,b.cm.comp$recruit3,type='l',lwd=2,lty=1,col='brown')
lines(b.cm.comp$stock.size,b.cm.comp$LCI.rec3,type='l',lwd=1,lty=2,col='brown')
lines(b.cm.comp$stock.size,b.cm.comp$UCI.rec3,type='l',lwd=1,lty=2,col='brown') 

points(x,y)
    
legend("top",lty=1,lwd=2,c('GCV','df=3'),bty='n',col=c('black','brown'))
 
mtext(side=2,line=3,outer=FALSE,r.name)     
mtext(side=1,line=2.2,outer=FALSE,s.name)

file.name <- paste(path.to.file,'sp_compare', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)
dev.off()

##### Residuals Plots #####

resid.matrix = matrix(NA,n,4)
smooth.resid.matrix = resid.matrix
ind = dat$wt==1 
sort.x = sort(x)

resid.matrix[,1] = b.cm$residuals[ind] 
resid.matrix[,2] = HS.raw.res      
resid.matrix[,3] = BH.raw.res      
resid.matrix[,4] = RK.raw.res

resid.scale = sqrt(mean(resid.matrix**2)) 

smooth.resid.matrix[,1] = predict(loess(resid.matrix[,1]~x,span=0.5),sort.x)/resid.scale
smooth.resid.matrix[,2] = predict(loess(resid.matrix[,2]~x,span=0.5),sort.x)/resid.scale
smooth.resid.matrix[,3] = predict(loess(resid.matrix[,3]~x,span=0.5),sort.x)/resid.scale  
smooth.resid.matrix[,4] = predict(loess(resid.matrix[,4]~x,span=0.5),sort.x)/resid.scale

ylim = c(min(smooth.resid.matrix),max(smooth.resid.matrix))
    
par(mar=c(4,4,0.5,1),las=0)
       
plot(sort.x,smooth.resid.matrix[,1],type='l',ylab='',xlab='',las=1,col='black',ylim=ylim,lwd=2)
lines(sort.x,smooth.resid.matrix[,2],col='green',lwd=2)                                       
lines(sort.x,smooth.resid.matrix[,3],col='blue',lwd=2)                                        
lines(sort.x,smooth.resid.matrix[,4],col='red',lwd=2)

abline(h=0,lty=2,lwd=2)
abline(h=-0.5,lty=2,lwd=2,col='grey')
abline(h=0.5,lty=2,lwd=2,col='grey')

legend("top",lty=1,col=c('black','green','blue','red'),lwd=2,
c('Scam','HS','BH','RK'),bty='n')
       
mtext(side=2,line=3,"Standardized residual pattern")
mtext(side=1,line=2.5,s.name)

file.name <- paste(path.to.file,'resid', sep = "")
savePlot(filename = file.name,type = "wmf", device = dev.cur(),restoreConsole = TRUE)

dev.off()

##### Fit table #####

fit.tab = matrix(NA,5,3)                  
fit.tab[1,1] = sum((b.cm$weight*b.cm$residuals)**2)/n 
fit.tab[2,1] = sum((b.cm3$weight*b.cm3$residuals)**2)/n 
fit.tab[3:5,1] = apply((resid.matrix[,2:4])**2,2,mean) 

fit.tab[1,2] = sum(b.cm$edf)
fit.tab[2,2] = sum(b.cm3$edf)
fit.tab[3,2] = 2        
fit.tab[4,2] = 2        
fit.tab[5,2] = 2

fit.tab[1,3] = mygcv(gcv.fit$par) 
fit.tab[2,3] = mygcv(log(sp3))  
fit.tab[3,3] = sum(HS.raw.res**2)/(n*(1 - 2/n)**2)        
fit.tab[4,3] = sum(BH.raw.res**2)/(n*(1 - 2/n)**2)        
fit.tab[5,3] = sum(RK.raw.res**2)/(n*(1 - 2/n)**2)

rownames(fit.tab) = c('Scam GCV','Scam df=3','HS','BH','RK')
colnames(fit.tab) = c('MSE','df','GCV')

ctext = paste('Fit statistics, n = ',n,sep='')
print(xtable(fit.tab,digits=c(7,3,1,3),caption=ctext,
  align=c('c','r','r','r')),type='html',
  file='fit_table.html',caption.placement="top")  
  
