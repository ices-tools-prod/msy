#' @title srmsymc boxplot
#' 
#' @export
#' 
#' @param x XXX
#' @param ... XXX
srmsymc_plotbox <- function(x, ...) 
  { 
  bp = boxplot(x,plot=FALSE,...,na.rm=TRUE)
  bp$stats = as.matrix(quantile(x,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE))
  bp$out = c(x[x<bp$stats[1]],x[x>bp$stats[5]])
  bp$out <- bp$out[!is.na(bp$out)]
  bp$group = rep(1,length(bp$out))
  bxp(bp,...)
}


#' @title Plot type 1
#' 
#' @export
#' 
#' @param p XXX
#' @param d XXX
#' @param xlim XXX
#' @param ylim XXX
plot_type1 <- function(p,d,xlim,ylim) {
  if(missing(xlim)) {
    xlim=c(0,max(p$variable))
  }
    if(missing(ylim)) {
    ylim=c(0,max(p$q95,d[,3]))
  }
  plot(xlim,ylim,type="n",axes=TRUE,xlab="",ylab="")
  lines(p$variable,p$q05,col="red",lty=3)
  lines(p$variable,p$q10,col="red",lty=2)
  lines(p$variable,p$q50,col="red",lty=1)
  lines(p$variable,p$q90,col="red",lty=2)
  lines(p$variable,p$q95,col="red",lty=3)
  lines(p$variable,p$mean,col="blue")
  points(d[,2],d[,3],type="b",cex=.7)
  # last year value
  i <- nrow(d)
  text(d[i,2], d[i,3], d[i,1],pos=4 )
  points(d[i,2], d[i,3],cex=.7,pch=19 )
}

#' @title Plot type 2
#' 
#' @export
#' 
#' @param sto XXX
#' @param det XXX
#' @param n XX
#' @param xlim XXX
#' @param ylim XXX
plot_type2 <- function(sto,det,n,xlim,ylim) {
  if(missing(xlim)) {
    xlim=c(0,max(sto$variable))
  }
  if(missing(ylim)) {
    ylim=c(0,max(sto$value))
  }
  plot(xlim,ylim,type="n",axes=TRUE,xlab="",ylab="")
  for (i in 1:n) lines(sto$variable[sto$iter == i],sto$value[sto$iter == i],col="red")
  lines(det$variable,det$value,col="blue")
}



#' @title XXX
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param path list
#' @param rby data.frame with stock summary data
#' @param n XXX
#' @param hair XXX
srmsymc_plotYield <- function(path="ricker",rby,n=100,hair=FALSE) {
  dat <- srmsymc_read_yield(path=path,doPlot=FALSE) 
  if(!hair) {
    plot_type1(dat$quantiles,rby[,c("year","fbar","catch")])
  } else {
    plot_type2(dat$stochastic,dat$deterministic,rby[,c("year","fbar","catch")],n=n)
  }
}

#' @title XXX
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param path list
#' @param rby data.frame with stock summary data
#' @param n XXX
#' @param hair XXX
srmsymc_plotSSB <- function(path="ricker",rby,n=100,hair=FALSE) {
  dat <- srmsymc_read_ssb(path=path,doPlot=FALSE) 
  if(!hair) {
    plot_type1(dat$quantiles,rby[,c("year","fbar","ssb")])
  } else {
    plot_type2(dat$stochastic,dat$deterministic,rby[,c("year","fbar","ssb")],n=100)
  }
}

#' @title XXX
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param path list
#' @param rby data.frame with stock summary data
#' @param n XXX
#' @param hair XXX
srmsymc_plotSSBR <- function(path="ricker",rby,n=30,hair=FALSE) {
  par <- srmsymc_read_par(path=path,longformat=FALSE)
  sto <- par$stochastic[,c("alpha","beta")]
  names(sto) <- c("a","b")
  n <- nrow(sto)
  det <- par$deterministic[,c("alpha","beta")]
  names(det) <- c("a","b")
  x <- c(1:(max(rby$ssb) * 1.2))
  eq <- switch(par$srno,ricker,bevholt,segreg)
  
  if(!hair) {
  dat <- data.frame(iter = rep(1:n,each=length(x)),
                    ssb=rep(x,n),
                    rec=NA)
  for (i in 1:n) dat$rec[dat$iter == i] <- exp(eq(sto[i,],x))
  p <- ddply(dat,c("ssb"),summarise,
                       p05=quantile(rec,0.05),
                       p10=quantile(rec,0.10),
                       p50=quantile(rec,0.50),
                       p90=quantile(rec,0.90),
                       p95=quantile(rec,0.95),
                       mean=mean(rec))
  plot(c(0,max(x)),c(0,max(rby$rec)),type="n",xlab="",ylab="")
  lines(p$ssb,p$p05,col="red",lty=3)
  lines(p$ssb,p$p10,col="red",lty=2)
  lines(p$ssb,p$p50,col="red",lty=1)
  lines(p$ssb,p$p90,col="red",lty=2)
  lines(p$ssb,p$p95,col="red",lty=3)
  lines(p$ssb,p$mean,col="blue")
  points(rby$ssb,rby$rec,cex=.7)
  } else {
  plot(c(0,max(x)),c(0,max(rby$rec)),type="n",xlab="",ylab="")
  for (i in 1:n) curve(exp(eq(sto[i,],x)), add=TRUE,col="red")
  curve(exp(eq(det,x)), add=TRUE, col="blue")
  }
}

#' @title XXX
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param path list
#' @param Fbar XXX
#' @param Fpa XXX
#' @param Flim XXX
srmsymc_plotPar <- function(path="ricker",Fbar,Fpa,Flim) {
  if(missing(Fbar)) Fbar <- NA
  if(missing(Fpa)) Fpa <- NA
  if(missing(Flim)) Flim <- NA
  par <- srmsymc_read_par(path=path,longformat=FALSE)
  cn <- c("f40","f35","fmsy","f01","fmax","fcrash")
  x <- par$stochastic[,cn]
  
  plot(0,0,axes=FALSE,type='n',xlim=xlim,ylim=c(-2.5,6.5),xlab=" ",ylab="")
  axis(2,(-2):6,c("fbartext","Fpa","Flim","Fcra","Fmax","F01","Fmsy","F35","F40"),cex.axis=1.3,las=1)
  axis(1)
  box()
  if(!is.na(Fbar)) srmsymc_plotbox(Fbar, add=TRUE, at=-2,horizontal=TRUE,axes=FALSE)
  if(!is.na(Fpa)) srmsymc_plotbox(Fpa, add=TRUE, at=-1,horizontal=TRUE,axes=FALSE)
  if(!is.na(Flim)) srmsymc_plotbox(Flim, add=TRUE, at=0,horizontal=TRUE,axes=FALSE)
  #if (!noredlines[[srtype]])
  #{
    srmsymc_plotbox(x$fcrash, add=TRUE, at=1,horizontal=TRUE,axes=FALSE)
    srmsymc_plotbox(x$fmax,   add=TRUE, at=2,horizontal=TRUE,axes=FALSE)
    srmsymc_plotbox(x$f01,    add=TRUE, at=3,horizontal=TRUE,axes=FALSE)
    srmsymc_plotbox(x$fmsy,   add=TRUE, at=4,horizontal=TRUE,axes=FALSE)
    srmsymc_plotbox(x$f35,    add=TRUE, at=5,horizontal=TRUE,axes=FALSE)
    srmsymc_plotbox(x$f40,    add=TRUE, at=6,horizontal=TRUE,axes=FALSE)
  #}
}

#' @title XXX
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param path list
#' @param rby data.frame with stock summary data
#' @param n XXX
#' @param Fbar XXX
#' @param Fpa XXX
#' @param Flim XXX
srmsymc_plotComposite <- function(path="ricker",rby,n=40,Fbar,Fpa,Flim) {
  #png(paste(outputfolder, stockname, "_Yield_",srsn[srtype],".png",sep=""),height=11.5,width=9,units="in",res=144)
  png("tmp.png",height=11.5,width=9,units="in",res=144)
  layout(t(matrix(c(1,1:7),2)), widths=c(5,5), heights=c(1,3,3.5,4))
  par(mai=c(0,0,0.8,0))
  plot.new()
  title(paste("stockname","srname[srtype]"))
  #srnhair = min(nhair,length(simdata[[srtype]]$fcrash))
  par(mai=c(0.5,1,0,0.5))
  
  srmsymc_plotPar(path=path)
  srmsymc_plotPar(path=path)
  srmsymc_plotYield(path=path,rby=rby)
  srmsymc_plotYield(path=path,rby=rby,n=n,hair=TRUE)
  srmsymc_plotSSB(path=path,rby=rby)
  srmsymc_plotSSB(path=path,rby=rby,n=n,hair=TRUE)
  dev.off()
}

#' @title XXX
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param rby XXX
#' @param n XXX

srmsymc_plotComposite2 <- function(rby,n=40) {
  #png(paste(outputfolder, stockname, "_SRR.png",sep=""),height=11.5,width=9,units="in",res=144)
  png("tmp2.png",height=11.5,width=9,units="in",res=144)
  layout(t(matrix(1:6,2)))
  srmsymc_plotSSBR("ricker",rby=rby)
  srmsymc_plotSSBR("ricker",rby=rby,hair=TRUE)
  srmsymc_plotSSBR("bevholt",rby=rby)
  srmsymc_plotSSBR("bevholt",rby=rby,hair=TRUE)
  srmsymc_plotSSBR("segreg",rby=rby)
  srmsymc_plotSSBR("segreg",rby=rby,hair=TRUE)

  dev.off()
}


#' @title XXX
#' 
#' @description XXX
#' 
#' @export 
#' 
#' @param path XXX
#' @param variable XXX
#' @param srweights XXX
#' @param stockname XXX

srmsymc_plotDistribution <- function(path,variable="fmsy",srweights=NA,stockname="The stock") {
  
  dat <- srmsymc_combineParameterFiles(srweights)
  
  Title <- paste(stockname,"- combined",variable,"distribution:",
                 ifelse(any(is.na(srweights)),"equally weighted","automatically weighted"))
                 
  
  f.breaks <- hist(c(dat$fmsy,dat$fcrash),50, plot=FALSE)$breaks
  mx = max(f.breaks)
  my <-  max(hist(c(dat$fmsy,dat$fcrash),50, plot=FALSE)$counts)
  
  y <- dat[,c(variable,"srno")]
  names(y) <- c("variable","srno")
  hist(y$variable, #rbind(simdata[[1]],simdata[[2]],simdata[[3]])$fmsy,
       f.breaks, 
       col=grey(0.75), 
       ylim = c(0,1.10*my),
       xlab="Fmsy", 
       main = Title)
  hist(y$variable[y$srno %in% c(1,2)],
       f.breaks,
       add =TRUE, 
       col=grey(0.5))
  hist(y$variable[y$srno %in% c(1)],
       f.breaks,
       add =TRUE, 
       col=grey(0.25))
  qs <- cbind(quantile(y$variable,c(0.05, 0.25, 0.50, 0.75, 0.95),na.rm=TRUE))
  
  lines(c(qs[1,1],qs[1,1]),c(0,my),col='red',lty=3)      #Confidence intervals
  lines(c(qs[2,1],qs[2,1]),c(0,my*1.05),col='red',lty=2)
  lines(c(qs[3,1],qs[3,1]),c(0,my),col='red')
  lines(c(qs[4,1],qs[4,1]),c(0,my*1.05),col='red',lty=2)
  lines(c(qs[5,1],qs[5,1]),c(0,my),col='red',lty=3)
  text(qs[1,1],my,"5%",pos=3,col="red")
  text(qs[2,1],my*1.05,"25%",pos=3,col="red")
  text(qs[3,1],my,"50%",pos=3,col="red")
  text(qs[4,1],my*1.05,"75%",pos=3,col="red")
  text(qs[5,1],my,"95%",pos=3,col="red")
  
  legend(0.75*mx,my,c("Ricker","Beverton-Holt","Hockeystick"), fill=grey(0.25*1:3))
}
