#' @title Calculate quantiles
#' 
#' @description xxx
#' 
#' @param d xx
#' @param d.det xx
#' @param variable xx
#' @export

calc.quantiles <- function(d, d.det, variable="variable") {
  q05 <- q10 <- q16 <- q50 <- q84 <- q90 <- q95 <- value <- NULL
  x <- ddply(d,c("variable"),summarise,q05=quantile(value,0.05),q10=quantile(value,0.10),q16=quantile(value,0.16),q50=quantile(value,0.50),q84=quantile(value,0.84),q90=quantile(value,0.90),q95=quantile(value,0.95))
  if(!missing(d.det)) x$mean <- d.det$value
  
  return(x)
}

#' @title Do a ggplot
#' 
#' @description xx
#' 
#' @export
#' 
#' @param x xx
#' @param d xx
#' @export 

do.ggplot <- function(x,d) {
  q05 <- q10 <- q16 <- q50 <- q84 <- q90 <- q95 <- value <- variable <- NULL
  gg.plot <- ggplot(x,aes(variable)) +  
    geom_ribbon(aes(ymin=q05,ymax=q95),fill='grey',alpha=1/2) +
    geom_ribbon(aes(ymin=q10,ymax=q90),fill='grey',alpha=1/2) +
    geom_ribbon(aes(ymin=q16,ymax=q84),fill='grey',alpha=1/2) +
    geom_line(data=d,aes(variable,value,group=iter),alpha = 0.05,col='red') +
    geom_line(aes(y=q50),col='red',lwd=1) +
    geom_line(aes(y=mean),col='blue',lwd=1)
  
  return(gg.plot)
}

#' @title Convert FLStock to rby
#' 
#' @description XXX
#' 
#' @export
#' 
#' @param stk FLStock object
fls2rby <- function(stk) {
  dms <- dims(stk)
  rage <- dms $ min
  if (rage == 0)
    {
    rby <- data.frame(rec = stock.n(stk)[1,drop=TRUE],
                    ssb = ssb(stk)[drop=TRUE],
                    year = with(dms, 1:year + minyear - 1),
                    catch=catch(stk)[drop=TRUE],
                    landings=landings(stk)[drop=TRUE],
                    fbar=fbar(stk)[drop=TRUE]) 
    } else {
      rby <- data.frame(rec = stock.n(stk)[1,-seq(rage),drop=TRUE],
                    ssb = ssb(stk)[1,seq(dms$year - rage),drop=TRUE],
                    year = with(dms, (rage+1):year + minyear - 1),
                    catch=catch(stk)[1,seq(dms$year - rage),drop=TRUE],
                    landings=landings(stk)[1,seq(dms$year - rage),drop=TRUE],
                    fbar = fbar(stk)[1,seq(dms$year - rage),drop=TRUE])
    }
  return(rby)
}