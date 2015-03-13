
#calc.quantiles <- function(d, d.det, variable="variable") {
#  q05 <- q10 <- q16 <- q50 <- q84 <- q90 <- q95 <- value <- NULL
#  x <- plyr::ddply(d,c("variable"),plyr::summarise,q05=quantile(value,0.05),q10=quantile(value,0.10),q16=quantile(value,0.16),q50=quantile(value,0.50),q84=quantile(value,0.84),q90=quantile(value,0.90),q95=quantile(value,0.95))
#  if(!missing(d.det)) x$mean <- d.det$value
#  
#  return(x)
#}


#do.ggplot <- function(x,d) {
#  q05 <- q10 <- q16 <- q50 <- q84 <- q90 <- q95 <- value <- 
#    aes <- variable <- iter <- NULL
#  gg.plot <- 
#    ggplot2::ggplot(x,aes(variable)) +  
#    ggplot2::geom_ribbon(aes(ymin=q05,ymax=q95),fill='grey',alpha=1/2) +
#    ggplot2::geom_ribbon(aes(ymin=q10,ymax=q90),fill='grey',alpha=1/2) +
#    ggplot2::geom_ribbon(aes(ymin=q16,ymax=q84),fill='grey',alpha=1/2) +
#    ggplot2::geom_line(data=d,aes(variable,value,group=iter),alpha = 0.05,col='red') +
#    ggplot2::geom_line(aes(y=q50),col='red',lwd=1) +
#    ggplot2::geom_line(aes(y=mean),col='blue',lwd=1)
#  
#  return(gg.plot)
#}

