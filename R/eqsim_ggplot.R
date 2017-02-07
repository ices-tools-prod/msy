#' @title Plots of the results from eqsim
#'
#' @description XXX
#'
#' @author Einar Hjorleifsson \email{einar.hjorleifsson@@gmail.com}
#'
#' @export
#'
#' @param sim An object returned from the function eqsim_run
#' @param Scale A value, the scaling on the yaxis
#' @param plotit Boolean, if TRUE (default) returns a plot

eqsim_ggplot <- function(sim, Scale=1, plotit=TRUE)
{

  # dummy
  Ftarget <- p05 <- p95 <- p50 <- variable <- value <- year <-
    Mean <- fbar <- rec <- ssb <- catch <- landings <- x <- y <- 0

  rby <- sim$rby
  for (i in c(2,3,5,6)) rby[,i] <- rby[,i]/Scale

  rbp <- sim$rbp
  for(i in 3:9) rbp[,i] <- rbp[,i]/Scale

  refs <- sim$Refs
  refs[3:6,] <- refs[3:6,]/Scale

  sim$Blim <- sim$Blim/Scale
  sim$Bpa <- sim$Bpa/Scale
  pProfile <- sim$pProfile

  i <- rbp$variable %in% "Recruitment"
  plotR <-
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) +
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) +
    ggplot2::geom_line(ggplot2::aes(y=Mean),linetype=2) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::labs(y = "",x="") +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,rec)) +
    ggplot2::coord_cartesian(ylim=c(0,rby$rec * 1.2),xlim=c(0,rby$fbar * 1.2))


  i <- rbp$variable %in% "Spawning stock biomass"
  plotSSB <-
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) +
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) +
    ggplot2::geom_hline(yintercept=sim$Blim,col="red",lwd=1) +
    ggplot2::annotate("text",x=0,y=sim$Blim,label="Blim",col="red",hjust=0,vjust=0) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,ssb)) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::coord_cartesian(ylim=c(0,rby$ssb * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(y = "",x="")

  i <- rbp$variable %in% "Catch"
  plotCatch <-
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) +
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) +
    ggplot2::geom_line(ggplot2::aes(y=Mean),linetype=2) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,catch)) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::coord_cartesian(ylim=c(0,rby$catch * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(y = "",x="")

  i <- rbp$variable %in% "Landings"
  plotLandings <-
    ggplot2::ggplot(rbp[i,],ggplot2::aes(Ftarget)) +
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95),fill="grey90") +
    ggplot2::geom_line(ggplot2::aes(y=p50)) +
    ggplot2::geom_line(ggplot2::aes(y=Mean),linetype=2) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[2,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[2,5],y=0,label="Fmsl",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=rby,ggplot2::aes(fbar,landings)) +
    ggplot2::facet_wrap(~ variable) +
    ggplot2::coord_cartesian(ylim=c(0,rby$landings * 1.2),xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(y = "",x="")


  d2 <- rby[,c("fbar","catch","landings")]
  names(d2) <- c("Ftarget","Catch","Landings")
  d2 <- reshape2::melt(d2,id.vars="Ftarget")
  d2$dummy <- "Yield"
  d2$Ftarget <- as.numeric(d2$Ftarget)

  i <- rbp$variable %in% c("Catch","Landings")
  d1 <- rbp[i,]
  d1$dummy <- "Yield"
  plotYield <-
    ggplot2::ggplot(d1,ggplot2::aes(Ftarget)) +
    ggplot2::theme_bw() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=p05,ymax=p95,fill=variable),alpha=0.15) +
    #geom_line(ggplot2::aes(y=p05,colour=variable)) +
    #geom_line(ggplot2::aes(y=p95,colour=variable)) +
    ggplot2::geom_line(ggplot2::aes(y=p50,colour=variable)) +
    ggplot2::geom_vline(xintercept=refs[1,1],col="red",lwd=1) +
    ggplot2::annotate("text",x=refs[1,1],y=0,label="F05",col="red",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[1,5],col="darkgreen",lwd=1) +
    ggplot2::annotate("text",x=refs[1,5],y=0,label="Fmsy",col="darkgreen",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_vline(xintercept=refs[2,5],col="blue",lwd=1,linetype=2) +
    ggplot2::annotate("text",x=refs[2,5],y=0,label="Fmsl",col="blue",hjust=0,vjust=0,angle=90) +
    ggplot2::geom_point(data=d2,ggplot2::aes(Ftarget,value,colour=variable)) +
    ggplot2::facet_wrap(~ dummy) +
    ggplot2::coord_cartesian(ylim=c(0,max(rby$catch) * 1.2),xlim=c(0,max(rby$fbar) * 1.2)) +
    ggplot2::labs(y = "",x="") +
    ggplot2::theme(legend.position="none") +
    ggplot2::scale_colour_manual(values=c("Catch"="darkgreen","Landings"="blue")) +
    ggplot2::scale_fill_manual(values=c("Catch"="darkgreen","Landings"="blue")) +
    ggplot2::annotate("text",x=max(rby$fbar) * 1.2,y=max(rby$catch) * 1.1,label="Catch",colour="darkgreen",hjust=1) +
    ggplot2::annotate("text",x=max(rby$fbar) * 1.2,y=max(rby$landings) * 1.1,label="Landings",colour="blue",hjust=1)



  pProfile$dummy <- "Probability plot"
  df <- data.frame(x=rep(max(rby$fbar),4),
                   y=c(0.80,0.75,0.70,0.65),
                   variable=c("p(SSB<Blim)","p(SSB<Bpa)","Fmsy","Fmsy - landings"))
  plotProbs <-
    ggplot2::ggplot(pProfile,ggplot2::aes(Ftarget,value,colour=variable)) +
    ggplot2::scale_colour_manual(values=c("pFmsyCatch"="darkgreen",
                                 "pFmsyLandings"="blue",
                                 "Blim"="red",
                                 "Bpa"="orange",
                                 "p(SSB<Blim)"="red",
                                 "p(SSB<Bpa)"="orange",
                                 "Fmsy"="darkgreen",
                                 "Fmsy - landings"="blue")) +
    ggplot2::theme_bw() +
    ggplot2::geom_line(lwd=1) +
    ggplot2::geom_text(data=df,ggplot2::aes(x,y,label=variable,colour=variable)) +
    ggplot2::geom_hline(yintercept=0.05,colour="black") +
    ggplot2::coord_cartesian(xlim=c(0,rby$fbar * 1.2)) +
    ggplot2::labs(x="",y="") +
    ggplot2::facet_wrap(~ dummy) +
    ggplot2::theme(legend.position="none")

  if(plotit) {
    vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
    print(plotSSB + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                          plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp = vplayout(1, 1))
    print(plotR + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                        plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(1,2))
    print(plotYield + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                            plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(2,1))
    print(plotProbs + ggplot2::theme(panel.margin = grid::unit(c(0,0,0,0), "cm"),
                            plot.margin = grid::unit(c(0,0.25,0,0), "cm")),
          vp=vplayout(2,2))
  } else {
    return(list(plotR=plotR,plotSSB=plotSSB,plotCatch=plotCatch,
              plotLandings=plotLandings,plotYield=plotYield,plotProbs=plotProbs))
  }
}
