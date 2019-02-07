#' @title Plots of the results from eqsim
#'
#' @description XXX
#'
#' @author Einar Hjorleifsson \email{einar.hjorleifsson@@gmail.com}
#'
#' @export
#'
#' @param sim An object returned from the function eqsim_run
#' @param ymax.multiplier A value that acts as a multiplier of the maximum observed
#' variable being plotted. E.g. 1.2 means that for each of the three panels a, b and c
#' the ymax is set to 1.2 of the maximum observed recruitment, spawning stock biomass
#' and yield (catch or landings, depending on user input.
#' @param catch Boolean, if TRUE (default) returns a plot based on catch. If false
#' returns a plot based on landings.

eqsim_plot <- function(sim, ymax.multiplier=1.2, catch=TRUE)
{

  # littleHelper function
  littleHelper <- function(rbp,dat,Flim) {
    ymax <- max(dat[,2]*ymax.multiplier)
    with(rbp,plot(Ftarget,p95,type="l",lty=4,ylab="", xlab="", ylim=c(0,ymax)))
    #title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
    #mtext(text = paste("c)", Variable), cex = 0.75, side = 3, line = 0.5)
    with(rbp,lines(Ftarget,p50, lty = 1))
    with(rbp,lines(Ftarget,p05, lty = 4))
    points(dat[,1],dat[,2],pch=21,cex=0.75,bg=1)
    abline(v=Flim,col="red")
    text(0.98*Flim,0,"F05",srt=90,pos=3,col="red",cex=0.75)
  }




  rby <- sim$rby
  #for (i in c(1,2,4,5)) rby[,i] <- rby[,i]/Scale

  rbp <- sim$rbp
  #for(i in 3:9) rbp[,i] <- rbp[,i]/Scale

  refs <- sim$Refs
  #refs[3:6,] <- refs[3:6,]/Scale

  pProfile <- sim$pProfile

  Flim <- sim$Refs[1,1]
  FCrash5 <- sim$Refs[1,6]
  FCrash50 <- sim$Refs[1,7]
  op <- par(mfrow = c(2, 2), mar = c(2.5, 4, 1.5, 1), oma = c(0, 0, 0, 0),
              cex.axis = 0.75, tcl = 0.25, mgp = c(0, 0.25, 0), las = 1)

  # A: Recruitment plot
  Variable <- "Recruitment"
  i <- rbp$variable == Variable
  littleHelper(rbp[i,],dat=rby[,c("fbar","rec")],Flim)
  title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
  mtext(text = paste(sim$id.sim," a) Recruits"), cex = 0.75, side = 3, line = 0.5)

  # B: SSB plot
  Variable <- "Spawning stock biomass"
  i <- rbp$variable == Variable
  littleHelper(rbp[i,],dat=rby[,c("fbar","ssb")],Flim)
  title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
  mtext(text = paste("b)", Variable), cex = 0.75, side = 3, line = 0.5)

  # C: Yield plot
  if (catch)  {
    # catch versus Fbar
    Variable <- "Catch"
    i <- rbp$variable == Variable
    littleHelper(rbp[i,],dat=rby[,c("fbar","catch")],Flim)
    # Einar added 29.1.2014
    with(rbp[i,],lines(Ftarget,Mean, lty = 1, lwd=4,col="red"))
    title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
    mtext(text = paste("c)", Variable), cex = 0.75, side = 3, line = 0.5)
    #add landings
    j <- rbp$variable == "Landings"
    with(rbp[j,],lines(Ftarget,p50, lty = 2))

    medianFmsy <- sim$Refs[1,4]
    catch_at_medianFmsy <- sim$Refs[3,4]
    meanFmsy <- sim$Refs[1,5]
    catch_at_meanFmsy <- sim$Refs[3,5]

    lines(rep(medianFmsy,2),c(0,catch_at_medianFmsy),col="brown")
    text(0.98*medianFmsy,0,"median Fmsy",srt=90,pos=4,col="brown",cex=0.75)
    lines(rep(meanFmsy,2),c(0,catch_at_meanFmsy),col="brown",lty=2)
    text(0.98*meanFmsy,0,"mean Fmsy",srt=90,pos=4,col="brown",cex=0.75)
    #lines(Fscan, catm, lty=1, col = 2)
    #lines(rep(Fscan[maxcatm], 2), c(0, y.max), lty = 1, col = 5)
  } else {
    # landings versus Fbar
    Variable <- "Landings"
    i <- rbp$variable == Variable
    littleHelper(rbp[i,],dat=rby[,c("fbar","landings")],Flim)
    # Einar added 29.1.2014
    with(rbp[i,],lines(Ftarget,Mean, lty = 1, lwd=4,col="red"))
    title(ylab=Variable, xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
    mtext(text = paste("c)", Variable), cex = 0.75, side = 3, line = 0.5)
    #add catch
    j <- rbp$variable == "Catch"
    with(rbp[j,],lines(Ftarget,p50, lty = 2))

    medianFmsy <- sim$Refs[2,4]
    landings_at_medianFmsy <- sim$Refs[4,4]
    meanFmsy <- sim$Refs[2,5]
    landings_at_meanFmsy <- sim$Refs[4,5]
    lines(rep(medianFmsy,2),c(0,landings_at_medianFmsy),col="brown")
    lines(rep(meanFmsy,2),c(0,landings_at_meanFmsy),col="brown",lty=2)
    #lines(Fscan, lanm, lty=1, col = 2)
    #lines(rep(Fscan[maxlanm], 2), c(0, y.max), lty = 1, col = 5)
  }

  # D: F versus SSB probability

  Variable = "Blim"
  i <- pProfile$variable == Variable
  with(pProfile[i,],plot(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="red"))

  title(ylab="Prob MSY, SSB<Bpa or Blim", xlab="F bar", cex.lab=0.75, line=2.5, cex.main=0.75)
  mtext(text = "d) Prob MSY and Risk to SSB", cex = 0.75, side = 3, line = 0.5)

  Variable = "Bpa"
  i <- pProfile$variable == Variable
  with(pProfile[i,],lines(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="darkgreen"))

  Variable = "pFmsyCatch"
  i <- pProfile$variable == Variable
  with(pProfile[i,],lines(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="cyan"))

  Variable = "pFmsyLandings"
  i <- pProfile$variable == Variable
  with(pProfile[i,],lines(Ftarget,value,type="l",lty=1,ylab="", xlab="",col="brown"))

  text(0.01,0.8, "SSB<Blim", cex = 0.75, col="red",pos=4)
  text(0.01,0.85, "SSB<Bpa", cex = 0.75, col="darkgreen",pos=4)
  text(0.01,0.9, "Prob of cFmsy", cex = 0.75, col="cyan",pos=4)
  text(0.01,0.95, "Prob of lFmsy", cex = 0.75, col="brown",pos=4)

  lines(c(Flim,Flim), c(0.0,0.05), lty = 2, col = "red")
  lines(c(0,Flim), c(0.05,0.05), lty = 2, col = "red")
  text(x = 0.1, y = 0.075, "5%", cex = 0.75, col = "red")


}  # end of eqsim_plot
