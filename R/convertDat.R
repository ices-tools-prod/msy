#' plotMSY
#'
#' @param datfilename XXX
#' @param srtype XXX
#' @return something
#' @author Tim Earl \email{timothy.earl@@cefas.co.uk}
#' @export
convertDat = function(datfilename, srtype)
{
##
#  datfilename = ".\\data\\cod\\srmsy.dat"
#  srtype = 1
##
  datfilename = tolower(datfilename)
  if (!file.exists(datfilename)) stop("file not found:", datfilename)
  dathead1 = scan(datfilename,"",comment.char='#',nlines=2,quiet=TRUE)
  dathead2 = scan(datfilename,0,comment.char='#',nlines=2,skip=2,quiet=TRUE)
  dathead2[5] = srtype
  datbody = scan(datfilename,list(0,0),comment.char='#',skip=5,quiet=TRUE)
  outfilename = paste(dirname(datfilename),dathead1[[2]],sep="/")
  if (!file.exists(outfilename)) stop("file not found:", outfilename)
  outtext = scan(outfilename,"",quiet=TRUE,sep='\n')
  dattext = scan(datfilename,"",quiet=TRUE,sep='\n')

#  write to out.dat and srmsymc.dat
  cat(outtext, sep='\n', file=".\\out.dat")
  cat(dattext[1], sep='\n', file=".\\srmsymc.dat")
  cat(dathead1[1],'out.dat\n', sep=' ', file=".\\srmsymc.dat", append=TRUE)
  cat(dattext[3], sep='\n', file=".\\srmsymc.dat", append=TRUE)
  cat(dathead2,'\n', sep=' ', file=".\\srmsymc.dat", append=TRUE)
  cat(dattext[-(1:4)], sep='\n', file=".\\srmsymc.dat", append=TRUE)

  return(list(Year = dathead2[1]:dathead2[2],Recruits = datbody[[1]],SSB=datbody[[2]],recruitage=dathead2[3],stock=dathead1[1]))
}
