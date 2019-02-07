#' plotMSY
#'
#' @param filename XXX
#' @param index XXX
#' @param expand XXX
#' @param silent XXX
#' @return something
#' @author Tim Earl \email{timothy.earl@@cefas.co.uk}
#' @export
read.Lowestoft = function(filename=NA,index=1:10,expand=TRUE,silent=FALSE)
{
  if (all(!expand,!silent)) cat("Expand is FALSE, this suppresses warnings about inconsistent array sizes\n")
  retval = list()
  if (is.na(filename)) filename = utils::choose.files("*.idx", "Choose index file",multi=FALSE)
  if (file.access(filename,4)==-1) stop("File \'", filename,"\' cannot be read",sep="")
  fieldcounts = utils::count.fields(filename)
  if (length(fieldcounts) < max(index)+2) stop("Not enough rows in index file - expected 13 lines")
  path = paste(rev(rev(strsplit(filename, "\\\\")[[1]])[-1]),collapse="\\")
  datfiles = scan(filename, "", nlines=11, skip=2, quiet=TRUE, flush=TRUE)
  if (path != "") datfiles = paste(path, datfiles, sep="\\")

  #Check file exists
  #read in files to datfiles[]
  filetitles = c("Landings","CAA","Catch WAA","Stock WAA","Natural mortality","Maturity","pF","pM","FOldest","F at Age","Tuning")

  if (all(unique(index)!=index))
  {
    warning("Dropping duplicate inputs to index")
    index = unique(index)
  }
  yearrange=NA
  agerange=NA
  rangename=""

  for (i in index)
  {
    if (i != 11) retval = c(retval, list(read.Stockfile(datfiles[i],expand,i,silent)))
    if (i == 11) retval = c(retval, list("Not implemented yet"))
    if (all(expand,i!=11) )  #remove i!=11 condition when tuning dat has been implemented
    {
      if (is.na(yearrange[1]))
      {
        yearrange = range(as.numeric(rownames(retval[[length(retval)]])))
        agerange = range(as.numeric(colnames(retval[[length(retval)]])))
        rangename = filetitles[i]
      } else {
        if (any(yearrange != range(as.numeric(rownames(retval[[length(retval)]])))))
           warning("\nYear ranges don't match up:\n", yearrange[1],"-",yearrange[2], " in ", rangename, "\n", paste(range(as.numeric(rownames(retval[[length(retval)]]))),collapse="-"), " in ",filetitles[i])
        if (any(agerange != range(as.numeric(colnames(retval[[length(retval)]])))))
           warning("\nAge ranges don't match up:\n", agerange[1],"-",agerange[2], " in ", rangename, "\n", paste(range(as.numeric(colnames(retval[[length(retval)]]))),collapse="-"), " in ",filetitles[i])
      }
    }

#    check years and ages match up
  }
  names(retval) = filetitles[index]
  return(retval)

}
