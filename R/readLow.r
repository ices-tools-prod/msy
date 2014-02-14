
#' plotMSY
#'
#' @param  filename XXX
#' @param expand XXX
#' @param index XXX
#' @param silent XXX
#' @return something
#' @author Tim Earl \email{timothy.earl@@cefas.co.uk}
#' @export
read.Stockfile = function(filename,expand=FALSE,index=NA,silent=FALSE)
{
  #Check that the file exist, can be read and is not empty
  if (file.access(filename,4)==-1) stop("File \'", filename,"\' cannot be read",sep="")
  #Create a te,mporary version with commas replaced by white space
  filecontents = scan(filename,"",sep="\n",quote=NULL,quiet=silent)
  filecontents = gsub(","," ",filecontents)
  cat(filecontents,file = "stockfile.tmp",sep="\n",quote=NULL)
  filename = "stockfile.tmp"
  
  
  fieldcounts = count.fields(filename)
  if (is.null(fieldcounts)) stop("File \'", filename,"\' appears to be empty",sep="")
  if (length(fieldcounts) < 5) stop("Header not complete - expected 5 lines, found",length(fieldcounts),".",sep="")  
  
  #Check that the header fits the file specification
  if (fieldcounts[2]<2) stop("Error on line 2. Expected 2 entries, found ",fieldcounts[2],".",sep="")
  if (fieldcounts[3]<2) stop("Error on line 3. Expected 2 entries, found ",fieldcounts[3],".",sep="")
  if (fieldcounts[4]<2) stop("Error on line 4. Expected 2 entries, found ",fieldcounts[4],".",sep="")
  if (fieldcounts[5]<1) stop("Error on line 5. Expected 1 entry, found ",  fieldcounts[5],".",sep="")  
  
  #read in header information
  #  Line 2
  linein = scan(filename, integer(0), nlines=1, skip=1,quiet=TRUE)
  if (all(linein[1] != 1,!silent)) warning("Expected sex value is 1, found ", linein[1], ".", sep="") 
  if (!is.na(index)) if (!index==linein[2]) stop("Index does not match expected value. Found ",linein[2]," expected ",index,".",sep="")
  if (!(linein[2] %in% 1:10)) stop("Index not valid. Found ",linein[2]," expected range 1-10.",sep="")
  sex = linein[1]
  index = linein[2]    

  #  Line 3
  linein = scan(filename, nlines=1, skip=2,quiet=TRUE)
  if(linein[1]>linein[2]) stop("First year after last year.")
  years = linein[1]:linein[2]  
  
  #  Line 4  
  linein = scan(filename, nlines=1, skip=3,quiet=TRUE)
  if(linein[1]>linein[2]) stop("First age greater than last age.")  
  ages = linein[1]:linein[2]  
  
  #  Line 5
  linein = scan(filename, nlines=1, skip=4,quiet=TRUE)
  if (!(linein[1] %in% c(1,2,3,5))) stop("Unexpected DFI value found (", linein[1], "). Expecting 1, 2, 3 or 5.",sep="")
  DFI = linein[1]  
  
  #Check enough data
  #  First enough rows
  if (((DFI==2)||(DFI==3))&&(length(fieldcounts) < 6)) stop("Error in file: ", filename,"Data row not found")
  if (((DFI==1)||(DFI==5))&&(length(fieldcounts) < 5+length(years))) stop("Not enough data rows found. Expected ",length(years)," found ",length(fieldcounts)-5)

  #  then enough entries in each row
  if ((DFI==1) && any(fieldcounts[6:(5+length(years))]<length(ages))) stop("Rows not long enough. Expected ",length(ages)," values in each row.")  
  if ((DFI==2) && any(fieldcounts[6]                  <length(ages))) stop("Rows not long enough. Expected ",length(ages)," values in each row.") 
  if ((DFI==3) && any(fieldcounts[6]                  <1))            stop("Rows not long enough. Expected 1 value in in each row.") 
  if ((DFI==5) && any(fieldcounts[6:(5+length(years))]<1))            stop("Rows not long enough. Expected 1 value in each row.") 

  #There's enough data!!! read it in
  if ((DFI==2)||(DFI==3)) rows = 1
  if ((DFI==1)||(DFI==5)) rows = length(years)
  if ((DFI==3)||(DFI==5)) cols = 1
  if ((DFI==1)||(DFI==2)) cols = length(ages)
  dat = matrix(unlist(scan(filename, as.list(rep(0,cols)), nlines=rows, skip=5,quiet=TRUE,flush=TRUE)),rows,cols)
  
  
  if (expand)
  {
   if(DFI==2) dat = matrix(dat,length(years),length(ages),byrow=TRUE)
   if(DFI==3) dat = matrix(dat,length(years),length(ages))
   if(DFI==5) dat = matrix(dat,length(years),length(ages))
  }

  if ((expand)||(DFI==1)||(DFI==5)) rownames(dat) = years
  if ((expand)||(DFI==1)||(DFI==2)) colnames(dat) = ages
  
  if (!expand) dat = drop(dat)
  
  unlink(filename)
  return(dat)
}       


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
  if (is.na(filename)) filename = choose.files("*.idx", "Choose index file",multi=FALSE)
  if (file.access(filename,4)==-1) stop("File \'", filename,"\' cannot be read",sep="")
  fieldcounts = count.fields(filename)
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


extractSenSum <- function(in.stock, y, optLand=TRUE)
{           
  #require(FLCore)                        
  d.flag <- (sum(discards(in.stock)) >0)
  ret <- vector("list",9)
  names(ret) <- c("Hist","sH","sD","wH","wD","wS","M","MT","dims")
  ret$Hist <- cbind(year=an(dimnames(in.stock@stock.n)$year), R= c(stock.n(in.stock)[1,]), ssb=c(ssb(in.stock)),
                    F=c(fbar(in.stock)), yield=c(catch(in.stock)))
  years <- rev(rev(dimnames(in.stock@stock.n)$year)[1:y])
  ages <- dimnames(in.stock@stock.n)$age
  my_sum <- function(x){r=cbind(mean=apply(x[,years,],1,mean),CV=apply(x[,years,],1,sd)/apply(x[,years,],1,mean));r[is.na(r)]<- 0; r}
  
  
  if (d.flag & optLand) {
    sC <- harvest(in.stock)[,years]
    sH <- harvest(in.stock)[,years] * (landings.n(in.stock)/catch.n(in.stock))[,years]
    sD <- harvest(in.stock)[,years] * (discards.n(in.stock)/catch.n(in.stock))[,years]    
    for (yy in years) sH[,yy] <- sH[,yy]/mean(sC[range(in.stock)["minfbar"]:range(in.stock)["maxfbar"],yy])
    for (yy in years) sD[,yy] <- sD[,yy]/mean(sC[range(in.stock)["minfbar"]:range(in.stock)["maxfbar"],yy])    
    ret$sH <- my_sum(sH)  
    ret$sD <- my_sum(sD)      
  } else {
    sH <- harvest(in.stock)[,years]
    for (yy in years) sH[,yy] <- sH[,yy]/mean(sH[range(in.stock)["minfbar"]:range(in.stock)["maxfbar"],yy])
    ret$sH <- my_sum(sH)
    ret$sD <- ret$sH*0  
  }
  

  ret$wH <- my_sum(landings.wt(in.stock))
  ret$wD <- my_sum(discards.wt(in.stock))
  ret$wD[is.na(ret$wD)] <- 0
  ret$wS <- my_sum(stock.wt(in.stock))  
  ret$M <- my_sum(m(in.stock)) 
  ret$MT <- my_sum(mat(in.stock)) 
  ret$MT[is.na(ret$MT)] <- 0 
  ret$dims <- list(age=ages,year=years,fbarage=range(in.stock)["minfbar"]:range(in.stock)["maxfbar"])
  return (ret)
}

writeSenSum <- function(stock.list,filename,stockname=filename)
{
  senfilename <- paste0(filename,".sen")
  sumfilename <- paste0(filename,".sum")
  cat(stockname, '\n', file=senfilename)
  cat(range(an(stock.list$dim$age)),max(an(stock.list$dim$year)),0,'\n', append=TRUE,file=senfilename)
  cat('1 1 0\n', append=TRUE,file=senfilename)
  cat(paste(paste0("N",stock.list$dims$age),0,0,collapse='\n'),'\n', append=TRUE,file=senfilename)
  for (x in names(stock.list)[2:8]) cat(paste(paste0(x,stock.list$dims$age),stock.list[[x]][,"mean"],stock.list[[x]][,"CV"],collapse='\n'),'\n', append=TRUE,file=senfilename)
  
  cat('Slightly stunted sum file\n12\n0\n0\n', file=sumfilename)
  cat(range(stock.list[[1]][,"year"]),'\n', file=sumfilename, append=TRUE)  
  cat("Recruitment\n", file=sumfilename, append=TRUE)
  cat("Dummy line 7\n", file=sumfilename, append=TRUE)
  cat("SSB\n", file=sumfilename, append=TRUE)
  for (i in 9:13) cat("dummy line", i,'\n', file=sumfilename, append=TRUE)
  cat("Yield\n", file=sumfilename, append=TRUE)
  for (i in 15:20) cat("dummy line", i,'\n', file=sumfilename, append=TRUE)
  cat(range(an(stock.list[["dims"]]$fbarage)),'\n', file=sumfilename, append=TRUE)  
  for (i in 22:27) cat("dummy line", i,'\n', file=sumfilename, append=TRUE)  
  cat(paste(stock.list$Hist[,"year"], stock.list$Hist[,"R"], stock.list$Hist[,"ssb"],0,0,stock.list$Hist[,"yield"],0,0,0,stock.list$Hist[,"F"],0,0,collapse='\n'), file=sumfilename, append=TRUE)  
  
}
