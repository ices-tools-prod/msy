# Code to run plotMSY
#
# see plotymsy.r and documentation for more details of how to use the program
#
# Please report any bugs to timothy.earl@cefas.co.uk to see whether these 
# can be removed in future versions. 
#
# Requires 
#          PlotMSY.r
#          convertSumSen.r
#          srmsymc2.exe
#          srmsymc.exe
# R2.14
# Possibly
#          readlow.r ... if pf pm are to be taken from the data files and not entered in the settings 
#                    
#
#Use setwd in R to set the directory to the file containing the .r files
#
#setwd("c:\\msy\\test\\")
#setwd("\\\\Lowfile3\\stock_assess$\\ICES_Working_Groups\\2012\\WGWIDE\\WHM_assessment\\MSY")
#
#Settings
#For full description of paramters, see "plotMSY instructions" document
#                        
#location of sen and sum file (both same name) ('.\\' indicates relative to working directory
senfile = ".\\data\\nsh\\nsh.sen"
titlename = "NS Herring"
fpa <- NA                    # Can add these to plots if known, otherwise set to NA
flim <- NA                   #
bpa <- 80
blim <- 150

index <- NA                   # Don't use Lowestoft format files to find pf and pm
pfpm <- c(0,0)                # Proportions of f and m occuring before spawning (pf and pm) specified manually)
nits <- 2000                   # Number of MCMC fits calculated and used for confidence interval, typically 100 for
                              #   investigatory analysis, 1,000 for final run to produce plots
nhair <- 100                  # Number of MCMC fits to plot as individual lines
varybiodata  <- TRUE          # Assume uncertainty on stock weights, maturity and natural mortality
srweights <- c(NA,NA,NA)      # Weights for the Ri, BH and HS S-R models. NA indicates automatic values
trimming  <- NA               # proportion of harmonic mean to plot / NA to produce diagnostic plot for trimming

#
#End of Settings                                                                    



source("plotMSY.r")
##Takes some time to run, use commented line below instead for a quicker test
stock = plotMSY(senfile, index, pfpm, srweights, trimming, nits, nhair, varybiodata, titlename, fpa, flim, bpa, blim, silent=TRUE, onlyYPR=FALSE)
#stock = plotMSY(senfile, index, pfpm, c(NA,NA,NA), 10, 10, varybiodata, titlename, fpa, flim, bpa, blim, silent=TRUE, onlyYPR=FALSE)



##1.5	0.122229664	0.156743017	0.721027319
