[![Build Status](https://travis-ci.org/wgmg/msy.svg?branch=master)](https://travis-ci.org/wgmg/msy)

msy
===

The `msy` R package is a collection of methods to estimate equilibrium reference points for fish stocks

Install
-------

To install the latest version of msy:

    # install.packages("devtools")
    library(devtools)
    install_github("flr/FLCore")
    install_github("wgmg/msy")
    library(msy)

Sometimes a more recent developmental version of msy can be obtained from:

    # install.packages("devtools")
    library(devtools)
    install_github("einarhjorleifsson/msy")
    library(msy)

Useage
------

FIT <- eqsr_fit(icesStocks$codNS,
                nsamp = 1000, 
                models = c("ricker", "smooth_hockey", "bevholt"))
eqsr_plot(FIT,n=2e4)
SIM <- eqsim_run(FIT, Fcv=0.25, Fphi=0.30,
                 Blim=70000,Bpa=150000,
                 Fscan = seq(0,1.2,len=40),
                 verbose=FALSE)
SIM$Refs
eqsim_plot(SIM,catch=TRUE)
eqsim_plot_range(SIM)

Contact
-------

You are welcome to:

* submit suggestions and bug-reports at: https://github.com/wgmg/msy/issues
* send a pull request on: https://github.com/wgmg/msy
