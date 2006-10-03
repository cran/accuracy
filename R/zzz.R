# for R


if (is.R()) {
  .onLoad <- function(lib, pkg) {
    setHook(packageEvent("Zelig", "attach"), ZeligHooks)
    if (sum(search()=="package:Zelig")>0) {
  	   ZeligHooks()
    }
    print(citation("accuracy"))
  return(TRUE)
  }

}

#for S-Plus

if (!is.R()) {
  # hack for R check
  myRequire <- require

  .First.lib <- function(lib.loc, section) {
  
    if (myRequire(pkgutils,quietly=TRUE)){
      print(citation("accuracy", lib.loc = lib.loc))
    } else  {
      
            cat("To cite the accuracy package in publications use:\n\n",
                paste("Altman, M., Gill, J. and M.P. McDonald (2003)",
               "Numerical Issues in Statistical Computing for the Social Scientist.",
               "John Wiley and Sons, New York.",
              "(Software version: R Package, accuracy, version 1.23)",
               "ISBN 0-471-23633-0.\n\n"))
    }
  return(TRUE)
  }

}
