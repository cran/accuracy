.onLoad <- function(lib, pkg) {
  if((version$major<1) ||(version$minor<6 && version$major==1)) {
    stop("This version for R 1.6 or later")
  }
  if (!frexpTest(FALSE)) {
  	warning("frexp: failed self-test")
  }
  if (!perturbTest()) {
  	warning("perturb: failed self-test")
  } 
  if (!secholTest()) {
  	warning("sechol: failed self-test")
  }
  setHook(packageEvent("perturb", "attach"), PerturbHooks)
  if (sum(search()=="package:perturb")>0) {
	PerturbHooks()
  }
  setHook(packageEvent("Zelig", "attach"), ZeligHooks)
  if (sum(search()=="package:Zelig")>0) {
  	ZeligHooks()
  }
  print(citation("accuracy"))
  return(TRUE)
}


PerturbHooks<-function (...) {
     cat("\n\n****ACCURACY: detaching perturb package because of conflicts with accuracy, but keeping reclassify\n(use help('reclassify',package='perturb') for help)\n\n")
     penv = as.environment("package:perturb")
     rec = get("reclassify",envir=penv)
     assign("reclassify",rec, envir=.GlobalEnv)
     detach("package:perturb")
}

