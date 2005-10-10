pzelig = function (z,
  ptb.R=50,
  ptb.ran.gen=NULL,
	ptb.s=NULL,
	explanatoryOnly=FALSE
  ) {

	if (is.null(z)) {
		warning("z is null")
		return(NULL)
	}
	if (is.null(z$zelig) || is.null(z$call) || is.null(z$model)) {
		warning("z is not a zelig object")
		return(NULL)
	}
  tm=terms(as.formula(as.list(z$call)$formula),data=z$model)
	res = names(attr(tm,"factors")[attr(tm,"response")])
	if (is.null(ptb.ran.gen)) {
    if (is.null(ptb.s)) {
         ptb.ran.gen = sapply(as.data.frame(z$model), PTBdefaultfn)
    } else {
     ptb.default.size=ptb.s
     ptb.s=NULL
     ptb.ran.gen = sapply(as.data.frame(z$model),
                  function(x)PTBdefaultfn(x,ptb.default.size))

    }

    # don't perturb responser terms
    if (explanatoryOnly) {
      for (i in which(names(z$model)==res)) {
       ptb.ran.gen[[i]] = PTBi
      } 
    }
  }
  
  zfun <- function (data,...) {
       zfcall = as.list(z$call)
       zfcall$data = data
       eval(as.call(zfcall))
  }
  
	pz=perturb(
          z$model,
          zfun,
          ptb.R=ptb.R,
          ptb.ran.gen=ptb.ran.gen,
          ptb.s=ptb.s
          )
  attr(pz,"origZelig")=z;
  return(pz)
}



setx.perturb<-function(obj,...) {
     setx(attr(obj,"origZelig"),...)
}


ZeligHooks<-function (...) {

   origsim=get("sim",envir=as.environment("package:Zelig"))
   sim.replacement=function (object, x, ...) {
    if  (inherits(object,"perturb")) {
       psim(object,x,...)
    } else {
      origsim(object,x,...)
    }
   }
   assignInNamespace("sim",sim.replacement,"Zelig")
   unlockBinding("sim",as.environment("package:Zelig"))
   assign("sim",sim.replacement, envir=as.environment("package:Zelig"))
   assign("sim",sim.replacement, envir=.GlobalEnv)
}


mergeSims<-function(x) {

  tsim=x[[1]]
  tqev=NULL
  tqpr=NULL
  tpar=NULL
  for (tmp in x) {
      tpar=rbind(tpar, tmp$par)
      tqev=rbind(tqev, tmp$qi$ev)
      tqpr=rbind(tqpr, tmp$qi$ev)
  }
  dimnames(tpar)=dimnames(tsim$par)
  tsim$par=tpar
  if (!is.null(tsim$qi$ev)) {
     tsim$qi$ev=tqev
  }
  if (!is.null(tsim$qi$pr)) {
         tsim$qi$pr=tqpr
  }
  return(tsim)
}

psim<-function(object,x=setx(object),...) {
   if (length(object)==0) {
	warning("zero length perturbation list")
	return(NULL)
   }

  if ( sum( names( object[[1]] )=="zelig" )==0
       && sum( methods( class=class( object[[1]] ) ) =="sim" ) ==0  )
  {
       	warning("perturbed model does not support sim method")
       	return(NULL)
  }

  res = sapply(object,function(tmp){sim(tmp,x,...)}, simplify=F)
  attr(res,"mergedSims")=mergeSims(res);
  class(res)="perturb.sim"
  return(res)
}

print.perturb.sim<-function(x,...) {
  cat("\n\n****",length(x), " COMBINED perturbation simulations","\n")
  return(print(attr(x,"mergedSims"),...))
}

plot.perturb.sim<-function(x,...) {
  cat("\n\n****",length(x), " COMBINED perturbation simulations","\n")
  return(plot(attr(x,"mergedSims"),...))
}

summary.perturb.sim<-function(object,...) {
  cat("\n\n****",length(object), " COMBINED perturbation simulations","\n")
  return(summary(attr(object,"mergedSims"),...))
}







