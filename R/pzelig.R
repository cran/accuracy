pzelig<-function(...) {
	sensitivityZelig(...)
}

sensitivityZelig = function (z,
  ptb.R=30,
  ptb.ran.gen=NULL,
        ptb.s=NULL,
        explanatoryOnly=FALSE,
	summarize=FALSE,
	simulate=FALSE,
	simArgs=NULL
  ) {
  if (simulate) {
	summarize=TRUE
  }
  if (!inherits(z,"strata") &&
           (is.null(z$zelig) || is.null(z$call) || is.null(z$model))) {
                warning("z is not a zelig object")
                return(NULL)
  }
  if (inherits(z,"strata")) {
          tmpz=z[[1]]
  } else {
          tmpz=z
  }
  perturbedData=eval(tmpz$call$data,parent.frame())[c(names(tmpz$model),tmpz$call$by)]
  tm=terms(as.formula(as.list(tmpz$call)$formula),data=perturbedData)
        res = names(attr(tm,"factors")[attr(tm,"response")])
  if (is.null(ptb.ran.gen)) {
    if (is.null(ptb.s)) {
         ptb.ran.gen = sapply(as.data.frame(perturbedData), PTBdefaultfn)
    } else {
     ptb.default.size=ptb.s
     ptb.s=NULL
     ptb.ran.gen = sapply(as.data.frame(perturbedData),
                  function(x)PTBdefaultfn(x,ptb.default.size))

    }

    # don't perturb response terms
    if (explanatoryOnly) {
      for (i in which(names(perturbedData)==res)) {
       ptb.ran.gen[[i]] = PTBi
      }
      if (inherits(z,"strata")) {
        ptb.ran.gen[[which(names(perturbedData)==tmpz$call$by)]] = PTBi
      }
    }
  }

  if ( (ptb.R<=0) || (length(ptb.R)>1) ) {
   	stop("ptb.R must be int > 0")
  }

  zfun <- function (data,...) {
       zfcall=tmpz$call
       zfcall$data = data
       eval(zfcall)
  }

  if (simulate) {
	if (is.null(simArgs)) {
		simArgs = list(x=setx(z))
	}
      pz=perturbAndSim(
          perturbedData,
          zfun,
          ptb.R=ptb.R,
          ptb.ran.gen=ptb.ran.gen,
          ptb.s=ptb.s,
	  summarize=summarize,
	  simulate=simulate,
	  simArgs=simArgs
	)
	
  } else {
      pz=perturb(
          perturbedData,
          zfun,
          ptb.R=ptb.R,
          ptb.ran.gen=ptb.ran.gen,
          ptb.s=ptb.s,
	  summarize=summarize
	)
  }
  attr(pz,"origZelig")=z;
  attr(pz,"origZeligEnv")=sys.parent();
  return(pz)
}

perturbAndSim<-function (data, statistic , ptb.R, ptb.ran.gen, ptb.s, summarize, simulate, simArgs) {
        retval = list(ptb.R);
	simlist = list(ptb.R);

        for (i in 1:ptb.R) {
                if (is.null(ptb.s)) {
                   retval[[i]]=
                    perturbHarness(data=data,ran.gen= ptb.ran.gen,statistic=statistic,nameData=TRUE)
                } else  {
                   retval[[i]]=
                    perturbHarness(data=data,ran.gen= ptb.ran.gen,statistic=statistic, ptb.s=ptb.s, nameData=TRUE)
                }
		if (simulate) {
			simArgs$object = retval[[i]]
			simlist[[i]]=do.call("sim",simArgs)
		}
                if (summarize) {
                        retval[[i]] = sensitivitySummaryIteration(retval[[i]])
                }
        }

        attr(retval, "R") = ptb.R
        if (summarize) {
          retval = sensitivitySummaryMerge(retval)
          class(retval)="sensitivity.summary"
        } else {
           attr(retval,"ran.gen")= ptb.ran.gen
           attr(retval, "s") = ptb.s
           attr(retval, "statistic") = statistic
           class(retval)="sensitivity"
        }
	if (simulate) {
  		attr(simlist,"mergedSims")=mergeSims(simlist);
  		class(simlist)="sensitivity.sim"
		retval$sim=simlist
	}
        return(retval)

}

psim<-function(object,x=setx(object),...) {
   if (length(object)==0) {
        warning("zero length perturbation list")
        return(NULL)
   }

  if (inherits(try(sim(object[[1]],setx(object),num=c(5,2)),silent=F ),"try-error"))
  {
        warning("perturbed model does not support sim method")
        return(NULL)
  }

  res = sapply(object,function(tmp){sim(tmp,x,...)}, simplify=F)
  attr(res,"mergedSims")=mergeSims(res);
  class(res)="sensitivity.sim"
  return(res)
}


setx.sensitivity<-function (obj, ...) {
    oz = attr(obj,"origZelig")
    oze = attr(obj,"origZeligEnv")
    eval(call("setx",oz),oze)
}

ZeligHooks<-function (...) {
   if (exists(".simHooked",envir=.GlobalEnv)) {
	return(TRUE)
   }
   origsim=get("sim",envir=as.environment("package:Zelig"))
   sim.replacement=function (object, x, ...) {
    if  (inherits(object,"sensitivity")) {
       psim(object,x,...)
    } else {
      origsim(object,x,...)
    }
   }
   assignInNamespace("sim",sim.replacement,"Zelig")
   unlockBinding("sim",as.environment("package:Zelig"))
   assign("sim",sim.replacement, envir=as.environment("package:Zelig"))
   assign("sim",sim.replacement, envir=.GlobalEnv)
   assign(".simHooked",TRUE,envir=.GlobalEnv)
}


mergeSims<-function(x) {
  tsim=x[[1]]
  if (inherits(x[[1]],"zelig.strata")) {
     simlevs=names(tsim)
  } else {
     simlevs=c(1)
     tsim=list(tsim)
  }
  for (lev in simlevs) {
     tqev=NULL
     tqpr=NULL
     tpar=NULL

     for (tmpl in x) {
        if (!inherits(x[[1]],"zelig.strata")){
           tmp=list(tmpl)
        } else {
           tmp=tmpl
        }
        tpar=rbind(tpar, tmp[[lev]]$par)
        tqev=rbind(tqev, tmp[[lev]]$qi$ev)
        tqpr=rbind(tqpr, tmp[[lev]]$qi$ev)
     }

     dimnames(tpar)=dimnames(tsim[[lev]]$par)
     tsim[[lev]]$par=tpar
     if (!is.null(tsim[[lev]]$qi$ev)) {
         tsim[[lev]]$qi$ev=tqev
     }
     if (!is.null(tsim[[lev]]$qi$pr)) {
         tsim[[lev]]$qi$pr=tqpr
     }
  }
  if (!inherits(x[[1]],"zelig.strata")){
     tsim=tsim[[1]]
  }
  return(tsim)
}

print.sensitivity.sim<-function(x,...) {
  print(summary(x),...)
}

print.sensitivity.sim.summary<-function(x,...) {
  cat("\n\n****",x$reps, " COMBINED perturbation simulations","\n")
  print(x$zsum,...)
}

HTML.sensitivity.sim.summary<-function(x,...) {
  HTML(
    paste("\n\n****",x$reps, " COMBINED perturbation simulations","\n",sep=""), ...)
  HTML(x$zsum,...)
}

plot.sensitivity.sim<-function(x,...) {
  cat("\n\n****",length(x), " COMBINED perturbation simulations","\n")
  return(plot(attr(x,"mergedSims"),...))
}

summary.sensitivity.sim<-function(object,...) {
  ret=list();
  ret$zsum=summary(attr(object,"mergedSims"),...)
  ret$reps=length(object)
  class(ret)="sensitivity.sim.summary"
  ret
}

