# 
# Sensitivity Analysis for Zelig
#
#
#    Copyright (C) 2004-6  Micah Altman
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#  Runs sensitivity() bases sensitivitivity analysi
# The core functions extract the calls and data from Zelig objects and
# then use sensitivity() on the constructed function. 
#
# Zelig sim is extended so that sim() can be called on the results, yielding
# a vector of indepedent simulation runs
#
# A separate summary function merges the resulting simulations into
# one zelig simulation object so that Zelig can summarize the results
#
# Also included: pretty printing functions


######################################################
#       
# pzelig/sensitivityZelig
#       
# This is the core function. It extracts the call and data
# information from the zelig object
#
# It then passes the call to perturbAndSim to do the actual
# call to sensitivity()
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################

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
  if (!inherits(z,"strata") && is.null(z$zelig)) {
                warning("z is not a zelig object")
                return(NULL)
  }
  if (inherits(z,"strata")) {
          tmpz=z[[1]]
  } else {
          tmpz=z
  }

  tm = tmpz$terms
  if ( class ( datanames <- try ( names( eval(tmpz$call$data,parent.frame())) , silent = T )) !=
        "try-error" ) {
	tmpenv = parent.frame()
  } else {
	tmpenv = attr(tm,".Environment")
	datanames = eval(tmpz$call$data, tmpenv )
  }
  depvars  = as.character(intersect (datanames, c(tmpz$call$by, call2names(tm[[2]]) )))
  expvars  = as.character(intersect (datanames, c(tmpz$call$by, call2names(tm[[3]]) )))
  framevars = c(depvars, expvars)
  perturbedData=eval(tmpz$call$data,tmpenv)[framevars]
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
      for (i in which(is.element(names(perturbedData),depvars))) {
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

######################################################
#       
# perturbAndSim
#       
# [Internal Function]
#       
# This is used to manage calling summary() and sim() after
# each sensitivity replication. Doing these in parallel, rather
# than running them at the end saves a lot of memory
#
# Parameters:
#
# - These get passed to sensitivity() as described in the docs
# 	+ data
# 	+ statistic
# 	+ ptb.ran.gen
# 	+ ptb.s
#
# - summarize -- do summaries?
# - simulate -- run sim() (implies summarize)
# - simArgs -- a list() of arguments to pass to sim()
# 
######################################################

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

######################################################
#       
# psim
#       
# A wrapper around Zelig's sim() that detects whether 
# a sensitivity() produced replication vector is used, and if so,
# applys sim() to each of the replications
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


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

######################################################
#       
# setx.sensitivity
#
# Zelig setx method for sensitivity
#       
# Parameters:
#
# See the R documentation file for setx
# 
######################################################


setx.sensitivity<-function (obj, ...) {
    oz = attr(obj,"origZelig")
    oze = attr(obj,"origZeligEnv")
    eval(call("setx",oz),oze)
}

######################################################
#       
# ZeligHooks
#       
#  [Internal]
#
# Hooks sim() so that it calls psim() instead, since 
# Zelig subclassing mechanism is not general enough to handle this case
#       
######################################################

ZeligHooks<-function (...) {
   if (!is.null(zeligOrigSim())) {
        return(TRUE)
   } 
   zeligOrigSim(new=Zelig::sim)
   sim.replacement=function (object, x, ...) {
    if  (inherits(object,"sensitivity")) {
       psim(object,x,...)
    } else {
       getFromNamespace("zeligOrigSim","accuracy")()(object,x,...)
    }
   }
   body(sim,envir=as.environment("package:Zelig"))=body(sim.replacement)
   environment(sim.replacement)=environment(zelig);
   assignInNamespace("sim",sim.replacement,"Zelig")
   assign("sim",sim.replacement, envir=.GlobalEnv)
}

######################################################
#       
# zeligOrigSim
#       
# [Internal]
#
# used by ZeligHooks to patch Zelig sim call
# zelig hook will update this replace this
######################################################

zeligOrigSim <- local({
           osim <- NULL
           function(new)
              if(!missing(new)) osim <<- new else osim 
        })


######################################################
#       
# mergeSims
#       
# [Internal]
#
# Merges the simulations of a vector of sensitivity analyses
# so that Zelig summary() methods can be applied
#       
# Parameters:
#
# x - a sensitivity.sim object
# 
######################################################

mergeSims<- function (x) 
{
    tsim = x[[1]]
    if (inherits(x[[1]], "zelig.strata")) {
        simlevs = names(tsim)
    }
    else {
        simlevs = c(1)
        tsim = list(tsim)
    }
    for (lev in simlevs) {
	tqi=NULL
        tpar = NULL
        for (tmpl in x) {
            if (!inherits(x[[1]], "zelig.strata")) {
                tmp = list(tmpl)
            }
            else {
                tmp = tmpl
            }
            tpar = rblist(tpar, tmp[[lev]]$par)
            tqi = rblist(tqi, tmp[[lev]]$qi)
        }
        dimnames(tpar) = dimnames(tsim[[lev]]$par)
        tsim[[lev]]$par = tpar
        tsim[[lev]]$qi = tqi
    }
    if (!inherits(x[[1]], "zelig.strata")) {
        tsim = tsim[[1]]
    }
    return(tsim)
}

######################################################
#       
# rblist
#       
# [Internal]
#
# Merges the rows of all elements of two lists
#       
# Parameters:
#
# l1,l2 - lists to merge
# 
######################################################


rblist<-function(l1,l2,deparse.level=1) {
	if (is.null(l1) || length(l1)==0) {
		return(l2)
	} 
	if (is.null(l2) || length(l2)==0) {
		return(l1)
	} 

	if(!is.list(l1) || !is.list(l2)) {
		return(rbind(l1,l2,deparse.level=0))
	}
	for (i in 1:length(l1)) {
		l1[[i]]=rbind(l1[[i]],l2[[i]], deparse.level=deparse.level)	
	}
	return(l1)
}
	

######################################################
#       
# call2names
#
# [Internal]
#       
# Recursively conversts a call to a list of names (symbols)
# This is a helper function used to assist in determining
# the original response variables in a formula
#       
# Parameters:
#
# x -- a call (from a formula)
# 
######################################################

call2names<-function(x) {
  if (class(x)=="name") {
     return(x)
  }
  if(class(x)=="call") {
     if (length(x)==1) {
        return(NULL)
     }
     cl = sapply(x[2:length(x)],class)
     return( c(
        as.list(x[which(cl=="name")+1]),
        sapply(x[which(cl=="call")+1],call2names),
         recursive=T))
  }
  return(NULL)
}

#
# Generic methods for display. See the corresponding 
# documentation for print(), summary(), and [R2HTML] HTML()
#
#

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

