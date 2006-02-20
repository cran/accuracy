# perturb.R
#
# Performs sensitivity analysis of non-linear and linear models
# through data pertubations
#
# Part of the Accuracy package. Available from www.r-project.org and
# www.hmdc.harvard.edu/numerical_issues/
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
# This is a set of functions that run multiple replications of
# the analysis, each with a different data perturbation
#
# The resulting analysis objects are stored in a vector, and summary()
# iterates over this vector, extracting the coefficients of the analysis
# into a results matrix


######################################################
#       
# perturb/sensitivity
#
# Replicates an analysis with randomly perturbed data
#
# Calls perturbHarness on each replicate
#
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################

"perturb"<-function(...) {
	sensitivity(...)
}

sensitivity<-function(data,statistic,..., ptb.R=50,
  ptb.ran.gen=NULL,
	ptb.s=NULL,
	summarize=FALSE
	) {

	ptb.R=as.integer(trunc(ptb.R))
	
	if (is.null(ptb.ran.gen)) {
    		if (is.null(ptb.s)) {
     			ptb.ran.gen = sapply(as.data.frame(data), PTBdefaultfn)
    		} else {
     			ptb.ran.gen = sapply(as.data.frame(data),
                  		function(x)PTBdefaultfn(x,ptb.s))
     			ptb.s=NULL
    		}
   	}

	if ( (ptb.R<=0) || (length(ptb.R)>1) ) {
		stop("ptb.R must be int > 0")
	}
	retval = list(ptb.R);
	for (i in 1:ptb.R) {
		if (is.null(ptb.s)) {
		   retval[[i]]=
		    perturbHarness(data=data,ran.gen= ptb.ran.gen,statistic=statistic,...)
		} else  {
		   retval[[i]]=
		    perturbHarness(data=data,ran.gen= ptb.ran.gen,statistic=statistic,...,
			ptb.s=ptb.s)
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
	return(retval)
}

######################################################
#          
# perturbHarness
#          
# [Internal]
#
# Performs one perturbation replication 
#       
# Parameters:
#
# data, ran.gen,statistic, ptb.s -- as in sensitivity() above
# nameData -- should data be sent to anlysis as a name  (by reference)
#		or by value. The former saves memory, but
#		some analysis methods do not handle it well.
#       
######################################################

perturbHarness<-function(data,ran.gen,statistic,..., ptb.s=NULL,nameData=F) {
	ndata = as.data.frame(data)
	if (length(ran.gen)==1) {
		ind=which(sapply(ndata,is.numeric))
		ran.gen=replicate(ncol(ndata),ran.gen)
		if (!is.null(ptb.s)) {
			ptb.s = replicate(ncol(ndata),ptb.s)
		}
	} else if (length(ran.gen)!=ncol(ndata)) {
		stop("ran.gen must be a single function, or a vector of functions of length ncol(df)")
	} else {
		ind = which(sapply(ran.gen,is.function))
	}

	for ( i in ind ) {
		if (is.null(ptb.s[[i]])) {
			ndata[[i]] = ran.gen[[i]](ndata[[i]])
		} else {
			ndata[[i]] = ran.gen[[i]](ndata[[i]],size=ptb.s[i])
		}
	}
	
        res=NULL
        if (class(statistic)=="function")    {
           ne = new.env(parent=environment(statistic))
           environment(statistic)=ne
           assign("ndata",ndata,envir=ne)
	   if (nameData) {
           	res = statistic(data=as.name("ndata"),...)
	   } else {
           	res = statistic(data=ndata,...)
	   }
        }  else  {
	   if (nameData) {
           	res= do.call(statistic,list(...,data=as.name("ndata")))
	   } else {
           	res= do.call(statistic,list(...,data=ndata))
	   }
        }

        return(res)
	
}

######################################################
#       
# perturbTest 
#       
# [Internal Function]
#
# self test of perturb, sanity checks
#       
# Parameters:
#
# silent - print debugging output
# 
######################################################


perturbTest<-function(silent=TRUE) {
	status=TRUE
	if(version$major<2) {
		data(longley, package="base")
	} else {
		data(longley, package="datasets")
	}

	for (i in c(PTBnc,PTBuc)) {
		pl = i(longley)
		if (sum(pl==longley) != 0) {
			status=FALSE
			if (!silent) {
				warning("not perturbing --")
				print(i)
			}
		}
		if (abs(mean(mean(pl-longley)))>.000001) {
			status=FALSE
			if (!silent) {
				warning("perturbations too big--")
				print(i)
			}
		}
	}
	for (i in c(PTBms,PTBn,PTBu)) {
		pl = i(longley)
		if (sum(pl==longley) != 0) {
			status=FALSE
			if (!silent) {
				warning("not perturbing --")
				print(i)
			}
		}
		if (abs(mean(mean(pl-longley)))>1) {
			status=FALSE
			if (!silent) {
				warning("perturbations too big--")
				print(i)
			}
		}
	}
	for (i in c(PTBns,PTBus)) {
		pl = i(longley)
		if (sum(pl==longley) != 0) {
			status=FALSE
			if (!silent) {
				warning("not perturbing --")
				print(i)
			}
		}
		if (abs(mean(mean(pl-longley)))>abs(mean(mean(longley)))) {
			status=FALSE
			if (!silent) {
				warning("perturbations too big--")
				print(i)
			}
		}
	}
	for (i in c(PTBnbr,PTBubr)) {
		t = runif(20)*2-1
		pl = i(t, lbound=-1, ubound=1)
		if (sum(t>=1) > 0 || sum(t<=-1)> 0) {
			status = FALSE
			if (!silent) {
				warning("bounded perturbations not bounded--")
				print(i)
			}
		}
	}

	#data test
	plongley = perturb(longley,lm,Employed~., ptb.R=10,
    	   ptb.ran.gen=c(PTBi, replicate(5,PTBus),PTBi), ptb.s=c(1,replicate(5,.001),1))
	sp=summary(plongley)
	coef= attr(sp,"coef.betas.m")
	if (nrow(coef)!=10 || ncol(coef)!=7) {
		status = FALSE
		if (!silent) {
			warning("perturb framework malfunction")
			print(i)
		}
	}

	return(status)
}


#
# Summary generic
#
# Iterates over replications and summarizes each, then merges results into
# a results matrix


summary.sensitivity<-function (object,...) {

	tmps = list()
	for (i in 1:length(object)) {
		tmps[[i]] = sensitivitySummaryIteration(object[[i]])
	}
	ret = sensitivitySummaryMerge(tmps)	
	attr(ret, "R") = attr(object, "R")
	class(ret)="sensitivity.summary";	
	return (ret);
}

######################################################
#          
# sensitivitySummaryIteration
#          
# [Internal]
#
# Extracts the coefficients and standard errors from the
# results of running summary() on an arbitrary analysis
#       
# Parameters:
#
# iteration - object produced by summary, must support
#		coef() method
#       
######################################################


sensitivitySummaryIteration<-function(iteration) {
	# generate individual summaries for list of replications
	tmp = try(summary(iteration),silent=T)
	if (is.null(cf<-coef(tmp)) || inherits( tmp,"try-error")){
	   if(is.null(cf<-coef(iteration))) {
		return(NULL)
	   }
	}
	coef.names=names(cf)
	if ( inherits( coef.formula<-try(tmp[["call"]],silent=T),"try-error")) {
	   if ( inherits( coef.formula<-try(tmp[["formula"]],silent=T),"try-error")) {
		coef.formula=NULL
	    }
	}
	if ( is.null(dim(cf)[2]) || (dim(cf)[2]<2) ) {
		coef.betas=cf
		coef.stderrs=NULL
	} else {
		coef.betas=cf[,1]
		coef.stderrs=cf[,2]
	}	
	
	tmp = coef.betas
 	attr(tmp,"coef.names") = coef.names
	attr(tmp,"coef.betas") = coef.betas
	attr(tmp,"coef.stderrs") = coef.stderrs
	attr(tmp,"coef.formula") = coef.formula
	return(tmp)	
}

######################################################
#          
# sensitivitySummaryMerge
#          
# [Internal]
#
# Merges summary results for each replicate into a results matrix
#       
# Parameters:
#
# A list of summary() results, where each element in the results
# is summary() for a replicate
#       
######################################################

sensitivitySummaryMerge<-function(s) {
	n = length(s)

	# check consistency and summarize
	coef.names = coef.names.m = attr(s[[1]],"coef.names")
	coef.betas = coef.betas.m = attr(s[[1]],"coef.betas")
	coef.stderrs = coef.stderrs.m = attr(s[[1]],"coef.stderrs")
	nomatch = FALSE
	if (n>1) {
	   for (i in 2:n) {
		if ( sum( attr(s[[i]],"coef.names") != coef.names)>0
		   || sum( length(attr(s[[i]],"coef.betas")) != length( coef.betas))>0
		   || sum( length(attr(s[[i]],"coef.stderrs")) != length( coef.stderrs))>0
		) {
			nomatch=TRUE
		}
		coef.names.m = rbind(coef.names.m, attr(s[[i]],"coef.names"))
		coef.betas.m = rbind(coef.betas.m, attr(s[[i]],"coef.betas"))
		coef.stderrs.m = rbind(coef.stderrs.m, attr(s[[i]],"coef.stderrs"))
	   }
	}
	if (nomatch) {
		warning("replications do not match (", i,")")
	}
	ret = list(coef.betas.m)
	row.names(coef.betas.m) = NULL
	row.names(coef.stderrs.m) = NULL

	attr(ret,"coef.names.m") = coef.names.m
	attr(ret,"coef.betas.m") = coef.betas.m
	attr(ret,"coef.stderrs.m") = coef.stderrs.m

	return(ret)
}
