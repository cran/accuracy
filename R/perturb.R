# perturb.R
#
# Performs sensitivity analysis of non-linear and linear models
# through data pertubations
#
# Part of the Accuracy package. Available from www.r-project.org and
# www.hmdc.harvard.edu/numerical_issues/
#
#    Copyright (C) 2004  Micah Altman
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



# perturb
#
# Replicates an analysis with randomly perturbed data
#
# data = data
# statistic = stat function to run
# ptb.R = number of replications
# ran.gen= function or vector of functions to apply to data to add noise
# ptb.s = size, or vector of sizes for noise generators

perturb<-function(data,statistic,..., ptb.R=1,ptb.ran.gen=PTBms,
	ptb.s=NULL) {
	ptb.R=as.integer(trunc(ptb.R))

	if ( (is.integer(ptb.R)==FALSE) || (ptb.R<=0) || (length(ptb.R)>1) ) {
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
	}

	# replicate swallows ... in R 1.9.0 ???!!??, worked in 1.6-1.8
	#retval = replicate(ptb.R,
	#	perturbHarness(data=data,ran.gen= ptb.ran.gen,statistic=statistic,...)
	#, simplify=FALSE);

	attr(retval,"ran.gen")= ptb.ran.gen
	attr(retval, "R") = ptb.R
	attr(retval, "s") = ptb.s
	attr(retval, "statistic") = statistic
	class(retval)="perturb"
	return(retval)
}

print.perturb<-function(x,...) {
	cat("Replications: \n",attr(x,"R"),"\n\n")
	cat("ran.gen: \n")
	print(attr(x,"ran.gen"))
	cat("s: \n")
	print(attr(x,"s"))
	cat("statistic: \n")
	print(attr(x,"statistic"))
}

perturbHarness<-function(data,ran.gen,statistic,..., ptb.s=NULL) {
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
		if (is.null(ptb.s)) {
			ndata[i] = ran.gen[[i]](ndata[i])
		} else {
			ndata[i] = ran.gen[[i]](ndata[i],size=ptb.s[i])
		}
	}
	
	stat = statistic(data=ndata,...)
	
	return(stat)
}


############################################################
# 
# Perturbations for vectors
#
# These functions take a vector and apply a mean-zero
# random perturbation to them.
#
# Fairly forgiving as to inputs. Will accept a list,
# matrix, scalar, or dataframe as well. However,
# applying a centered perturbation to a dataframe or
# matrix centers it on the whole dataframe or matrix, 
# rather than centering it on the vector. This is probably
# not what you want.
#
#
# Perturb.unif -- uniform perturbations 
# Perturb.meps -- perturbations on the order of storage
#		  roundoff error
# Perturb.norm -- normal perturbations
# 
#
# Inputs: 
#	x - vector or matrix
#	s - size of perturbation
#	centered - boolean, center disturbance to guarantee
#		(up to numeric tolerance) mean 0 for uniform
#		and normal disturbances, when used by itself
#	scaled - make disturbance size scaled to observation
#		size
#	Note: Using both scaled=TRUE and centered=TRUE, 
#	 	cannot guarantee mean 0-- the centering
# 		occurs before the scaling of the noise.
#		Once rescaled, the noise is no longer
#		guaranteed to be mean zero.

#
# Returns:
#	 perturbed vector, matrix or dataframe
#
##############################################################

PTBunif<-function(x ,size=1 , centered=FALSE, scaled=FALSE) {
	delta=runif(length(as.matrix(x)))*size; 
	
	if (scaled && centered) {
		warning("Using scaled and centering together is rarely what you really want to do.")
	}

	# centering to guarantee mean 0
	if (centered) {		
		delta=delta-mean(delta)
	}

	# scaled disturbance
	if (scaled) {
		delta = delta * x
	}	
	
	return(x+delta)
}

PTBnorm<-function(x , size=1, centered=FALSE, scaled=FALSE) {
	delta = rnorm( length(as.matrix(x)), mean=0, sd= 1) * size
	

	if (scaled && centered) {
		warning("Using scaled and centering together is rarely what you really want to do.")
	}

	# centering to guarantee mean 0
	if (centered) {
		delta = (delta - mean(delta))/sqrt(var(delta))
	}

	if (scaled) {
		delta = delta * x
	}	
	
	return(x+delta)
}

PTBmeps<-function(x, centered=FALSE, scaled=FALSE, size=1) {	
	# scaled=FALSE leads to larger roundoff 
	if (!scaled) {
		warning("scaled=FALSE probably not what you want for PTBmeps")
	}

	n = length(as.matrix(x))
	
	# note centering here does not guarantee zero mean!
	if (centered) {
		delta=integer(n)
		delta[1:floor((n+1)/2)]=.Machine$double.eps *size
		delta[ceiling((n+1)/2):n]=-1*.Machine$double.neg.eps *size
		delta=sample(delta,n)
	} else {
		delta= runif(n)*2-1
		delta= floor(delta)*.Machine$double.neg.eps +
			ceiling(delta)*.Machine$double.eps
	}

	if (scaled) {
		delta = delta*x
	}
	return(x+delta)
}

# Wrapper functions for bounded perturbations

PTBmsb<-function(x,size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="trunc",scaled=TRUE,
		ran.gen=PTBmeps)
}

PTBmsbr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="resample", scaled=TRUE,
		ran.gen=PTBmeps)
}

PTBubr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="resample",
		ran.gen=PTBunif)
}

PTBubrr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="relresample",
		ran.gen=PTBunif)
}

PTBnbr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="resample",
		ran.gen=PTBunif)
}

PTBnbrr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="relresample",
		ran.gen=PTBunif)
}

PTBusbr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="resample",scaled=TRUE,
		ran.gen=PTBunif)
}

PTBusbrr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="relresample",scaled=TRUE,
		ran.gen=PTBunif)
}

PTBnsbr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="resample", scaled=TRUE,
		ran.gen=PTBunif)
}

PTBnsbrr<-function(x, size=1, lbound=0,ubound=1) {
	ptbBndHarness(x,lbound,ubound,size=size,mode="relresample",scaled=TRUE,
		ran.gen=PTBunif)
}

# bounded perturbation -- takes set of samples and ensure results
# are within variable bounds

ptbBndHarness<-function(x,lbound,ubound,s=NULL,mode="trunc",
	ran.gen=PTBunif, maxiter=5000,...) {

	# check args
	if (mode!="trunc" && mode !="resample" && mode != "relresample") {
		stop("unknown mode")
	}

	# initialize bounds vectors
	if (length(lbound)==1) {
		lbound=replicate(length(x),lbound)
	} else if (length(lbound)!=length(x)) {
		stop("lbound must be scalar or same length as x")
	}
	if (length(ubound)==1) {
		ubound=replicate(length(x),ubound)
	} else if (length(ubound)!=length(x)) {
		stop("ubound must be scalar or same length as x")
	}

	# set sizes for relresample mode
	if (mode=="relresample") {
		tmps=pmin(x-lbound,ubound-x)
		s = pmin(s,tmps)
	}


	# initial sample	
	if (is.null(s)) {
		retv = ran.gen(x,...)
	} else {
		retv = ran.gen(x, s=s,...)
	}
	
	if ( mode=="resample" || mode=="relresample") {
		badi = which( retv<lbound | retv>ubound )
		iter=0
		while (length(badi)>0 && iter<=maxiter) {
			if (is.null(s)) {
				retv[badi] = ran.gen(x[badi],...)
			} else {
				retv[badi] = ran.gen(x[badi], s=s,...)
			}
			badi = which( retv<lbound | retv>ubound )
			iter=iter+1;	
		}
		if (iter>maxiter) {
			warning("resampling: maximum iterations exceeded, truncating to bounds")
		}
	} 

	# this is 'trunc' mode, also catches exceeded iterations on
	# resampling, and numerical issues in comparisons

	retv=pmax(lbound,retv)
	retv=pmin(ubound,retv)

	return(retv)
	
}


###############################################################
#
# These are wrappers around ptb.{meps,unif,norm}
# to make it easier to use the ptb vector
# functions from the perturb() harness.
#
###############################################################

PTBi<-function(x,size=1) {
	return(x)
}

PTBus<-function(x ,size=1) {
	return(PTBunif(x,size=size,scaled=TRUE))
}

PTBuc<-function(x ,size=1) {
	return(PTBunif(x,centered=TRUE,size=size))
}

PTBu<-function(x ,size=1) {
	return(PTBunif(x,centered=FALSE,scaled=FALSE,size=size))
}


PTBns<-function(x ,size=1) {
	return(PTBnorm(x,scaled=TRUE,size=size))
}

PTBnc<-function(x ,size=1) {
	return(PTBnorm(x,centered=TRUE,size=size))
}

PTBn<-function(x ,size=1) {
	return(PTBnorm(x,centered=FALSE,scaled=FALSE,size=1))
}

PTBms<-function(x,size=1) {
	return(PTBmeps(x,centered=FALSE,scaled=TRUE,size=size))
}


#
# Summary Functions
#

summary.perturb<-function (object,...) {
	n = length(object)
	s = vector(n,mode="list")

	# generate individual summaries for list of replications
	for (i in 1:n) {
		tmp = summary(object[[i]])
		if (inherits(tmp,"summary.lm")) {
			coef.names = dimnames(tmp[["coefficients"]])[[1]]
			coef.betas = tmp[["coefficients"]][,1]
			coef.stderrs = tmp[["coefficients"]][,2]
			coef.formula = tmp[["formula"]]
		} else if (inherits(tmp,"summary.mle")) {
			coef.names = dimnames(attr(tmp,"coef"))[[1]]
			coef.betas = attr(tmp,"coef")[,1]
			coef.stderrs = attr(tmp,"coef")[,2]
			coef.formula = attr(tmp,"call")
		} else if (inherits(tmp,"summary.nls")) {
			coef.names = dimnames(tmp[["parameters"]])[[1]]
			coef.betas = tmp[["parameters"]][,1]
			coef.stderrs = tmp[["parameters"]][,2]
			coef.formula = tmp[["formula"]]
		} else if (inherits(tmp,"summary.glm")) {
			coef.names = dimnames(tmp[["coefficients"]])[[1]]
			coef.betas = tmp[["coefficients"]][,1]
			coef.stderrs = tmp[["coefficients"]][,2]
			coef.formula = tmp[["terms"]]
		} else {
			coef.names=names(coef(tmp))
			coef.betas=coef(tmp)
			coef.stderrs=NULL
			coef.formula=NULL
			if (is.null(coef.betas)) {
				stop("Don't know how to summarize replications of type ",
					class(tmp) )
			}
		}
	
 		attr(tmp,"coef.names") = coef.names
		attr(tmp,"coef.betas") = coef.betas
		attr(tmp,"coef.stderrs") = coef.stderrs
		attr(tmp,"coef.formula") = coef.formula
		s[[i]]=tmp
	}
	
	# check consistency and summarize
	coef.names = attr(s[[1]],"coef.names")
	coef.names.m = coef.names
	coef.betas = attr(s[[1]],"coef.betas")
	coef.betas.m = coef.betas
	coef.stderrs= attr(s[[1]],"coef.stderrs")
	coef.stderrs.m = coef.stderrs;	
	coef.formula = attr(s[[1]],"coef.formula")
	if (n>1) {
	   for (i in 2:n) {
		if ( sum( attr(s[[i]],"coef.names") != coef.names)>0
		   || sum( attr(s[[i]],"coef.formula") != coef.formula)>0
		   || sum( length(attr(s[[i]],"coef.betas")) != length( coef.betas))>0
		   || sum( length(attr(s[[i]],"coef.stderrs")) != length( coef.stderrs))>0
		) {
			warning("replications does not match (", i,")")
		}
		coef.names.m = rbind(coef.names.m, attr(s[[i]],"coef.names"))
		coef.betas.m = rbind(coef.betas.m, attr(s[[i]],"coef.betas"))
		coef.stderrs.m = rbind(coef.stderrs.m, attr(s[[i]],"coef.stderrs"))
	   }
	}


	row.names(coef.betas.m) = NULL
	row.names(coef.stderrs.m) = NULL
	attr(s,"coef.names.m") = coef.names.m
	attr(s,"coef.betas.m") = coef.betas.m
	attr(s,"coef.stderrs.m") = coef.stderrs.m
	attr(s,"coef.formula") = coef.formula


	attr(s, "ran.gen")= attr(object, "ran.gen")
	attr(s, "R") = attr(object, "R")
	attr(s, "s") = attr(object, "s")
	attr(s,"statistic") = attr(object, "statistic")

	class(s)="perturbS";	
	return(s)
}

print.perturbS<-function(x,...) {

	cat("Replications: \n",attr(x,"R"),"\n\n")
	cat("ran.gen: \n")
	print(attr(x,"ran.gen"))
	cat("s: \n")
	print(attr(x,"s"))
	cat("statistic: \n")
	print(attr(x,"statistic"))

	cat("formula: \n","\n\n")
	print(attr(x,"coef.formula"))

	cat("betas:\n\n")
	print(summary(attr(x,"coef.betas.m")))
	cat("stderrs:\n\n")
	print(summary(attr(x,"coef.stderrs.m")))
	
}


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
