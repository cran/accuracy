# Perturbation functions
#
# These functions perturb a vector of observations in various ways
# 
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

# These are used in conjunction sith sensitivity() to perform a sensitivity analysis.o
# This provides several types of perturbation functions.
#
# - perturbations for ordered discrete variables
#	Observations are reclassified via a reclassification matrix. 
#	The reclassification matrix is of a form such that the probability of
# 	reclassification is monotonically decreasing with the number of levels of
#	reclassification
#
# - perturbations of unordered discrete variables
#	Reclassification via a reclass matrix. Default is completely random
#
# - perturbations for continuous variables
#	Mean-zero noise of various forms is added to the variable
#
# - perturbations for bounded continuous variables
#	Mean-zero noise is added, rsults are truncated/resampled to bounds
#


##################################################################
#
# Perturbations functions for vectors of discrete observations
#
# There are really three core functions here:
#	PTBfactor -- perturb a factor using a reclassification matrix
#	reclass.mat.diag -- create a diagonal reclass matrix
#	reclass.mat.random -- create a uniform reclass matrix
#
#################################################################

######################################################
#       
# PTBfactor
#       
# This perturbs a factor, its the workhorse method
# It uses a classification matrix generated by reclass.mat.* or by reclassify()
# in the perturb package
#
# This provided discrete (permutation) perturbations
# it is essentially a wrapper for PTBfactor
#       
# Parameters:
#
# See the R documentation file for PTBdiscrete
# 
######################################################


PTBfactor<-function(x,r) {
  x=as.factor(x)

  if (!inherits(r,"reclassify") && !is.matrix(r)) {
		stop("r must be a reclassification object, or cumulative probability matrix")
	}
	if (inherits(r,"reclassify")) {
		cumprob = r$cum.reclass.prob
	} else {
		cumprob = r
	}

	ret= sapply( as.data.frame(t(cumprob[x,]>runif(length(x)))),
				function(tmp){min(which(tmp))}
		)

  ret=factor(ret,labels=levels(x),levels=1:length(levels(x)))
	return(ret)
}

######################################################
#
# reclass.mat.diag
#
# This returns a cumulative probability matrix useful for reclassifying
# discrete variables with PTBfactor
#
# n is the number of factor levels
# q is the probability of remaining at the same level
#
# The transition probability of changing from i->j is 0 if |i-j|>1
#
# See the R documentation file for PTBdiscrete
# 
######################################################

reclass.mat.diag<-function(n,q) {
     if (n==0) {
        return(integer())
     }
     if (n==1) {
        return(1)
     }
  p=1-q
  # mostly a diagnal with 1's in the upper
	d = diag(1-p/2,nrow=n,ncol=n)
	d[upper.tri(d)]=1

  # fixup corner cases
	d[n,n-1]=p
	d[1,1]=1-p
	d[1,2]=1
	d[n,n]=1

	for (i in 2:n-1) {
		d[i,i-1]=(p/2)
	}
	return(d)
}

######################################################
#
# reclass.mat.diag
#
# This returns a cumulative probability matrix useful for reclassifying
# discrete variables with PTBfactor
#
# n is the number of factor levels
# q is the probability of remaining at the same level
#
# The transition probability of changing from i->j is uniform
#
# See the R documentation file for PTBdiscrete
# 
######################################################

reclass.mat.random<-function(n,q) {
     if (n==0) {
        return(integer())
     }
     if (n==1) {
        return(1)
     }
  p=1-q
  d = matrix(data=0, nrow=n, ncol=n)
	d[upper.tri(d)]=1
	d[n,n]=1

	inc = p/(n-1)
  for (i in 1:n) {
     for (j in 1:n-1) {
         if (j<i) {
            d[i,j] = j*inc
         } else if (j==i) {
            d[i,j] = (1-p)+(j-1)*inc
         } else {
            d[i,j] = 1 - ((n-j)*inc)
         }
      }
  }
	return(d)
}
######################################################
#       
# PTBdiscrete
#       
# This is the main public perturbation function for factors and
# other discrete observations.
#
# This provided discrete (permutation) perturbations
# it is essentially a wrapper for PTBfactor
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


PTBdiscrete<-function(x,r=NULL,q=.99) {
  if (is.null(r)) {
     r=reclass.mat.default(x,q)
  }

  if (is.integer(x)) {
     return(as.integer(as.character(PTBfactor(as.factor(x),r))))
  } else if (is.numeric(x)) {
     return(as.numeric(as.character(PTBfactor(as.factor(x),r))))
  } else if (is.character(x)) {
    return(PTBcharacter(x,r));
  } else if (is.logical(x)) {
    return(PTBlogical(x,r))
  } else if(is.factor(x)) {
    return(PTBfactor(x,r))
  } else {
    warning("coercing x to a factor")
    return(PTBfactor(as.factor(x),r))
  }
}


######################################################
#       
# reclass.mat.default
#       
# Produces a default reclassification matrix dependent on type
#
# Parameters:
#
# See the R documentation file 
# 
######################################################

reclass.mat.default<-function(x,q=.99) {

     if (is.logical(x)) {
        r=reclass.mat.random(2,q)
     } else if(is.numeric(x)) {
       r=reclass.mat.diag(length(levels(as.factor(x))),q)
     } else if(is.ordered(x)) {
       r=reclass.mat.diag(length(levels(x)),q)
     } else {
        r=reclass.mat.random(length(levels(as.factor(x))),q)
     }
     return(r)
}

######################################################
#       
# PTBdefaultfn
#       
# Heuristic to supply  a default perturbation function given a vector of data
#
# Parameters:
#
# See the R documentation file 
# 
######################################################


PTBdefaultfn<-function(vec,q=.99) {

    if (is.numeric(vec) && !is.discrete(vec)) {
       eval(substitute(function(x){PTBus(x,s)},list(s={1-q})))
    } else {
       tmpr=reclass.mat.default(vec,q)
       eval(substitute(function(x){PTBdiscrete(x,s)},list(s=tmpr)))
    }
}


######################################################
#       
# PTBdefault
#       
# Heuristic to supply  a default perturbation function given a vector of data
#
# Parameters:
#
# See the R documentation file 
# 
######################################################


PTBdefault<-function(vec,q=.99) {

    if (is.numeric(vec) && !is.discrete(vec)) {
       return(PTBus(vec,1-q))
    } else {
       tmpr=reclass.mat.default(vec,q)
       return(PTBdiscrete(vec,tmpr))
    }
}

######################################################
#       
# is.discrete
#       
# Attempts to evaluate whether vector contains discrete values
#
# Parameters:
#
# x - a vector
# 
######################################################

is.discrete<-function(x) {
	if (is.data.frame(x) || is.list(x)) {
		return(sapply(x,is.discrete))
	}
	x=as.vector(x);
	if (is.integer(x)) {
		return (TRUE);
	}  else if (is.factor(x)) {
		return (TRUE);
  } else if (is.logical(x)) {
    return (TRUE);
  } else if (is.character(x)) {
    return(TRUE);
	} else  if (sum(as.integer(x)!=x)==0) {
		return (TRUE);
	}
	return (FALSE);
}

######################################################
#       
# PTBlogical
#       
# This perturbs a logical variable using a classification matrix
#
# This provided discrete (permutation) perturbations
# it is essentially a wrapper for PTBfactor
#       
# Parameters:
#
# See the R documentation file for PTBdiscrete
# 
######################################################


PTBlogical<-function(x,r) {

    as.logical(PTBdiscrete(as.factor(x),r))
}

######################################################
#       
# PTBcharacter
#       
# This perturbs a character variable using a classification matrix
#
# This provided discrete (permutation) perturbations
# it is essentially a wrapper for PTBfactor
#       
# Parameters:
#
# See the R documentation file for PTBdiscrete
# 
######################################################


PTBcharacter<-function(x,r) {
   as.character(PTBdiscrete(as.factor(x),r))
}


##################################################################
#
# Perturbations functions for continuous observations
#
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
# The workhorses are:
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
# See the R documentation file for more information
#
##############################################################

PTBunif<-function(x ,size=1 , centered=FALSE, scaled=FALSE) {
	delta=(runif(length(as.matrix(x)))-0.5)*size; 
	
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
############################################################
#
# Multiple Perturbations for Continuos vectors
#
# These functions take a vector and apply repeated uniform
# or normal perturbations to them.
#
#
# PTBmultiunif -- multiple uniform perturbations, similar to PTBunif
# PTBmutinorm -- multiple uniform perturbations, similar to PTBnorm
#
#
# Inputs:
#
#        Identical to PTBunif, with the addition of
#        reps = number of repeater perturbations
#
# Returns:
#        perturbed vector, matrix or dataframe
#
##############################################################


PTBmultinorm<-function(x , size=1 , centered=FALSE, scaled=FALSE, reps=1) {
       if (reps<1) {
	stop("Reps less than 1");
       }
       delta=vector(length=length(as.matrix(x)));
       for (i in 1:trunc(reps)) {
	         delta=delta+rnorm( length(as.matrix(x)), mean=0, sd= 1) * size;
       }
       delta = delta/reps;

        if (scaled && centered) {
                warning("Using scaled and centering together is rarely what you really want to do.")
        }

        # centering to guarantee mean 0
        if (centered) {
                delta = (delta - mean(delta))/sqrt(var(delta))
        }

        # scaled disturbance
        if (scaled) {
                delta = delta * x
        }

        return(x+delta)
}

PTBmultiunif<-function(x , size=1 , centered=FALSE, scaled=FALSE, reps=1) {
       if (reps<1) {
	stop("Reps less than 1");
       }
      delta=vector(length=length(as.matrix(x)));
       for (i in 1:trunc(reps)) {
	  delta=delta+((runif(length(as.matrix(x)))-0.5)*size);
       }
      delta = delta/reps;

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

############################################################
#
# Multi-Perturbation Generator Function
#
#
#
# PTBmu.gen (rep=1)
# 
# Returns a PTMmultiXXXX function with a pre-configured number of reps
# 
# PTBmu.gen (rep=1)
#
# Returns a PTMmultiXXXX function with a pre-configured number of reps
#
##############################################################

PTBmu.gen<-function(reps=1) {
     
     eval(substitute(function(x, size=1) 
        {return(PTBmultiunif(x,reps=reps,size=size))}
        ,list(reps=reps)))
}
        
PTBmn.gen<-function(reps=1) {
     eval(substitute(function(x, size=1) 
        {return(PTBmultinorm(x,reps=reps,size=size))}
        ,list(reps=reps)))              
}       
##################################################################
#
# perturbations for bounded continuous variables
#
#	Mean-zero noise is added, rsults are truncated/resampled to bounds
#
# The workhorse is ptbBndHarness(), the rest are simply wrappers.
#
##################################################################


######################################################
#       
# ptbBndHarness
#
# [Internal]
#       
# bounded perturbation -- takes set of samples and returns
# perturbed results that  are within variable bounds
#
# Parameters:
#
# x - vector
# lbound - lower bound
# ubound - upper bound
#
# mode - bounding mode. 
#	- trunc -- simply truncates perturbed vector to the bounds
#	- resample -- retries perturbations for observations that
#		exceed bounds
#	- rel.resample -- noise is rescaled when observations are close to
#		boundaries (see citations in documentation)
#
# maxiter - number of tries to resample, if necessary. 
#	If maxiter is exceeded, then 'trunc' is used on
# 	remaining observations
# s - size of perturbation (to be passed to ran.gen)
# ran.gen - random generation function
# 
# See the R documentation of PTBmsb for more details
#
######################################################

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

######################################################
#       
# PTB*b*
#       
# Wrapper functions around ptbBndHarness to produce
# various forms of bounded perturbations
#       
# See the R documentation of PTBmsb for more details
# on the perturbation forms
#               
# Parameters:           
#                               
# x - vector            
# lbound - lower bound          
# ubound - upper bound  
# size - size of perturbation (to be passed to ran.gen)
#
# See the R documentation of PTBmsb for more details
#
######################################################

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

