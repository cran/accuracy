# LRE
#
# Tools for comparing model results
#
# Part of the Accuracy package. Available from www.r-project.org and
# www.hmdc.harvard.edu/numerical_issues/
#
#    Copyright (C) 2006  Micah Altman
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


######################################################
#       
# LRE
#       
# Calculates the Log Relative Error between two values.
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################



"LRE"<-
function (x, correct, use.LAE = TRUE, cutoff=16)
{
    zeros = which(correct == 0)
    nonzeros = which(correct != 0)
    nas = which(x == NA)
    res = double(length(x))
    res[nas] = NA
    res[nonzeros] = -1 * log10(abs(x[nonzeros] - correct[nonzeros])/abs(correct[nonzeros]))
    if (use.LAE) {
        res[zeros] = -1 * log10(abs(x[zeros]))
    }
    else {
        res[zeros] = NA
    }
    res[which(res>cutoff)]=cutoff

    dim(res)=dim(correct)
    names(res)=names(correct)
    dimnames(res)=dimnames(correct)
    return(res)
}

######################################################
#       
# modelsCompare
#       
# Compare coefficients from multiple models and display discrepancies
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################

modelsCompare<-function(
    model, model2=NULL,
    param.function=modelBetas,
    similarity.function=LRE,
    summary.function=min
) {

  # check arguments
  if (is.null(model2)) {
	models=model
  } else {
  	models = list(model,model2)
  }

  baseline=param.function(models[[1]])
  
  tmpComparisons=vector(mode="list", length=length(models)-1)
  for (i in 1:length(tmpComparisons)) {
    tmpComparisons[[i]]=similarity.function(param.function(models[[i+1]]), baseline )
  }

  
  ret= apply(simplify(tmpComparisons), 1, summary.function)
  if (!is.null(dim(baseline))) {
	   dim(ret)=dim(baseline)
  }
  if (!is.null(names(baseline))) {
  	names(ret)=names(baseline)
  }
  if (!is.null(dimnames(baseline))) {
  	dimnames(ret)=dimnames(baseline)
  }
  return(ret)

}

######################################################
#       
# LRE
#       
# Check if the results of two models agree to a specific number of digits
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################


modelsAgree<-function(model, model2=NULL, digits=3, ...) {
  if (is.null(model2)) {
	models=model
  } else {
  	models = list(model,model2)
  }
  min(modelsCompare(models,...))>=digits
}

######################################################
#       
# betas
#       
# attempt to extract betas
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################

modelBetas<-function(model) {
	tmp=sensitivitySummaryIteration(model)
	return(attr(tmp,"coef.betas"))
}


######################################################
#       
# coefSummary
#       
# convenience function, extracts coefficients from summary() of model
#       
# Parameters:
#
# See the R documentation file for details of each argument and return value
# 
######################################################

modelSummary<-function(model) {
	tmp=sensitivitySummaryIteration(model)
	return(tmp)
}

######################################################
#       
# modelCompareSelfTest 
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


modelsCompareSelfTest<-function(silent=TRUE) {
	status=TRUE

	if (LRE(10,10)!=16 || LRE(10,1)>0 || !isTRUE(all.equal(LRE(1.000001,1),6))
       || !isTRUE(all.equal(LRE(.01,0),2)) ) {
		status = FALSE
		if (!silent) {
			warning("LRE malfunction")
		}
	}
	
	data(longley)
	longley2=as.data.frame(longley+2)
	
	l1 = lm(longley)
	l2 = lm(longley2)
	
	if (!modelsAgree(l1,l1)) {
	  status = FALSE
		if (!silent) {
			warning("Identical models should agree")
		}
	}
	
	if (modelsAgree(l1,l2)) {
	  status = FALSE
		if (!silent) {
			warning("Different models should disagree")
		}
	}
  
  if (!modelsAgree(l1,l2, digits=1)) {
	  status = FALSE
		if (!silent) {
			warning("Doesn't match rough agreement")
		}
	}

	return(status)
}
