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
    zeros = which(x == 0)
    nonzeros = which(x != 0)
    nas = which(x == NA)
    res = double(length(x))
    res[nas] = NA
    res[nonzeros] = -1 * log10(abs(x[nonzeros] - correct[nonzeros])/abs(correct[nonzeros]))
    if (use.LAE) {
        res[zeros] = -1 * log10(abs(x[zeros] - abs(correct[zeros])))
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
    models=list(),
    param.function=function(model)coef(summary(model)),
    similarity.function=LRE,
    summary.function=min
) {

  # check arguments

  baseline=param.function(models[[1]])
  tmpComparisons=sapply( models[2:length(models)],
      function(tmpCmp){similarity.function(param.function(tmpCmp), baseline )},
      simplify=TRUE)
  apply(tmpComparisons, 1, summary.function)
  ret= apply(tmpComparisons, 1, summary.function)
  names(ret)=names(baseline)
  dim(ret)=dim(baseline)
  dimnames(ret)=dimnames(baseline)
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


modelsAgree<-function(model, model2, digits=3, ...) {
  models = list(model,model2)
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
