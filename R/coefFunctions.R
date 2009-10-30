#
# Internal helper function to extract and manipulate
# model coefficients
#


######################################################
#          
# coefarrayFlatten
#          
# [Internal]
#
# Attemts to extract betas from coef arrays
#       
# Parameters:
#
#	x- 3d array of coefficients, as returned by lme
#       
######################################################

coefarrayFlatten<-function(x, silent=T) {
  if (!inherits(x,"array")) {
    return(x)
  } else {
    if (length(dim(x))!=3) {
        if (!silent) {
          warning("can't reshape")
        }
        return(x)
    }
  }
  tmpft=aperm(x,c(1,3,2))
  dft=dim(tmpft)
  tmpft=matrix(tmpft,ncol=dim(tmpft)[3])
  dimnames(tmpft)[[2]]=dimnames(x)[[2]]
  dimnames(tmpft)[[1]]=as.vector(
    outer(dimnames(x)[[1]],dimnames(x)[[3]],function(...)paste(...,sep=":")))
  return(tmpft)
}


######################################################
#          
# sensitivitySummaryIteratiion
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

  tmp=forceDispatch("summary",iteration)	
  cf<-NULL
  if (!is.null(tmp)) {
		cf<-forceDispatch("coef",tmp)
  } 
  if (is.null(cf)) {
		cf<-forceDispatch("coef",iteration)		
  }
  if (is.null(cf)) {
		 return(NULL)
  }
  

	coef.names=names(cf)
  if (is.null(coef.names)){
    coef.names=rownames(cf)
  }
	if ( inherits( coef.formula<-tryR(tmp[["call"]],silent=T),"try-error")) {
	   if ( inherits( coef.formula<-tryR(tmp[["formula"]],silent=T),"try-error")) {
		coef.formula=NULL
	    }
	}

	if (is.array(cf) && length(dim(cf)==3)) {
		cf=coefarrayFlatten(cf)
	}

	if ( is.null(dim(cf)[2]) || (dim(cf)[2]<2) ) {
		coef.betas=cf
		coef.stderrs=NULL
	} else {
		coef.betas=cf[,1]
		coef.stderrs=cf[,2]
	}	
	
	tmp = cf
 	attr(tmp,"coef.names") = coef.names
	attr(tmp,"coef.betas") = coef.betas
	attr(tmp,"coef.stderrs") = coef.stderrs
	attr(tmp,"coef.formula") = coef.formula
	return(tmp)	
}


######################################################                            
#          
# ForceDispatch
#          
# [Internal]
#
# Forces the dispatch of a method
#       
# Parameters:
#
#       
######################################################


forceDispatch<-function(f="summary",object=NULL) {
      res=tryR(do.call(f,list(object)),silent=T)
      if (!inherits(res,"try-error")) {
           return(res)
      }
      for ( cl in class(object) ) {
	  if (is.R()) {
          	ftmp<-tryR(methods:getFunction(paste(f,cl,sep=".")),silent=T)
	  } else {
          	ftmp<-tryR(getFunction(paste(f,cl,sep=".")),silent=T)
	  }
          if (!inherits(ftmp,"try-error")) {
              break;
          }
          if (is.R()) {
             ftmp<- tryR(getS3method(f=f,class=cl),silent=T)
             if (!inherits(ftmp,"try-error")) {
              break;
             }
             ftmp<- tryR(methods:getMethod(f=f,signature=cl),silent=T)
             if (!inherits(ftmp,"try-error")) {
              break;
             }
          } else {
             ftmp<- tryR(getMethod(f=f,sig=cl),silent=T)
             if (!inherits(ftmp,"try-error")) {
              break;
             }
          }
       }
       if (inherits(ftmp,"try-error")) {
          return(NULL)
       }
       res = tryR(ftmp(object),silent=T)
       if (inherits(res,"try-error")){
          return(NULL)
       }
       return(res)
}


