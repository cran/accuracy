# globaltests
#
# Tests for global optimality of non-linear and maximum-likelihood solutions
#
# Part of the Accuracy package. Available from www.r-project.org and
# www.hmdc.harvard.edu/numerical_issues/
#
#    Copyright (C) 2004  Micah Altman , Michael McDonald
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


dehaan<-function(llTest, llMax, pval=.05 ) {
  # returns TRUE if "llMax" likelihood is greater than
  # the (1-pval) confidence interval for the true optimum
  # derived from "llTest" likelihoods.

  if ( !is.numeric(llTest) ||  !is.numeric(llMax) ) {
        stop("llTest must be numeric")
  }

  if (length(llTest)<2) {
        stop("Size of llTest must be greater than 1")
  }

  if (length(pval) !=1 || pval>=1 || pval<0)  {
        stop("pval must be scalar: 0< pval <1")
  }

  x = sort(llTest)
  n = length(llTest)


  # De Hann (1981) proves this test statistic works for some k(n) where
  # k(n)/n->0 as n->infinity.  A likely candidate is k=sqrt(n)

  k = n^(1/2)
  alpha = k/2

  lp = x[n]+ ( (x[n])-(x[n-1]) ) / (pval ^ (-1/alpha) -1)

  if (lp>llMax)  {
    return(FALSE)
  } else  {
    return(TRUE)
  }
}

dehaanSelftest<-function(silent = TRUE) {
        
        # Tests of deHaan function

        ret = TRUE

        # This test data should always return a FALSE.

        l1= c(1,2,3,4,5)
        max1 = 4; 

        if (dehaan(c(1,2,3,4,5),4) ||
	    !dehaan(c(1,2,3,4,5),1000) )  {
                ret = FALSE
                 if (!silent) {
                        warning("failed selftest ")
                }
        }

	if (  dehaan(1:100,1) ||
	      dehaan(1:100,101) || 
	      dehaan(1:100,99) || 
	      dehaan(1:20,21,pval=.9999) ||
              !dehaan(1:100,101,pval=.001)|| 
	      !dehaan(1:20,21) 
 	) {
           ret = FALSE
           if (!silent) {
                  warning("failed selftest")
           }
        } 

        return(ret)
}


starr<-function(betas, tol=.0001, dmethod="euclidean") {

  if (!is.matrix(betas)) {
    stop("betas must be matrix")
  }
  if (length(tol)!=1 || (tol<0) || !is.numeric(tol)) {
    stop("tolmust be scalar >=0")
  }

        return(ret)
}


starr<-function(betas, tol=.0001, dmethod="euclidean") {

  if (!is.matrix(betas)) {
    stop("betas must be matrix")
  }
  if (length(tol)!=1 || (tol<0) || !is.numeric(tol)) {
    stop("tolmust be scalar >=0")
  }

  # find unique (within tolerance) optima, count them

  n = nrow(betas)
  if (n==1)  {
	return(1)
  }
 
  dm = as.matrix(dist(betas,method=dmethod))
  optc = integer(n)
  for (i in 1:n) {
        optc[i] = sum(dm[i,]<tol)
  }

  sortInd = sort(optc,decreasing=TRUE,index.return=TRUE)$ix
  for (i in sortInd) {
    if (optc[i]>0) {
    tmp=optc[i]
        optc[dm[i,]<tol]=0; 
        optc[i] = tmp
    }
  }
  nopt = sum(optc>0); 

  ndouble = sum(optc==2)
  nsingle = sum(optc==1)

  rv = nsingle/n + 2*ndouble/(n*(n-1))
  
  return(rv)
                              
}

starrSelftest<-function (silent=TRUE) {

  # BOD test
  #
  # For BOD example, see help("BOD")

  ret = TRUE
  

  x=rbind(c(1,1,1), c(1,2,1), c(1,1.1,1), c(1,2,1), c(3,4,5))
  if ( (starr(x)!=.7) || 
       (starr(rbind(1)) !=1) ||
       (starr(rbind(1,1,2,2)) !=1/3) ||
       (starr(rbind(1,1.00001,2,2)) !=1/3) ||
       (starr (rbind(1,1,1,1,1,1,1,1,1,1,1,1)) !=0) 
     ) 
   {
    	ret=FALSE
    	if (!silent) {
       		warning("Starr test for global optimum failed selftest")
    	}
   }
  return(ret)

}

starrRun<-function(start,optfunc,...) {
  r= NULL
  dstart = as.data.frame(start)
  for (i in 1:nrow(dstart)) {
    r = rbind(r, coef(optfunc(..., t(dstart)[,i])))
  }
}
