# html
#
# presentation generic routines for sensitivity()
#
# Part of the Accuracy package. Available from www.r-project.org and
# www.hmdc.harvard.edu/numerical_issues/
#
#    Copyright (C) 2004-6 Micah Altman
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
# See the documentation for plot.default, print.default, anova for arguments, etc.
# 
# Also see HTML.default in R2HTML, for the HTML pretty printing routined


#
# HTML() methods
#

###
# Workaround for codetools false warning
###

HTML<-try(R2HTML::HTML,silent=TRUE)
if (inherits(HTML,"try-error")){
	HTML<-print
}


HTML.sensitivity.summary<-function(x,quiet=TRUE,...) {
    printGoodRuns(x,printMethod=HTML,...)
    printBetaSummary(x,printMethod=HTML,...)
    if (!quiet) {
      printErrSummary(x,printMethod=HTML,...)
    }
}

HTML.sensitivity.anova<-function(x,...) {      
    printGoodRuns(x,printMethod=HTML,...)
    printAnovaSummary(x,printMethod=HTML,...)
}

#
# print() methods 
# 


print.sensitivity.summary<-function(x,quiet=TRUE,...) {
    printGoodRuns(x,...)
    printBetaSummary(x,...)
    if (!quiet) {
      printErrSummary(x,...)
    }
}

print.sensitivity<-function(x,quiet=TRUE,...) {
  print("perturbation list")
	if (!quiet) {
     cat("ran.gen: \n")
	   print(attr(x,"ran.gen"))
	   cat("s: \n")
	   print(attr(x,"s"))
	   cat("statistic: \n")
	   print(attr(x,"statistic"))
   }
  cat("Class: \n", class(attr(x,"baseline")),"\n\n")
  cat("Replications: \n",attr(x,"R"),"\n\n")
}

print.sensitivity.anova<-function(x,quiet=TRUE,...) {
   printGoodRuns(x,...)
   printAnovaSummary(x,...)
}

#
# plot() methods 
# 

plot.sensitivity.summary<-function(x,...) {
   cf=as.data.frame(attr(x,"coef.betas.m"))
   if (is.R()){
      op=par(no.readonly=T)
   } else {
      op=par()
   }
   par(mfrow=c(1,length(cf)))
   for (i in 1:length(cf)) {
       boxplot(cf[i],...)
       title(main=names(cf[i]))
   }
   par(op)
   return(invisible())
}

plot.sensitivity<-function(x,ask=dev.interactive(),...) {
  if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
  }
  count=1
  for (tmp in x) {
     cat("\n\n****Plotting perturbation #", count,"\n")
     count=count+1
     plot(tmp,...)
  }
  return(invisible())
}

plot.sensitivity.anova<-function(x,...) {
  plot.sensitivity(x,...)
  return(invisible())
}


#
# anova() methods

anova.sensitivity<-function(object,...) {
   # Applies  anova to each replication, summarizes  resulting matrix of 
   # anova coeficients

        n = length(object)
        s = vector(n,mode="list")

        # generate individual summaries for list of replications
        for (i in 1:n) {

      	    tmp = try(anova(object[[i]]),silent=FALSE)
                if (!inherits(tmp,"anova")) {
		               warning("Could not perform anova on perturbation #", i);
	                 	next;
                } 
                coef.names = attr(tmp,"row.names")
	              coef.df = tmp[[1]]
  	            coef.sumsq=tmp[[2]]
  	            coef.meansq=tmp[[3]]
	              coef.f=tmp[[4]]
	              coef.pr=tmp[[5]]
               
                attr(tmp,"coef.names") = coef.names
                attr(tmp,"coef.df") = coef.df
                attr(tmp,"coef.sumsq") = coef.sumsq    
                attr(tmp,"coef.meansq") = coef.meansq
                attr(tmp,"coef.f") = coef.f
                attr(tmp,"coef.pr") = coef.pr          
                s[[i]]=tmp
        }
        s=s[which(sapply(s,is.null)==FALSE)]
        # check consistency and summarize
        coef.names = attr(s[[1]],"coef.names")
        coef.names.m = coef.names
        coef.df = attr(s[[1]],"coef.df")
        coef.df.m = coef.df
        coef.sumsq= attr(s[[1]],"coef.sumsq")
        coef.sumsq.m = coef.sumsq;
        coef.meansq= attr(s[[1]],"coef.meansq")
        coef.meansq.m = coef.meansq;
        coef.f= attr(s[[1]],"coef.f")
        coef.f.m = coef.f;
        coef.pr= attr(s[[1]],"coef.pr")
        coef.pr.m = coef.pr;

        if (n>1) {
           for (i in 2:n) {
                if ( sum( attr(s[[i]],"coef.names") != coef.names)>0
                   || sum( length(attr(s[[i]],"coef.df")) != length( coef.df))>0
                   || sum( length(attr(s[[i]],"coef.sumsq")) != length( coef.sumsq))>0
                   || sum( length(attr(s[[i]],"coef.meansq")) != length( coef.meansq))>0
                   || sum( length(attr(s[[i]],"coef.f")) != length( coef.f))>0
                   || sum( length(attr(s[[i]],"coef.pr")) != length( coef.pr))>0
                ) {
                        warning("replications do not match (", i,")")
                }
                coef.sumsq.m = rbind(coef.sumsq.m, attr(s[[i]],"coef.sumsq"))
                coef.meansq.m = rbind(coef.meansq.m, attr(s[[i]],"coef.meansq"))
                coef.f.m = rbind(coef.f.m, attr(s[[i]],"coef.f"))
                coef.pr.m = rbind(coef.pr.m, attr(s[[i]],"coef.pr"))
           }
        }


        rownames(coef.names.m) = NULL
        attr(s,"coef.names.m") = coef.names.m

        attr(s,"coef.df.m") = rename.matrix(coef.df.m,"DF")
        attr(s,"coef.sumsq.m") = rename.matrix(coef.sumsq.m,coef.names.m)
        attr(s,"coef.meansq.m") = rename.matrix(coef.meansq.m,coef.names.m)
        attr(s,"coef.pr.m") = rename.matrix(coef.pr.m,coef.names.m)
        attr(s,"coef.f.m") = rename.matrix(coef.f.m,coef.names.m)

        attr(s, "ran.gen")= attr(object, "ran.gen")
        attr(s, "R") = length(s)
        attr(s, "origR") = attr(object, "R")

        attr(s, "s") = attr(object, "s")
        attr(s,"statistic") = attr(object, "statistic")

        class(s)="sensitivity.anova";
        return(s)
}

#
# Internal Helper functions for print methods
#

######################################################
#       
# rename.matrix
#       
# [Internal] 
#
# Renames a matrix or dataframe , data-frame style, based on a name vector
#
# Returns: data frame
#
#       
# Parameters:
#
# m - matrix or data frame 
# nm - name  vector
# 
######################################################

rename.matrix<-function(m,nm){
     rownames(m)=NULL
     m=as.data.frame(matrix(m,ncol=length(nm)))
     names(m)=nm
     return(m)
}

######################################################
#       
# q95
#       
# [Internal] 
#
# Returns 95% quantiles of a matrix or data frame
#
# Returns: data frame
#       
# Parameters:
#
# m - matrix or data frame 
# nm - name  vector
# 
######################################################

q95<-function(x) {
   sapply(as.data.frame(x),function(x){quantile(x,na.rm=TRUE,probs=c(0.025,.975))})
}


######################################################
#       
# hsummary
#       
# [Internal] 
#
# Returns a html friendly descriptive statistics summary of a matrix
#
# Returns: summary matrix
#       
# Parameters:
#
# df - matrix or data frame 
# 
######################################################

hSummary<-function(df) {
     df=as.data.frame(df)
     resm=NULL;
     for (i in df) {
         i=na.omit(i)
         if (length(i)==0) {
          res =    c(NA,NA,NA,NA,NA,NA)
         }   else {
           res = c(mean(i),sd(i), min(i),
                 quantile(as.matrix(i),
                  probs=c(0.025,.975)),
                  max(i))
        }
        res=matrix(res,ncol=length(res))
        resm=rbind(resm,res)
     }
     colnames(resm) = c("mean","stdev","min", "2.5%","97.5%","max")
     rownames(resm)=colnames(df)
     return(resm)
}


###########################################
#
# Internal subroutines used in printing
#
###########################################

printBetaSummary<-function(x,printMethod=print,...) {
	hs = hSummary(attr(x,"coef.betas.m"))
        baseline= attr(x,"baselineSummary")
   if(!is.null(baseline)) {
		hs[,2] = attr(baseline,"coef.betas")
		colnames(hs)[2]="(Orig. Est.)"
		colnames(hs)[1]="Perturb Est."
		if (!is.null(attr(baseline,"coef.stderrs"))) {
			hs[,3] = attr(baseline,"coef.stderrs")
			colnames(hs)[3]="(Orig. Stderr)"
			hs=hs[,1:5]
			flags=character(length=length(dim(hs)[1]))
			for (i in 1:dim(hs)[1]) {
			   if ( ((hs[i,2]-1.96*hs[i,3]) > hs[i,4]) ||
               ((hs[i,2]+1.96*hs[i,3]) < hs[i,5])
            ) {
               flags[i]="*"
         }
         else {
               flags[i]=""
         }
      }
      hs=cbind(as.data.frame(hs),flags)
      names(hs)[6]="[Unstable]"
		}
  }
	printMethod(hs,...)
}

printErrSummary<-function(x,printMethod=print,...) {
    if (!is.null(attr(x,"coef.stderrs.m"))) {
        printMethod("Sensitivity of stderrs to perturbations:",...)
	printMethod("",...);
        printMethod (hSummary(attr(x,"coef.stderrs.m")),...)
     }
}

printGoodRuns<-function(x,printMethod=print,catMethod=cat, ...) {
     printMethod(paste("Sensitivity of coefficients over ", attr(x,"R"), " perturbations:"),...)
     if (attr(x,"R")!=attr(x,"origR")) {
        failures=attr(x,"origR")-attr(x,"R")
        printMethod(paste("( ", failures, " failures )"),...)
     }   
}

printAnovaSummary<-function(x,printMethod=print,...) {
        printMethod("",...)
        printMethod("sumsq:",...)
        printMethod(hSummary(attr(x,"coef.sumsq.m")), ...)
        printMethod("meansq:",...)
        printMethod(hSummary(attr(x,"coef.meansq.m")),...)
        printMethod("pr:",...)
        printMethod(hSummary(attr(x,"coef.pr.m")),...)
}
