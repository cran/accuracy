\name{modelsCompare}
\alias{LRE}
\alias{modelsCompare}
\alias{modelsAgree}
\alias{modelSummary}
\alias{modelBetas}
\title{Functions for comparing results: LRE, compare, agrees}
\description{
	Use these functions to compare a set of results from an analysis, to check if they 
	agree to a certain number of significant digits.
}
\usage{
	# compare numeric values against correct values
	LRE(x,correct, use.LAE=TRUE, cutoff=16)

	# compare the results from two possibly identical models
	modelsAgree(model, model2=NULL, digits=3, ...)
 	modelsCompare(model, model2=NULL, param.function=modelBetas,
				 similarity.function=LRE, summary.function=min)
                                                        
	# extract betas or entire coefficient summaries
	modelBetas(model)
	modelSummary(model)
}

\arguments{
  \item{x}{ vector of computed values}
  \item{correct}{ vector of correct values}
  \item{use.LAE}{use log absolute error if $c_i$==0. If false, returns NA where $c_i$==0}
  \item{cutoff}{In LRE results, convert INF to cutoff. Since perfect agreement yields an infinite LRE otherwise. Default value of 16 is based on the number of significant digits for doubles.}
  \item{model}{first model for comparison, or if model2=NULL, supply a list  of models for comparison}
  \item{model2}{second model for comparison, if }
  \item{digits}{number of digits to use for comparison}
  \item{param.function}{function used to extract model parameters for comparison}
  \item{similarity.function}{function used to compute similarity between sets of model parameters}
  \item{summary.function}{function used to summarize differences}
  \item{\dots}{parameters to pass from modelsAgree to modelsCompare}
} 

\value{ 
        Returns a new vector of log-relative-errors (log absolute error where $c_i$ ==0).
	The resulting values are roughly interpretable as the number of significant digits
	of agreement between c and x. Larger numbers indicate that x is a relatively
	more accurate approximation of c. Numbers less than or equal to zero indicate
	complete disagreement between x and c.
}

\details{
	modelsAgree is a convenienence function that calls modelsCompare to summarize similarities
between models results. 
}

\references{
Altman, M., J. Gill and M. P. McDonald.  2003.  \emph{Numerical Issues in Statistical
Computing for the Social Scientist}.  John Wiley \& Sons.
\url{http://www.hmdc.harvard.edu/numerical_issues/}
}

\author{
Micah Altman
\email{Micah\_Altman@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah_altman/}, Michael McDonald
}


\examples{

# simple LRE examples
LRE(1.001,1) # roughly 3 significant digits agreement
LRE(1,1) # complete agreement
LRE(20,1) # complete disagreement

# compare two models

if (is.R()) {
 hasNmle=require(nlme,quietly=TRUE) 
}  else {
  hasNmle=require(nlme3,quietly=TRUE)
}

if (hasNmle) { 
  O2<-Orthodont
  O3<-Orthodont
  
  O2[,2]<-O2[,2]+.1                                    
  O3[,2]<-O2[,2]+.11 
  if (is.R()) {                                   
    fm1<- lmList(Orthodont)
    fm2<- lmList(O2)
    fm3<- lmList(O2)
  } else {
    # workaround bug in lmlist in S-Plus
    fm1<- eval(substitute(lmList(Orthodont)))
    fm2<- eval(substitute(lmList(O2)))
    fm3<- eval(substitute(lmList(O2)))
  }

  # do the three models agree? 
  print(modelsAgree(list(fm1,fm2,fm2)))

  # show details of diagreement between first 2 models
  print(modelsCompare(fm1,fm2))

  #compare betas at 2 significant digits
  print(modelsAgree(fm1,fm2,digits=2,param.function=modelBetas))
  #compare betas at 1 significant digit
  print(modelsAgree(fm1,fm2,digits=1,param.function=modelBetas))
}

\dontshow{
  # unit test
if (is.R()){
  selfTest=get("modelsCompareSelfTest",envir=environment(LRE))
} else {
  selfTest=modelsCompareSelfTest
}
if (!selfTest()) {
  	warning("failed self-test")
}
}

}
\seealso{ See Also as \code{\link{ttst}}}

\keyword{misc}
\keyword{debugging}
