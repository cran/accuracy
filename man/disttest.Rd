\name{Distribution Test Data, and log relative error function}
\alias{chisqtst}
\alias{ftst}
\alias{gammatst}
\alias{normtst}
\alias{ttst}
\alias{LRE}
\title{ Benchmark data to test the accuracy of statistical distribution functions, function to compute log relative error }
\description{
	This is benchmark data, used to test the accuracy of statistical distribution functions. Also included is a function for computing log-relative error, a measure of accuracy.
}
\usage{
	data(ttst)
	data(ftst)
	data(gammatst)
	data(normtst)
	data(chisqtst)
	LRE(x,correct, use.LAE=TRUE)
}
\arguments{
  \item{x}{ vector of computed values}
  \item{correct}{ vector of correct values}
  \item{use.LAE}{use log absolute error if $c_i$==0. If false, returns NA where $c_i$==0}
} 

\value{ 
        Returns a new vector of log-relative-errors (log absolute error where $c_i$ ==0).
	The resulting values are roughly interpretable as the number of significant digits
	of agreement between c and x. Larger numbers indicate that x is a relatively
	more accurate approximation of c. Numbers less than or equal to zero indicate
	complete disagreement between x and c.
}

\details{
These tests values are tab-delimited ascii datasets, with a variable
header line. Typical variables are:

        P - probablilty value
        DF - degrees of freedom
        INVDIST - inv the inverse function value for P:
                i.e.: invt(P,DF) = INVT
        INVDIST - inv the inverse function value for P:
                i.e.: invt(P,DF) = INVT
        PINVDIST - value for cumulative distribrution of
                        INVDIST
                i.e.: cumt(INVT,DF) = PINVT

Note that because of limits in numerical precision, cumt(invt(P,DF)) is not always equal to P

Of course, not all possible values of P and DF can be listed.
The test values for P  were created from systematic and random samples in 
{ [1e-12,1-1e-12],0,1 }, and [0,1E5] respectively.

We use Kneusel's ELV program (1989) to calculate the values the cumulative and 
inverse distributions. The results are claimed by Knsusel to be accurate to 6 digits. We
checked these results using Brown's (1998) DCDFLIB library. The results
of both calculations agreed to the 6 digits supplied by ELV, but in a small
number of cases, DCDFSTAT's calculation at the 7th digit indicated that
the sixth digit would change if rounded. 

Note that both ELV and DCDFLIB can generate many more distributions than
are included here. Other resources are described in the references
}

\references{
Altman, M., J. Gill and M. P. McDonald.  2003.  \emph{Numerical Issues in Statistical
Computing for the Social Scientist}.  John Wiley \& Sons.
\url{http://www.hmdc.harvard.edu/numerical_issues/}

Altman, Micah and Michael McDonald, 2001. "Choosing Reliable Statistical Software." 
\emph{PS: Political Science \& Politics} 24(3): 681-8

Barry W. Brown, James Lovato, Kathy Russell, 1998,   DCDFLIB,
\url{ftp://odin.mda.uth.tmc.edu in pub/unix/dcdflib.c-1.0-tar.Z}

Kneusel, L. 1989. \emph{Computergesteutzte Berechnung statistischer Verteilungen.} Olde
nbourg, Meunchen-Wien.
\url{http://www.stat.uni-muenchen.de/~knuesel/elv/elv.html}
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

#
# how accurate are student's t-test functions?
#

data(ttst)
# compute t quantiles using benchmark data
tqt = qt(ttst$p,ttst$df)

# compute log-relative-error (LRE) of  qt() results, compared to 
# correct answers

lrq = LRE(tqt, ttst$invt);

# if there are entries with LRE's of < 5, there may be
# significant inaccuracies in the qt() function

table(trunc(lrq))

# now repeat process, for pt()

tpt = pt(ttst$invt,ttst$df)
lrp= LRE(tpt, ttst$pinvt);
table(trunc(lrp))


}
\keyword{misc}
\keyword{debugging}
